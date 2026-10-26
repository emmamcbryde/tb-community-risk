"""Live Streamlit protocol/render audit, with cold/warm deterministic timings.

Run from the repository root. Generated evidence stays in abm/output/milestone3f.
This complements AppTest when no interactive browser surface is available.
"""
from __future__ import annotations

import argparse
from io import BytesIO
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter, sleep

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--streamlit-target", type=Path)
    parser.add_argument("--port", type=int, default=8520)
    parser.add_argument("--frozen", action="store_true")
    parser.add_argument("--apptest", action="store_true")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--test-filter")
    args = parser.parse_args()
    if args.streamlit_target:
        sys.path.insert(0, str(args.streamlit_target.resolve()))
    if args.apptest:
        import pytest
        pytest_args = [
            "tests/test_explicit_tbi_ui_isolation.py",
            "tests/test_explicit_tbi_streamlit_app.py",
            "tests/test_frozen_reference_loader.py",
            "-q", "--basetemp=abm/output/milestone3f/pytest-supported",
        ]
        if args.test_filter:
            pytest_args.extend(["-k", args.test_filter])
        raise SystemExit(pytest.main(pytest_args))
    if args.profile:
        import cProfile
        import pstats
        from app.explicit_tbi_workflow import default_workflow_state, run_explicit_workflow_analysis
        destination = ROOT / "abm/output/milestone3f/deterministic-profile.prof"
        profiler = cProfile.Profile()
        result = profiler.runcall(run_explicit_workflow_analysis, default_workflow_state())
        profiler.dump_stats(str(destination))
        pstats.Stats(profiler).strip_dirs().sort_stats("tottime").print_stats(20)
        print(json.dumps(result["timings"], indent=2))
        return
    import pyarrow.ipc as ipc
    import streamlit
    import websocket
    from streamlit.proto.BackMsg_pb2 import BackMsg
    from streamlit.proto.ForwardMsg_pb2 import ForwardMsg

    output = ROOT / "abm/output/milestone3f"
    output.mkdir(parents=True, exist_ok=True)
    name = ("frozen" if args.frozen else "isolated") + "-" + streamlit.__version__
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(sys.path)
    entry = "streamlit_app.py" if args.frozen else "explicit_tbi_ui/app.py"
    server_log = (output / f"{name}-server.log").open("w", encoding="utf-8")
    process = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", entry,
         "--server.headless", "true", "--server.address", "localhost",
         "--server.port", str(args.port), "--browser.gatherUsageStats", "false",
         "--global.developmentMode", "false"],
        cwd=ROOT, env=env, stdout=server_log, stderr=subprocess.STDOUT,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
    )
    ws = None
    try:
        deadline = perf_counter() + 60
        while perf_counter() < deadline:
            try:
                ws = websocket.create_connection(
                    f"ws://localhost:{args.port}/_stcore/stream", timeout=180,
                    subprotocols=["streamlit"], origin=f"http://localhost:{args.port}",
                )
                break
            except (OSError, websocket.WebSocketException):
                if process.poll() is not None:
                    raise RuntimeError(f"Server exited; inspect {server_log.name}")
                sleep(0.25)
        if ws is None:
            raise RuntimeError("Streamlit server did not start")

        def render(*, workflow_id=None, page=0, click_id=None, route=None):
            request = BackMsg()
            request.rerun_script.query_string = ""
            if route:
                request.rerun_script.page_name = route
            if workflow_id:
                widget = request.rerun_script.widget_states.widgets.add()
                widget.id = workflow_id
                widget.int_value = page
            if click_id:
                widget = request.rerun_script.widget_states.widgets.add()
                widget.id = click_id
                widget.trigger_value = True
            start = perf_counter()
            ws.send_binary(request.SerializeToString())
            record = {"serverPages": [], "navigationPages": [], "headings": [], "radios": {}, "buttons": {}, "json": [], "tables": [], "errors": [], "captions": [], "numberInputs": [], "pageNotFound": False}
            while True:
                message = ForwardMsg()
                message.ParseFromString(ws.recv())
                kind = message.WhichOneof("type")
                if kind == "new_session":
                    record["serverPages"] = [p.page_name for p in message.new_session.app_pages]
                elif kind == "navigation":
                    record["navigationPages"] = [p.page_name for p in message.navigation.app_pages]
                elif kind == "page_not_found":
                    record["pageNotFound"] = True
                elif kind == "delta" and message.delta.WhichOneof("type") == "new_element":
                    element = message.delta.new_element
                    element_type = element.WhichOneof("type")
                    if element_type == "heading":
                        record["headings"].append(element.heading.body)
                    elif element_type == "radio":
                        record["radios"][element.radio.label] = {"id": element.radio.id, "options": list(element.radio.options)}
                    elif element_type == "button":
                        record["buttons"][element.button.label] = element.button.id
                    elif element_type == "json":
                        record["json"].append(json.loads(element.json.body))
                    elif element_type == "arrow_data_frame":
                        frame = ipc.open_stream(BytesIO(element.arrow_data_frame.data)).read_all().to_pandas()
                        record["tables"].append(frame.to_dict(orient="records"))
                    elif element_type == "exception":
                        record["errors"].append(element.exception.message)
                    elif element_type == "markdown":
                        record["captions"].append(element.markdown.body)
                    elif element_type == "number_input":
                        record["numberInputs"].append({"label": element.number_input.label, "format": element.number_input.format, "default": element.number_input.default})
                elif kind == "script_finished":
                    break
            record["wallSeconds"] = perf_counter() - start
            assert not record["errors"], record["errors"]
            return record

        report = {"streamlitVersion": streamlit.__version__, "entry": entry}
        report["initial"] = initial = render()
        print(f"Initial render: {initial['wallSeconds']:.3f}s; server pages: {initial['serverPages']}; navigation: {initial['navigationPages']}", flush=True)
        if args.frozen:
            expected = ["Set up", "Run Analysis", "Results", "Health Economics", "Explore Decisions", "Evidence & Assumptions", "Technical settings", "Dynamic transmission model", "Model comparison"]
            assert initial["navigationPages"] == expected, initial["navigationPages"]
            report["results"] = render(route="Results")
            assert not report["results"]["pageNotFound"]
            assert report["results"]["tables"], "Frozen Results must render data tables"
        else:
            expected = ["Set up", "Run analysis", "Results", "Evidence and technical information"]
            assert initial["serverPages"] == ["app"], initial["serverPages"]
            assert initial["radios"]["Workflow"]["options"] == expected
            assert "Input known TBI prevalence" in initial["headings"]
            current = render(workflow_id=initial["radios"]["Workflow"]["id"], page=1)
            for index, temperature in enumerate(("cold", "warm"), start=1):
                run = render(workflow_id=current["radios"]["Workflow"]["id"], page=1, click_id=current["buttons"]["Run explicit TBI analysis"])
                evidence = render(workflow_id=run["radios"]["Workflow"]["id"], page=3)
                diagnostics = next(row for row in evidence["json"] if "analysisTimingsSeconds" in row)
                assert diagnostics["analysisDispatchCount"] == index, diagnostics
                report[temperature] = {"clickWallSeconds": run["wallSeconds"], **diagnostics}
                print(f"{temperature}: {json.dumps(report[temperature])}", flush=True)
                results = render(workflow_id=evidence["radios"]["Workflow"]["id"], page=2)
                report[temperature + "Results"] = results
                current = render(workflow_id=results["radios"]["Workflow"]["id"], page=1)
            report["legacyRouteAttempt"] = route = render(route="Economics")
            assert route["pageNotFound"]
            assert not any("Economics" in heading for heading in route["headings"])
        destination = output / f"{name}-audit.json"
        destination.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"PASS: {destination}", flush=True)
    finally:
        if ws:
            ws.close()
        process.terminate()
        process.wait(timeout=15)
        server_log.close()


if __name__ == "__main__":
    main()
