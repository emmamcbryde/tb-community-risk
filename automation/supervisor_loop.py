from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from openai import OpenAI


SUPERVISOR_INSTRUCTIONS = """
You are supervising a Codex coding agent working on a Python repository.

Your job:
1. Give Codex one small, concrete task at a time.
2. After each Codex result, inspect Codex output, git diff, git status, and test output.
3. Continue only when the next step is obvious and low-risk.
4. Pause for the human when there is an unclear requirement, a scientific/modelling judgement,
   an architectural decision, a destructive change, a change outside scope, or repeated test failure.
5. Prefer minimal, reversible changes.
6. Do not ask Codex to access secrets, API keys, credentials, or files outside the repository.

Respond using exactly one of these formats:

CODEX_TASK:
<one precise instruction for Codex>

ASK_HUMAN:
<short question or decision needed>

DONE:
<short summary of what was completed and what remains>
"""


def tail(text: str, limit: int) -> str:
    if not text:
        return ""

    text = text.replace("\r\n", "\n")

    if len(text) <= limit:
        return text

    return f"...[truncated to last {limit} characters]...\n{text[-limit:]}"


def run_command(
    command,
    cwd: Path,
    timeout: int = 1200,
    shell: bool = False,
) -> tuple[int, str, str]:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        timeout=timeout,
        shell=shell,
    )
    return result.returncode, result.stdout, result.stderr


def get_git_status(repo: Path) -> str:
    _, stdout, stderr = run_command(["git", "status", "--short"], cwd=repo)
    return stdout.strip() or stderr.strip()


def ask_supervisor(
    client: OpenAI,
    model: str,
    message: str,
    previous_response_id: str | None,
) -> tuple[str, str]:
    kwargs = {
        "model": model,
        "instructions": SUPERVISOR_INSTRUCTIONS,
        "input": message,
        "store": True,
    }

    if previous_response_id is not None:
        kwargs["previous_response_id"] = previous_response_id

    response = client.responses.create(**kwargs)
    return response.output_text.strip(), response.id


def run_codex(
    repo: Path,
    task: str,
    sandbox: str,
    test_command: str,
    skip_tests: bool,
) -> str:
    codex_executable = shutil.which("codex") or shutil.which("codex.cmd")

    if codex_executable is None:
        raise RuntimeError(
            "Could not find the 'codex' command on PATH. "
            "Try running 'where.exe codex' in PowerShell."
        )

    codex_command = [
        codex_executable,
        "exec",
        "--sandbox",
        sandbox,
        task,
    ]

    codex_code, codex_stdout, codex_stderr = run_command(
        codex_command,
        cwd=repo,
        timeout=1800,
    )

    _, git_status, git_status_err = run_command(
        ["git", "status", "--short"],
        cwd=repo,
    )

    _, diff_stat, diff_stat_err = run_command(
        ["git", "diff", "--stat"],
        cwd=repo,
    )

    _, diff_full, diff_full_err = run_command(
        ["git", "diff"],
        cwd=repo,
    )

    if skip_tests:
        test_code = 0
        test_stdout = "Tests skipped by --skip-tests."
        test_stderr = ""
    else:
        test_code, test_stdout, test_stderr = run_command(
            test_command,
            cwd=repo,
            timeout=1200,
            shell=True,
        )

    return f"""
CODEX_COMMAND:
{' '.join(codex_command[:-1])} <task>

CODEX_EXIT_CODE:
{codex_code}

CODEX_STDOUT:
{tail(codex_stdout, 6000)}

CODEX_STDERR:
{tail(codex_stderr, 6000)}

GIT_STATUS:
{tail(git_status or git_status_err, 4000)}

GIT_DIFF_STAT:
{tail(diff_stat or diff_stat_err, 4000)}

GIT_DIFF:
{tail(diff_full or diff_full_err, 14000)}

TEST_COMMAND:
{test_command if not skip_tests else 'SKIPPED'}

TEST_EXIT_CODE:
{test_code}

TEST_STDOUT:
{tail(test_stdout, 8000)}

TEST_STDERR:
{tail(test_stderr, 8000)}
""".strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an OpenAI supervisor loop that delegates coding tasks to Codex."
    )

    parser.add_argument(
        "repo",
        nargs="?",
        default=".",
        help="Path to the Git repository. Default: current directory.",
    )

    parser.add_argument(
        "--model",
        default=os.environ.get("OPENAI_SUPERVISOR_MODEL", "gpt-5.5"),
        help="OpenAI model to use for the supervisor.",
    )

    parser.add_argument(
        "--sandbox",
        choices=["read-only", "workspace-write"],
        default="workspace-write",
        help="Codex sandbox mode. Default: workspace-write.",
    )

    parser.add_argument(
        "--test-command",
        default=f'"{sys.executable}" -m pytest -q',
        help="Command to run after each Codex turn.",
    )

    parser.add_argument(
        "--max-turns",
        type=int,
        default=4,
        help="Maximum number of Codex turns before stopping.",
    )

    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Allow starting even if the Git working tree has uncommitted changes.",
    )

    parser.add_argument(
        "--skip-tests",
        action="store_true",
        help="Skip the test command after Codex runs.",
    )

    parser.add_argument(
        "--yes",
        action="store_true",
        help="Skip the initial confirmation prompt.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    repo = Path(args.repo).expanduser().resolve()

    if not repo.exists():
        raise SystemExit(f"Repository path does not exist: {repo}")

    if not (repo / ".git").exists():
        raise SystemExit(f"Not a Git repository: {repo}")

    if not os.environ.get("OPENAI_API_KEY"):
        raise SystemExit(
            "OPENAI_API_KEY is not set in this terminal. "
            "Set it before running this script."
        )

    if (shutil.which("codex") or shutil.which("codex.cmd")) is None:
        raise SystemExit(
            "Could not find the 'codex' command. "
            "Install Codex CLI or check that it is on PATH."
        )

    git_status = get_git_status(repo)

    if git_status and not args.allow_dirty:
        print("Git working tree is not clean:")
        print(git_status)
        raise SystemExit(
            "\nCommit or stash these changes first, or rerun with --allow-dirty."
        )

    print(f"Repository: {repo}")
    print(f"Supervisor model: {args.model}")
    print(f"Codex sandbox: {args.sandbox}")
    print(f"Test command: {'SKIPPED' if args.skip_tests else args.test_command}")
    print(f"Maximum Codex turns: {args.max_turns}")

    if not args.yes:
        answer = input(
            "\nType YES to allow this script to run Codex in the selected sandbox: "
        ).strip()

        if answer != "YES":
            raise SystemExit("Stopped before running Codex.")

    goal = input("\nGoal for this coding loop:\n> ").strip()

    if not goal:
        raise SystemExit("No goal provided.")

    client = OpenAI()
    previous_response_id: str | None = None

    message = f"""
Repository path:
{repo}

Human goal:
{goal}

Start by deciding the first small Codex task.
""".strip()

    for turn_number in range(1, args.max_turns + 1):
        print(f"\n=== SUPERVISOR TURN {turn_number} ===")

        supervisor_reply, previous_response_id = ask_supervisor(
            client=client,
            model=args.model,
            message=message,
            previous_response_id=previous_response_id,
        )

        print(supervisor_reply)

        if supervisor_reply.startswith("ASK_HUMAN:"):
            guidance = input("\nYour guidance, or type quit to stop:\n> ").strip()

            if guidance.lower() in {"quit", "exit", "stop"}:
                print("Stopped by human.")
                break

            message = f"""
Human guidance:
{guidance}

Continue supervising the Codex workflow.
""".strip()
            continue

        if supervisor_reply.startswith("DONE:"):
            print("\nSupervisor marked the task done.")
            break

        if supervisor_reply.startswith("CODEX_TASK:"):
            task = supervisor_reply.split("CODEX_TASK:", 1)[1].strip()

            if not task:
                message = "The previous CODEX_TASK was empty. Please provide a concrete task."
                continue

            print("\n=== RUNNING CODEX TASK ===")
            print(task)

            codex_report = run_codex(
                repo=repo,
                task=task,
                sandbox=args.sandbox,
                test_command=args.test_command,
                skip_tests=args.skip_tests,
            )

            print("\n=== CODEX REPORT ===")
            print(codex_report)

            message = f"""
Codex has completed the task. Review the result and decide the next action.

Task given to Codex:
{task}

Result:
{codex_report}
""".strip()
            continue

        message = f"""
Your previous response did not use one of the required formats.

Previous response:
{supervisor_reply}

Please respond with exactly one of:
CODEX_TASK:
ASK_HUMAN:
DONE:
""".strip()

    else:
        print("\nReached maximum Codex turns.")

    print("\n=== FINAL GIT STATUS ===")
    final_status = get_git_status(repo)
    print(final_status or "Clean working tree.")


if __name__ == "__main__":
    main()
