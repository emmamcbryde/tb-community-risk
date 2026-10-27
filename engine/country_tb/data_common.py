"""Offline data-layer primitives. No model, UI, or network imports."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal, InvalidOperation
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
from uuid import uuid4

LOADER_VERSION = "country_tb_data_loader_1.0.0"
PROFILE_VERSION = "1.1.0"
REPO_ROOT = Path(__file__).resolve().parents[2]


class DataValidationError(ValueError):
    pass


def decimal(value, *, nonnegative=True):
    if isinstance(value, bool) or value is None:
        raise DataValidationError("Expected a finite number, not missing/boolean.")
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise DataValidationError(f"Invalid number: {value!r}") from exc
    if not result.is_finite() or (nonnegative and result < 0):
        raise DataValidationError(f"Non-finite or negative value: {value!r}")
    return result


def numeric_text(value):
    result = decimal(value, nonnegative=False)
    if result == 0:
        return "0"
    text = format(result, "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def canonical(payload):
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")


def scientific_payload(payload):
    """Presentation is explicit. Only identified set-valued arrays are sorted."""
    unordered = {"rows", "source_files", "warnings", "candidates", "inactive_alternatives", "published_both_checks"}
    def visit(value, key=""):
        if isinstance(value, dict):
            return {k: visit(v, k) for k, v in value.items() if k not in {"display", "profile_hash", "receipt"}}
        if isinstance(value, list):
            values = [visit(v) for v in value]
            return sorted(values, key=canonical) if key in unordered else values
        return value
    return visit(payload)


def profile_hash(payload):
    return digest(canonical(scientific_payload(payload)))


def seal(payload):
    result = deepcopy(payload)
    result["profile_hash"] = profile_hash(result)
    return result


def file_bytes(path, *, normalize_lf=False):
    value = Path(path).read_bytes()
    return value.replace(b"\r\n", b"\n") if normalize_lf else value


def csv_rows(payload, *, compressed=False):
    raw = gzip.decompress(payload) if compressed else payload
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""))
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise DataValidationError("Missing or duplicate CSV columns.")
    rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise DataValidationError("Malformed CSV row.")
    return tuple(reader.fieldnames), rows


def csv_bytes(rows, columns):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def require_hash(value):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise DataValidationError("Invalid SHA-256.")
    return value


def contained_file(directory, filename):
    directory = Path(directory).resolve()
    path = (directory / filename).resolve()
    if path.parent != directory or not path.is_file():
        raise DataValidationError(f"Missing or unsafe snapshot file: {filename}")
    return path


def publish_immutable(directory, snapshot_id, files):
    """Atomic publish of a NEW snapshot directory; identical reimport is a no-op.

    Never appends, replaces an existing snapshot or deletes source files.
    Caller validates data before invoking. Temporary drafts are retained on error.
    """
    directory = Path(directory)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", snapshot_id):
        raise DataValidationError("Unsafe snapshot identifier.")
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / snapshot_id
    if any(Path(name).name != name for name in files):
        raise DataValidationError("Unsafe snapshot filename.")
    if target.exists():
        existing = {p.name: p.read_bytes() for p in target.iterdir() if p.is_file()}
        if existing != files:
            raise DataValidationError("Immutable snapshot already exists with different content.")
        return target
    # Normal inherited repository permissions, including on Windows. mkdtemp's
    # private ACL would survive rename and make a published dataset unreadable
    # to other workspace processes.
    draft = directory / (".draft-" + uuid4().hex)
    draft.mkdir(exist_ok=False)
    for name, payload in sorted(files.items()):
        with (draft / name).open("xb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
    try:
        draft.rename(target)
    except FileExistsError:
        if {p.name: p.read_bytes() for p in target.iterdir() if p.is_file()} != files:
            raise DataValidationError("Concurrent immutable import conflict.")
    return target
