#!/usr/bin/env python3
"""Check that every WDL command block writes the required execution fields."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


REQUIRED_FIELDS = (
    "stage=",
    "start_time=",
    "completion_time=",
    "dimensions=",
    "outputs=",
)
COMMAND_BLOCK = re.compile(r"command\s*<<<(.*?)>>>", re.DOTALL)


def wdl_files(paths: list[Path]) -> list[Path]:
    """Return explicit WDL files and WDL files below the requested directories."""
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            files.extend(sorted(path.rglob("*.wdl")))
        elif path.suffix == ".wdl":
            files.append(path)
        else:
            raise ValueError(f"Expected a WDL file or directory: {path}")
    return files


def missing_logging_fields(path: Path) -> list[str]:
    """Return errors for WDL command blocks that omit required log fields."""
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    for index, block in enumerate(COMMAND_BLOCK.findall(text), start=1):
        missing = [field for field in REQUIRED_FIELDS if field not in block]
        if missing:
            errors.append(
                f"{path}: command block {index} is missing "
                f"{', '.join(missing)}"
            )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Require standard logging fields in WDL command blocks."
    )
    parser.add_argument("paths", nargs="+", type=Path)
    arguments = parser.parse_args()

    try:
        files = wdl_files(arguments.paths)
    except ValueError as error:
        print(error, file=sys.stderr)
        return 2

    errors = [error for path in files for error in missing_logging_fields(path)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1

    print("All WDL command blocks contain required logging fields.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
