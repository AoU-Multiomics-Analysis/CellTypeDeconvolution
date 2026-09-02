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
COMMAND_START = re.compile(r"\bcommand\s*(<<<|\{)")


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


def matching_brace(text: str, start: int) -> int | None:
    """Return the closing brace for a brace-form command block."""
    depth = 1
    quote: str | None = None
    index = start

    while index < len(text):
        character = text[index]
        if quote is not None:
            if quote == '"' and character == "\\":
                index += 2
                continue
            if character == quote:
                quote = None
        elif character in ("'", '"'):
            quote = character
        elif character == "#":
            newline = text.find("\n", index)
            index = len(text) if newline == -1 else newline
            continue
        elif character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return index
        index += 1

    return None


def command_blocks(text: str) -> list[str]:
    """Extract heredoc and brace-form WDL command blocks without truncation."""
    blocks: list[str] = []
    search_start = 0

    while match := COMMAND_START.search(text, search_start):
        delimiter = match.group(1)
        block_start = match.end()
        if delimiter == "<<<":
            block_end = text.find(">>>", block_start)
        else:
            block_end = matching_brace(text, block_start)
        if block_end is None or block_end == -1:
            raise ValueError(
                f"unterminated {delimiter}-form command block at offset "
                f"{match.start()}"
            )
        blocks.append(text[block_start:block_end])
        search_start = block_end + (3 if delimiter == "<<<" else 1)

    return blocks


def missing_logging_fields(path: Path) -> list[str]:
    """Return errors for WDL command blocks that omit required log fields."""
    text = path.read_text(encoding="utf-8")
    try:
        blocks = command_blocks(text)
    except ValueError as error:
        return [f"{path}: {error}"]

    errors: list[str] = []
    for index, block in enumerate(blocks, start=1):
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
