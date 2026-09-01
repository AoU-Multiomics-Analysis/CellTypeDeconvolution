#!/usr/bin/env python3
"""Update the generated GitHub Actions badge block in README.md."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path


README = Path("README.md")
WORKFLOW_DIR = Path(".github/workflows")
START = "<!-- workflow-badges:start -->"
END = "<!-- workflow-badges:end -->"


def workflow_name(path: Path) -> str:
    for line in path.read_text().splitlines():
        match = re.match(r"^name:\s*(.+?)\s*$", line)
        if match:
            return match.group(1).strip("\"'")
    return path.stem.replace("-", " ").title()


def repository_slug() -> str:
    if os.environ.get("GITHUB_REPOSITORY"):
        return os.environ["GITHUB_REPOSITORY"]

    remote = subprocess.check_output(
        ["git", "remote", "get-url", "origin"],
        text=True,
    ).strip()

    if remote.endswith(".git"):
        remote = remote[:-4]

    if remote.startswith("git@github.com:"):
        return remote.removeprefix("git@github.com:")

    match = re.match(r"https://github.com/(.+)$", remote)
    if match:
        return match.group(1)

    raise RuntimeError(f"Cannot infer GitHub repository from origin remote: {remote}")


def badge_block() -> str:
    repo = repository_slug()
    workflows = sorted(WORKFLOW_DIR.glob("*.yml")) + sorted(WORKFLOW_DIR.glob("*.yaml"))

    badges = []
    for workflow in workflows:
        name = workflow_name(workflow)
        filename = workflow.name
        url = f"https://github.com/{repo}/actions/workflows/{filename}"
        badges.append(f"[![{name}]({url}/badge.svg)]({url})")

    return "\n".join([START, *badges, END])


def main() -> None:
    readme = README.read_text()
    replacement = badge_block()

    if START in readme and END in readme:
        pattern = re.compile(f"{re.escape(START)}.*?{re.escape(END)}", re.DOTALL)
        updated = pattern.sub(replacement, readme)
    else:
        updated = re.sub(r"^(# .+\n)", rf"\1\n{replacement}\n", readme, count=1)

    if updated != readme:
        README.write_text(updated)


if __name__ == "__main__":
    main()
