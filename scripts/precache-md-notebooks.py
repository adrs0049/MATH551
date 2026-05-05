#!/usr/bin/env python3
"""
Pre-execute .md notebooks (MyST markdown files with {code-cell} blocks) and
populate the mystmd execution cache (_build/execute/) with their outputs.

This avoids mystmd's flaky kernel manager during the HTML build:
its executionSemaphore is hardcoded to cpus()-1, --execute-parallel is
silently ignored, and one hung kernel deadlocks the whole queue. By
seeding the cache here (using nbclient directly, one notebook at a time,
with a per-notebook timeout), the HTML-build job finds every notebook
already cached and never needs to spin up a kernel.

Usage:
    precache-md-notebooks.py [--timeout SECONDS] [<file.md> ...]

If no files are given, every .md file with a {code-cell} block under the
project source tree (excluding _build/, node_modules/, myst-theme/) is
processed. Skips files whose code cells are already in the mystmd cache.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError


CACHE_DIR = Path("_build/execute")

# Matches a fenced ``` {code-cell} <lang> ``` block, including any options
# block that follows the opening fence. We only need the code body.
CODE_CELL_RE = re.compile(
    r"^```\{code-cell\}[^\n]*\n(.*?)^```\s*$",
    re.DOTALL | re.MULTILINE,
)

# Frontmatter is the first --- ... --- block in the file.
FRONTMATTER_RE = re.compile(r"\A---\s*\n(.*?)\n---\s*\n", re.DOTALL)


def parse_md(path: Path) -> tuple[str, list[str]]:
    """Return (kernel_name, [code-cell sources]) for a .md file."""
    text = path.read_text()

    # Pull kernel name out of the frontmatter (default: python3).
    kernel_name = "python3"
    fm_match = FRONTMATTER_RE.match(text)
    if fm_match:
        fm = fm_match.group(1)
        # Cheap YAML scan — kernelspec.name on its own line, two-space indent.
        m = re.search(r"^\s*kernelspec\s*:\s*\n((?:[ \t]+.*\n?)+)", fm, re.MULTILINE)
        if m:
            n = re.search(r"^\s+name\s*:\s*(\S+)\s*$", m.group(1), re.MULTILINE)
            if n:
                kernel_name = n.group(1)

    # Strip leading-cell-options block from each code cell's body.
    sources: list[str] = []
    for m in CODE_CELL_RE.finditer(text):
        body = m.group(1)
        # mystmd treats a leading `:::` -less options block (lines starting
        # with `:`) as cell metadata. Drop those + the blank line after.
        lines = body.split("\n")
        i = 0
        while i < len(lines) and lines[i].startswith(":"):
            i += 1
        if i < len(lines) and lines[i] == "":
            i += 1
        sources.append("\n".join(lines[i:]).rstrip("\n"))

    return kernel_name, sources


def build_cache_key(kernel_name: str, sources: list[str]) -> str:
    """Replicate mystmd's buildCacheKey for the unmodified-cells case."""
    items = [
        {"kind": "block", "content": s, "raisesException": False}
        for s in sources
    ]
    h = hashlib.md5()
    h.update(kernel_name.encode())
    h.update(json.dumps(items, separators=(",", ":"), ensure_ascii=False).encode())
    return h.hexdigest()


def execute_in_memory(kernel_name: str, sources: list[str], timeout: int) -> list[list]:
    """Execute the cells with nbclient and return per-cell outputs."""
    nb = nbformat.v4.new_notebook()
    nb.metadata["kernelspec"] = {"name": kernel_name, "display_name": kernel_name}
    nb.cells = [nbformat.v4.new_code_cell(src) for src in sources]

    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name=kernel_name,
        allow_errors=False,
    )
    client.execute()
    return [cell.get("outputs", []) for cell in nb.cells]


def has_code_cells(path: Path) -> bool:
    return bool(CODE_CELL_RE.search(path.read_text()))


def discover_md_files() -> list[Path]:
    skip_prefixes = ("_build/", "node_modules/", "myst-theme/", ".git/")
    out: list[Path] = []
    for p in Path(".").rglob("*.md"):
        s = str(p)
        if any(s.startswith(prefix) for prefix in skip_prefixes):
            continue
        if has_code_cells(p):
            out.append(p)
    return sorted(out)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=300, help="seconds per notebook")
    parser.add_argument("files", nargs="*", help="explicit .md files to process")
    args = parser.parse_args()

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    files = [Path(f) for f in args.files] if args.files else discover_md_files()
    if not files:
        print("No .md files with {code-cell} blocks found.")
        return 0

    print(f"Pre-executing {len(files)} .md notebook(s) into {CACHE_DIR}/")
    print("=" * 60)

    failed: list[tuple[Path, str]] = []
    for md in files:
        kernel, sources = parse_md(md)
        if not sources:
            print(f"  SKIP {md} (no code cells)")
            continue

        key = build_cache_key(kernel, sources)
        cache_file = CACHE_DIR / f"{key}.json"
        if cache_file.exists():
            print(f"  HIT  {md} (already in cache)")
            continue

        print(f"  RUN  {md} ({len(sources)} cells, kernel={kernel}) ... ", end="", flush=True)
        try:
            outputs = execute_in_memory(kernel, sources, args.timeout)
        except CellExecutionError as e:
            print(f"FAILED (cell error)")
            failed.append((md, str(e).splitlines()[0]))
            continue
        except Exception as e:  # nbclient timeouts, kernel crashes, etc.
            print(f"FAILED ({type(e).__name__}: {e})")
            failed.append((md, f"{type(e).__name__}: {e}"))
            continue

        cache_file.write_text(json.dumps(outputs, separators=(",", ":")))
        print(f"OK -> {key[:12]}.json")

    if failed:
        print()
        print(f"::warning::{len(failed)} notebook(s) failed to pre-execute; mystmd will retry live:")
        for md, err in failed:
            print(f"  {md}: {err}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
