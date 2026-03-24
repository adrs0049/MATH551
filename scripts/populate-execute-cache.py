#!/usr/bin/env python3
"""
Populate the MyST execution cache (_build/execute/) from pre-executed .ipynb files.

This allows `jupyter-book build --execute --html` to skip re-execution of
notebooks that were already executed (e.g., in a separate CI job).

The cache key is an MD5 hash of the kernel name + JSON of all code cell
contents (matching mystmd's buildCacheKey function).
"""
import hashlib
import json
import os
import sys
from pathlib import Path


def build_cache_key(kernel_name: str, cells: list[dict]) -> str:
    """Reproduce mystmd's buildCacheKey function."""
    hashable_items = []
    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        hashable_items.append({
            "kind": "block",
            "content": source,
            "raisesException": False,
        })

    h = hashlib.md5()
    h.update(kernel_name.encode())
    h.update(json.dumps(hashable_items).encode())
    return h.digest().hex()


def extract_outputs(cells: list[dict]) -> list:
    """Extract outputs from cells in the format MyST expects."""
    results = []
    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        results.append(cell.get("outputs", []))
    return results


def populate_cache(notebook_path: str, cache_dir: str) -> None:
    """Read an executed notebook and write its outputs to the MyST cache."""
    with open(notebook_path) as f:
        nb = json.load(f)

    # Get kernel name
    kernel_name = nb.get("metadata", {}).get("kernelspec", {}).get("name", "python3")
    cells = nb.get("cells", [])

    # Check if notebook has any code cells with outputs
    code_cells = [c for c in cells if c["cell_type"] == "code"]
    if not code_cells:
        return

    has_outputs = any(c.get("outputs") for c in code_cells)
    if not has_outputs:
        print(f"  SKIP {notebook_path} (no outputs)")
        return

    # Build cache key and extract outputs
    cache_key = build_cache_key(kernel_name, cells)
    outputs = extract_outputs(cells)

    # Write to cache
    cache_file = os.path.join(cache_dir, f"{cache_key}.json")
    with open(cache_file, "w") as f:
        json.dump(outputs, f)

    print(f"  OK   {notebook_path} -> {cache_key}.json")


def main():
    cache_dir = "_build/execute"
    os.makedirs(cache_dir, exist_ok=True)

    # Find notebooks to process
    if len(sys.argv) > 1:
        notebooks = sys.argv[1:]
    else:
        notebooks = sorted(Path("notebooks").glob("*.ipynb"))

    print(f"Populating MyST execution cache for {len(notebooks)} notebooks")
    print(f"Cache directory: {cache_dir}")
    print("=" * 50)

    for nb in notebooks:
        populate_cache(str(nb), cache_dir)


if __name__ == "__main__":
    main()
