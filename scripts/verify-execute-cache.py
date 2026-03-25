#!/usr/bin/env python3
"""Verify that the MyST execution cache contains entries for all notebooks."""
import hashlib
import json
import os
import sys


def build_cache_key(kernel_name, cells):
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
    h.update(json.dumps(hashable_items, separators=(",", ":"), ensure_ascii=False).encode())
    return h.hexdigest()


def main():
    cache_dir = "_build/execute"
    nb_dir = "notebooks"

    if not os.path.isdir(cache_dir):
        print(f"ERROR: Cache directory {cache_dir} does not exist")
        sys.exit(1)

    cache_files = set(os.listdir(cache_dir))
    print(f"Cache directory: {cache_dir} ({len(cache_files)} files)")
    print("=" * 60)

    matched = 0
    missing = 0

    for nb_file in sorted(os.listdir(nb_dir)):
        if not nb_file.endswith(".ipynb"):
            continue
        path = os.path.join(nb_dir, nb_file)
        with open(path) as f:
            nb = json.load(f)

        kernel = nb.get("metadata", {}).get("kernelspec", {}).get("name", "python3")
        cells = nb.get("cells", [])
        cache_key = build_cache_key(kernel, cells)
        found = f"{cache_key}.json" in cache_files

        status = "OK" if found else "MISSING"
        print(f"  {status:7s} {nb_file} -> {cache_key}.json")

        if found:
            matched += 1
        else:
            missing += 1

    print("=" * 60)
    print(f"Results: {matched} cached, {missing} missing")

    if missing > 0:
        print("WARNING: Some notebooks not cached (may be re-executed during build)")


if __name__ == "__main__":
    main()
