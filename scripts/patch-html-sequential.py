#!/usr/bin/env python3
"""
Patch mystmd/jupyter-book buildHtml to render pages sequentially.

The default Promise.all + p-limit(5) causes the Remix server to accumulate
memory, which OOMs on 16GB GitHub runners. This replaces the parallel fetch
with a sequential for-loop.
"""
import sys


def patch(cjs_path):
    with open(cjs_path, 'r') as f:
        content = f.read()

    if 'PATCHED_SEQUENTIAL_HTML' in content:
        print(f"  OK: already patched")
        return True

    # Find the start of the Promise.all block
    marker = 'await Promise.all(routes.map(async (route) => limitConnections(async () => {'
    start = content.find(marker)
    if start == -1:
        print(f"  WARNING: could not find Promise.all block in {cjs_path}")
        return False

    # Walk back to include the leading whitespace/await
    line_start = content.rfind('\n', 0, start) + 1

    # Find the closing })));
    search_from = start + len(marker)
    end_marker = '})));'
    end = content.find(end_marker, search_from)
    if end == -1:
        print(f"  WARNING: could not find closing }}))) in {cjs_path}")
        return False
    end += len(end_marker)

    # Extract the inner block (the route processing logic)
    inner_start = start + len(marker) + 1  # skip the opening {
    inner = content[inner_start:end - len(end_marker)]

    # Detect the resolve variable name (resolve9, resolve10, etc.)
    import re
    resolve_match = re.search(r'new Promise\((resolve\d+)\)', inner)
    resolve_var = resolve_match.group(1) if resolve_match else 'resolve'

    # Build the replacement: sequential for-loop
    replacement = (
        '  // PATCHED_SEQUENTIAL_HTML: render pages one at a time to reduce memory\n'
        '  for (const route of routes) {\n'
        + inner.replace('return;', 'continue;')
        + '  }'
    )

    content = content[:line_start] + replacement + content[end:]

    with open(cjs_path, 'w') as f:
        f.write(content)

    print(f"  PATCHED: sequential HTML rendering enabled")
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: patch-html-sequential.py <path-to-cjs>")
        sys.exit(1)
    success = patch(sys.argv[1])
    sys.exit(0 if success else 1)
