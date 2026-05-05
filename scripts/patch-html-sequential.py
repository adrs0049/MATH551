#!/usr/bin/env python3
"""
Patch mystmd/jupyter-book buildHtml to render pages sequentially.

The default Promise.all + p-limit(5) hits the Remix server with 5 concurrent
requests, which raises peak memory on CI runners. Sequential rendering keeps
peak lower at a small wall-clock cost.

(Earlier versions of this patch also restarted the Remix server every N
pages to free accumulated memory. That stopped working with mystmd 1.9.0:
Remix v1.17 takes ~30s to rebuild on restart, longer than the fetch retry
window, so post-restart fetches fail and the post-fetch copySync step never
runs, producing an artifact with no _build/html/build/_assets/ directory.)
"""
import sys


def patch(cjs_path):
    with open(cjs_path, 'r') as f:
        content = f.read()

    if 'PATCHED_SEQUENTIAL_HTML' in content:
        print(f"  OK: already patched")
        return True

    marker = 'await Promise.all(routes.map(async (route) => limitConnections(async () => {'
    start = content.find(marker)
    if start == -1:
        print(f"  WARNING: could not find Promise.all block in {cjs_path}")
        return False

    line_start = content.rfind('\n', 0, start) + 1

    end_marker = '})));'
    end = content.find(end_marker, start)
    if end == -1:
        print(f"  WARNING: could not find closing }}))) in {cjs_path}")
        return False
    end += len(end_marker)

    inner_start = start + len(marker) + 1
    inner = content[inner_start:end - len(end_marker)]

    replacement = f'''  // PATCHED_SEQUENTIAL_HTML: sequential rendering to lower peak memory
  for (const route of routes) {{
{inner.replace('return;', 'continue;')}  }}'''

    content = content[:line_start] + replacement + content[end:]

    with open(cjs_path, 'w') as f:
        f.write(content)

    print(f"  PATCHED: sequential HTML rendering (no server restart)")
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: patch-html-sequential.py <path-to-cjs>")
        sys.exit(1)
    success = patch(sys.argv[1])
    sys.exit(0 if success else 1)
