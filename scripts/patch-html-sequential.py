#!/usr/bin/env python3
"""
Patch mystmd/jupyter-book buildHtml to render pages sequentially and
restart the Remix server every BATCH_SIZE pages to free memory.

The default Promise.all + p-limit(5) causes the Remix server to accumulate
memory for all route modules, which OOMs on 16GB GitHub runners.
"""
import sys

BATCH_SIZE = 15  # Restart server every N pages


def patch(cjs_path):
    with open(cjs_path, 'r') as f:
        content = f.read()

    if 'PATCHED_SEQUENTIAL_HTML' in content:
        print(f"  OK: already patched")
        return True

    # Find the Promise.all block
    marker = 'await Promise.all(routes.map(async (route) => limitConnections(async () => {'
    start = content.find(marker)
    if start == -1:
        print(f"  WARNING: could not find Promise.all block in {cjs_path}")
        return False

    # Find line start (include leading whitespace)
    line_start = content.rfind('\n', 0, start) + 1

    # Find the closing })));
    end_marker = '})));'
    end = content.find(end_marker, start)
    if end == -1:
        print(f"  WARNING: could not find closing }}))) in {cjs_path}")
        return False
    end += len(end_marker)

    # Extract the inner block (route processing logic)
    inner_start = start + len(marker) + 1
    inner = content[inner_start:end - len(end_marker)]

    # Build replacement: sequential for-loop with server restart every BATCH_SIZE pages
    replacement = f'''  // PATCHED_SEQUENTIAL_HTML: sequential rendering with server restarts
  // to reduce memory on CI runners (batch size: {BATCH_SIZE})
  const BATCH_SIZE = {BATCH_SIZE};
  for (let i = 0; i < routes.length; i += BATCH_SIZE) {{
    const batch = routes.slice(i, i + BATCH_SIZE);
    for (const route of batch) {{
{inner.replace('return;', 'continue;')}  }}
    // Restart server to free accumulated memory
    if (i + BATCH_SIZE < routes.length) {{
      session.log.info(`\\u267B\\uFE0F  Restarting server after ${{i + BATCH_SIZE}}/${{routes.length}} pages to free memory`);
      appServer.stop();
      const newServer = await startServer(session, {{ ...opts, buildStatic: true, baseurl }});
      if (!newServer) {{
        session.log.error('Failed to restart server');
        break;
      }}
      Object.assign(appServer, newServer);
    }}
  }}'''

    content = content[:line_start] + replacement + content[end:]

    with open(cjs_path, 'w') as f:
        f.write(content)

    print(f"  PATCHED: sequential HTML rendering with server restart every {BATCH_SIZE} pages")
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: patch-html-sequential.py <path-to-cjs>")
        sys.exit(1)
    success = patch(sys.argv[1])
    sys.exit(0 if success else 1)
