#!/usr/bin/env python3
"""
Patch mystmd's buildHtml to render pages sequentially instead of in parallel.

The default Promise.all + p-limit(5) causes the Remix server to accumulate
memory for all pages simultaneously, which OOMs on 16GB GitHub runners.

This patch replaces the parallel fetch with a sequential for-loop.
"""
import sys

def patch(cjs_path):
    with open(cjs_path, 'r') as f:
        content = f.read()

    if 'PATCHED_SEQUENTIAL_HTML' in content:
        print(f"  OK: already patched")
        return True

    old = (
        '  await Promise.all(routes.map(async (route) => limitConnections(async () => {\n'
        '    const resp = await fetchWithRetry(session, route.url);\n'
        '    if (!resp.ok) {\n'
        '      session.log.error(`Error fetching ${route.url}`);\n'
        '      return;\n'
        '    }\n'
        '    if (route.binary && resp.body) {\n'
        '      await new Promise((resolve10) => {\n'
        '        const filename = import_node_path58.default.join(htmlDir, route.path);\n'
        '        if (!import_fs_extra.default.existsSync(filename))\n'
        '          import_fs_extra.default.mkdirSync(import_node_path58.default.dirname(filename), { recursive: true });\n'
        '        const fileWriteStream = import_fs_extra.default.createWriteStream(filename);\n'
        '        resp.body.pipe(fileWriteStream);\n'
        '        fileWriteStream.on("finish", resolve10);\n'
        '      });\n'
        '    } else {\n'
        '      const content3 = await resp.text();\n'
        '      writeFileToFolder(import_node_path58.default.join(htmlDir, route.path), content3);\n'
        '    }\n'
        '  })));'
    )

    new = (
        '  // PATCHED_SEQUENTIAL_HTML: render pages one at a time to reduce memory\n'
        '  for (const route of routes) {\n'
        '    const resp = await fetchWithRetry(session, route.url);\n'
        '    if (!resp.ok) {\n'
        '      session.log.error(`Error fetching ${route.url}`);\n'
        '      continue;\n'
        '    }\n'
        '    if (route.binary && resp.body) {\n'
        '      await new Promise((resolve10) => {\n'
        '        const filename = import_node_path58.default.join(htmlDir, route.path);\n'
        '        if (!import_fs_extra.default.existsSync(filename))\n'
        '          import_fs_extra.default.mkdirSync(import_node_path58.default.dirname(filename), { recursive: true });\n'
        '        const fileWriteStream = import_fs_extra.default.createWriteStream(filename);\n'
        '        resp.body.pipe(fileWriteStream);\n'
        '        fileWriteStream.on("finish", resolve10);\n'
        '      });\n'
        '    } else {\n'
        '      const content3 = await resp.text();\n'
        '      writeFileToFolder(import_node_path58.default.join(htmlDir, route.path), content3);\n'
        '    }\n'
        '  }'
    )

    if old not in content:
        print(f"  WARNING: could not find target block in {cjs_path}")
        return False

    content = content.replace(old, new, 1)

    with open(cjs_path, 'w') as f:
        f.write(content)

    print(f"  PATCHED: sequential HTML rendering enabled")
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: patch-html-sequential.py <path-to-myst.cjs>")
        sys.exit(1)
    success = patch(sys.argv[1])
    sys.exit(0 if success else 1)
