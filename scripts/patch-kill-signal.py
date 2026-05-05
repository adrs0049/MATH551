#!/usr/bin/env python3
"""
Patch mystmd's killProcessTree to use SIGKILL instead of SIGTERM.

WHY: appServer.stop() calls killProcessTree on the Remix dev server. With
SIGTERM, Remix runs its exit handler which deletes public/build/. The next
line in mystmd is copySync(theme/public, _build/html) — it ENOENTs on the
just-deleted files. SIGKILL prevents the exit handler from running, so
public/build/_assets/* stay intact and the copy succeeds with hashes that
match what mystmd already wrote into the rendered HTML.

WHEN TO REMOVE: when myst-theme is updated to a Remix version that doesn't
delete its build output on shutdown, or when mystmd is fixed to copy theme
assets BEFORE stopping the dev server.
"""
import re
import sys


# Different bundles rename the parameter (process8 vs process9) and the
# resolve callback (resolve10 vs resolveN). Match with regex so we patch
# both jupyter-book.cjs and myst.cjs without hardcoding the names.
PATTERN = re.compile(
    r'function killProcessTree\((\w+)\) \{\n'
    r'  return new Promise\(\((\w+)\) => \{\n'
    r'    if \(!\1\.pid\) \{\n'
    r'      \1\.kill\("SIGTERM"\);\n'
    r'      \2\(\);\n'
    r'      return;\n'
    r'    \}\n'
    r'    \(0, import_tree_kill\.default\)\(\1\.pid, "SIGTERM", \(err\) => \{\n'
    r'      if \(err && \1\.pid\) \{\n'
    r'        \(0, import_tree_kill\.default\)\(\1\.pid, "SIGKILL", \(\) => \2\(\)\);\n'
    r'      \} else \{\n'
    r'        \2\(\);\n'
    r'      \}\n'
    r'    \}\);\n'
    r'  \}\);\n'
    r'\}'
)


def patch(cjs_path):
    with open(cjs_path, 'r') as f:
        content = f.read()

    if 'PATCHED_KILL_SIGNAL' in content:
        print(f"  OK: already patched")
        return True

    m = PATTERN.search(content)
    if not m:
        print(f"  WARNING: killProcessTree pattern not found in {cjs_path}")
        return False

    proc, resolve = m.group(1), m.group(2)
    new_body = (
        f'function killProcessTree({proc}) {{\n'
        f'  // PATCHED_KILL_SIGNAL: SIGKILL prevents Remix dev\'s exit handler from\n'
        f'  // deleting public/build/ before mystmd\'s post-fetch copySync reads it.\n'
        f'  return new Promise(({resolve}) => {{\n'
        f'    if (!{proc}.pid) {{\n'
        f'      {proc}.kill("SIGKILL");\n'
        f'      {resolve}();\n'
        f'      return;\n'
        f'    }}\n'
        f'    (0, import_tree_kill.default)({proc}.pid, "SIGKILL", () => {resolve}());\n'
        f'  }});\n'
        f'}}'
    )

    content = content[:m.start()] + new_body + content[m.end():]
    with open(cjs_path, 'w') as f:
        f.write(content)

    print(f"  PATCHED: killProcessTree uses SIGKILL")
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: patch-kill-signal.py <path-to-cjs>")
        sys.exit(1)
    success = patch(sys.argv[1])
    sys.exit(0 if success else 1)
