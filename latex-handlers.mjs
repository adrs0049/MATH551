/**
 * MyST plugin: document transform that rewrites proof nodes with unsupported
 * kinds (algorithm, property) so they render in LaTeX.
 *
 * Also converts `details` nodes to render their children inline for LaTeX,
 * since there is no LaTeX equivalent of collapsible dropdowns.
 *
 * Register in myst.yml:
 *   project:
 *     plugins:
 *       - latex-handlers.mjs
 */

function visit(node, type, callback) {
  if (!node) return;
  if (node.type === type) callback(node);
  if (node.children) {
    for (const child of node.children) {
      visit(child, type, callback);
    }
  }
}

function visitAll(node, callback) {
  if (!node) return;
  callback(node);
  if (node.children) {
    for (const child of node.children) {
      visitAll(child, callback);
    }
  }
}

const plugin = {
  name: 'latex-proof-fix-plugin',
  transforms: [
    {
      name: 'fix-unsupported-proof-kinds',
      doc: 'Map unsupported prf: kinds (algorithm, property) to supported environments by rewriting the kind field.',
      stage: 'document',
      plugin: (opts, utils) => (tree) => {
        // Convert details/summary, grid, and bibliography nodes to plain content.
        // details nodes come from :class: dropdown or :::{dropdown} directives.
        // In LaTeX there's no collapsible content, so we just render inline.
        // grid nodes come from card/grid layouts with no LaTeX equivalent.
        // bibliography nodes are handled by the template's \bibliography command.
        const flattenTypes = new Set(['details', 'grid', 'card', 'cardTitle', 'header', 'footer', 'bibliography']);
        const skipTypes = new Set(['summary']);

        // Recursively flatten unsupported node types.
        // Repeat until stable since flattening can expose more nodes to flatten.
        function flattenUnsupported(node) {
          if (!node.children) return;
          // First recurse into children
          for (const child of node.children) {
            flattenUnsupported(child);
          }
          // Then flatten this node's children
          let changed = true;
          while (changed) {
            changed = false;
            const newChildren = [];
            for (const child of node.children) {
              if (flattenTypes.has(child.type)) {
                changed = true;
                if (child.children) {
                  for (const innerChild of child.children) {
                    if (!skipTypes.has(innerChild.type)) {
                      newChildren.push(innerChild);
                    }
                  }
                }
              } else {
                newChildren.push(child);
              }
            }
            node.children = newChildren;
          }
        }
        flattenUnsupported(tree);
      },
    },
  ],
};

export default plugin;
