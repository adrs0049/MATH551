#!/usr/bin/env node
/**
 * Compute MyST execution cache keys for notebooks using the same logic as mystmd.
 * Compare with the populate-execute-cache.py output to find mismatches.
 */
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const cacheDir = '_build/execute';
const nbDir = 'notebooks';

function buildCacheKey(kernelName, cells) {
    const hashableItems = [];
    for (const cell of cells) {
        if (cell.cell_type !== 'code') continue;
        const source = Array.isArray(cell.source) ? cell.source.join('') : cell.source;
        hashableItems.push({
            kind: 'block',
            content: source,
            raisesException: false
        });
    }
    const hash = crypto.createHash('md5')
        .update(kernelName)
        .update(JSON.stringify(hashableItems));
    return hash.digest('hex');
}

const notebooks = fs.readdirSync(nbDir)
    .filter(f => f.endsWith('.ipynb'))
    .sort();

const cacheFiles = new Set(fs.existsSync(cacheDir) ? fs.readdirSync(cacheDir) : []);

console.log(`Notebooks: ${notebooks.length}, Cache files: ${cacheFiles.size}`);
console.log('='.repeat(70));

let matched = 0, missing = 0;
for (const nbFile of notebooks) {
    const nb = JSON.parse(fs.readFileSync(path.join(nbDir, nbFile), 'utf8'));
    const kernelName = nb.metadata?.kernelspec?.name || 'python3';
    const cacheKey = buildCacheKey(kernelName, nb.cells);
    const found = cacheFiles.has(`${cacheKey}.json`);
    const status = found ? 'OK' : 'MISSING';
    console.log(`  ${status.padEnd(7)} ${nbFile} -> ${cacheKey}.json`);
    if (found) matched++; else missing++;
}

console.log('='.repeat(70));
console.log(`Results: ${matched} cached, ${missing} missing`);
