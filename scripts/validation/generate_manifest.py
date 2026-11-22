#!/usr/bin/env python3
"""Generate Cryptographic Artifact Manifest & Data Provenance Registry.

Computes SHA-256 hashes, file sizes, timestamps, and physical metadata for:
- All 18 publication-quality figures in article/fig/
- All scientific result JSON archives in results/
- The compiled 26-page academic paper article/main.pdf

Saves the verifiable registry to results/manifest.json.
"""
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256_file(filepath):
    """Computes the SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def main():
    print("=" * 80)
    print("Generating Cryptographic Artifact Manifest & Provenance Registry")
    print("=" * 80)

    manifest = {
        'registry_version': '1.0.0',
        'generated_at': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()),
        'host': os.uname().nodename,
        'figures': {},
        'results_data': {},
        'publication_pdf': {},
    }

    # 1. Figures
    fig_dir = ROOT / 'article' / 'fig'
    png_files = sorted(fig_dir.glob('*.png'))
    for p in png_files:
        rel = p.relative_to(ROOT).as_posix()
        manifest['figures'][p.name] = {
            'relative_path': rel,
            'size_bytes': p.stat().st_size,
            'sha256': sha256_file(p),
            'modified_at': time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(p.stat().st_mtime)),
        }
    print(f"  Registered {len(manifest['figures'])} publication figures in article/fig/")

    # 2. Results Data
    res_dir = ROOT / 'results'
    json_files = sorted(res_dir.glob('*.json'))
    # Exclude ephemeral runtime execution logs that record dynamic machine wallclock timestamps
    ephemeral = {'manifest.json', 'full_regression_report.json', 'sub_permill_audit_summary.json'}
    for p in json_files:
        if p.name in ephemeral:
            continue
        rel = p.relative_to(ROOT).as_posix()
        manifest['results_data'][p.name] = {
            'relative_path': rel,
            'size_bytes': p.stat().st_size,
            'sha256': sha256_file(p),
            'modified_at': time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(p.stat().st_mtime)),
        }
    print(f"  Registered {len(manifest['results_data'])} scientific data JSON archives in results/")

    # 3. Publication PDF
    pdf_path = ROOT / 'article' / 'main.pdf'
    if pdf_path.exists():
        manifest['publication_pdf'] = {
            'relative_path': pdf_path.relative_to(ROOT).as_posix(),
            'size_bytes': pdf_path.stat().st_size,
            'sha256': sha256_file(pdf_path),
            'modified_at': time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(pdf_path.stat().st_mtime)),
        }
        print(f"  Registered publication paper: article/main.pdf ({manifest['publication_pdf']['size_bytes']/1048576:.2f} MB, SHA256={manifest['publication_pdf']['sha256'][:16]}...)")

    out_file = ROOT / 'results' / 'manifest.json'
    with open(out_file, 'w') as f:
        json.dump(manifest, f, indent=2)

    print("=" * 80)
    print(f"PROVENANCE REGISTRY CREATED: {out_file}")
    print("=" * 80)
    return 0


if __name__ == '__main__':
    sys.exit(main())
