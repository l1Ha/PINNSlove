#!/usr/bin/env python3
"""Verify Cryptographic Artifact Manifest & Provenance Registry.

Checks SHA-256 integrity for all registered files in results/manifest.json.
Returns 0 if all checksums match, 1 otherwise.
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256_file(filepath):
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def main():
    manifest_file = ROOT / 'results' / 'manifest.json'
    if not manifest_file.exists():
        print(f"Error: {manifest_file} not found.")
        return 1

    with open(manifest_file) as f:
        manifest = json.load(f)

    print("=" * 80)
    print(f"Verifying Integrity of Artifacts against {manifest_file.name}")
    print("=" * 80)

    total_files = 0
    mismatches = []

    # Check figures
    for name, info in manifest.get('figures', {}).items():
        total_files += 1
        path = ROOT / info['relative_path']
        if not path.exists():
            mismatches.append(f"Missing file: {info['relative_path']}")
            continue
        h = sha256_file(path)
        if h != info['sha256']:
            mismatches.append(f"SHA-256 mismatch in {info['relative_path']}")

    # Check results
    for name, info in manifest.get('results_data', {}).items():
        total_files += 1
        path = ROOT / info['relative_path']
        if not path.exists():
            mismatches.append(f"Missing file: {info['relative_path']}")
            continue
        h = sha256_file(path)
        if h != info['sha256']:
            mismatches.append(f"SHA-256 mismatch in {info['relative_path']}")

    # Check PDF
    pdf_info = manifest.get('publication_pdf')
    if pdf_info:
        total_files += 1
        path = ROOT / pdf_info['relative_path']
        if not path.exists():
            mismatches.append(f"Missing file: {pdf_info['relative_path']}")
        else:
            h = sha256_file(path)
            if h != pdf_info['sha256']:
                mismatches.append(f"SHA-256 mismatch in {pdf_info['relative_path']}")

    if mismatches:
        print(f"FAILED: {len(mismatches)} files had integrity errors:")
        for m in mismatches:
            print(f"  [ERROR] {m}")
        return 1

    print(f"SUCCESS: All {total_files} artifacts verified with matching SHA-256 checksums! (100% OK)")
    print("=" * 80)
    return 0


if __name__ == '__main__':
    sys.exit(main())
