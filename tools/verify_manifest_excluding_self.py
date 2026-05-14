#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, sys
root=Path(__file__).resolve().parents[1]
manifest=root/'FILE_MANIFEST.csv'
missing=[]; bad=[]; checked=0
with manifest.open(newline='',encoding='utf-8') as f:
    for row in csv.DictReader(f):
        rel=row['path']
        if rel in {'FILE_MANIFEST.csv','FILE_MANIFEST.json'}: continue
        p=root/rel
        if not p.exists(): missing.append(rel); continue
        h=hashlib.sha256(p.read_bytes()).hexdigest(); checked+=1
        if h!=row['sha256']: bad.append((rel,row['sha256'],h))
if missing or bad:
    print('Manifest verification FAILED')
    for x in missing[:50]: print('MISSING',x)
    for x in bad[:50]: print('MISMATCH',x[0],x[1],x[2])
    sys.exit(1)
print(f'Manifest verification PASS: verified {checked} entries; manifest files excluded by policy.')
