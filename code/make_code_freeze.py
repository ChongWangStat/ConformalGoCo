"""Write checksums of the current tracked reproducibility files (not Git history)."""
from pathlib import Path
import hashlib, subprocess
ROOT=Path(__file__).resolve().parents[1]
paths=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
lines=['# Current reproducibility file checksums','',
       'SHA-256 hashes cover the current checkout. This manifest excludes itself and temporary build output.','',
       'The complete FunMap score file and large STRING view are distributed separately through the Zenodo archive.','',
       '| File | Bytes | SHA-256 |','|---|---:|---|']
for name in sorted(paths):
    p=ROOT/name
    if not name or not p.is_file() or name=='CODE_FREEZE.md' or name.startswith('code/out/') or name.endswith('_sync_payload.b64'): continue
    digest=hashlib.sha256(p.read_bytes()).hexdigest()
    lines.append(f'| `{name}` | {p.stat().st_size} | `{digest}` |')
(ROOT/'CODE_FREEZE.md').write_text('\n'.join(lines)+'\n')
print('Wrote checksums for',len(lines)-8,'files')
