"""Write portable checksums of current tracked files after Git clean filters."""
from pathlib import Path
import hashlib, subprocess
ROOT=Path(__file__).resolve().parents[1]
paths=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
lines=['# Current reproducibility file checksums','',
       'SHA-256 hashes cover Git-normalized file contents (text line endings may differ in a Windows checkout). This manifest excludes itself and temporary build output.','',
       'Large inputs intentionally absent from Git are not covered here; see README.md for their availability and archive scope.','',
       '| File | Bytes | SHA-256 |','|---|---:|---|']
for name in sorted(paths):
    p=ROOT/name
    if not name or not p.is_file() or name=='CODE_FREEZE.md' or name.startswith('code/out/') or name.endswith('_sync_payload.b64'): continue
    blob_id=subprocess.check_output(['git','hash-object','-w','--path='+name,'--',name],cwd=ROOT,stderr=subprocess.PIPE).strip()
    content=subprocess.check_output(['git','cat-file','blob',blob_id.decode()],cwd=ROOT)
    digest=hashlib.sha256(content).hexdigest()
    lines.append(f'| `{name}` | {len(content)} | `{digest}` |')
(ROOT/'CODE_FREEZE.md').write_text('\n'.join(lines)+'\n')
print('Wrote checksums for',len(lines)-8,'files')
