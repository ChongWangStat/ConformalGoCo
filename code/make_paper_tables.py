"""Regenerate the current manuscript tables and verify their displayed numbers.

The numerical calculations remain in _table_builder.py and the original
special-purpose builders. table_specs.json supplies the authors' current
captions/notes and the displayed values to verify, never calculation inputs.
"""
from pathlib import Path
import json, os, re, subprocess, sys
ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / 'paper' / 'tables'


def displayed_rows(text):
    body = text[text.index(r'\midrule'):text.index(r'\bottomrule')]
    rows = []
    for line in body.splitlines():
        if '&' not in line or r'\multicolumn' in line or r'\cmidrule' in line:
            continue
        row = []
        for cell in line.split('&'):
            s = cell.strip()
            for a,b in [(r'\\',''),(r'\textbf',''),('$',''),('{',''),('}',''),(r'\%','')]:
                s=s.replace(a,b)
            s=s.strip()
            if re.fullmatch(r'[+-]?[0-9]+(?:\.[0-9]+)?',s): row.append(float(s))
        if row: rows.append(row)
    return rows


def replace_command(text, name, replacement):
    start=text.index('\\'+name+'{');i=start+len(name)+2;depth=1;j=i
    while depth:
        if text[j] in '{}': depth += 1 if text[j]=='{' else -1
        j+=1
    return text[:start]+replacement+text[j:]


def main():
    TABLES.mkdir(parents=True,exist_ok=True)
    spec=json.loads((ROOT/'paper'/'table_specs.json').read_text())
    static={p.name:p.read_text() for p in TABLES.glob('*.tex') if p.name in ('Table1_constants.tex','Table2_methods_compared.tex')}
    for script in ['make_tables.py','make_supp_tables.py','make_block_diagnostics.py','make_case_study.py','make_table_s9.py','_table_builder.py']:
        subprocess.run([sys.executable,str(ROOT/'code'/script)],cwd=ROOT,check=True)
    final=dict(static); report={}; failures=[]
    for source,item in spec.items():
        p=TABLES/source
        if not p.exists(): raise FileNotFoundError('Required table not generated: '+str(p))
        text=p.read_text(); actual=displayed_rows(text); expected=item['expected_rows']
        ok=actual==expected
        report[item['output']]={'displayed_numeric_rows':len(actual),'matches_manuscript':ok}
        if not ok:
            failures.append(item['output'])
            report[item['output']]['actual']=actual; report[item['output']]['expected']=expected
        text=replace_command(text,'caption',item['caption'])
        text=replace_command(text,'label',r'\label{'+item['label']+'}')
        note=re.search(r'\\par\\vspace\{.*?\\end\{minipage\}',text,re.S)
        if note: text=text[:note.start()]+item['note']+text[note.end():]
        elif item['note']:
            at=text.rfind(r'\end{table')
            text=text[:at]+item['note']+'\n'+text[at:]
        final[item['output']]=text
    out=ROOT/'code'/'out';out.mkdir(parents=True,exist_ok=True)
    (out/'table_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    if failures: raise RuntimeError('Displayed numbers differ from manuscript: '+', '.join(failures))
    # Only current, manuscript-used tables remain in the published directory.
    for p in TABLES.iterdir():
        if p.is_file() and p.name!='key_numbers.json': p.unlink()
    for name,text in final.items(): (TABLES/name).write_text(text)
    print('PASS: all 13 calculated tables match the manuscript; 2 fixed specification tables retained.')

if __name__=='__main__': main()
