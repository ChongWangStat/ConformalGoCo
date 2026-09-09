"""Body word count for the GoCo manuscript, per section, excluding floats, displays, captions and endmatter."""
import re, os
os.chdir(r"C:\Users\chwang\Box\Research\26 09 ConformalGOCo\GoCo_BIB_revision_v2")

BS = chr(92)

def detex(t):
    t = re.sub(r'(?m)^\s*%.*$', '', t)                       # comment lines
    for env in ['table', r'table\*', 'figure', r'figure\*', 'equation', r'equation\*',
                'align', r'align\*', 'tabular', 'threeparttable']:
        t = re.sub(r'\\begin\{' + env + r'\}.*?\\end\{' + env + r'\}', ' ', t, flags=re.S)
    t = re.sub(r'\\(input|includegraphics|label|ref|eqref|cite|citep|caption)\s*\{[^{}]*\}', ' X ', t)
    t = re.sub(r'\$[^$]*\$', ' X ', t)                        # inline math -> one token
    t = re.sub(r'\\[a-zA-Z]+\*?', ' ', t)                     # remaining macros
    t = re.sub(r'[{}' + BS + r'&~^_]', ' ', t)
    return [w for w in t.split() if re.search(r'[A-Za-z]', w)]

main = open('GoCo_BIB_v2.tex', encoding='utf-8').read()
meth = open('sections/methods_theory_v2.tex', encoding='utf-8').read()

start = main.find(BS + 'section{Introduction}')
end = main.find(BS + 'section*{Key Points}')
if end < 0: end = main.find(BS + 'section*{Data availability}')
body = main[start:end].replace(BS + 'input{sections/methods_theory_v2}', meth)

print('BODY (Introduction -> Conclusions):', len(detex(body)), 'words')
marks = [(m.start(), m.group(1)) for m in re.finditer(r'\\section\{([^}]*)\}', body)] + [(len(body), 'END')]
for (a, name), (b, _) in zip(marks, marks[1:]):
    seg = body[a:b]
    subs = [(m.start(), m.group(1)) for m in re.finditer(r'\\subsection\{([^}]*)\}', seg)]
    print(f'  {name:32s} {len(detex(seg)):5d}')
    for k, (sa, sname) in enumerate(subs):
        sb = subs[k + 1][0] if k + 1 < len(subs) else len(seg)
        print(f'      {sname[:56]:58s} {len(detex(seg[sa:sb])):5d}')

ab = main[main.find(BS + 'abstract{'): main.find(BS + 'keywords{')]
kp = main[main.find(BS + 'section*{Key Points}'): start]
print('\nabstract:', len(detex(ab)), '| key points:', len(detex(kp)))
print('methods file alone:', len(detex(meth)))
