"""Verify local reconstructions, not pre-committed manuscript tables or figures."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/"code/out/table_generation.json"
if not p.exists():raise FileNotFoundError("Run python code/make_paper_tables.py first")
report=json.loads(p.read_text())
assert len(report)>=13 and all(x["bytes"]>120 for x in report.values())
for panel in ("funmap","string"):
    d=pd.read_csv(ROOT/"results/sf"/f"{panel}_v4_100.csv")
    counts=d.groupby(["arm","alpha","delta"]).repeat.nunique()
    assert (counts==100).all(),(panel,counts)
    assert not d.duplicated(["repeat","arm","alpha","delta"]).any()
    assert len(pd.read_csv(ROOT/"results/sf"/f"{panel}_generic_dense.csv"))>0
expected={"GoCo_Workflow.pdf","Figure2_ties.pdf","Figure3_primary_delta050.pdf",
          "Figure4_alpha_sweep_delta050.pdf","FigureS1_primary_delta010.pdf"}
generated=ROOT/"paper/figures"
assert expected.issubset({p.name for p in generated.glob("*.pdf")})
assert all((generated/x).stat().st_size>1000 for x in expected)
assert not (ROOT/"paper/table_specs.json").exists()
assert not (ROOT/"paper/submission").exists()
result={"generated_analytical_tables":len(report),"neighborhood_splits":"100 per arm/target/delta",
        "generated_original_figures":sorted(expected),
        "funmap_string_from_score_rerun":"requires missing frozen per-call prediction scores"}
(ROOT/"code/out/release_verification.json").write_text(json.dumps(result,indent=2)+"\n")
print("PASS: validated regenerated numerical tables and original figures")
