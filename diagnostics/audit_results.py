"""Independent checks of computed estimates and archived UN-DID summaries.

Run after replicate_merit.do, using Python 3 (standard library only).
The archived summaries are evidence to audit, never inputs to estimation.
"""
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def mean(values):
    values = list(values)
    return sum(values) / len(values)


def jackknife(y, w):
    total = sum(w)
    w = [v / total for v in w]
    att = sum(a * b for a, b in zip(y, w))
    deleted = [(att - a * b) / (1 - b) for a, b in zip(y, w)]
    se = math.sqrt((len(y) - 1) / len(y) * sum((v - att) ** 2 for v in deleted))
    return att, se


data_hash = hashlib.sha256((ROOT / "merit.dta").read_bytes()).hexdigest()
assert data_hash == "1509b32bf680bf34783c5f27d58027e67931c85eead8f58c235c004b8887abdc"
raw = read(OUT / "merit_input.csv")
assert len(raw) == 42161
assert len({r["state"] for r in raw}) == 51
assert len({r["state"] for r in raw if int(r["gvar"]) > 0}) == 10
assert {int(r["year"]) for r in raw} == set(range(1989, 2001))
cells = defaultdict(list)
for r in raw:
    cells[int(r["state"]), int(r["year"])].append(float(r["coll"]))
diffs = []
errors = []
for path in (OUT / "undid_stages").glob("filled_diff_df_*.csv"):
    for r in read(path):
        if r["RI"] != "0":
            continue
        s = int(r["silo_name"])
        post, pre = map(int, r["diff_times"].split(";"))
        delta = mean(cells[s, post]) - mean(cells[s, pre])
        errors.append(abs(delta - float(r["diff_estimate"])))
        assert int(r["n"]) == len(cells[s, post]) + len(cells[s, pre])
        assert int(r["n_t"]) == len(cells[s, post])
        diffs.append(r)
assert len(diffs) == 1394
assert max(errors) < 1e-12

specifications = []
for r in read(OUT / "undid_specifications.csv"):
    group_key = "gt" if r["aggregation"] == "simple" else "gvar"
    col = "diff_estimate_covariates" if r["covariates"] == "true" else "diff_estimate"
    grouped = defaultdict(list)
    for d in diffs:
        grouped[d[group_key]].append(d)
    y, w = [], []
    for group in grouped.values():
        means = []
        for treatment in (1, 0):
            subset = [d for d in group if int(d["treat"]) == treatment]
            weights = [float(d["n"]) if r["weighting"] in ("both", "diff") else 1 for d in subset]
            means.append(sum(float(d[col]) * wt for d, wt in zip(subset, weights)) / sum(weights))
        y.append(means[0] - means[1])
        w.append(sum(float(d["n_t"]) for d in group if int(d["treat"]) == 1)
                 if r["weighting"] in ("both", "att") else 1)
    att, se = jackknife(y, w)
    assert abs(att - float(r["att"])) < 1e-10
    assert abs(se - float(r["subgroup_jackknife_se"])) < 1e-10
    specifications.append({**r, "independent_att": att, "independent_se": se})

historical = []
components = []
for agg, col, key in (("gt", "ATT_gt", "gt"), ("g", "ATT_g", "gvar")):
    path = ROOT / "diagnostics/historical" / f"UNDID_results_{agg}_new.csv"
    rows = [r for r in read(path) if r[agg]]  # discard blank/annotation rows
    y = [float(r[col]) for r in rows]
    att, se = jackknife(y, [1] * len(y))
    stored_att = float(rows[0]["agg_ATT"])
    assert abs(att - stored_att) < 1e-9
    historical.append({"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n": len(y), "stored_att": stored_att, "component_mean": att,
        "stored_se": float(rows[0]["jackknife_SE"]), "component_jackknife_se": se,
        "stored_se_times_sqrt50": float(rows[0]["jackknife_SE"]) * math.sqrt(50)})
    for r in rows:
        subset = [d for d in diffs if d[key] == r[agg]]
        actual = mean(float(d["diff_estimate"]) for d in subset if d["treat"] == "1") - mean(
            float(d["diff_estimate"]) for d in subset if d["treat"] == "0")
        components.append({"aggregation": agg, "subgroup": r[agg], "archived": float(r[col]),
                           "fresh_unweighted": actual, "difference": float(r[col]) - actual})

didint = []
for agg, col in (("simple", "att_gt"), ("group", "att_cohort")):
    rows = read(OUT / f"didint_{agg}_subgroups.csv")
    att, se = jackknife([float(r[col]) for r in rows], [float(r["weights"]) for r in rows])
    assert abs(att - float(rows[0]["agg_att"])) < 1e-12
    assert abs(se - float(rows[0]["se_agg_att"]) * math.sqrt((len(rows)-1)/len(rows))) < 1e-12
    didint.append({"aggregation": agg, "att": att, "historical_subgroup_jackknife_se": se})

targets = [[.0466, .0113, .0464, .0133, .0464, .0102],
           [.0458, .0133, .0339, .0211, .0458, .0084]]
previous_targets = [[.0485, .0110, .0464, .0133, .0464, .0102],
                    [.0459, .0188, .0339, .0211, .0458, .0084]]
columns = ["undid_att", "undid_se", "csdid_att", "csdid_se", "didint_att", "didint_se"]
checks = []
previous_checks = []
csdid_source = ROOT / "diagnostics/verified_csdid_20260922.txt"
csdid_log = csdid_source.read_text()
for row in read(ROOT / "diagnostics/verified_csdid_20260922.csv"):
    assert row["data_sha256"] == data_hash
    for name in ("att", "se"):
        match = re.search(r"^RESULT_CSDID_" + row["aggregation"].upper() + "_" +
                          name.upper() + r"=\s*([0-9.]+)", csdid_log, re.MULTILINE)
        assert match and float(match.group(1)) == float(row["csdid_" + name])
computed_rows = read(OUT / "panel_c_results.csv")
assert [r["aggregation"] for r in computed_rows] == ["simple", "group"]
for row, target in zip(computed_rows, targets):
    for name, expected in zip(columns, target):
        value = float(row[name])
        checks.append({"aggregation": row["aggregation"], "column": name, "computed": value,
                       "target": expected, "matches_four_decimals": round(value, 4) == expected})

for row, target in zip(computed_rows, previous_targets):
    for name, expected in zip(columns, target):
        previous_checks.append({"aggregation": row["aggregation"], "column": name,
            "computed": float(row[name]), "target": expected,
            "matches_four_decimals": round(float(row[name]), 4) == expected})

# Ensure the paper's committed numerical snapshot agrees with computed output.
saved_rows = read(ROOT / "results/panel_c_results.csv")
assert [r["aggregation"] for r in saved_rows] == ["simple", "group"]
for computed, saved in zip(computed_rows, saved_rows):
    for name in columns:
        assert abs(float(computed[name]) - float(saved[name])) < 1e-9

report = {"data_sha256": data_hash, "raw_difference_cells_checked": len(diffs), "max_raw_difference_error": max(errors),
          "csdid_source_log_sha256": hashlib.sha256(csdid_source.read_bytes()).hexdigest(),
          "run_status": (OUT / "replication_status.txt").read_text().strip(),
          "independent_undid_specifications": specifications, "historical_summaries": historical,
          "didint_independent_checks": didint, "table_checks": checks,
          "baseline": "revised paper, 22 September 2026",
          "earlier_draft_checks": previous_checks,
          "earlier_draft_matches": sum(c["matches_four_decimals"] for c in previous_checks),
          "matches": sum(c["matches_four_decimals"] for c in checks), "total": len(checks)}
(OUT / "replication_audit.json").write_text(json.dumps(report, indent=2) + "\n")
with (OUT / "historical_undid_component_comparison.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(components[0]))
    writer.writeheader()
    writer.writerows(components)
assert report['matches'] == report['total'] == 12
print(f"Independent checks passed; revised baseline: {report['matches']}/{report['total']} match; "
      f"earlier draft: {report['earlier_draft_matches']}/12 match.")
