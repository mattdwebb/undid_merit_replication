# LaTeX for the revised MERIT results

The revision uses the computed estimates in
[`../results/panel_c_results.csv`](../results/panel_c_results.csv), dated
22 September 2026.

- [`merit_results.tex`](merit_results.tex) replaces the passage beginning
  "Lastly, we estimate average treatment effects" and the complete MERIT table,
  ending immediately before the Conclusion section. Panels A and B are retained
  without numerical changes.
- [`merit_discussion.tex`](merit_discussion.tex) contains only the replacement
  discussion.
- [`panel_c.tex`](panel_c.tex) contains only the Panel C heading and rows, for
  insertion into the existing seven-column `tabular` environment. Retain the
  closing horizontal rules, `end{tabular}`, label, and `end{table}`.

The discussion is part of the correction: UN-DID uses no covariates, whereas
CSDID and DID-INT include `male`, `black`, and `asian`. The UN-DID and DID-INT
SEs delete final-stage subgroups; the CSDID SEs delete states. Replacing only
the numbers while retaining the previous claims that all specifications
include covariates and all SEs cluster by state would be inaccurate.

The fragments rely on the paper's existing citation keys and labels. They
are intended for the current paper, not as standalone LaTeX documents.
The appendix's other CCC specifications are outside this numerical revision.
