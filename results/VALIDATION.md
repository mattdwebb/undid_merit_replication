# Validation of the revised table

Verified locally on 22 September 2026 using Julia 1.11.7, DiDInt.jl 0.7.5,
Undid.jl 0.5.2, Stata 18, `csdid` 1.81, and `csdidjack` 0.5.2.

`do replicate_merit.do reuse_csdid` completed with exit code 0. It regenerated
the Julia estimates from `merit.dta` and reused the completed CSDID calculations
preserved in `diagnostics/verified_csdid_20260922.csv` and its companion log.
This was not a second full CSDID jackknife recomputation.

`python diagnostics/audit_results.py` passed:

- 1,394 UN-DID difference cells and their sample counts agreed with raw data.
- All 16 UN-DID specification/aggregation results agreed with independent
  reconstruction.
- The historical DID-INT subgroup jackknife agreed with an independent
  deletion calculation.
- CSDID values agreed with the preserved completed calculation log.
- All 12 entries agreed with the revised four-decimal paper baseline and the
  numerical snapshot in `panel_c_results.csv`.
- Eight of 12 entries agreed with the earlier draft; its four UN-DID entries
  are superseded, not retrospectively validated.

The replacement LaTeX was compiled within the full paper using pdfLaTeX and
BibTeX. The result has 81 US Letter pages, no unresolved references or
citations, and no new overfull-box warnings. The affected pages were rendered
for visual inspection. Existing unrelated layout warnings remain.

See the main README for covariate and standard-error definitions. This
validation covers Panel C; Panels A and B and the other appendix CCC
specifications were not re-estimated in this revision.
