# Synthetic KING oracle

`oracle.vcf` and `panel.tsv` encode the eight-marker, two-sample numerical
example from `docs/user-story-king-robust-plugin.md` (acceptance B).
All sequence and VRS identifiers are **synthetic, syntactically valid test
identifiers**, not computed annotations of biological alleles. Do not use this
panel for biological inference. No download or reference database is required.

Expected A/B counts: M=8, H_A=4, H_B=6, HH=3, O=1.
Between-family kinship=0; within-family estimate=0.1; IBS0 fraction=0.125.
Self comparisons have kinship=0.5. Tests generate missing-call, QC, and invalid
input variants of these fixtures in temporary directories.
