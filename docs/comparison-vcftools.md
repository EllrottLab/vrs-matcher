# How does vrs-matcher compare to VCFtools?

Bioinformatics response · 2026-09-30

Based on repository commit `9382c731bad1335fcdc47c31e71fd44353bc30b6`.

Implementation update: the [KING-robust plugin](how-to-king-robust.md) now adds
kinship estimation using an optional called-genotype index. The comparison below
describes the original `identity` algorithm and carried-allele index; its
callability limitations do not apply to the new explicit genotype table.
Biological validation and any performance superiority claim remain pending.

## Response to the analyst/engineer

**VCFtools is a general VCF analysis toolkit; vrs-matcher is a specialized sample-similarity index.** There is overlap in comparing genetic observations, but they answer different questions and their scores are not interchangeable.

VCFtools provides filtering, summary statistics, population-genetic calculations, and file comparisons. Its individual-discordance workflow is useful for asking, “How often do these two callsets disagree for the mapped individual?” [VCFtools manual](https://vcftools.github.io/man_latest.html)

This repo asks, “Which indexed samples carry similar sets of VRS alleles?” It ingests already annotated VCFs into SQLite, compares pairs, ranks a query against other samples, and returns shared allele identifiers. That makes it useful for exploring duplicate submissions, sample swaps, and repeated processing of a specimen. It does not yet establish a validated identity threshold or demonstrate a performance advantage over VCFtools. [Loader](../src/vrs_matcher/loader.py), [matcher](../src/vrs_matcher/matcher.py), [identity guide](../docs/how-to-sample-identity-confirmation.md)

## Practical differences

| Question | vrs-matcher | VCFtools |
|---|---|---|
| What is compared? | Carried non-reference `VRS_Allele_IDs`, with genotype-state metadata. | VCF sites and genotypes in the selected comparison workflow. |
| What must I prepare? | VRS annotation upstream; unique sample names across independently indexed observations. | Compatible, consistently sorted VCFs; individual mapping when names differ. |
| Can I search for a sample match? | Built-in pairwise comparison and one-versus-all ranking by Jaccard. | Individual discordance compares supplied individual correspondences; a cohort search needs orchestration. |
| What persists? | A reusable SQLite sample–allele index. | Analysis outputs from the requested file operation. |
| How can I customize the matcher? | Matching plugins and a Python API for restricting candidate VRS IDs. | Configure the selected VCF analysis and its filters. |

Repository behavior: [CLI](../src/vrs_matcher/cli.py), [database](../src/vrs_matcher/db.py), [plugins](../docs/plugins.md). VCFtools comparison requirements and options: [manual](https://vcftools.github.io/man_latest.html).

“VCFtools” also includes Perl utilities: `vcf-compare` reports overlap and non-reference discordance and offers sequence comparison for indels. It would be inaccurate to describe the entire suite as a position-only comparator. The proposed benchmark names the C++ `vcftools --diff-indv-discordance` operation explicitly. [Perl utilities documentation](https://vcftools.github.io/perl_module.html)

## Why the scores differ

The default matcher computes:

- **Jaccard:** shared carried VRS IDs divided by the union of carried IDs.
- **Weighted concordance:** over shared IDs only, matching zygosity scores 1 and differing zygosity scores 0.5. No shared IDs gives 0. This is not quality-weighted or allele-frequency-weighted, and stored phasing does not affect the score.

These definitions come directly from [matcher.py](../src/vrs_matcher/matcher.py). In contrast, the VCFtools individual-discordance implementation counts genotype disagreement among compared calls at common sites, including reference genotypes; diploid allele order does not affect equality. [VCFtools comparison implementation](https://github.com/vcftools/vcftools/blob/master/src/cpp/variant_file_diff.cpp)

For example, consider three fully called biallelic loci with matching REF/ALT in both files:

| Locus | Sample A | Sample B |
|---|---|---|
| 1 | `0/0` | `0/0` |
| 2 | `0/1` | `0/1` |
| 3 | `1/1` | `0/1` |

Both samples carry the same two ALT IDs, so Jaccard is **1.0**. Weighted concordance is **0.75**. Exact genotype concordance is **2/3**, corresponding to VCFtools discordance **1/3**. These are three valid answers to different questions. A shared-allele score can also look strong when very few alleles are shared; report counts alongside scores.

## What VRS does—and what this implementation assumes

The intended benefit is that independently prepared records can match through a consistently generated allele identifier. However, this loader consumes strings from `VRS_Allele_IDs`; it neither generates nor validates their biological equivalence. Annotation, reference compatibility, normalization, and ALT-to-ID ordering are upstream responsibilities. Do not assume automatic cross-assembly matching, equivalence between a complex allele and decomposed events, or general structural-variant support merely because the index accepts IDs. [Loader](../src/vrs_matcher/loader.py), [annotation fixture](../tests/integration/conftest.py)

Several implementation details affect biological interpretation:

- Reference calls and fully missing calls produce no allele rows. An absent allele in the index cannot distinguish reference, missing coverage, filtering, or absent annotation. Missingness and different assay footprints can therefore reduce overlap.
- PASS/unfiltered records are accepted; defaults are GQ ≥20 and DP ≥0 **when those values exist**. Missing quality values are not automatically excluded.
- Partial calls such as `./1` contribute an ALT with `HET` state. Multiallelic calls use a record-level zygosity category, not a complete allele-dosage comparison.
- Two empty sets have Jaccard 1.0: this is a mathematical convention, not evidence of identical specimens.
- Sample names are global database keys. `source_dataset` does not prevent repeated names from merging observations or replacing overlapping allele rows.

Evidence: [loader](../src/vrs_matcher/loader.py), [matcher](../src/vrs_matcher/matcher.py), [database](../src/vrs_matcher/db.py).

## How I would choose

For conventional VCF QC and concordance between known sample correspondences, use the relevant VCFtools operation. For repeated searches against a stored cohort using VRS allele identity, evaluate vrs-matcher. Its current identity algorithm is neither kinship estimation nor haplotype comparison.

The existing 1000 Genomes test checks that mean within-super-population Jaccard exceeds between-super-population Jaccard. That is a useful integration check, but population similarity is not proof of sample identity. We should test labeled replicates and difficult nonmatches, harmonize callable sites, and measure annotation/indexing costs before claiming better accuracy or speed. The [benchmark ADR](adr-benchmark-vcftools.md) specifies that evaluation. [Existing integration test](../tests/integration/test_1kg_cluster.py)
