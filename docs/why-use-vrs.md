# Why VRS identifiers are used alongside VCF coordinates

VRS-Matcher uses **VRS identifiers to identify alleles across indexed sample
observations**, while retaining VCF information to interpret calls and trace
results to their source. This serves a bioinformatics workflow in which samples
may arrive in separate files, from different processing runs, or under different
observation names. The comparison key should describe the allele, rather than
its position within a particular file.

VRS is an addition to the VCF workflow. It does not replace genotype calls,
reference provenance, quality filtering, or a consistently defined comparison
panel. In this repository, annotation happens upstream; indexing consumes the
annotated VCF, and matching reads the resulting SQLite database.

## Why a chromosome and position are insufficient

A location is not an allele. At one position, `A>C` and `A>G` are different
alleles, even though their chromosome and position match. Chromosome names also
need reference context: `chr1:100` does not identify which reference sequence
was used. Conversely, two files can use different names for the same sequence.

A reference-qualified, consistently normalized `(CHROM, POS, REF, ALT)` tuple
can be a valid comparison key. VRS offers a standardized computed identifier
for exchanging that identity between pipelines and data stores. A VRS Allele
represents a state at a location on a reference sequence; its identifier is
computed from the identifying content, rather than assigned by a central
registry. Consistent normalization and identifier generation are necessary for
interoperability. See the [GA4GH computed identifier specification](https://vrs.ga4gh.org/en/stable/conventions/computed_identifiers.html).

This is particularly relevant to indels in repetitive sequence, where equivalent
sequence changes can have different textual representations. VRS normalization
addresses representation ambiguity, but comparing arbitrary ID strings cannot
repair inconsistent upstream annotation. The underlying representation and
normalization are described in the [VRS specification paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8929418/).

For this project, the practical benefit is a reusable allele key for sample
comparison and shared-allele reporting. It is not evidence that VRS-based
matching is intrinsically more accurate or faster than a correctly harmonized
VCF-based analysis.

## What each piece of information contributes

| Information | Meaning and use here | What it does not establish |
|---|---|---|
| `CHROM`, `POS` | Input locus; retained in the identity index and used to locate KING panel records | Allele identity or reference sequence identity by themselves |
| `REF`, `ALT` | Record allele definitions; checked against the KING panel | That two differently represented complex events are equivalent |
| `ga4gh:SQ.…` sequence identifier | Declared reference sequence identity in the KING panel and VCF contig header | Independent verification that the input was actually called against that sequence |
| `ga4gh:VA.…` allele identifier | Shared allele key for identity matching and KING marker lookup | A sample's genotype, call quality, donor identity, or relationship |
| `GT`, `GQ`, `DP`, `FILTER` | Observed genotype and available evidence used during ingestion | Allele identity independent of the record's allele definitions |
| Reference/annotation metadata and file hashes | Reproducibility and compatibility checks for genotype-enabled ingestion | Biological correctness of the source calls |

A VRS ID is also distinct from the VCF `ID` column: the loader reads the
`INFO/VRS_Allele_IDs` annotation, not a record's rsID or other local identifier.
The identifier is not a human-readable coordinate and cannot be decoded back
into the source VCF record. Preserve the input files and annotation provenance.

## Where and how the repository uses VRS

### 1. Upstream: create allele annotations

Prepare VCFs against a known reference and use a consistent VRS annotation
implementation and normalization policy. The expected annotation contains one
identifier per ALT, in ALT order:

```text
##INFO=<ID=VRS_Allele_IDs,Number=A,Type=String,Description="VRS ALT allele IDs">
```

For example, a record with `REF=A`, `ALT=C,G`, and `GT=0/2` carries the second
ALT's VRS identifier. A `1/2` genotype carries both ALT identifiers. These are
illustrative allele mappings, not real variant annotations. Do not insert a
reference-allele identifier at the start of this list: that would change its
alignment with the loader's ALT indexing.

The [integration fixture](../tests/integration/conftest.py) illustrates upstream
annotation with `ga4gh.vrs` and a SeqRepo data proxy. The production loader does
not generate VRS objects, normalize variants, fetch reference sequence, or
recompute identifier digests. Matching does not resolve identifiers through an
external service.

### 2. Identity ingestion: index the carried ALT alleles

The [default loader](../src/vrs_matcher/loader.py) maps each sample's positive
GT allele indexes to the corresponding annotated ALT IDs, applies FILTER and
available GQ/DP thresholds, and emits carried-allele rows. The database key is
`(sample_id, vrs_id)` in `sample_allele`.

The row also retains the original GT string, record-level zygosity, chromosome,
one-based VCF position, quality fields, and source dataset label. These provide
context for the observation; chromosome and position are not the equality key
used by the matcher. This table does not preserve REF/ALT or a full VRS object,
and repeated rows for the same sample/ID can replace earlier row metadata.
It is an index, not a lossless VCF archive.

The default loader trusts the annotation strings and their ordering. Records
without usable annotations are skipped; it does not independently establish
biological equivalence or enforce the KING panel's identifier checks.
Incorrect IDs can therefore create false matches, and inconsistent annotation
can hide real matches.

### 3. Identity matching: compare allele sets

The [identity matcher](../src/vrs_matcher/matcher.py) retrieves VRS ID sets and
associated genotype states. It computes Jaccard similarity from their
intersection and union, and genotype-state concordance over shared IDs.
Shared-variant output reports those IDs. Candidate restrictions in the Python
API also use VRS IDs.

This supports comparisons of independently indexed observations, including
replicate submissions and reprocessed specimens. However, the result depends
on which alleles were measured and retained. A high overlap is evidence for
sample similarity, not a validated donor-identity decision. The
[identity guide](how-to-sample-identity-confirmation.md) describes that use case.

A crucial limitation is that this index contains **carried ALT alleles**:
`0/0` produces no allele row, and neither does a fully missing call. Absence
cannot distinguish reference from unobserved or filtered data. Some partially
called genotypes can contribute their called ALT. This representation cannot
supply the jointly called genotype counts needed for KING.

### 4. KING ingestion: use VRS plus explicit locus and genotype information

The [genotype loader](../src/vrs_matcher/genotypes.py) requires a declared panel
of human autosomal, biallelic SNVs. Each panel marker has an ALT VRS ID **and**
a reference sequence identifier, one-based VCF position, REF, and ALT.

For records on panel contigs, ingestion checks the declared contig sequence
digest. It locates panel markers using `(sequence_id, pos)`; eligible SNV
records must then agree with the panel's REF/ALT orientation and ALT VRS ID.
The supported `1`/`chr1` aliases are accepted only with matching declared
sequence identity. Identifier syntax and matching reference/annotation header
values are checked, but the underlying VRS digest is not recomputed.

Both representations are necessary here: the ID is the persistent marker key,
while the explicit locus and allele checks establish that the incoming genotype
is interpreted against the intended panel definition. The loader does not
silently swap REF/ALT or recode an incompatible record.

`genotype_marker` stores the panel definitions. `genotype_call` stores
`(sample_id, vrs_id, dosage)`, with dosage 0, 1, or 2 relative to the panel ALT.
A dosage-0 row uses the **ALT marker ID** as its key; it does not assert that the
sample carries that ALT allele. Missing, excluded, or absent calls have no row
and are never imputed as reference.

The [KING matcher](../src/vrs_matcher/king.py) joins observations by these marker
IDs in memory and calculates kinship from jointly called dosages. The estimator
uses genotype counts, not properties of the identifier strings. A consistently
harmonized VCF-based calculation on the same calls should therefore be compared
for numerical agreement; VRS itself does not change the KING equations.

See the [KING guide](how-to-king-robust.md) for input contracts and the
[index description](index-description.md) for the schema and read/write lifecycle.

## Boundaries that matter biologically

- **No automatic cross-assembly matching.** Reference sequence identity is part
  of the allele representation. This repository performs no liftover or
  cross-reference equivalence mapping. Harmonize inputs upstream; KING requires
  a compatible declared reference and panel.
- **No automatic complex-variant equivalence.** A complex allele and several
  decomposed variants are not made interchangeable by string equality. The
  current KING path excludes indels and multiallelic records; acceptance of an
  ID by the identity loader is not general structural-variant support.
- **No removal of assay or QC effects.** Capture territory, coverage, caller
  behavior, missingness, and filtering still affect observed overlap and
  jointly called markers. Comparable IDs do not guarantee comparable evidence.
- **No genotype or haplotype inference from IDs.** An allele key carries neither
  sample dosage nor phase. The identity index keeps GT metadata, while the
  KING call table retains dosage without phase.
- **No substitute for provenance.** Keep the actual annotation software/version,
  reference snapshot, normalization policy, panel, and source VCFs. Merely
  supplying matching provenance labels cannot validate their truth.

## Choosing this approach

Use the VRS index when repeatedly comparing observations across files or
processing runs, and when shared allele identifiers are useful beyond one VCF.
For a single, already harmonized cohort, coordinate-and-allele-based analysis
can be entirely appropriate; VRS annotation adds preparation work and reference
data dependencies.

For fair evaluation, separate annotation, indexing, and repeated query costs,
and compare identical marker sets, calls, and QC policies. The
[VCFtools comparison](comparison-vcftools.md) and
[benchmark ADR](adr-benchmark-vcftools.md) explain those workload distinctions.
