# User story: GREGoR phenotypes for interpreting genetic matches

Status: Proposed; phenotype ingestion and reporting are not implemented.

## User story

As a rare-disease bioinformatician, I want to connect GREGoR phenotype assertions
to participants represented by indexed genomic observations, so that I can
review shared alleles and estimated relationships in clinical context, identify
cases worth investigating together, and revisit hypotheses as phenotyping
improves.

The scientific question is: **Do these genetic observations warrant joint
investigation given the participants' reported features and family context?**
Phenotype similarity alone should neither establish sample identity nor turn a
shared allele into a causal finding.

## Basis in GREGoR

GREGoR's data model connects families, participants, and molecular data, with
versioned releases intended to support joint analysis across workspaces.
The linked repository identifies its JSON as the authoritative model.
Sources: [GREGoR data model overview](https://gregorconsortium.org/data-model)
and [model repository](https://github.com/UW-GAC/gregor_data_models).

The authoritative JSON reviewed on 2026-09-30 declares version **1.12**.
Its `phenotype` table links assertions through `participant_id`. Required
assertion fields are `term_id`, `ontology`, and `presence`; presence distinguishes
`Present`, `Absent`, and `Unknown`. Supported vocabularies include HPO, MONDO,
OMIM, ORPHANET, SNOMED, and ICD10, with HPO required for phenotype description.
Optional fields include `additional_details`, `onset_age_range`,
`additional_modifiers`, and `syndromic`. `phenotype_id` is generated during
AnVIL deposition from participant, term, and presence, rather than supplied in
the uploaded TSV. Family identifiers, overall affected status, and age at last
observation belong to `participant`, not `phenotype`.
Source: [GREGoR JSON model](https://raw.githubusercontent.com/UW-GAC/gregor_data_models/main/GREGoR_data_model.json).

The requirements below are proposed behavior for VRS-Matcher, not claims that
GREGoR or this repository already provides these analyses. Implementation must
pin the actual imported model version and source snapshot; the linked `main`
branch can change.

## Scientific use cases

### Compare potentially informative cases across families

After finding a shared candidate VRS allele, an analyst should be able to
inspect the participants' recorded positive features, explicit negative
findings, and onset information side by side. Similar features in apparently
unrelated participants could motivate case review or functional follow-up.
Discordant features could motivate reconsidering the hypothesis or checking
phenotyping completeness.

The report must retain genetic evidence separately from phenotype evidence.
Common alleles can be shared by chance, and similar presentations can have
different causes. Gene-level convergence involving different alleles would
require an additional variant-to-gene annotation source; the phenotype table
and current allele matcher cannot establish it on their own.

### Interpret candidate relationships within families

Use the [KING results](how-to-king-robust.md), declared family links, and
phenotype assertions to select relatives for segregation review. An analyst
could inspect whether a feature is reported in relatives who carry a candidate
allele and whether relevant negative findings were actually recorded.

This supports hypothesis generation about variable expression or age-dependent
presentation. It does not estimate penetrance from a small, ascertained family.
A young relative without a recorded feature is not automatically an unaffected
control. KING establishes neither the direction of a relationship nor disease
status; candidate segregation requires explicit genotype evidence at the
candidate locus, beyond the background KING SNP panel.

### Review possible duplicate observations and metadata mismatches

When identity matching suggests two observations may represent the same donor,
show their participant mappings and phenotype source snapshots. A discrepancy
can prompt review of sample attribution, record linkage, or assessment timing.
It must not automatically relabel a sample: clinical records can differ because
of age, incomplete assessment, or updates.

Conversely, similar phenotypes must not cause genetically distinct observations
to be merged. Repeated sequencing of one participant must not count as multiple
independent cases supporting a disease hypothesis.

### Build reproducible phenotype-defined review cohorts

Allow an analyst to select participants with an explicitly present term and
inspect their available genetic observations. Initially, use exact ontology
and term matches, with explicit inclusion/exclusion rules. Report participants
with unknown status or no assertion separately from those with an explicit
negative assertion.

This enables targeted review without claiming a validated phenotype similarity
score. Parent/child ontology expansion and semantic ranking can be added later
with a pinned ontology release, a specified algorithm, and separate validation.
An exact-match result must disclose that it can miss differently granular terms.

### Revisit cases after phenotype updates

Permit a new phenotype snapshot to be associated with existing genetic results
without rebuilding allele or genotype calls. Show which assertions changed and
which snapshot informed a review. This supports reanalysis when new features
emerge or previously uncertain findings are clarified; it does not imply that
the underlying genotype evidence changed.

## Minimum proposed data and reporting contract

1. **Link observations explicitly.** Import a mapping from indexed `sample_id`
   to source-scoped `participant_id`. Multiple observations may map to one
   participant; one observation must not silently map to multiple participants.
   Do not assume matching strings prove a linkage. Preserve unmapped observations
   as unmapped, rather than attaching guessed phenotypes.
2. **Keep participant assertions separate from genotype rows.** Store phenotype
   snapshots and the mapping alongside the index, without copying assertions
   into each allele row. Import participant/family context only where required
   by the selected use case. The existing index has no participant registry;
   this linkage is new work. See the [index description](index-description.md).
3. **Preserve assertion meaning and provenance.** Retain original fields, source
   identity, model version, snapshot/file hash, and import time. Accept both
   deposition TSVs without generated IDs and exports with IDs. Keep conflicting
   assertions visible with provenance; never silently choose one or turn a
   missing row into an explicit absence. Import time is not assessment time.
4. **Report at the correct level.** Return observation IDs, linked participants,
   genetic scores and call counts, phenotype assertions, and mapping/data
   completeness. Count unique participants and, where available, families when
   reporting independent case support. Keep ontology namespaces distinct.
5. **Preserve genetic behavior.** Enrichment must leave identity and KING scores,
   rankings, and marker selection unchanged. Any explicit phenotype cohort
   filter must be recorded in the output. No combined diagnostic score is
   proposed. [VRS identifiers](why-use-vrs.md) remain allele keys, not phenotype
   or participant identifiers.
6. **Use authorized data boundaries.** Import and export only phenotype records
   available for the intended analysis under the source dataset's permissions
   and consent restrictions. Public schema availability does not make
   participant-level clinical data public.

## Acceptance scenarios

| Scenario | Expected result |
|---|---|
| Two sequencing observations map to one participant | Both show the participant's assertions; case counts include that participant once. |
| An observation has no mapping, or a mapping is ambiguous | Missing linkage is explicit; ambiguous linkage is rejected for publication until resolved. No guessed phenotype is attached. |
| One participant has a present term, another an absent assertion, a third unknown, and a fourth no row | All four states remain distinguishable; only the first enters an exact positive-term cohort. |
| Two sources disagree about the same participant's feature | Both assertions and their provenance remain reviewable; no silent overwrite or majority vote. |
| A source uses the same local participant ID as another source | Records remain separate unless an explicit reconciliation mapping connects them. |
| A deposition TSV lacks `phenotype_id` | Import succeeds using a deterministic local assertion identity and preserves the source fields; malformed required fields fail validation. |
| A query encounters different ontologies or a broader/narrower term | No unconfigured equivalence is inferred; exact-match limitations are reported. |
| A candidate allele is shared by relatives and unrelated cases | Output separates genetic sharing, family context, and phenotype evidence; it does not label the allele causal. |
| A phenotype snapshot changes while genetic data stay fixed | The new report identifies the snapshot change; genetic scores and call counts are identical. The earlier report remains reproducible. |
| A genetically similar pair has discordant phenotypes | The discrepancy is available for human review; neither sample identity nor participant mapping is automatically changed. |

Use synthetic participants and assertions for automated acceptance tests.
Scientific evaluation should then ask analysts to review authorized, curated
case sets for useful retrieval, missed cases due to term granularity, and
incorrect joins. Include relatives, repeated observations, sparse phenotyping,
and explicit negative findings. Do not use phenotype agreement as ground truth
for genetic identity or causal variant validation.

## Definition of done

A reviewer can reproduce a phenotype-enriched genetic comparison from pinned
inputs, trace every assertion to its source and participant mapping, and
distinguish genetic evidence from clinical context. The acceptance scenarios
pass, and genotype-only workflows retain their existing results. Automated
diagnosis, causal classification, semantic phenotype ranking, and full GREGoR
molecular-table ingestion remain outside this initial story.
