# Examples

Sample data and plugins for trying out `vrs-matcher`.

| File                        | Notes                                                                     |
|-----------------------------|---------------------------------------------------------------------------|
| `example-cohort.vcf`        | Small VRS-annotated VCF (6 alleles, plain text) you can open and read.    |
| `example-cohort.vcf.gz`     | The same VCF, bgzip-compressed, as most genomics tooling expects it.      |
| `example-cohort.vcf.gz.tbi` | The tabix index for the `.gz` file.                                       |
| `vrs-matcher.ipynb`         | Notebook walking through the matching workflow.                           |
| `plugins/`                  | Example matcher plugins (see [`../docs/plugins.md`](../docs/plugins.md)). |

`load-samples` reads either the plain `.vcf` or the `.vcf.gz` (both load the same 6 records)
