# Examples

Sample data and plugins for trying out `vrs-matcher`.

## Files

| File | What it is |
| --- | --- |
| `example-cohort.vcf` | Small VRS-annotated VCF (6 alleles, plain text) you can open and read. |
| `example-cohort.vcf.gz` | The same VCF, **bgzip**-compressed, as most genomics tooling expects it. |
| `example-cohort.vcf.gz.tbi` | The **tabix** index for the `.gz` file. |
| `vrs-matcher.ipynb` | Notebook walking through the matching workflow. |
| `plugins/` | Example matcher plugins (see [`../docs/plugins.md`](../docs/plugins.md)). |

`load-samples` reads either the plain `.vcf` or the `.vcf.gz`; both load the
same 6 records. The `.db` files you may see here are local scratch databases
built by `load-samples`. They are git-ignored, not part of the example set.

## Plain VCF vs. bgzip + tabix

- **`.vcf` (plain text)** is fine for small examples. You can `cat` it, and
  `load-samples` reads it directly.
- **`.vcf.gz` (bgzip)** is what most genomics tooling expects. It must be
  compressed with `bgzip`, not plain `gzip`. `bgzip` produces a
  block-compressed file (BGZF) that tools can seek into; a plain `gzip` file
  looks similar but `tabix` rejects it.
- **`.tbi` (tabix index)** lets tools jump to a genomic region without reading
  the whole file. Build it with `tabix` from a bgzipped VCF.

The full chain is `bgzip` the VCF, then `tabix` the result, leaving you with
`.vcf.gz` and `.vcf.gz.tbi`.

### Regenerating the compressed files

If you edit `example-cohort.vcf`, rebuild the `.gz` and `.tbi` so they stay in
sync:

```bash
bgzip -f -k example-cohort.vcf       # -k keeps the plain .vcf
tabix -f example-cohort.vcf.gz       # writes example-cohort.vcf.gz.tbi
```

### A header note (the `##contig` warning)

A VCF header should declare every contig its records use, e.g.:

```text
##contig=<ID=chr1,length=248956422>
```

If that line is missing, htslib-based tools (including `load-samples`) still
work but print a harmless warning like
`Contig 'chr1' is not defined in the header`. The bundled VCF declares its
contig, so you won't see this. Add the line if you hit the warning with your
own data.
