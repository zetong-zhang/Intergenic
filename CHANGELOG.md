# Changelog

All notable project changes are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [0.0.1] - 2026-11-13

### Added

- FASTA, FASTQ, GenBank, and EMBL sequence input.
- GFF and BED annotation input, including per-sequence fallback to embedded
  GenBank/EMBL annotations.
- Linear and circular intergenic-region extraction.
- GFF, BED, and FASTA output with normalized `IGS000001` identifiers.
- Normalized origin-spanning coordinates and 30 bp adjacent-feature reports.
- Automated smoke and comprehensive test suites.
