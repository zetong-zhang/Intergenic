# Software design

Intergenic separates sequence parsing, annotation selection, interval
extraction, and formatting.

## Parsing

`scaffold` reads one sequence record at a time. FASTA headers are held as a
pending line rather than requiring random access, so multi-record FASTA input
does not require loading and indexing the whole file.

GenBank/EMBL feature tables are parsed while scanning to the sequence section.
External GFF/BED loading is handled separately by `feature_loader`.

## Annotation selection

`feature_loader` records every external sequence name independently of type and
score filtering. This allows `resolve()` to make a sequence-level decision:
external annotation wins when the name exists; embedded annotation is used
otherwise.

## Interval extraction

`scaffold::extract` builds a boolean mask for each sequence, collects unmasked
runs, merges runs touching both sequence ends for circular records, and applies
the minimum length after merging. Coordinates remain signed so a wrapped run can
be represented without splitting it.

## Output

GFF, BED, and FASTA serializers consume the same `region` objects. This ensures
that IDs, adjacency metadata, region length, and origin-spanning behavior remain
consistent across formats.

## Extension points

- Add annotation formats in `feature_loader`.
- Add sequence formats in `scaffold::load`.
- Add output formats in the output serializer section of `Intergenic.cpp`.
- Extend comprehensive tests for every new coordinate or fallback rule.
