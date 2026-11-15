# File formats and coordinate conventions

## Sequence input

Intergenic accepts:

- FASTA;
- FASTQ;
- GenBank;
- EMBL.

Sequence names are the first token after `>`/`@` for FASTA/FASTQ, the `LOCUS`
name for GenBank, and the `ID` name for EMBL.

Topology is read from `circular` or `linear` in the record header. Records
without an explicit topology use the value of `-c`.

## External annotation input

GFF uses sequence ID, type, start, end, score, and strand columns. GFF
coordinates are 1-based and inclusive.

BED requires at least three columns: sequence, start, and end. BED coordinates
are 0-based, end-exclusive. Optional BED score and strand columns are honored.

Target feature selection is controlled by `--types` and `--score`.

## Embedded annotation input

GenBank `FEATURES` and EMBL `FT` records are parsed from the sequence file when
the external annotation does not provide the sequence.

Single coordinates, `complement(...)`, `join(...)`, and `order(...)` are
supported. For a multi-segment location, the first through last coordinate
defines the masked span.

## Internal coordinates

Intergenic stores coordinates as 0-based, end-exclusive intervals. A wrapped
interval has `start > end` and represents:

```text
[start, sequence_length) followed by [0, end)
```

## GFF output

One row is emitted per IGR:

```text
host  INTERGENIC_v0.0.2  intergenic  start  end  .  .  .  ID=IGS000001;Adjacent=...
```

GFF coordinates are 1-based and inclusive. An origin-spanning region remains a
single row and intentionally has `start > end`.

## BED output

Five columns are emitted:

```text
host  start  end  IGS000001  Adjacent=...
```

Coordinates are 0-based and end-exclusive. Origin-spanning regions remain a
single row with `start > end`.

## FASTA output

The header format is:

```text
>host:start..end(strand) ID=IGS000001 Adjacent=...
```

Coordinates use the GFF 1-based, inclusive convention. The sequence for an
origin-spanning region is joined across the origin.

## IDs

IDs start at `IGS000001`. At least six digits are used, and the field expands
automatically when more than 999999 regions are emitted.

## Adjacent target features

A target feature is adjacent when its endpoint is within 30 bp of either IGR
boundary. All matching feature coordinates are listed in stable order:

```text
host:start..end(strand),host:start..end(strand)
```

Coordinates are 1-based and inclusive. `Adjacent=.` indicates no adjacent
target feature.

## Motif scan CSV

When motif scanning is enabled, Intergenic writes a separate CSV file:

```text
id,host,igr_start,igr_end,motif,strand,start,end,score,sequence
```

- `id`: ID of the containing IGR.
- `host`: sequence name.
- `igr_start`, `igr_end`: 1-based inclusive IGR coordinates.
- `motif`: MEME motif name.
- `strand`: `+` or `-`.
- `start`, `end`: 1-based inclusive genomic coordinates of the hit.
- `score`: sum of log2 likelihood ratios against the MEME background.
- `sequence`: matched bases in motif orientation.

The report contains a header even when no hit passes the threshold. Circular
hits may have `start > end`.
