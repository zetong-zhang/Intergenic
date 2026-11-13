# Intergenic

`Intergenic` extracts intergenic regions from FASTA, FASTQ, GenBank, and EMBL
sequences. Annotations may come from an external GFF/BED file or from embedded
GenBank/EMBL records. The output format is selected from the output file
extension:

- `.gff` or `.gff3`: GFF3
- `.bed`: BED
- `.fa`, `.fasta`, or `.fna`: FASTA

When no output file is given, GFF3 is written to standard output.

## Build

```sh
make
```

On Windows with MinGW:

```powershell
mingw32-make
```

Gzip input is supported when zlib is available:

```sh
make USE_ZLIB=1
```

## Usage

```sh
./intergenic -i example.fa -a example.gff -o intergenic.gff
./intergenic -i example.fa -a example.bed -f BED -o intergenic.bed
./intergenic -i example.fa -a example.gff -o intergenic.fa
```

Use `-l` to set the minimum IGR length, `-t` to choose annotation feature
types, `-s` to set a score threshold, and `-c` to treat records without an
explicit topology as circular.

GFF coordinates are 1-based and inclusive. BED and internal coordinates are
0-based and end-exclusive. Every IGR receives one normalized ID such as
`IGS000001`; IDs are used in the GFF attribute, the BED name column, and as the
second field of the FASTA header. The first FASTA header field is the coordinate
form `<host>:<start>..<end>(<strand>)`, for example
`>NC_000913.3:5021..5233(.) ID=IGS000001`.

A circular IGR crossing the origin remains a single record. Its GFF/BED
coordinates intentionally have `start > end`, while its FASTA sequence is
joined across the origin.

Annotation precedence is resolved independently for every sequence. If the
external GFF/BED file contains the sequence name, that sequence uses only the
external annotations. Otherwise, GenBank `FEATURES` or EMBL `FT` annotations
are used as the fallback. FASTA/FASTQ records without an external entry have no
annotations.

For multi-segment locations such as `join(...)` or `order(...)`, the first
through last coordinates define the masked interval.

Every IGR output includes `Adjacent=<coordinates>` for target features whose
endpoint distance to either IGR boundary is at most 30 bp. Each item uses
`host:start..end(strand)`, items are separated by commas, and `.` means no
adjacent target feature. The field is a GFF attribute, the fifth BED column,
and the third FASTA header field.

## Test

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tests/smoke.ps1
```

## License

Intergenic is distributed under the GNU General Public License, version 3.
See [LICENSE](LICENSE) for the project and third-party notices.
