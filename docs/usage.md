# Usage

## Requirements

- A C++17 compiler;
- GNU Make or MinGW Make;
- PowerShell for the provided test scripts;
- zlib headers and libraries only for `.gz` input.

## Build

Linux and macOS:

```sh
make
make test
```

Windows with MinGW:

```powershell
mingw32-make
mingw32-make test
```

Enable gzip input:

```sh
make clean all USE_ZLIB=1
```

The executable is `intergenic` on Unix and `intergenic.exe` on Windows.

## Quick start

```sh
./intergenic -i example.fa -a example.gff -o intergenic.gff
```

Write BED or FASTA using the output extension:

```sh
./intergenic -i example.fa -a example.gff -o intergenic.bed
./intergenic -i example.fa -a example.gff -o intergenic.fa
```

Without `-o`, GFF is written to stdout.

## Options

| Option | Behavior |
| --- | --- |
| `-a, --anno FILE` | External GFF/BED annotation. Use `-` for stdin. |
| `-i, --input FILE` | FASTA, FASTQ, GenBank, or EMBL sequence input. Use `-` for stdin. |
| `-f, --format GFF\|BED` | External annotation format. Default: `GFF`. |
| `-o, --output FILE` | Output file. Format is selected by extension; default: stdout/GFF. |
| `-c, --circ` | Treat records without explicit topology as circular. |
| `-l, --minlen N` | Minimum IGR length. Default: `90`. |
| `-t, --types LIST` | Comma-separated target feature types. Default: `CDS`. |
| `-s, --score X` | Inclusive minimum target-feature score. Default: `0`. |
| `-q, --quiet` | Suppress informational stderr output. |
| `-h, --help` | Print help and exit. |
| `-v, --version` | Print version and exit. |

## Annotation precedence

Annotation lookup is performed independently for every sequence:

1. If the external GFF/BED contains the sequence name, only its matching
   features are used.
2. Otherwise, GenBank `FEATURES` or EMBL `FT` features are used.
3. FASTA/FASTQ records without an external entry have no target features and
   are treated as entirely unmasked.

If a sequence name is present externally but all of its features are removed by
`-t` or `-s`, Intergenic does not fall back to embedded annotations.

## stdin

Only one input may use stdin at a time.

```sh
cat genome.fa | ./intergenic -i - -a genes.gff -o intergenic.fa
cat genes.gff | ./intergenic -i genome.fa -a - -o intergenic.gff
```

Using both `-i -` and `-a -` is rejected.

## Output selection

- `.gff`, `.gff3`, and unknown extensions: GFF
- `.bed`: BED
- `.fa`, `.fasta`, `.fna`: FASTA

See [file-formats.md](file-formats.md) for the exact column and coordinate
contract.
