# Intergenic

`Intergenic` extracts intergenic regions from FASTA, FASTQ, GenBank, and EMBL
sequences. Annotations may come from an external GFF/BED file or from embedded
GenBank/EMBL records. The output format is selected from the output file
extension:

- `.gff` or `.gff3`: GFF3
- `.bed`: BED
- `.fa`, `.fasta`, or `.fna`: FASTA

When no output file is given, GFF3 is written to standard output.

## Documentation

- [Complete usage and CLI reference](docs/usage.md)
- [Input/output and coordinate contract](docs/file-formats.md)

## Build

Requirements: a C++17 compiler. CMake 3.16 or newer is recommended for
cross-platform builds; CMake Presets require 3.21 or newer. GNU Make or MinGW
Make remains available as a lightweight alternative. zlib development files
are optional and only needed for `.gz` input.

The build is configured for Linux, macOS, and Windows with GCC, Clang,
AppleClang, and MSVC-compatible CMake generators.

### CMake

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

On multi-configuration generators such as Visual Studio, omit
`-DCMAKE_BUILD_TYPE` if preferred; `--config Release` selects the build type.
CMake Presets are available:

```sh
cmake --preset default
cmake --build --preset default
```

Enable gzip input with CMake:

```sh
cmake -S . -B build-zlib -DINTERGENIC_USE_ZLIB=ON
cmake --build build-zlib --config Release
```

Install with:

```sh
cmake --install build --config Release --prefix /your/install/prefix
```

### Make

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
./intergenic -i genome.fa -a annotation.bed -f BED -o intergenic.bed
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

Continuous integration covers Linux, macOS, and Windows in
[.github/workflows/ci.yml](.github/workflows/ci.yml).

## License

Intergenic is distributed under the GNU General Public License, version 3.
See [LICENSE](LICENSE) for the project and third-party notices and
[COPYING](COPYING) for the complete license text.
