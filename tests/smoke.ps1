param(
    [string]$Binary = ""
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($Binary)) {
    $Binary = Join-Path $root "intergenic.exe"
}
$candidate = $Binary
if (-not (Test-Path -LiteralPath $candidate) -and
        [string]::IsNullOrEmpty([IO.Path]::GetExtension($candidate))) {
    $candidate += ".exe"
}
if (-not (Test-Path -LiteralPath $candidate)) {
    throw "Intergenic executable not found: $Binary"
}
$Binary = [IO.Path]::GetFullPath($candidate)

$temp = Join-Path ([IO.Path]::GetTempPath()) ("intergenic-smoke-" + [guid]::NewGuid())
New-Item -ItemType Directory -Path $temp | Out-Null

function Assert-Equal([object]$Actual, [object]$Expected, [string]$Message) {
    if ($Actual -ne $Expected) {
        throw "$Message (expected: $Expected, actual: $Actual)"
    }
}

function Invoke-Intergenic([string[]]$Arguments) {
    & $Binary -q @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Intergenic failed with exit code $LASTEXITCODE"
    }
}

function Read-FastaRecords([string]$Path) {
    $text = Get-Content -LiteralPath $Path -Raw
    $records = @()
    foreach ($match in [regex]::Matches($text, "(?m)^>(.*)\r?\n([^>]*)")) {
        $records += [pscustomobject]@{
            Header = $match.Groups[1].Value
            Sequence = ($match.Groups[2].Value -replace "\s", "")
        }
    }
    return $records
}

try {
    $sequence = "AAAACCCCGGGGTTTTACGT"
    $gff = Join-Path $temp "features.gff"
    Set-Content -LiteralPath $gff -Encoding ascii -Value @(
        "##gff-version 3"
        "chr`tTEST`tCDS`t3`t8`t.`t+`t0`tID=g1"
        "chr`tTEST`tCDS`t16`t18`t0.5`t-`t0`tID=g2"
    )

    $circular = Join-Path $temp "circular.fa"
    Set-Content -LiteralPath $circular -Encoding ascii -Value @(
        ">chr circular"
        $sequence
    )

    $gffout = Join-Path $temp "circular.gff"
    Invoke-Intergenic @("-i", $circular, "-a", $gff, "-l", "3", "-o", $gffout)
    $rows = @(Get-Content -LiteralPath $gffout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $rows.Count 2 "circular GFF row count"
    Assert-Equal ($rows[0].Split("`t")[3,4] -join ",") "19,2" "origin-crossing coordinates"
    Assert-Equal ($rows[0].Split("`t")[8]) "ID=IGS000001;Adjacent=chr:3..8(+),chr:16..18(-)" "origin-crossing ID"
    Assert-Equal ($rows[1].Split("`t")[3,4] -join ",") "9,15" "internal IGR"
    Assert-Equal ($rows[1].Split("`t")[8]) "ID=IGS000002;Adjacent=chr:3..8(+),chr:16..18(-)" "internal IGR ID"

    $bedout = Join-Path $temp "circular.bed"
    Invoke-Intergenic @("-i", $circular, "-a", $gff, "-l", "3", "-o", $bedout)
    $bed = @(Get-Content -LiteralPath $bedout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $bed.Count 2 "circular BED row count"
    Assert-Equal (($bed[0].Split("`t"))[1,2] -join ",") "18,2" "circular BED coordinates"
    Assert-Equal $bed[0].Split("`t")[3] "IGS000001" "circular BED ID"
    Assert-Equal $bed[0].Split("`t")[4] "Adjacent=chr:3..8(+),chr:16..18(-)" "circular BED adjacency"

    $fastaout = Join-Path $temp "circular-igr.fa"
    Invoke-Intergenic @("-i", $circular, "-a", $gff, "-l", "3", "-o", $fastaout)
    $records = @(Read-FastaRecords $fastaout)
    Assert-Equal $records.Count 2 "circular FASTA record count"
    $headerFields = @($records[0].Header -split "\s+")
    Assert-Equal $headerFields[0] "chr:19..2(.)" "circular FASTA coordinates"
    Assert-Equal $headerFields[1] "ID=IGS000001" "circular FASTA ID"
    Assert-Equal $headerFields[2] "Adjacent=chr:3..8(+),chr:16..18(-)" "circular FASTA adjacency"
    Assert-Equal $records[0].Sequence "GTAA" "origin-crossing sequence"
    Assert-Equal $records[1].Sequence "GGGGTTT" "internal IGR sequence"

    $linear = Join-Path $temp "linear.fa"
    Set-Content -LiteralPath $linear -Encoding ascii -Value @(
        ">chr linear"
        $sequence
    )
    $linearout = Join-Path $temp "linear.gff"
    Invoke-Intergenic @("-i", $linear, "-a", $gff, "-l", "3", "-o", $linearout)
    $linearRows = @(Get-Content -LiteralPath $linearout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $linearRows.Count 1 "linear GFF row count"
    Assert-Equal ($linearRows[0].Split("`t")[3,4] -join ",") "9,15" "linear IGR"

    $default = Join-Path $temp "default.fa"
    Set-Content -LiteralPath $default -Encoding ascii -Value @(
        ">chr"
        $sequence
    )
    $defaultout = Join-Path $temp "default.gff"
    Invoke-Intergenic @("-i", $default, "-a", $gff, "-l", "3", "-c", "-o", $defaultout)
    $defaultRows = @(Get-Content -LiteralPath $defaultout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $defaultRows.Count 2 "circular CLI row count"

    $genbank = Join-Path $temp "record.gbk"
    Set-Content -LiteralPath $genbank -Encoding ascii -Value @(
        "LOCUS       chr 20 bp DNA linear"
        "ORIGIN"
        "        1 aaaaccccgg ggttttacgt"
        "//"
    )
    $genbankout = Join-Path $temp "genbank.gff"
    Invoke-Intergenic @("-i", $genbank, "-a", $gff, "-l", "3", "-o", $genbankout)
    $genbankRows = @(Get-Content -LiteralPath $genbankout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $genbankRows.Count 1 "GenBank GFF row count"
    Assert-Equal ($genbankRows[0].Split("`t")[3,4] -join ",") "9,15" "GenBank IGR"

    $embedded = Join-Path $temp "embedded.gbk"
    Set-Content -LiteralPath $embedded -Encoding ascii -Value @(
        "LOCUS       embedded 20 bp DNA circular"
        "FEATURES             Location/Qualifiers"
        "     source          1..20"
        "     CDS             3..8"
        "     CDS             complement(16..18)"
        "ORIGIN"
        "        1 aaaaccccgg ggttttacgt"
        "//"
    )
    $embeddedout = Join-Path $temp "embedded.gff"
    Invoke-Intergenic @("-i", $embedded, "-l", "3", "-o", $embeddedout)
    $embeddedRows = @(Get-Content -LiteralPath $embeddedout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $embeddedRows.Count 2 "GenBank embedded annotation row count"
    Assert-Equal ($embeddedRows[0].Split("`t")[3,4] -join ",") "19,2" "GenBank embedded crossing"

    $overrideGff = Join-Path $temp "override.gff"
    Set-Content -LiteralPath $overrideGff -Encoding ascii -Value @(
        "chr1`tTEST`tCDS`t3`t8`t.`t+`t0`tID=g1"
    )
    $mixed = Join-Path $temp "mixed.gbk"
    Set-Content -LiteralPath $mixed -Encoding ascii -Value @(
        "LOCUS       chr1 20 bp DNA circular"
        "FEATURES             Location/Qualifiers"
        "     CDS             3..8"
        "     CDS             16..18"
        "ORIGIN"
        "        1 aaaaccccgg ggttttacgt"
        "//"
        "LOCUS       chr2 20 bp DNA circular"
        "FEATURES             Location/Qualifiers"
        "     CDS             3..8"
        "     CDS             16..18"
        "ORIGIN"
        "        1 aaaaccccgg ggttttacgt"
        "//"
    )
    $mixedout = Join-Path $temp "mixed.gff"
    Invoke-Intergenic @("-i", $mixed, "-a", $overrideGff, "-l", "3", "-o", $mixedout)
    $mixedRows = @(Get-Content -LiteralPath $mixedout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $mixedRows.Count 3 "mixed annotation row count"
    Assert-Equal @($mixedRows | Where-Object { $_.StartsWith("chr1`t") }).Count 1 "GFF overrides chr1"
    Assert-Equal @($mixedRows | Where-Object { $_.StartsWith("chr2`t") }).Count 2 "embedded annotates chr2"
    Assert-Equal ($mixedRows[0].Split("`t")[3,4] -join ",") "9,2" "GFF override replaces embedded CDS"

    $overrideBed = Join-Path $temp "override.bed"
    Set-Content -LiteralPath $overrideBed -Encoding ascii -Value @(
        "chr1`t2`t8"
    )
    $mixedBedOut = Join-Path $temp "mixed-bed.gff"
    Invoke-Intergenic @("-i", $mixed, "-a", $overrideBed, "-f", "BED", "-l", "3", "-o", $mixedBedOut)
    $mixedBedRows = @(Get-Content -LiteralPath $mixedBedOut | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $mixedBedRows.Count 3 "mixed BED annotation row count"
    Assert-Equal ($mixedBedRows[0].Split("`t")[3,4] -join ",") "9,2" "BED override replaces embedded CDS"

    $embl = Join-Path $temp "record.embl"
    Set-Content -LiteralPath $embl -Encoding ascii -Value @(
        "ID   chr; SV 1; linear; DNA; STD; UNC; 20 BP."
        "SQ   Sequence 20 BP;"
        "     aaaaccccgg ggttttacgt"
        "//"
    )
    $emblout = Join-Path $temp "embl.bed"
    Invoke-Intergenic @("-i", $embl, "-a", $gff, "-l", "3", "-o", $emblout)
    $emblRows = @(Get-Content -LiteralPath $emblout | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $emblRows.Count 1 "EMBL BED row count"
    Assert-Equal (($emblRows[0].Split("`t"))[1,2] -join ",") "8,15" "EMBL IGR"

    $embeddedEmbl = Join-Path $temp "embedded.embl"
    Set-Content -LiteralPath $embeddedEmbl -Encoding ascii -Value @(
        "ID   embedded; SV 1; circular; DNA; STD; UNC; 20 BP."
        "FT   source          1..20"
        "FT   CDS             3..8"
        "FT   CDS             complement(16..18)"
        "SQ   Sequence 20 BP;"
        "     aaaaccccgg ggttttacgt"
        "//"
    )
    $embeddedEmblOut = Join-Path $temp "embedded-embl.gff"
    Invoke-Intergenic @("-i", $embeddedEmbl, "-l", "3", "-o", $embeddedEmblOut)
    $embeddedEmblRows = @(Get-Content -LiteralPath $embeddedEmblOut | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $embeddedEmblRows.Count 2 "EMBL embedded annotation row count"
    Assert-Equal ($embeddedEmblRows[0].Split("`t")[3,4] -join ",") "19,2" "EMBL embedded crossing"

    $multi = Join-Path $temp "multi.fa"
    Set-Content -LiteralPath $multi -Encoding ascii -Value @(
        ">chr1"
        "ACGTACGT"
        ">chr2"
        "TTTTGGGG"
    )
    $multiout = Join-Path $temp "multi-igr.fa"
    Invoke-Intergenic @("-i", $multi, "-l", "1", "-o", $multiout)
    $multiRecords = @(Read-FastaRecords $multiout)
    Assert-Equal $multiRecords.Count 2 "multi-FASTA record count"
    Assert-Equal ($multiRecords[0].Header -split "\s+")[1] "ID=IGS000001" "first FASTA ID"
    Assert-Equal ($multiRecords[1].Header -split "\s+")[1] "ID=IGS000002" "second FASTA ID"
    Assert-Equal $multiRecords[0].Sequence "ACGTACGT" "first FASTA sequence"
    Assert-Equal $multiRecords[1].Sequence "TTTTGGGG" "second FASTA sequence"

    $fastq = Join-Path $temp "single.fastq"
    Set-Content -LiteralPath $fastq -Encoding ascii -Value @(
        "@chr linear"
        $sequence
        "+"
        ("I" * $sequence.Length)
    )
    $fastqout = Join-Path $temp "fastq.fa"
    Invoke-Intergenic @("-i", $fastq, "-l", "1", "-o", $fastqout)
    $fastqRecords = @(Read-FastaRecords $fastqout)
    Assert-Equal $fastqRecords.Count 1 "FASTQ record count"
    Assert-Equal ($fastqRecords[0].Header -split "\s+")[1] "ID=IGS000001" "FASTQ FASTA ID"
    Assert-Equal ($fastqRecords[0].Header -split "\s+")[2] "Adjacent=." "FASTQ empty adjacency"
    Assert-Equal $fastqRecords[0].Sequence $sequence "FASTQ sequence"

    $distanceGff = Join-Path $temp "adjacent-distance.gff"
    Set-Content -LiteralPath $distanceGff -Encoding ascii -Value @(
        "chr`tTEST`tCDS`t1`t100`t.`t+`t0`tID=outer"
        "chr`tTEST`tCDS`t1`t70`t.`t+`t0`tID=near"
        "chr`tTEST`tCDS`t1`t69`t.`t+`t0`tID=far"
    )
    $distanceFasta = Join-Path $temp "adjacent-distance.fa"
    Set-Content -LiteralPath $distanceFasta -Encoding ascii -Value @(
        ">chr linear"
        ("A" * 200)
    )
    $distanceOut = Join-Path $temp "adjacent-distance-out.gff"
    Invoke-Intergenic @("-i", $distanceFasta, "-a", $distanceGff, "-l", "90", "-o", $distanceOut)
    $distanceRows = @(Get-Content -LiteralPath $distanceOut | Where-Object { $_ -notmatch "^#" })
    Assert-Equal $distanceRows.Count 1 "adjacent distance row count"
    Assert-Equal $distanceRows[0].Split("`t")[8] "ID=IGS000001;Adjacent=chr:1..70(+),chr:1..100(+)" "30 bp adjacency"

    Write-Host "All smoke tests passed."
}
finally {
    $resolvedTemp = [IO.Path]::GetFullPath($temp)
    $tempRoot = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
    if ($resolvedTemp.StartsWith($tempRoot, [StringComparison]::OrdinalIgnoreCase) -and
            (Split-Path -Leaf $resolvedTemp).StartsWith("intergenic-smoke-")) {
        Remove-Item -LiteralPath $resolvedTemp -Recurse -Force
    }
}
