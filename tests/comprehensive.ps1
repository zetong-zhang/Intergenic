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

$temp = Join-Path ([IO.Path]::GetTempPath()) (
    "intergenic-comprehensive-" + [guid]::NewGuid()
)
New-Item -ItemType Directory -Path $temp | Out-Null
$script:runCounter = 0
Get-Random -SetSeed 20261114 | Out-Null

function Assert-Equal([object]$Actual, [object]$Expected, [string]$Message) {
    if ($Actual -ne $Expected) {
        throw "$Message (expected: $Expected, actual: $Actual)"
    }
}

function Assert-True([bool]$Condition, [string]$Message) {
    if (-not $Condition) {
        throw $Message
    }
}

function Assert-Contains([string]$Text, [string]$Expected, [string]$Message) {
    if (-not $Text.Contains($Expected)) {
        throw "$Message (missing: $Expected)"
    }
}

function Invoke-Intergenic {
    param(
        [string[]]$Arguments = @(),
        [AllowNull()][string[]]$InputLines = $null
    )

    $script:runCounter ++
    $start = [Diagnostics.ProcessStartInfo]::new()
    $start.FileName = $Binary
    $start.UseShellExecute = $false
    $start.CreateNoWindow = $true
    $start.RedirectStandardOutput = $true
    $start.RedirectStandardError = $true
    $start.RedirectStandardInput = $null -ne $InputLines
    $start.Arguments = (@($Arguments | ForEach-Object {
        '"' + ($_ -replace '"', '\"') + '"'
    }) -join " ")
    $process = [Diagnostics.Process]::new()
    $process.StartInfo = $start
    $process.Start() | Out-Null
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    if ($null -ne $InputLines) {
        $process.StandardInput.Write(($InputLines -join "`r`n") + "`r`n")
        $process.StandardInput.Close()
    }
    $process.WaitForExit()
    return [pscustomobject]@{
        ExitCode = $process.ExitCode
        Stdout = $stdoutTask.Result
        Stderr = $stderrTask.Result
    }
}

function Assert-Success([string[]]$Arguments, [AllowNull()][string[]]$InputLines = $null) {
    $result = Invoke-Intergenic -Arguments $Arguments -InputLines $InputLines
    if ($result.ExitCode -ne 0) {
        throw "Command failed: $($Arguments -join ' ')`n$($result.Stderr)"
    }
    return $result
}

function Assert-Failure(
    [string[]]$Arguments,
    [string]$ExpectedMessage,
    [string]$Label,
    [AllowNull()][string[]]$InputLines = $null
) {
    $result = Invoke-Intergenic -Arguments $Arguments -InputLines $InputLines
    Assert-True ($result.ExitCode -ne 0) "Expected failure: $Label"
    Assert-Contains ($result.Stdout + $result.Stderr) $ExpectedMessage "$Label error"
}

function Write-Fixture([string]$Name, [string[]]$Lines) {
    $path = Join-Path $temp $Name
    if ($null -eq $Lines) { $Lines = @() }
    [IO.File]::WriteAllLines($path, $Lines, [Text.Encoding]::ASCII)
    return $path
}

function Read-Gff([string]$Path) {
    $rows = @()
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line.StartsWith("#") -or [string]::IsNullOrWhiteSpace($line)) {
            continue
        }
        $fields = $line.Split("`t")
        Assert-Equal $fields.Count 9 "GFF column count"
        $id = if ($fields[8] -match "ID=([^;]+)") { $Matches[1] } else { "" }
        $adjacent = if ($fields[8] -match "Adjacent=([^;]*)") {
            $Matches[1]
        } else { "" }
        $rows += [pscustomobject]@{
            Host = $fields[0]
            Source = $fields[1]
            Type = $fields[2]
            Start = [int]$fields[3]
            End = [int]$fields[4]
            Strand = $fields[6]
            Attributes = $fields[8]
            ID = $id
            Adjacent = $adjacent
        }
    }
    return $rows
}

function Read-Bed([string]$Path) {
    $rows = @()
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line.StartsWith("#") -or [string]::IsNullOrWhiteSpace($line)) {
            continue
        }
        $fields = $line.Split("`t")
        Assert-Equal $fields.Count 5 "BED column count"
        $rows += [pscustomobject]@{
            Host = $fields[0]
            Start = [int]$fields[1]
            End = [int]$fields[2]
            ID = $fields[3]
            Adjacent = $fields[4].Substring("Adjacent=".Length)
        }
    }
    return $rows
}

function Read-Fasta([string]$Path) {
    $text = Get-Content -LiteralPath $Path -Raw
    $records = @()
    foreach ($match in [regex]::Matches($text, "(?m)^>(.*)\r?\n([^>]*)")) {
        $records += [pscustomobject]@{
            Header = $match.Groups[1].Value
            Sequence = ($match.Groups[2].Value -replace "\s", "")
            Fields = @($match.Groups[1].Value -split "\s+")
        }
    }
    return $records
}

function Get-ReferenceRuns(
    [int]$Length,
    [object[]]$Features,
    [bool]$Circular,
    [int]$MinLength
) {
    $masked = New-Object "bool[]" $Length
    foreach ($item in $Features) {
        for ($i = $item.Start; $i -lt $item.End; $i ++) {
            $masked[$i] = $true
        }
    }

    $runs = @()
    $i = 0
    while ($i -lt $Length) {
        while ($i -lt $Length -and $masked[$i]) { $i ++ }
        if ($i -eq $Length) { break }
        $start = $i
        while ($i -lt $Length -and -not $masked[$i]) { $i ++ }
        $runs += [pscustomobject]@{ Start = $start; End = $i }
    }

    if ($Circular -and $runs.Count -ge 2 -and $runs[0].Start -eq 0 -and
            $runs[-1].End -eq $Length) {
        $runs[0].Start = $runs[-1].Start
        $runs = @($runs | Select-Object -First ($runs.Count-1))
    }

    $expected = @()
    foreach ($run in $runs) {
        $size = if ($run.Start -le $run.End) {
            $run.End-$run.Start
        } else {
            $run.End+$Length-$run.Start
        }
        if ($size -ge $MinLength) {
            $expected += $run
        }
    }
    return $expected
}

function Merge-Intervals([object[]]$Intervals) {
    if ($Intervals.Count -eq 0) { return @() }
    $sorted = @($Intervals | Sort-Object Start, End)
    $merged = @([pscustomobject]@{
        Start = $sorted[0].Start
        End = $sorted[0].End
    })
    foreach ($item in $sorted | Select-Object -Skip 1) {
        if ($item.Start -le $merged[-1].End) {
            if ($item.End -gt $merged[-1].End) {
                $merged[-1].End = $item.End
            }
        } else {
            $merged += [pscustomobject]@{
                Start = $item.Start
                End = $item.End
            }
        }
    }
    return $merged
}

function Test-Overlap([object[]]$Intervals, [int]$Start, [int]$End) {
    if ($Start -ge $End -or $Intervals.Count -eq 0) { return $false }
    $low = 0
    $high = $Intervals.Count-1
    $candidate = -1
    while ($low -le $high) {
        $middle = [int](($low+$high)/2)
        if ($Intervals[$middle].Start -lt $End) {
            $candidate = $middle
            $low = $middle+1
        } else {
            $high = $middle-1
        }
    }
    return $candidate -ge 0 -and $Intervals[$candidate].End -gt $Start
}

try {
    # CLI behavior and exit codes.
    $version = Assert-Success @("--version")
    Assert-Equal $version.Stdout.Trim() "INTERGENIC_v0.0.1" "Version output"
    $help = Assert-Success @("--help")
    Assert-Contains $help.Stdout "INTERGENIC: AN INTERGENIC REGION EXTRACTION TOOLKIT" "Help banner"
    Assert-Contains $help.Stdout "--minlen" "Help minlen"
    $noArgs = Assert-Success @()
    Assert-Contains $noArgs.Stdout "Usage:" "No-argument help"
    Assert-Failure @("--bogus") "Error:" "Unknown option"

    $sequence = "AAAACCCCGGGGTTTTACGT"
    $basicGff = Write-Fixture "basic.gff" @(
        "##gff-version 3"
        "chr`tTEST`tCDS`t3`t8`t.`t+`t0`tID=g1"
        "chr`tTEST`tCDS`t16`t18`t0.5`t-`t0`tID=g2"
    )
    $circularFa = Write-Fixture "circular.fa" @(">chr circular", $sequence)

    # Missing files, invalid values and conflicting stdin.
    Assert-Failure @("-i", (Join-Path $temp "missing.fa")) "cannot open genome file" "Missing genome"
    Assert-Failure @("-i", $circularFa, "-l", "-1") "cannot be negative" "Negative minlen"
    Assert-Failure @("-i", $circularFa, "-a", $basicGff, "-f", "XYZ") "unsupported annotation" "Bad annotation format"
    Assert-Failure @("-i", $circularFa, "-a", (Join-Path $temp "missing.gff")) "cannot open annotation file" "Missing annotation"
    Assert-Failure @("-i", $circularFa, "-o", $temp) "cannot open output file" "Output is directory"
    Assert-Failure @("-i", "-", "-a", "-") "cannot both be read from stdin" "Conflicting stdin" -InputLines @("") | Out-Null

    # Default minimum length is exactly 90.
    $fa89 = Write-Fixture "len89.fa" @(">chr linear", ("A"*89))
    $fa90 = Write-Fixture "len90.fa" @(">chr linear", ("A"*90))
    $out89 = Join-Path $temp "len89.gff"
    $out90 = Join-Path $temp "len90.gff"
    Assert-Success @("-i", $fa89, "-o", $out89) | Out-Null
    Assert-Success @("-i", $fa90, "-o", $out90) | Out-Null
    Assert-Equal @(Read-Gff $out89).Count 0 "89 bp below default"
    $rows90 = @(Read-Gff $out90)
    Assert-Equal $rows90.Count 1 "90 bp at default"
    Assert-Equal "$($rows90[0].Start),$($rows90[0].End)" "1,90" "Default region bounds"

    # Output format inference, normalized fields and stdout fallback.
    $fa100 = Write-Fixture "len100.fa" @(">chr linear", ("A"*100))
    $outGff3 = Join-Path $temp "result.gff3"
    $outBed = Join-Path $temp "result.BED"
    $outFasta = Join-Path $temp "result.FASTA"
    $outUnknown = Join-Path $temp "result.data"
    Assert-Success @("-i", $fa100, "-o", $outGff3) | Out-Null
    Assert-Success @("-i", $fa100, "-o", $outBed) | Out-Null
    Assert-Success @("-i", $fa100, "-o", $outFasta) | Out-Null
    Assert-Success @("-i", $fa100, "-o", $outUnknown) | Out-Null
    $gff100 = @(Read-Gff $outGff3)
    Assert-Equal $gff100.Count 1 "gff3 extension"
    Assert-Equal $gff100[0].Source "INTERGENIC_v0.0.1" "GFF source"
    Assert-Equal $gff100[0].ID "IGS000001" "GFF ID"
    Assert-Equal $gff100[0].Adjacent "." "No adjacent features"
    Assert-Equal @(Read-Bed $outBed).Count 1 "uppercase BED extension"
    Assert-Equal @(Read-Fasta $outFasta).Count 1 "FASTA extension"
    Assert-Equal @(Read-Gff $outUnknown).Count 1 "Unknown extension defaults to GFF"
    $stdoutGff = Assert-Success @("-i", $fa100)
    Assert-Contains $stdoutGff.Stdout "`tINTERGENIC_v0.0.1`tintergenic`t" "stdout GFF"

    # Multiple scaffolds keep one ID sequence across records.
    $multiFa = Write-Fixture "multi.fa" @(
        ">chr1"
        "AAAA"
        ">chr2"
        "CCCC"
    )
    $multiOut = Join-Path $temp "multi-out.fa"
    Assert-Success @("-i", $multiFa, "-l", "1", "-o", $multiOut) | Out-Null
    $multiRecords = @(Read-Fasta $multiOut)
    Assert-Equal $multiRecords.Count 2 "Multiple scaffold count"
    Assert-Equal $multiRecords[0].Fields[1] "ID=IGS000001" "First scaffold ID"
    Assert-Equal $multiRecords[1].Fields[1] "ID=IGS000002" "Second scaffold ID"

    # Annotation type and score filtering.
    $filterGff = Write-Fixture "filters.gff" @(
        "chr`tTEST`tCDS`t10`t20`t0.4`t+`t0`tID=low"
        "chr`tTEST`tCDS`t40`t50`t0.6`t-`t0`tID=high"
        ("chr`tTEST`t" + "tRNA" + "`t70`t80`t1.0`t+`t0`tID=rna")
    )
    $filterOut = Join-Path $temp "scores.gff"
    Assert-Success @("-i", $fa100, "-a", $filterGff, "-s", "0.5", "-l", "30", "-o", $filterOut) | Out-Null
    $scoreRows = @(Read-Gff $filterOut)
    Assert-Equal $scoreRows.Count 2 "Score filter row count"
    Assert-Equal "$($scoreRows[0].Start),$($scoreRows[0].End)" "1,39" "Low score excluded"
    Assert-Equal "$($scoreRows[1].Start),$($scoreRows[1].End)" "51,100" "High score retained"

    $typeOut = Join-Path $temp "types.gff"
    Assert-Success @("-i", $fa100, "-a", $filterGff, "-t", "CDS,tRNA", "-s", "0.5", "-l", "30", "-o", $typeOut) | Out-Null
    $typeRows = @(Read-Gff $typeOut)
    $typeBounds = @($typeRows | ForEach-Object { "$($_.Start),$($_.End)" }) -join ";"
    Assert-Equal $typeRows.Count 1 "Type filter row count: $typeBounds"
    Assert-Equal "$($typeRows[0].Start),$($typeRows[0].End)" "1,39" "Type filter result"

    $equalScoreOut = Join-Path $temp "equal-score.gff"
    Assert-Success @("-i", $fa100, "-a", $filterGff, "-s", "0.6", "-l", "30", "-o", $equalScoreOut) | Out-Null
    Assert-Equal @(Read-Gff $equalScoreOut).Count 2 "Score threshold is inclusive"

    # Genome and annotation stdin.
    $stdinGffOut = Join-Path $temp "stdin-annotation.gff"
    $annotationLines = @(Get-Content -LiteralPath $basicGff)
    Assert-Success @("-i", $circularFa, "-a", "-", "-l", "3", "-o", $stdinGffOut) $annotationLines | Out-Null
    Assert-Equal @(Read-Gff $stdinGffOut).Count 2 "Annotation from stdin"

    $stdinFaOut = Join-Path $temp "stdin-genome.gff"
    $genomeLines = @(Get-Content -LiteralPath $circularFa)
    Assert-Success @("-i", "-", "-a", $basicGff, "-l", "3", "-o", $stdinFaOut) $genomeLines | Out-Null
    Assert-Equal @(Read-Gff $stdinFaOut).Count 2 "Genome from stdin"

    # GenBank multi-segment locations and external-over-embedded precedence.
    $joined = Write-Fixture "joined.gbk" @(
        "LOCUS       chr 20 bp DNA circular"
        "FEATURES             Location/Qualifiers"
        "     CDS             complement(join(3..8,16..18))"
        "ORIGIN"
        "        1 aaaaccccgg ggttttacgt"
        "//"
    )
    $joinedOut = Join-Path $temp "joined.gff"
    Assert-Success @("-i", $joined, "-l", "3", "-o", $joinedOut) | Out-Null
    $joinedRows = @(Read-Gff $joinedOut)
    Assert-Equal $joinedRows.Count 1 "Join feature count"
    Assert-Equal "$($joinedRows[0].Start),$($joinedRows[0].End)" "19,2" "Join bounds"

    $externalNoCds = Write-Fixture "external-no-cds.gff" @(
        ("chr`tTEST`t" + "tRNA" + "`t3`t8`t.`t+`t0`tID=t1")
    )
    $overrideOut = Join-Path $temp "external-empty.gff"
    Assert-Success @("-i", $joined, "-a", $externalNoCds, "-l", "3", "-o", $overrideOut) | Out-Null
    $overrideRows = @(Read-Gff $overrideOut)
    Assert-Equal $overrideRows.Count 1 "Provided host disables embedded fallback"
    Assert-Equal "$($overrideRows[0].Start),$($overrideRows[0].End)" "1,20" "External empty mask wins"

    # Malformed inputs must fail cleanly.
    $badGff = Write-Fixture "bad.gff" @("chr`tTEST`tCDS`t1`t2")
    Assert-Failure @("-i", $fa100, "-a", $badGff) "GFF file format" "Malformed GFF"
    $badBed = Write-Fixture "bad.bed" @("chr`t1")
    Assert-Failure @("-i", $fa100, "-a", $badBed, "-f", "BED") "BED file format" "Malformed BED"
    $badFastq = Write-Fixture "bad.fastq" @("@chr", $sequence, "+", ("I"*19))
    Assert-Failure @("-i", $badFastq) "quality lengths differ" "Malformed FASTQ"
    $badSequence = Write-Fixture "bad.txt" @("not a sequence")
    Assert-Failure @("-i", $badSequence) "unsupported input sequence format" "Unknown genome format"
    $emptyHeader = Write-Fixture "empty-header.fa" @(">", "ACGT")
    Assert-Failure @("-i", $emptyHeader) "invalid sequence header" "Empty FASTA name"
    $noOrigin = Write-Fixture "no-origin.gbk" @("LOCUS       chr 20 bp DNA linear", "//")
    Assert-Failure @("-i", $noOrigin) "no sequence section" "GenBank without sequence"

    # Random linear/circular extraction cases against an independent model.
    for ($case = 0; $case -lt 20; $case ++) {
        $length = Get-Random -Minimum 20 -Maximum 121
        $circular = (Get-Random -Minimum 0 -Maximum 2) -eq 1
        $minlen = Get-Random -Minimum 1 -Maximum 31
        $featureCount = Get-Random -Minimum 0 -Maximum 8
        $features = @()
        $featureLines = @()
        for ($i = 0; $i -lt $featureCount; $i ++) {
            $start = Get-Random -Minimum 0 -Maximum $length
            $width = Get-Random -Minimum 1 -Maximum 21
            $end = [Math]::Min($length, $start+$width)
            if ($end -le $start) { continue }
            $features += [pscustomobject]@{ Start = $start; End = $end }
            $strand = @("+", "-")[$i % 2]
            $featureLines += "seq$case`tTEST`tCDS`t$($start+1)`t$end`t.`t$strand`t0`tID=f$i"
        }
        $topology = if ($circular) { "circular" } else { "linear" }
        $caseFa = Write-Fixture "fuzz-$case.fa" @(">seq$case $topology", ("A"*$length))
        $caseGff = Write-Fixture "fuzz-$case.gff" $featureLines
        $caseOut = Join-Path $temp "fuzz-$case-out.gff"
        Assert-Success @("-i", $caseFa, "-a", $caseGff, "-l", "$minlen", "-o", $caseOut) | Out-Null
        $actual = @(Read-Gff $caseOut | ForEach-Object {
            "$($_.Start-1),$($_.End)"
        }) -join ";"
        $expectedRuns = @(Get-ReferenceRuns $length $features $circular $minlen)
        $expected = @($expectedRuns | ForEach-Object {
            "$($_.Start),$($_.End)"
        }) -join ";"
        Assert-Equal $actual $expected "Random extraction case $case"
    }

    # Real project example: coordinate validity, no annotation overlap and
    # consistency between GFF, BED and FASTA output.
    $exampleFa = Join-Path $root "example.fa"
    $exampleGff = Join-Path $root "example.gff"
    if ((Test-Path -LiteralPath $exampleFa) -and (Test-Path -LiteralPath $exampleGff)) {
        $realGff = Join-Path $temp "real.gff"
        $realBed = Join-Path $temp "real.bed"
        $realFasta = Join-Path $temp "real.fa"
        $sw = [Diagnostics.Stopwatch]::StartNew()
        Assert-Success @("-i", $exampleFa, "-a", $exampleGff, "-o", $realGff) | Out-Null
        Assert-Success @("-i", $exampleFa, "-a", $exampleGff, "-o", $realBed) | Out-Null
        Assert-Success @("-i", $exampleFa, "-a", $exampleGff, "-o", $realFasta) | Out-Null
        $sw.Stop()

        $genomeRecords = @(Read-Fasta $exampleFa)
        Assert-Equal $genomeRecords.Count 1 "Example scaffold count"
        $genomeLength = $genomeRecords[0].Sequence.Length
        $realRows = @(Read-Gff $realGff)
        $realBeds = @(Read-Bed $realBed)
        $realFastas = @(Read-Fasta $realFasta)
        Assert-True ($realRows.Count -gt 0) "Real example produced no IGRs"
        Assert-Equal $realBeds.Count $realRows.Count "Real BED/GFF count"
        Assert-Equal $realFastas.Count $realRows.Count "Real FASTA/GFF count"

        $ids = New-Object "System.Collections.Generic.HashSet[string]"
        $annotationIntervals = @()
        foreach ($line in Get-Content -LiteralPath $exampleGff) {
            if ($line.StartsWith("#") -or [string]::IsNullOrWhiteSpace($line)) {
                continue
            }
            $fields = $line.Split("`t")
            if ($fields.Count -eq 9 -and $fields[2] -eq "CDS") {
                $start = [Math]::Max(0, [int]$fields[3]-1)
                $end = [Math]::Min($genomeLength, [int]$fields[4])
                if ($start -lt $end) {
                    $annotationIntervals += [pscustomobject]@{
                        Start = $start
                        End = $end
                    }
                }
            }
        }
        $merged = @(Merge-Intervals $annotationIntervals)

        for ($i = 0; $i -lt $realRows.Count; $i ++) {
            $row = $realRows[$i]
            $bed = $realBeds[$i]
            $fasta = $realFastas[$i]
            Assert-Equal $row.Source "INTERGENIC_v0.0.1" "Real GFF source $i"
            Assert-True ($row.ID -match "^IGS\d{6,}$") "Real GFF ID $i"
            Assert-True ($ids.Add($row.ID)) "Duplicate real GFF ID $($row.ID)"
            Assert-True ($row.Start -ge 1 -and $row.Start -le $genomeLength) "GFF start $i"
            Assert-True ($row.End -ge 1 -and $row.End -le $genomeLength) "GFF end $i"
            Assert-Equal $bed.Host $row.Host "Real BED host $i"
            Assert-Equal $bed.Start ($row.Start-1) "Real BED start $i"
            Assert-Equal $bed.End $row.End "Real BED end $i"
            Assert-Equal $bed.ID $row.ID "Real BED ID $i"
            Assert-Equal $bed.Adjacent $row.Adjacent "Real BED adjacency $i"
            Assert-Equal $fasta.Fields[0] "$($row.Host):$($row.Start)..$($row.End)(.)" "Real FASTA coordinate $i"
            Assert-Equal $fasta.Fields[1] "ID=$($row.ID)" "Real FASTA ID $i"
            Assert-Equal $fasta.Fields[2] "Adjacent=$($row.Adjacent)" "Real FASTA adjacency $i"

            $start = $row.Start-1
            $end = $row.End
            $regionLength = if ($start -le $end) {
                $end-$start
            } else {
                $end+$genomeLength-$start
            }
            Assert-True ($regionLength -ge 90) "Real IGR shorter than 90 bp at $i"
            Assert-Equal $fasta.Sequence.Length $regionLength "Real FASTA sequence length $i"

            if ($start -le $end) {
                Assert-True (-not (Test-Overlap $merged $start $end)) "Real IGR overlaps CDS $i"
            } else {
                Assert-True (-not (Test-Overlap $merged $start $genomeLength)) "Real wrapped IGR overlaps CDS $i"
                Assert-True (-not (Test-Overlap $merged 0 $end)) "Real wrapped IGR overlaps CDS $i"
            }
        }
        Write-Host "Real example: $($realRows.Count) IGRs in $([Math]::Round($sw.Elapsed.TotalSeconds, 2)) s."
    }

    Write-Host "All comprehensive tests passed."
}
finally {
    $resolvedTemp = [IO.Path]::GetFullPath($temp)
    $tempRoot = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
    if ($resolvedTemp.StartsWith($tempRoot, [StringComparison]::OrdinalIgnoreCase) -and
            (Split-Path -Leaf $resolvedTemp).StartsWith("intergenic-comprehensive-")) {
        Remove-Item -LiteralPath $resolvedTemp -Recurse -Force
    }
}
