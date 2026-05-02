$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $PSScriptRoot
$quarto = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe"
$outputDir = Join-Path $projectRoot "Output\manuscript"

New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$env:R_LIBS_USER = "C:\Users\User\AppData\Local\R\win-library\4.3"
$env:LOCALAPPDATA = Join-Path $projectRoot ".localappdata"
New-Item -ItemType Directory -Force -Path $env:LOCALAPPDATA | Out-Null

function Render-Output {
  param(
    [string]$InputFile,
    [string]$OutputFile,
    [string]$Format
  )

  $expectedPath = Join-Path $outputDir $OutputFile
  if (Test-Path $expectedPath) {
    Remove-Item -LiteralPath $expectedPath -Force
  }

  & $quarto render $InputFile `
    --to $Format `
    --output $OutputFile `
    --output-dir $outputDir `
    --execute-dir $projectRoot `
    --no-clean

  $nestedPath = Join-Path $outputDir (Join-Path "Manuscript" $OutputFile)
  $rootOutputPath = Join-Path (Join-Path $projectRoot "Output") $OutputFile

  if (Test-Path $nestedPath) {
    Move-Item -LiteralPath $nestedPath -Destination $expectedPath -Force
  }

  if (Test-Path $rootOutputPath) {
    Move-Item -LiteralPath $rootOutputPath -Destination $expectedPath -Force
  }

  if ($Format -eq "apaquarto-pdf") {
    $texFile = [System.IO.Path]::ChangeExtension($OutputFile, ".tex")
    $expectedTexPath = Join-Path $outputDir $texFile
    if (Test-Path $expectedTexPath) {
      Remove-Item -LiteralPath $expectedTexPath -Force
    }
    $texCandidates = @(
      (Join-Path $projectRoot $texFile),
      (Join-Path (Join-Path $projectRoot "Output") $texFile),
      (Join-Path (Join-Path $projectRoot "Manuscript") $texFile),
      (Join-Path $outputDir (Join-Path "Manuscript" $texFile))
    )

    foreach ($candidate in $texCandidates) {
      if (Test-Path $candidate) {
        Move-Item -LiteralPath $candidate -Destination $expectedTexPath -Force
      }
    }
  }
}

Render-Output "Manuscript\vmac_selective_manuscript.qmd" "vmac_selective_manuscript.docx" "apaquarto-docx"
Render-Output "Manuscript\vmac_selective_manuscript_anom.qmd" "vmac_selective_manuscript_anom.docx" "apaquarto-docx"

Render-Output "Manuscript\vmac_selective_manuscript.qmd" "vmac_selective_manuscript.pdf" "apaquarto-pdf"
Render-Output "Manuscript\vmac_selective_manuscript_anom.qmd" "vmac_selective_manuscript_anom.pdf" "apaquarto-pdf"

$cleanupPaths = @(
  (Join-Path $outputDir "Manuscript"),
  (Join-Path $outputDir "Output"),
  (Join-Path $outputDir "pre_registation")
)

foreach ($path in $cleanupPaths) {
  if (Test-Path $path) {
    Remove-Item -LiteralPath $path -Recurse -Force
  }
}

if (Test-Path $env:LOCALAPPDATA) {
  Remove-Item -LiteralPath $env:LOCALAPPDATA -Recurse -Force
}
