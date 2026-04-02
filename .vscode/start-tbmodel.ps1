$ErrorActionPreference = "Stop"

$projectRoot = "C:\Users\emmas\Documents\GITHUB\tb-community-risk"
$condaExe    = "C:\Users\emmas\anaconda3\Scripts\conda.exe"
$envName     = "tbmodel"

if (-not (Test-Path $projectRoot)) {
    throw "Project folder not found: $projectRoot"
}

if (-not (Test-Path $condaExe)) {
    throw "Conda not found at: $condaExe"
}

Set-Location $projectRoot

# Load conda into this PowerShell session
(& $condaExe "shell.powershell" "hook") | Out-String | Invoke-Expression

# Activate the project environment
conda activate $envName

Write-Host ""
Write-Host "Ready: $envName"
Write-Host "Folder: $projectRoot"
Write-Host ""