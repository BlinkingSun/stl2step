<#
package-engine.ps1 — per-release ENGINE ASSETS the SolidOut in-app updater consumes (Windows half).

  .github\scripts\package-engine.ps1 <tag> [-BuildDir <dir>] [-Out <dir>] [-Vcpkg <dir>]

Produces, for the tag's engine:

(1) stl2step-engine-<tag>-windows-x64.zip containing:
    - stl2step.exe
    - *.dll closure in the SAME directory (Windows DLL search order; no subdirs).
    The closure is derived exactly as scripts\update-engine.ps1 does (dumpbin /dependents
    BFS, resolve against build dir then vcpkg bin then VC++ redist). No hand list.
    No Authenticode signing (zip is hash-verified by the updater; unlike macOS there is
    no notarization step).

(2) engine-manifest.json — windows-x64 name/sha256/size/minWindows filled; macos-arm64
    preserved if an existing manifest is already in -Out (from the macOS packaging step).

(3) SHA256SUMS-engine.txt — windows line plus any macOS line from the merged manifest.

Does not upload. Does not mutate the build tree (copies the binary out).
#>
[CmdletBinding()]
param(
    [Parameter(Position = 0, Mandatory = $true)][string]$Tag,
    [string]$BuildDir = '',
    [string]$Out = '',
    [string]$Vcpkg = ''
)

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'

$E_USAGE = 2

function Say([string]$m) { [Console]::Out.WriteLine("package-engine: $m") }
function Die([string]$m) { [Console]::Error.WriteLine("package-engine: $m"); exit 1 }

function Invoke-Native {
    param([Parameter(Mandatory = $true)][scriptblock]$Block, [switch]$PassOutput)
    $saved = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $out = & $Block 2>&1
        $code = $LASTEXITCODE
    } finally { $ErrorActionPreference = $saved }
    if ($PassOutput) { return [pscustomobject]@{ Code = $code; Out = $out } }
    return $code
}

$Root = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path))
$BuildDir = if ($BuildDir) { $BuildDir } else { Join-Path $Root 'build-wip' }
$Out = if ($Out) { $Out } else { Join-Path $Root 'dist\engine' }

$BuiltDir = $null
$VcpkgBin = $null

$SystemDlls = @(
    'kernel32.dll', 'kernelbase.dll', 'ntdll.dll', 'user32.dll', 'gdi32.dll', 'gdi32full.dll',
    'advapi32.dll', 'ole32.dll', 'oleaut32.dll', 'shell32.dll', 'shlwapi.dll', 'ws2_32.dll',
    'crypt32.dll', 'bcrypt.dll', 'ncrypt.dll', 'secur32.dll', 'comdlg32.dll', 'comctl32.dll',
    'dbghelp.dll', 'version.dll', 'winmm.dll', 'wintrust.dll', 'setupapi.dll', 'iphlpapi.dll',
    'userenv.dll', 'psapi.dll', 'rpcrt4.dll', 'msvcrt.dll', 'ucrtbase.dll', 'combase.dll',
    'd3d11.dll', 'dxgi.dll', 'd3dcompiler_47.dll', 'opengl32.dll', 'glu32.dll', 'gdiplus.dll',
    'imm32.dll', 'uxtheme.dll', 'dwmapi.dll', 'propsys.dll', 'powrprof.dll', 'cfgmgr32.dll',
    'normaliz.dll', 'mswsock.dll', 'dnsapi.dll', 'winhttp.dll', 'wininet.dll', 'urlmon.dll',
    'oleacc.dll', 'msimg32.dll', 'avrt.dll', 'mfplat.dll', 'winspool.drv'
)
function Test-SystemDll([string]$name) {
    $n = $name.ToLowerInvariant()
    if (($n -like 'api-ms-*') -or ($n -like 'ext-ms-*') -or ($SystemDlls -contains $n)) { return $true }
    foreach ($d in @($BuiltDir, $VcpkgBin)) {
        if ($d -and (Test-Path (Join-Path $d $name))) { return $false }
    }
    return (Test-Path (Join-Path "$env:SystemRoot\system32" $name))
}
function Test-VcRuntimeDll([string]$name) {
    return $name -match '(?i)^(msvcp140|vcruntime140|concrt140|vccorlib140)'
}

function Find-VsRoot {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (Test-Path $vswhere) {
        $p = (Invoke-Native -PassOutput { & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath }).Out
        if ($p) { return ("$($p | Select-Object -First 1)").Trim() }
    }
    foreach ($c in @(
            'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools',
            'C:\Program Files\Microsoft Visual Studio\2022\BuildTools',
            'C:\Program Files\Microsoft Visual Studio\2022\Community',
            'C:\Program Files\Microsoft Visual Studio\2022\Professional')) {
        if (Test-Path $c) { return $c }
    }
    return $null
}

function Get-OcctVersion([string]$vcpkgRoot) {
    $info = Join-Path $vcpkgRoot 'installed\vcpkg\info'
    if (-not (Test-Path $info)) { return $null }
    $f = Get-ChildItem $info -Filter 'opencascade_*' -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $f) { return $null }
    if ($f.Name -match '^opencascade_([0-9]+(\.[0-9]+)*)') { return $Matches[1] }
    return $null
}
function Select-VcpkgRoot {
    $candidates = @()
    if ($Vcpkg) { $candidates += $Vcpkg }
    if ($env:VCPKG_INSTALLATION_ROOT) { $candidates += $env:VCPKG_INSTALLATION_ROOT }
    if ($env:VCPKG_ROOT) { $candidates += $env:VCPKG_ROOT }
    if ($candidates.Count -eq 0) { $candidates += @('D:\vcpkg', 'C:\vcpkg') }
    $best = $null; $bestVer = $null
    foreach ($c in $candidates) {
        $toolchain = Join-Path $c 'scripts\buildsystems\vcpkg.cmake'
        $info = Join-Path $c 'installed\vcpkg\info'
        if (-not (Test-Path $toolchain) -and -not (Test-Path $info)) { continue }
        $v = Get-OcctVersion $c
        if (-not $v) { continue }
        if (-not $bestVer -or ([version]$v -gt [version]$bestVer)) { $best = $c; $bestVer = $v }
    }
    return $best
}

$VsRoot = Find-VsRoot
if (-not $VsRoot) { Die 'no Visual Studio C++ toolset found' }
$DumpBin = Get-ChildItem (Join-Path $VsRoot 'VC\Tools\MSVC\*\bin\Hostx64\x64\dumpbin.exe') -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending | Select-Object -First 1
if (-not $DumpBin) { Die 'no dumpbin.exe in the VS install' }

$BuiltExe = $null
foreach ($c in @(
        (Join-Path $BuildDir 'Release\stl2step.exe'),
        (Join-Path $BuildDir 'stl2step.exe'),
        (Join-Path $BuildDir 'bin\Release\stl2step.exe'))) {
    if (Test-Path $c) { $BuiltExe = $c; break }
}
if (-not $BuiltExe) { Die "no stl2step.exe under $BuildDir (pass -BuildDir at a CMake engine build tree)" }
Say "build binary $BuiltExe"

$VcpkgRoot = Select-VcpkgRoot
if (-not $VcpkgRoot) { Die 'no usable vcpkg root with an opencascade port (set VCPKG_ROOT or pass -Vcpkg)' }
$VcpkgBin = Join-Path $VcpkgRoot 'installed\x64-windows\bin'
Say "vcpkg root $VcpkgRoot"

$BuiltDir = Split-Path -Parent $BuiltExe
$RedistDirs = Get-ChildItem (Join-Path $VsRoot 'VC\Redist\MSVC\*\x64\Microsoft.VC*.CRT') -Directory -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending

function Get-Dependents([string]$path) {
    $out = (Invoke-Native -PassOutput { & $DumpBin.FullName /nologo /dependents $path }).Out
    $names = @()
    foreach ($line in $out) {
        if ("$line" -match '^\s{4}(\S+\.[Dd][Ll][Ll])\s*$') { $names += $Matches[1] }
    }
    return $names
}
function Resolve-Dll([string]$name) {
    foreach ($d in @($BuiltDir, $VcpkgBin)) {
        $p = Join-Path $d $name
        if (Test-Path $p) { return $p }
    }
    if (Test-VcRuntimeDll $name) {
        foreach ($d in $RedistDirs) {
            $p = Join-Path $d.FullName $name
            if (Test-Path $p) { return $p }
        }
    }
    return $null
}

$closure = [ordered]@{}
$queue = New-Object System.Collections.Generic.Queue[string]
$queue.Enqueue($BuiltExe)
$seen = New-Object System.Collections.Generic.HashSet[string]
[void]$seen.Add($BuiltExe.ToLowerInvariant())
$unresolved = @()
while ($queue.Count -gt 0) {
    $cur = $queue.Dequeue()
    foreach ($name in (Get-Dependents $cur)) {
        $key = $name.ToLowerInvariant()
        if ($closure.Contains($key)) { continue }
        if (-not (Test-VcRuntimeDll $name) -and (Test-SystemDll $name)) { continue }
        $src = Resolve-Dll $name
        if (-not $src) { $unresolved += "$name (needed by $(Split-Path -Leaf $cur))"; continue }
        $closure[$key] = $src
        if (-not $seen.Contains($src.ToLowerInvariant())) {
            [void]$seen.Add($src.ToLowerInvariant())
            $queue.Enqueue($src)
        }
    }
}
if ($unresolved.Count -gt 0) {
    foreach ($u in $unresolved) { Say "UNRESOLVED $u" }
    Die "$($unresolved.Count) closure member(s) could not be resolved"
}
$dllCount = $closure.Count
Say "closure: $dllCount DLLs (same derivation as update-engine.ps1)"

$ZipName = "stl2step-engine-$Tag-windows-x64.zip"
New-Item -ItemType Directory -Path $Out -Force | Out-Null
$Out = (Resolve-Path $Out).Path
$Work = Join-Path ([System.IO.Path]::GetTempPath()) ("package-engine-" + [guid]::NewGuid().ToString('N').Substring(0, 8))
$Stage = Join-Path $Work 'payload'
New-Item -ItemType Directory -Path $Stage -Force | Out-Null
try {
    Copy-Item -LiteralPath $BuiltExe -Destination (Join-Path $Stage 'stl2step.exe') -Force
    foreach ($key in $closure.Keys) {
        $src = $closure[$key]
        Copy-Item -LiteralPath $src -Destination (Join-Path $Stage (Split-Path -Leaf $src)) -Force
    }

    $savedPath = $env:PATH
    $savedEap = $ErrorActionPreference
    $env:PATH = "$env:SystemRoot\system32;$env:SystemRoot"
    $ErrorActionPreference = 'Continue'
    try {
        $r = Invoke-Native -PassOutput { & (Join-Path $Stage 'stl2step.exe') --version }
        $verOut = $r.Out; $verCode = $r.Code
    } finally { $env:PATH = $savedPath; $ErrorActionPreference = $savedEap }
    $versionLine = (($verOut | ForEach-Object { "$_" }) -join "`n").Trim().Split("`n")[0].Trim()
    if ($verCode -ne 0 -or -not $versionLine) { Die 'staged stl2step --version produced no output' }
    if ($versionLine -notmatch '^stl2step\s+(.+)$') { Die "could not parse engine version from '$versionLine'" }
    $engineVersion = $Matches[1]
    Say "version $versionLine"

    $ZipPath = Join-Path $Out $ZipName
    if (Test-Path $ZipPath) { Remove-Item $ZipPath -Force }
    Compress-Archive -Path (Join-Path $Stage '*') -DestinationPath $ZipPath -CompressionLevel Optimal
    $size = (Get-Item $ZipPath).Length
    $sha256 = (Get-FileHash -LiteralPath $ZipPath -Algorithm SHA256).Hash.ToLowerInvariant()
    Say "zip $ZipName size=$size sha256=$($sha256.Substring(0, 12))…"

    $manifestPath = Join-Path $Out 'engine-manifest.json'
    $sumsPath = Join-Path $Out 'SHA256SUMS-engine.txt'
    $doc = $null
    if (Test-Path $manifestPath) {
        $doc = Get-Content -Raw $manifestPath | ConvertFrom-Json
    } else {
        $doc = [pscustomobject]@{
            schema         = 1
            tag            = $Tag
            engineVersion  = $engineVersion
            assets         = [pscustomobject]@{
                'macos-arm64' = [pscustomobject]@{
                    name   = "stl2step-engine-$Tag-macos-arm64.tar.gz"
                    sha256 = $null
                    size   = $null
                    note   = 'filled by the macOS packaging step'
                }
                'windows-x64' = [pscustomobject]@{}
            }
        }
    }
    if (-not $doc.assets.'windows-x64') {
        $doc.assets | Add-Member -NotePropertyName 'windows-x64' -NotePropertyValue ([pscustomobject]@{}) -Force
    }
    $doc.tag = $Tag
    $doc.engineVersion = $engineVersion
    $doc.assets.'windows-x64' = [pscustomobject]@{
        name       = $ZipName
        sha256     = $sha256
        size       = [int64]$size
        minWindows = '10'
    }
    $doc | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $manifestPath -Encoding UTF8

    $sumLines = @("$sha256  $ZipName")
    if ($doc.assets.'macos-arm64'.sha256 -and $doc.assets.'macos-arm64'.name) {
        $sumLines = @("$($doc.assets.'macos-arm64'.sha256)  $($doc.assets.'macos-arm64'.name)") + $sumLines
    }
    [System.IO.File]::WriteAllText($sumsPath, ($sumLines -join "`n") + "`n", (New-Object System.Text.UTF8Encoding($false)))
    Say "wrote $manifestPath"
    Say "wrote $sumsPath"
    Say "ok tag=$Tag version=$versionLine zip=$ZipName dlls=$dllCount (no Authenticode)"
} finally {
    Remove-Item $Work -Recurse -Force -ErrorAction SilentlyContinue
}
