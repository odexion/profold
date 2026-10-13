# The profold installer for Windows.
#
#   powershell -ExecutionPolicy ByPass -c "irm https://github.com/odexion/profold/releases/latest/download/install.ps1 | iex"
#
# Installs for the current user, with no administrator rights:
#   %LOCALAPPDATA%\Programs\profold\profold.exe, a Start menu shortcut,
#   and `profold` on the user PATH.
#
# Settings, all optional:
#   $env:PROFOLD_VERSION = "0.2.0"     a particular release rather than the latest
#   $env:PROFOLD_INSTALL_DIR = "D:\x"  where profold.exe goes
#   $env:PROFOLD_NO_MODIFY_PATH = "1"  leave PATH alone
#   $env:PROFOLD_DOWNLOAD_BASE = URL   fetch the archives from a mirror

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
$Repo = "odexion/profold"

# Windows PowerShell 5.1 may not offer TLS 1.2 by default.
try {
    [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
} catch {}

# Written as character codes: `irm | iex` may not read the script as UTF-8.
$Heavy = [string][char]0x2501
$Tip = [string][char]0x257A
$Tick = [string][char]0x221A
$Dot = [string][char]0x00B7

function Step([string]$Text) {
    Write-Host "==> " -ForegroundColor Cyan -NoNewline
    Write-Host $Text -ForegroundColor White
}
function Ok([string]$Text) {
    Write-Host "  $Tick " -ForegroundColor Green -NoNewline
    Write-Host $Text
}
function Warn([string]$Text) {
    Write-Host "  ! " -ForegroundColor Yellow -NoNewline
    Write-Host $Text
}
function Fail([string]$Text) {
    Write-Host "error: " -ForegroundColor Red -NoNewline
    Write-Host $Text
    throw "profold was not installed"
}

function Get-Target {
    $arch = $null
    try {
        $arch = [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture.ToString()
    } catch {
        # Older .NET: an x64 PowerShell on ARM reports AMD64 in the first.
        $arch = if ($env:PROCESSOR_ARCHITEW6432) { $env:PROCESSOR_ARCHITEW6432 } else { $env:PROCESSOR_ARCHITECTURE }
    }
    switch -Regex ($arch) {
        "^(Arm64|ARM64)$" { return "aarch64-pc-windows-msvc" }
        "^(X64|AMD64)$" { return "x86_64-pc-windows-msvc" }
        default { Fail "unsupported Windows architecture: $arch" }
    }
}

function Format-Size([double]$Bytes) {
    $units = "B", "KiB", "MiB", "GiB"
    $i = 0
    while ($Bytes -ge 1024 -and $i -lt 3) { $Bytes /= 1024; $i++ }
    if ($i -eq 0) { return "{0:N0} {1}" -f $Bytes, $units[$i] }
    return "{0:N1} {1}" -f $Bytes, $units[$i]
}

function Show-Bar([long]$Have, [long]$Total, [double]$Seconds) {
    $width = 30
    $speed = if ($Seconds -gt 0) { "  " + (Format-Size ($Have / $Seconds)) + "/s" } else { "" }
    Write-Host "`r  Downloading  " -NoNewline
    if ($Total -gt 0) {
        $f = [Math]::Min(1.0, $Have / $Total)
        $n = [int][Math]::Floor($f * $width)
        Write-Host ($Heavy * $n) -ForegroundColor Magenta -NoNewline
        if ($n -lt $width) { Write-Host ($Tip + ($Heavy * ($width - $n - 1))) -ForegroundColor DarkGray -NoNewline }
        Write-Host ("  {0,3:N0}%  {1} / {2}{3}    " -f ($f * 100), (Format-Size $Have), (Format-Size $Total), $speed) -NoNewline
    } else {
        Write-Host ("{0}{1}    " -f (Format-Size $Have), $speed) -NoNewline
    }
}

function Get-File([string]$Url, [string]$Path, [switch]$Quiet) {
    $request = [System.Net.HttpWebRequest]::Create($Url)
    $request.UserAgent = "profold-installer"
    $response = $request.GetResponse()
    $total = $response.ContentLength
    $in = $response.GetResponseStream()
    $out = [System.IO.File]::Create($Path)
    $clock = [System.Diagnostics.Stopwatch]::StartNew()
    $last = -1000
    $have = 0
    try {
        $buffer = New-Object byte[] 262144
        while (($n = $in.Read($buffer, 0, $buffer.Length)) -gt 0) {
            $out.Write($buffer, 0, $n)
            $have += $n
            if (-not $Quiet -and $clock.ElapsedMilliseconds - $last -ge 100) {
                $last = $clock.ElapsedMilliseconds
                Show-Bar $have $total $clock.Elapsed.TotalSeconds
            }
        }
    } finally {
        $out.Close()
        $in.Close()
        $response.Close()
    }
    if (-not $Quiet) {
        Show-Bar $have $have $clock.Elapsed.TotalSeconds
        Write-Host ""
    }
}

function Add-ToPath([string]$Dir) {
    $current = [Environment]::GetEnvironmentVariable("Path", "User")
    $parts = if ($current) { $current.Split(";") | Where-Object { $_ } } else { @() }
    if ($parts -contains $Dir) { return }
    if ($env:PROFOLD_NO_MODIFY_PATH) {
        Warn "$Dir is not on your PATH; add it to run ``profold`` from a terminal"
        return
    }
    [Environment]::SetEnvironmentVariable("Path", (@($parts) + $Dir) -join ";", "User")
    $env:Path = "$env:Path;$Dir"
    Ok "added $Dir to your PATH (open a new terminal to use ``profold``)"
}

function Install-Profold {
    Write-Host ""
    Write-Host "  profold" -ForegroundColor Magenta -NoNewline
    Write-Host " installer" -ForegroundColor White -NoNewline
    Write-Host "  a molecular viewer $Dot github.com/$Repo" -ForegroundColor DarkGray
    Write-Host ""

    if (-not (Get-Command tar.exe -ErrorAction SilentlyContinue)) {
        Fail "this installer needs tar.exe, part of Windows 10 1803 and later"
    }

    $target = Get-Target
    $version = $env:PROFOLD_VERSION
    if ($version) {
        $version = $version.TrimStart("v")
        $base = "https://github.com/$Repo/releases/download/v$version"
    } else {
        $base = "https://github.com/$Repo/releases/latest/download"
    }
    if ($env:PROFOLD_DOWNLOAD_BASE) { $base = $env:PROFOLD_DOWNLOAD_BASE }
    $archive = "profold-$target.tar.gz"
    $label = if ($version) { "profold $version" } else { "profold" }
    Step "Installing $label for $target"

    $tmp = Join-Path ([System.IO.Path]::GetTempPath()) ("profold-" + [Guid]::NewGuid())
    New-Item -ItemType Directory -Path $tmp | Out-Null
    try {
        $file = Join-Path $tmp $archive
        try {
            Get-File "$base/$archive" $file
            Get-File "$base/$archive.sha256" "$file.sha256" -Quiet
        } catch {
            Fail "could not download $base/$archive`n  Is there a release for ${target}? See https://github.com/$Repo/releases"
        }
        $want = ((Get-Content "$file.sha256" -Raw).Trim() -split "\s+")[0].ToLower()
        $got = (Get-FileHash $file -Algorithm SHA256).Hash.ToLower()
        if ($want -ne $got) { Fail "the download is damaged (checksum mismatch); try again" }
        Ok "checksum verified"

        $unpacked = Join-Path $tmp "unpacked"
        New-Item -ItemType Directory -Path $unpacked | Out-Null
        & tar.exe -xzf $file -C $unpacked
        if ($LASTEXITCODE -ne 0) { Fail "could not unpack $archive" }
        $new = Join-Path $unpacked "profold.exe"
        if (-not (Test-Path $new)) { Fail "the download has no profold.exe in it" }

        $dir = if ($env:PROFOLD_INSTALL_DIR) { $env:PROFOLD_INSTALL_DIR } else { Join-Path $env:LOCALAPPDATA "Programs\profold" }
        New-Item -ItemType Directory -Force -Path $dir | Out-Null
        $exe = Join-Path $dir "profold.exe"
        # A running profold.exe cannot be overwritten, but it can be moved.
        $old = Join-Path $dir ".profold.exe.old"
        Remove-Item $old -Force -ErrorAction SilentlyContinue
        if (Test-Path $exe) { Move-Item $exe $old -Force }
        Move-Item $new $exe
        Remove-Item $old -Force -ErrorAction SilentlyContinue
        Ok "installed $exe"

        $programs = [Environment]::GetFolderPath("Programs")
        $shell = New-Object -ComObject WScript.Shell
        $link = $shell.CreateShortcut((Join-Path $programs "profold.lnk"))
        $link.TargetPath = $exe
        $link.WorkingDirectory = $dir
        $link.Description = "View proteins, molecules, crystals and simulations"
        $link.Save()
        Ok "added profold to the Start menu"

        Add-ToPath $dir
    } finally {
        Remove-Item $tmp -Recurse -Force -ErrorAction SilentlyContinue
    }

    Write-Host ""
    Write-Host "  profold is installed!" -ForegroundColor Green -NoNewline
    Write-Host " To start, open it from the Start menu, or run profold."
    Write-Host "  It keeps itself up to date: new versions arrive with a `"Restart to update`" button." -ForegroundColor DarkGray
    Write-Host ""
}

Install-Profold
