[CmdletBinding()]
param(
  [Parameter(Mandatory = $true)][string]$RepRoot,
  [Parameter(Mandatory = $true)][string]$EnvRoot
)
$ErrorActionPreference = "Stop"
$root = [IO.Path]::GetFullPath($RepRoot)
$envRoot = [IO.Path]::GetFullPath($EnvRoot)
$adapter = Join-Path $root "adapter"
$od = Join-Path $root "methods\susiex"
if (-not (Test-Path -LiteralPath $adapter -PathType Container)) { throw "adapter output is missing" }
if (Test-Path -LiteralPath $od) { throw "refusing to overwrite SuSiEx output: $od" }
New-Item -ItemType Directory -Path $od | Out-Null

function Convert-ToWsl([string]$Path) {
  $p = [IO.Path]::GetFullPath($Path)
  if ($p -notmatch '^([A-Za-z]):\\(.*)$') { throw "cannot convert path to WSL: $p" }
  "/mnt/$($Matches[1].ToLower())/$($Matches[2].Replace('\','/'))"
}

$producer = @(Import-Csv -Delimiter "`t" -LiteralPath (Join-Path $root "status\producer_status.tsv"))
if ($producer.Count -ne 1 -or $producer[0].status -ne "success") { throw "producer status is invalid" }
$repName = ([int]$producer[0].rep).ToString("000")
$name = "$($producer[0].case)_rep$repName"
$aw = Convert-ToWsl $adapter
$ow = Convert-ToWsl $od
$ew = Convert-ToWsl $envRoot
$bin = "$ew/vendor/SuSiEx/bin_static/SuSiEx"
$plink = "$ew/vendor/SuSiEx/utilities/plink"
$cmd = "$bin --sst_file=$aw/susiex/sst_EUR.txt,$aw/susiex/sst_EAS.txt,$aw/susiex/sst_AFR.txt --n_gwas=1000000,150000,150000 --ld_file=$aw/susiex/LD_EUR,$aw/susiex/LD_EAS,$aw/susiex/LD_AFR --out_dir=$ow --out_name=$name --chr=22 --bp=50000001,50000500 --chr_col=1,1,1 --snp_col=2,2,2 --bp_col=3,3,3 --a1_col=4,4,4 --a2_col=5,5,5 --eff_col=6,6,6 --se_col=7,7,7 --pval_col=9,9,9 --plink=$plink --n_sig=5"
$snp = Join-Path $od "$name.snp"
$summary = Join-Path $od "$name.summary"
$cs = Join-Path $od "$name.cs"
$component = Join-Path $od "component_pip.tsv"
function Write-NativeStatus([string]$Status, [string]$Terminal,
                            [string]$MetricStatus, [int]$Ncs) {
  $summaryHash = if (Test-Path -LiteralPath $summary -PathType Leaf) { (Get-FileHash -Algorithm SHA256 -LiteralPath $summary).Hash.ToLower() } else { "" }
  $csHash = if (Test-Path -LiteralPath $cs -PathType Leaf) { (Get-FileHash -Algorithm SHA256 -LiteralPath $cs).Hash.ToLower() } else { "" }
  $snpHash = if (Test-Path -LiteralPath $snp -PathType Leaf) { (Get-FileHash -Algorithm SHA256 -LiteralPath $snp).Hash.ToLower() } else { "" }
  $componentPath = if (Test-Path -LiteralPath $component -PathType Leaf) { [IO.Path]::GetFullPath($component) } else { "" }
  $componentHash = if ($componentPath) { (Get-FileHash -Algorithm SHA256 -LiteralPath $component).Hash.ToLower() } else { "" }
  @([pscustomobject]@{method="susiex";status=$Status;native_terminal=$Terminal;
    elapsed_sec=$sw.Elapsed.TotalSeconds;input_md5=$producer[0].input_md5;
    version="binary";commit="1db2f5838a55f9af8595626a0a01f88940bb502f";
    V_policy="native_auto_default";metric_status=$MetricStatus;n_cs=$Ncs;
    summary_path=[IO.Path]::GetFullPath($summary);summary_sha256=$summaryHash;
    cs_path=[IO.Path]::GetFullPath($cs);cs_sha256=$csHash;
    snp_path=[IO.Path]::GetFullPath($snp);snp_sha256=$snpHash;
    component_pip_path=$componentPath;component_pip_sha256=$componentHash;call=$cmd}) |
    Export-Csv -Delimiter "`t" -NoTypeInformation -LiteralPath (Join-Path $od "status.tsv")
}
$sw = [Diagnostics.Stopwatch]::StartNew()
& wsl.exe -d Ubuntu -- bash -lc $cmd
$ec = $LASTEXITCODE
$sw.Stop()
if ($ec -ne 0) {
  Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
  throw "official SuSiEx failed with exit code $ec"
}
foreach ($f in @($snp, $summary, $cs)) {
  if (-not (Test-Path -LiteralPath $f -PathType Leaf) -or (Get-Item -LiteralPath $f).Length -eq 0) {
    Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
    throw "expected SuSiEx native output is missing: $f"
  }
}
$summaryLines = @(Get-Content -LiteralPath $summary | ForEach-Object { $_.Trim() } |
  Where-Object { $_ -and $_ -notmatch '^\s*#' })
$csLines = @(Get-Content -LiteralPath $cs | ForEach-Object { $_.Trim() } |
  Where-Object { $_ -and $_ -notmatch '^\s*#' })
$summaryMarkers = @($summaryLines | Where-Object { $_ -eq "FAIL" -or $_ -eq "NULL" })
$csMarkers = @($csLines | Where-Object { $_ -eq "FAIL" -or $_ -eq "NULL" })
$hasFail = "FAIL" -in $summaryMarkers -or "FAIL" -in $csMarkers
$summaryNull = "NULL" -in $summaryMarkers
$csNull = "NULL" -in $csMarkers
if ($hasFail -or $summaryNull -ne $csNull -or
    ($summaryNull -and (@($summaryLines | Where-Object { $_ -ne "NULL" }).Count -gt 0 -or
                        @($csLines | Where-Object { $_ -ne "NULL" }).Count -gt 0)) -or
    (-not $summaryNull -and ($summaryLines.Count -eq 0 -or $csLines.Count -eq 0))) {
  Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
  throw "SuSiEx returned an invalid native terminal marker combination"
}
$rows = @(Import-Csv -Delimiter "`t" -LiteralPath $snp)
if ($rows.Count -ne 500) {
  Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
  throw "SuSiEx .snp must have exactly 500 rows"
}
for ($i = 0; $i -lt 500; $i++) {
  $v = "v{0:D4}" -f ($i + 1)
  if ([string]$rows[$i].SNP -cne $v) {
    Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
    throw "SuSiEx SNP order is not canonical"
  }
}
$pc = @($rows[0].PSObject.Properties.Name | Where-Object { $_ -match '^PIP\(CS[0-9]+\)$' })
if ($summaryNull -and $pc.Count -ne 0) {
  Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
  throw "SuSiEx NULL terminal must not contain component PIP columns"
}
if (-not $summaryNull -and $pc.Count -lt 1) {
  Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
  throw "SuSiEx non-NULL terminal returned no component PIP columns"
}
if ($summaryNull) {
  Write-NativeStatus "success" "success_no_cs_NULL" "unavailable_no_cs_NULL" 0
  exit 0
}
$pipOut = @("variant`tvariant_row`tcomponent`tpip`trank")
foreach ($c in $pc) {
  $vals = New-Object double[] 500
  for ($i = 0; $i -lt 500; $i++) {
    [double]$x = 0
    if (-not [double]::TryParse([string]$rows[$i].PSObject.Properties[$c].Value,
      [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture,
      [ref]$x) -or [double]::IsNaN($x) -or [double]::IsInfinity($x) -or
      $x -lt 0 -or $x -gt 1) {
      Write-NativeStatus "failed" "failure" "unavailable_native_failure" 0
      throw "SuSiEx component PIP is invalid"
    }
    $vals[$i] = $x
  }
  for ($i = 0; $i -lt 500; $i++) {
    $rank = 1 + @($vals | Where-Object { $_ -gt $vals[$i] }).Count
    $variant = "v{0:D4}" -f ($i + 1)
    $pipOut += "$variant`t$($i+1)`t$c`t$($vals[$i].ToString('G17',[Globalization.CultureInfo]::InvariantCulture))`t$rank"
  }
}
$pipOut | Set-Content -LiteralPath (Join-Path $od "component_pip.tsv")
Write-NativeStatus "success" "success_with_cs" "unsupported_no_posterior_moments" $pc.Count
