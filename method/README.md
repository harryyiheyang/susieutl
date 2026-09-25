# Method workflow

This directory is the entry point for continuing the analysis in a local
chat. From the repository root, run
`Rscript method/seur_sim1000/seur_sim1000_make.R` to generate outputs under
local `tmp/` paths. The workflow order is:

1. `seur_sim1000_make.R`
2. `seur_sim1000_utl.R` (UTL6)
3. `seur_sim1000_adapter.R`
4. `seur_sim1000_single.R`
5. `seur_sim1000_mesusie.R` (MESuSiE)
6. `seur_sim1000_multisusie.py` (MultiSuSiE)
7. `seur_sim1000_susiex.ps1` (SuSiEx)
8. Run the frozen 20-case MOM fit, then collect its outputs and run the LFSR
   step using `mom_sa_random20.R` and `mom_sa_random20_lfsr.R`.

All data and generated results—including RDS, CSV, TSV, pickle, binary, and
log files—belong only in local `tmp/` and are not committed to Git.

External inputs and software are not vendored here:

- XMAP example input: `example_data.RData`.
- Installed R packages: SuSiEUTL, digest, susieR, and MESuSiE 1.0.
- Python with numpy and the official MultiSuSiE implementation.
- PowerShell, WSL with Ubuntu, SuSiEx, and plink.

Previously checked upstream commits were MultiSuSiE `10351eb…` and SuSiEx
`1db2f583…`; their source trees are not vendored.
