# Research Log — NHANES Nutrition + Cognition + Depression

## 2026-07-24 — Bootstrap → Inner Loop 1 Complete

### What Happened
- Literature survey: 21 papers, 4 critical gaps identified
- Data pipeline: bypassed nhanesA timeouts with local XPT download
- Fixed CFQ variable names (CFDCST1/2/3 and CFDCSR, not CFD_WL_IMM/DEL)
- Fixed DPQ columns (010-090 not 010-019)
- Fixed BIOPRO column typo (LBDSCRSI not LBDSCRPSI)
- Fixed ALQ columns
- Fixed RCS (ns() instead of rcs() for svyglm compatibility)
- Fixed svyglm predict return value handling

### Key Results
- **H1**: CONUT > E-DII > PNI ≈ GNRI (CONUT strongest, β=0.052, p=0.001)
- **H2**: Depression mediation refuted for DII (-6.8%, NS)
- **H3**: Mixed — CONUT outperforms, PNI/GNRI don't
- **H4**: Inconclusive (power limitation with N=284 food insecure)

### Files Changed
- Fixed: 00_config.R, 01_load_and_derive.R, 02_analysis.R, 03_figures.R
- Added: 01a_load_local_xpt.R (local loader)
- Added: literature/ (21 papers), to_human/ (progress report)

### Next
- Inner loop 2: Test depression mediation for PNI/CONUT/GNRI
- Inner loop 3: Multi-index composite score
- Inner loop 4: RCS for CONUT
