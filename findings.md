# Findings — NHANES Nutrition + Cognition + Depression

**Last Updated**: 2026-07-24 (loop 2)
**Status**: Inner loops 1-2 complete — mediation + composite done

## Current Understanding

### H1: All 4 indices predict cognitive function — PARTIALLY SUPPORTED
Only CONUT and E-DII independently significant. PNI/GNRI show trends.

### H2: Depression mediates nutrition-cognition — DEPENDS ON INDEX
**Critical finding**: Mediation differs by index type:
- **E-DII**: a-path NS (p=0.127), mediated -6.8% (refuted)
- **CONUT**: a-path NS (p=0.307), mediated 4.4% (refuted)
- **PNI**: a-path p=0.0004, mediated **9.2%** ✓
- **GNRI**: a-path p=0.0006, mediated **9.1%** ✓

This makes biological sense: PNI and GNRI both incorporate albumin (systemic health marker), and better systemic nutrition → less depression → better cognition. DII (dietary inflammation) and CONUT (no albumin weighting) lack this mood pathway.

Consistent with published AHEI (10.5%) and Korean NHRS (17%) mediation studies.

### H3: Blood-based indices vs dietary recall — INCONCLUSIVE, BUT COMPOSITE WINS
- **Composite (All 4)**: β=0.105, p=0.001, ΔR²=0.0073 ← BEST
- **Composite (CONUT+E-DII)**: β=0.091, p=0.001, ΔR²=0.0073
- CONUT alone: β=0.052, p=0.001, ΔR²=0.0041

Combining dietary + blood-based indices captures more nutritional-cognitive variance than any single index. **Novel finding.**

### H4: Food security interaction — INCONCLUSIVE (low power N=284)

### Key Subgroup Patterns
- E-DII effect modified by SES and metabolic health
- Effect strongest in 60-69, higher income, non-diabetic
- Effect null in low income (p=0.95) and diabetic (p=0.94)

### Sensitivity
- Robust to excluding stroke, cancer, depression
- Attenuated at age ≥65

## Completed Arms
### GI Cancer + PNI/CONUT/GNRI (manuscript complete)
Dual-cohort NHANES III + Continuous NHANES (1988-2019), 313 GI cancer patients.
PNI HR 0.78, GNRI HR 0.67. Manuscript targeting *Clinical Nutrition*.

### Cognition + Nutrition (inner loops 1-2 complete)
**Main findings**:
1. CONUT strongest single index predictor of cognition
2. Depression mediates PNI/GNRI→cognition (9.1-9.2%), NOT DII/CONUT
3. Multi-index composite (All 4) outperforms every individual index (β=0.105)

## Open Questions
1. CKD stratification — need eGFR data merged
2. RCS for CONUT and PNI (currently only DII done)
3. Write up as manuscript targeting nutrition/cognition journal?
4. Is the SES modification real or recall artifact?

## Lessons and Constraints
- SEM mediation: lavaan bootstrap with 1000 reps takes ~2 min per model
- Composite scores: compute z-scores outside mutate to avoid select() scoping issues
