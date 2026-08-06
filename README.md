# GreedyComBat

GreedyComBat uses a two-step outcome-association-based procedure involving cross-validated univariate filtering followed by greedy subset selection to make use of the high-dimensional covariate sets available. Covariate selection may subsequently be used for harmonisation via ComBatFamily (https://github.com/andy1764/ComBatFamily/tree/main).

## Code

- `GreedyComBat.R`: GreedyComBat covariate-selection functions.

### Additional code
- `CovariateCuration.R`: construction of the ADNI covariate matrix used for the manuscript.
- `ADNIAnalysis.R`: analysis code used to generate Figures 1-6; these provide examples of how GreedyComBat functions may be called.
