# LanguageEndangerment
Custom R code for ordinal probit regression with correction for autocorrelation

The repository is for Appendix A in "Global predictors of language endangerment and the future of the world’s linguistic diversity"

data.Rdata and Codes.R include all the data and R codes used in all the analyses reported in the paper.

Codes.R has seperate blocks for each analysis. 

The analyses mainly use functions in R package "ordinalNet" for ordinal probit regression.

Custom R function 'autoord' and its related functions modify functions in "ordinalNet" to correct for autocorrelation in ordinal probit regression.

Details of these functions are in Codes.R.


Note:
Since publishing this code, the authors have been made aware of an error in Dow (2008) (see this preprint by Alexander Koplenig for details https://osf.io/preprints/osf/xbu2k_v1). Reanalysis using an alternative method (semiparametric eigenvector spatial fltering) suggests that this issue does not have any substantive impact on the conclusions of the study.
