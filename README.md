# missing-value-imputation-methods

this is an accompanying simulation for another [repo](https://github.com/bshashikadze/maternaldiabetes-offspring-liver-omics-paper)


Proteomics dataset is generated with Data-Independent-Acquisition (DIA)

[three imputation methods were compared]()

1. imputation with random forest
2. imputation with KNN
3. left censored imputation (similar as in Perseus dafault and DEP R package)

Additionally, all three imputation were applied row-wise (proteins are rows) and column-wise (proteins are columns)


This script is NOT intended for general comparison of imputation methods (for this see e.g. [1-3], also [R package for proteomics data imputation](https://cran.rstudio.com/web/packages/imp4p/index.html)), but to estimate which imputation will work better in a specific proteomics dataset. For this, missing values were introduced randomly in a matrix (without missingness), finally the Normalized Root Mean Square Error (NRMSE) was calculated to compare imputed value with a real value. NRMSE close to 0 indicates better model. 

$\sqrt{\frac{\sum_{i=1}^{n} (Si-Oi)^2}{n}}\$

1.	Jin, L., Y. Bi, C. Hu, J. Qu, S. Shen, X. Wang, and Y. Tian, A comparative study of evaluating missing value imputation methods in label-free proteomics. Sci Rep, 2021. 11(1): p. 1760.
2.	Liao, S.G., Y. Lin, D.D. Kang, D. Chandra, J. Bon, N. Kaminski, F.C. Sciurba, and G.C. Tseng, Missing value imputation in high-dimensional phenomic data: imputable or not, and how? BMC Bioinformatics, 2014. 15: p. 346.
3.	Stekhoven, D.J. and P. Bühlmann, MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics, 2011. 28(1): p. 112-118.

