# A simulation study: Detection of missing data mechanisms with graphical models and their imputation

Knowing the missing mechanisms in a given dataset is central to handling missing data, minimizing biases that can arise due to missingness and also adapt the imputation method. 
To detect these mechanisms a graphical representation of these is tested in this simulation study. Graphical models provide a good tool for comprehending, encoding
and communicating assumptions about the missingness process. To improve the traditional graphical models which assume the data to be fully observed an alternative
algorithm, named mvpc-algorithm is also tested, based on the theory of Mohan et. al. (Mohan, Karthika ; Pearl, Judea: Graphical models for processing missing data. In:
Journal of the American Statistical Association 116 (2021), Nr. 534, S. 1023–1037) and the code base of Tu (Tu, Ruibo ; Zhang, Cheng ; Ackermann, Paul ; Mohan, Karthika ; Kjellström,
Hedvig ; Zhang, Kun: Causal discovery in the presence of missing data. In: The22nd International Conference on Artificial Intelligence and Statistics PMLR, 2019, S.1762–1770) (see also code folder "mvpc").
Based on the detection of the missingness mechanisms, the imputation should be adapted for MNAR variables. Whereas for MCAR/MAR variables multiple imputation with
Amelia was used. For the MNAR variables a weighted knn imputation is carried out. The results in the simulation study show that, although MNAR imputations are biased
when using traditional multiple imputation methods which include model assumptions that are violated through this missing mechanism, the imputations
are quite similar to those of the MCAR/MAR variables when using Amelia. On the other side the imputations suffered from changing the method to the weighted knn imputation.

The main part of the R-code for this simulation study can be found in the "src" folder. The simulated data sets, containing not only the data itself but also results from the models can be reproduced by this code but also downloaded from the "data" folder.



