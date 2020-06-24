=====
Usage
=====


Introduction
#############

It has been shown that both DNA methylation and RNA transcription are linked 
to chronological age and age related diseases. Several estimators have 
been developed to predict human aging from DNA level and RNA level. Most of the 
human transcriptional age predictor are based on microarray data and limited 
to only a few tissues. To date, transcriptional studies on aging using 
RNASeq data from different human tissues is limited. The aim of this package 
is to provide a tool for across-tissue and tissue-specific transcriptional age 
calculation based on Genotype-Tissue Expression (GTEx) RNASeq data 
[1]_.


Description of RNASeq age calculator
#######################################

We utilized the GTEx data to construct our across-tissue and tissue-specific 
transcriptional age calculator. GTEx is a public available genetic database 
for studying tissue specific gene expression and regulation. GTEx V6 release 
contains gene expression data at gene, exon, and transcript level of 9,662 
samples from 30 different tissues. To avoid the influence of tumor on gene 
expression, the 102 tumor samples from GTEx V6 release are dropped and 
the remaining 9,560 samples were used in the subsequent analysis. To 
facilitate integrated analysis and direct comparison of multiple datasets, 
we utilized recount2 [2]_ version of GTEx data, where 
all samples were processed with the same analytical pipeline. FPKM values 
were calculated for each individual sample using `getRPKM` function in 
Bioconductor package `recount <http://bioconductor.org/packages/release/bioc/html/recount.html>`_.

For the tissue-specific RNASeq age calculator, elastic net 
[3]_ algorithm was used to train the predictors for each 
individual tissue. Chronological age was response variable whereas logarithm 
transformed FPKM of genes were predictors. The across-tissue calculator was 
constructed by first performing differential expression analysis on the 
RNASeq count data for each individual tissue. To identify genes consistently 
differentially expressed across tissues, we adapted the binomial test 
discussed in de Magalhaes et al. [4]_ to find the genes with the 
largest number of age-related signals. A detailed explanation can be found 
in our paper.

The package is implemented as follows. For each tissue, signature and 
sample type (see below for the descriptions), we pre-trained the calculator 
using elastic net based on the GTEx samples. We saved the pre-trained model 
coefficients as internal data in the package. The package takes gene 
expression data as input and then match the input genes to the genes in the 
internal data. This matching process is automatic so that the users just 
need to provide gene expression data without having to pull out the 
internal coefficients.


Usage of RNASeq age calculator
#################################

To use racpy in a project::

    from racpy import RNAAgeCalc

Then construct an `RNAAgeCalc` object (here we use "brain" as an example)::

    rac_obj = RNAAgeCalc(tissue = "brain")

Next we use an example of FPKM data to make prediction::

    from racpy import fpkm
    res = rac_obj.predict_age(fpkm)
    print(res)

Here we explain the options in `RNAAgeCalc`.

tissue
************

`tissue` is a string indicates which tissue the gene expression data is
obtained from. Users are expected to provide one of the following tissues.
If the tissue argument is not provided or the provided tissue is not in this 
list, the age predictor trained on all tissues will be used to calculate 
RNA age.

* adipose_tissue    
* adrenal_gland    
* blood    
* blood_vessel    
* brain    
* breast    
* colon    
* esophagus    
* heart    
* liver    
* lung    
* muscle    
* nerve    
* ovary    
* pancreas    
* pituitary     
* prostate    
* salivary_gland    
* skin    
* small_intestine     
* spleen      
* stomach        
* testis       
* thyroid       
* uterus       
* vagina       


exprtype
************

`exprtype` is either "count" or "FPKM". If `exprtype` is count, the 
expression data will be converted to FPKM by the internal function and 
the calculator will be applied on FPKM data. When calculating FPKM, by default 
gene length is obtained from the package's internal database. The internal 
gene length information was obtained from recount2. However, users are able 
to provide their own gene length information by using `genelength` argument 
in `predict_age` function (see below).


idtype
**********
`idtype` is a string which indicates the gene id type in `exprdata`. Default 
is "symbol". The following id types are supported.  

* symbol    
* ensembl.gene    
* entrezgene   
* refseq   


stype
***********
`stype` is a string which specifies which version of pre-trained calculators 
to be used. Two versions are provided. If `stype="all"`, the calculator 
trained on samples from all races (American Indian/Alaska Native, Asian, 
Black/African American, and Caucasian) will be used. If `stype="Caucasian"`, 
the calculator trained on Caucasian samples only will be used. We found that 
RNA Age signatures could be different in different races (see our paper for 
details). Thus we provide both the universal calculator and race specific 
calculator. The race specific calculator for American Indian/Alaska Native, 
Asian, or Black/African American are not provided due to the small sample 
size in GTEx data.


signature
************
`signature` is a string which indicate the age signature to use when 
calculating RNA age. This argument is not required. 

In the case that this argument is not provided, if `tissue` argument is also
provided and the tissue is in the list above, the tissue specific age
signature given by our DESeq2 analysis result on GTEx data will be used. 
Otherwise, the across tissue signature "GTExAge" will be used. 

In the case that this argument is provided, it should be one of the following 
signatures. 

* DESeq2. DESeq2 signature was obtained by performing differential expression 
  analysis on each tissue and select the top differential expressed genes.  
* Pearson. Pearson signature represents the genes highly correlated with 
  chronological age by Pearson correlation.    
* Dev. Dev signature contains genes with large variation in expression across 
  samples. We adapted the gene selection strategy discussed in [5]_, which is 
  a gene must have at least a :math:`t_1`-fold difference in expression between 
  any two samples in the training set and at least one sample have expression 
  level > :math:`t_2` FPKM to be included in the prediction models. :math:`t_1` 
  and :math:`t_2` (typically 5 or 10) are thresholds to control the degree of 
  deviance of the genes. We used :math:`t_1 = t_2 = 10` for most tissues. 
  For some tissues with large sample size, in order to maximize the prediction 
  accuracy while maintaining low computation cost, we increased :math:`t_1` and 
  :math:`t_2` such that the number of genes retained in the model is between 
  2,000 and 7,000.    
* deMagalhaes. deMagalhaes signature contains the 73 age-related genes by [4]_.    
* GenAge. GenAge signature contains the 307 age-related genes in the Ageing 
  Gene Database [6]_.    
* GTExAge. GTExAge signature represents the genes consistently differentially 
  expressed across tissues discussed in our paper.   
* Peters. Peters signature contains the 1,497 genes differentially expressed 
  with age discussed in [7]_.    
* all. "all" represents all the genes used when constructing the RNAAge 
  calculator. 

If the genes in `exprdata` do not cover all the genes in the signature, 
imputation will be made automatically by the `KNNImputer` function in 
`missingpy <https://pypi.org/project/missingpy/>`__.

Below are the options for the `predict_age` function.

exprdata
**********

`exprdata` a pandas DataFrame which contains gene expression data
with each row represents a gene and each column represents a sample. Users are 
expected to use the argument "exprtype" to specify raw count or FPKM. The index 
of "exprdata" should be gene ids and columns names of "exprdata" should be sample ids.
Here is an example of FPKM expression data::
    
    from racpy import fpkm
    fpkm.head()


genelength
*************

`genelength` is a pandas Series, DataFrame, numpy array, or list which contains gene 
length in bp. The size of `genelength` should be equal to the number of rows in `exprdata`. 
This argument is optional. When using `exprtype = "FPKM"`, `genelength` argument is ignored. 
When using `exprtype = "count"`, the raw count will be converted to FPKM. If `genelength` 
is provided, the function will convert raw count to FPKM based on the user-supplied gene 
length. Otherwise, gene length is obtained from the internal database.


chronage
***********
`chronage` is a pandas DataFrame which contains the chronological age of each
sample. This argument is optional. 

If provided, it should be a DataFrame with 1st column sample id and 2nd column 
chronological age. The sample order in `chronage` doesn't have to be in the 
same order as in `exprdata`. However, the samples in `chronage` and `exprdata` 
should be the same. If some samples' chronological age are not available, 
users are expected to set the chronological age in `chronage` to NaN. If 
`chronage` contains more than 2 columns, only the first 2 columns will be 
considered. If more than 30 samples' chronological age are available, age 
acceleration residual will be calculated. Age acceleration residual is 
defined as the residual of linear regression with RNASeq age as dependent 
variable and chronological age as independent variable.

If this argument is not provided, the age acceleration residual will not be
calculated. 


Example
#############

This example is just for illustration purpose. It does not represent any real data::

    import pandas as pd
    from racpy import RNAAgeCalc
    from racpy import fpkm
    # construct a gene expression data
    fpkm_large = pd.concat([fpkm, fpkm+1, fpkm+2, fpkm+3], axis = 1)
    fpkm_large = pd.concat([fpkm_large, fpkm_large, fpkm_large, fpkm_large], axis = 1)
    fpkm_large.columns = ["sample"+str(item+1) for item in range(32)]
    # construct the samples' chronological age
    chronage2 = pd.DataFrame()
    chronage2["sampleid"] = fpkm_large.columns
    chronage2["age"] = range(31, 63)

    rac_obj2 = RNAAgeCalc(tissue = "brain")
    res2 = rac_obj2.predict_age(exprdata=fpkm_large, chronage=chronage2)
    print(res2)



Visualization
################

We suggest visualizing the results by plotting RNAAge vs chronological age. 
This can be done by calling `makeplot` function and passing in the DataFrame 
returned by `predict_age` function::

    import matplotlib.pyplot as plt
    from racpy import makeplot
    makeplot(res2)
    plt.show()


References
#############

.. [1] Lonsdale, John, et al. "The genotype-tissue expression (GTEx) project." Nature genetics 45.6 (2013): 580.
.. [2] Collado-Torres, Leonardo, et al. "Reproducible RNA-seq analysis using recount2." Nature biotechnology 35.4 (2017): 319-321.
.. [3] Zou, Hui, and Trevor Hastie. "Regularization and variable selection via the elastic net." Journal of the royal statistical society: series B (statistical methodology) 67.2 (2005): 301-320.
.. [4] De Magalhães, João Pedro, João Curado, and George M. Church. "Meta-analysis of age-related gene expression profiles identifies common signatures of aging." Bioinformatics 25.7 (2009): 875-881.
.. [5] Fleischer, Jason G., et al. "Predicting age from the transcriptome of human dermal fibroblasts." Genome biology 19.1 (2018): 221.
.. [6] de Magalhaes, Joao Pedro, and Olivier Toussaint. "GenAge: a genomic and proteomic network map of human ageing." FEBS letters 571.1-3 (2004): 243-247.
.. [7] Peters, Marjolein J., et al. "The transcriptional landscape of age in human peripheral blood." Nature communications 6.1 (2015): 1-14.