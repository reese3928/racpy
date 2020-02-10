import itertools
import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
import mygene
from missingpy import KNNImputer
from . import exceptions as e


class RNAAgeCalc(object):
    '''Calculate RNA age.

    Attributes:
    - tissue: a string indicates which tissue the gene expression data
    is obtained from. Users are expected to provide one of the these tissues:

    adipose_tissue, adrenal_gland, blood, blood_vessel, brain, breast, colon,
    esophagus, heart, liver, lung, muscle, nerve, ovary, pancreas, pituitary,
    prostate, salivary_gland, skin, small_intestine, spleen, stomach, testis,
    thyroid, uterus, vagina.

    If the tissue argument is not provide or the provided tissue is not in
    this list, then the age predictor trained on all tissues will be used to
    calculate RNA age.

    - exprtype: either "count" or "FPKM". Default is "FPKM". For RPKM data,
    please use exprtype="FPKM".

    - idtype: a string which indicates the gene id type in "exprdata". It
    should be one of "symbol", "ensembl.gene", "entrezgene" or "refseq".
    Default is "symbol".

    - stype: a string which specifies which version of pre-trained
    calculators to be used. It should be either "all" or "Caucasian". "all"
    means samples from all races (American Indian/Alaska Native, Asian,
    Black/African American, and Caucasian) are used to obtain the pre-trained
    calculator. "Caucasian" means only the Caucasian samples are used to build
    up the pre-trained calculator. Default is "all".

    - signature: a string which indicate the age signature to use when
    calculating RNA age. This argument is not required.

    In the case that this argument is not provided, if `tissue` argument is
    also provided and the tissue is in the list above, the tissue specific age
    signature given by our DESeq2 analysis result on GTEx data will be used.
    Otherwise, the across tissue signature "GTExAge" will be used.

    In the case that this argument is provided, it should be one of the
    following signatures. A detailed description of the meaning of these
    signatures is given in the package tutorial.

    DESeq, Pearson, Dev, deMagalhaes, GenAge, GTExAge, Peters, all
    '''

    def __init__(self, tissue=None, exprtype="FPKM", idtype="symbol",
                 stype="all", signature=None):
        # check input
        self._check_exprtype(exprtype)
        self._exprtype = exprtype

        self._check_idtype(idtype)
        self._idtype = idtype

        self._check_stype(stype)
        self._stype = stype
        # end check input

        # process tissue and signature argument
        if tissue is None:
            print("tissue is not provided, using the RNA age predictor "
                  "trained by all tissues automatically.")
            tissue = "all_tissues"
            if signature is None:
                print("signature is not provided, using GTExAge signature "
                      "automatically.")
                signature = "GTExAge"
            else:
                assert signature in self._sig_set(), \
                    "signature should be one of DESeq, Pearson, Dev, " \
                    "deMagalhaes, GenAge, GTExAge, Peters, all."
                if signature == "DESeq":
                    print("DESeq signature is currently not available for all "
                          "tissues, using Pearson signature automatically.")
                    signature = "Pearson"
        else:
            assert isinstance(tissue, str), "tissue should be a string."
            if tissue.lower() in self._tissue_set():
                tissue = tissue.lower()
                if signature is None:
                    print("signature is not provided, using DESeq signature "
                          "automatically.")
                    signature = "DESeq"
                else:
                    assert signature in self._sig_set(), \
                        "signature should be one of DESeq, Pearson, Dev, " \
                        "deMagalhaes, GenAge, GTExAge, Peters, all."
            else:
                print("the provided tissue was not found, using the RNA age "
                      "predictor trained by all tissues automatically.")
                tissue = "all_tissues"
                if signature is None:
                    print("signature is not provided, using GTExAge "
                          "signature automatically.")
                    signature = "GTExAge"
                else:
                    assert signature in self._sig_set(), \
                        "signature should be one of DESeq, Pearson, Dev, " \
                        "deMagalhaes, GenAge, GTExAge, Peters, all."
                    if signature == "DESeq":
                        print("DESeq signature is currently not available for "
                              "all tissues, using Pearson signature "
                              "automatically.")
                        signature = "Pearson"
        self._tissue = tissue
        self._signature = signature
        # end: process tissue and signature argument

    @property
    def exprtype(self):
        return self._exprtype

    @exprtype.setter
    def exprtype(self, exprtype):
        self._check_exprtype(exprtype)
        self._exprtype = exprtype

    @staticmethod
    def _check_exprtype(exprtype):
        assert isinstance(exprtype, str), "exprtype should be a string."
        assert exprtype in ["count", "FPKM"], "exprtype should be count " \
                                              "or FPKM."

    @property
    def idtype(self):
        return self._idtype

    @idtype.setter
    def idtype(self, idtype):
        self._check_idtype(idtype)
        self._idtype = idtype

    @staticmethod
    def _check_idtype(idtype):
        assert isinstance(idtype, str), "idtype should be a string."
        assert idtype in ["symbol", "ensembl.gene", "entrezgene", "refseq"], \
            "idtype should be one of symbol, ensembl.gene, entrezgene" \
            " or refseq."

    @property
    def stype(self):
        return self._stype

    @stype.setter
    def stype(self, stype):
        self._check_stype(stype)
        self._stype = stype

    @staticmethod
    def _check_stype(stype):
        assert isinstance(stype, str), "stype should be a string."
        stype = stype.lower()
        assert stype in ["all", "caucasian"], \
            "stype should be either all or Caucasian."

    @property
    def tissue(self):
        return self._tissue

    @tissue.setter
    def tissue(self, tissue):
        raise e.SetterError("Cannot set tissue, please construct a new "
                            "RNAAgeCalc object and use the tissue argument"
                            " to set tissue. ")

    @property
    def signature(self):
        return self._signature

    @signature.setter
    def signature(self, signature):
        raise e.SetterError("Cannot set signature, please construct a new "
                            "RNAAgeCalc object and use the signature argument"
                            " to set signature. ")

    def _processlength(self, rawcount, genelength):
        if genelength is None:
            location = os.path.dirname(os.path.realpath(__file__))
            lengthDB = pd.read_csv(os.path.join(location, "internal_data",
                                                "lengthDB.csv"), index_col=0)
            if self.idtype != "ensembl.gene":
                # convert the gene id to ensembl to match the gene id in gene
                # length database
                mg = mygene.MyGeneInfo()
                genes = list(rawcount.index)
                print("converting the gene id to ensembl to match the gene id "
                      "in gene length database. This may take some time.")
                temp = mg.querymany(genes, scopes='symbol',
                                    fields='ensembl.gene',
                                    species='human', returnall=True,
                                    as_dataframe=True)["out"]
                temp1 = temp.loc[~temp["ensembl.gene"].isna(), "ensembl.gene"]
                temp2 = temp.loc[~temp["ensembl"].isna(), "ensembl"]
                index_expand = []
                value_expand = []
                for index, value in temp2.iteritems():
                    index_expand.append([index]*len(value))
                    value_expand.append(list(map(lambda x: x["gene"], value)))
                temp3 = pd.Series(
                    list(itertools.chain.from_iterable(value_expand)),
                    index=list(itertools.chain.from_iterable(index_expand)))
                temp1 = temp1.append(temp3)

                temp1 = temp1[temp1.isin(lengthDB.index)]
                temp1 = temp1[~temp1.index.duplicated(keep="first")]
                temp1 = temp1.drop_duplicates(keep=False)

                genename = temp1[rawcount.index]
                genename[genename.isnull()] = "unknown"
                genelength = lengthDB.loc[genename]
            else:
                genelength = lengthDB.loc[rawcount.index]
        else:
            bool1 = isinstance(genelength, pd.DataFrame) or \
                    isinstance(genelength, pd.Series) or \
                    isinstance(genelength, list) or \
                    isinstance(genelength, np.ndarray)
            assert bool1, "genelength should be a pandas DataFrame, Series, "\
                          "numpy array or list."

            if isinstance(genelength, pd.DataFrame):
                assert genelength.shape[0] == rawcount.shape[0], \
                    "The number of rows in genelength should be the same as " \
                    "that of rawcount."
                if genelength.shape[1] > 1:
                    print("genelength contains more than 1 column, only the"
                          " 1st column will be used.")
                    genelength = pd.DataFrame(genelength.iloc[:, 0])
                assert genelength.applymap(np.isreal).all().bool()
            elif isinstance(genelength, pd.Series):
                assert len(genelength) == rawcount.shape[0], \
                    "The length of genelength should be the same as " \
                    "number of rows in rawcount."
                genelength = pd.DataFrame(genelength)
                assert genelength.applymap(np.isreal).all().bool()
            elif isinstance(genelength, list):
                assert all(isinstance(item, float) or isinstance(item, int)
                           for item in genelength)
                assert len(genelength) == rawcount.shape[0], \
                    "The length of genelength should be the same as " \
                    "number of rows in rawcount."
                genelength = pd.DataFrame(genelength)
            else:
                if genelength.shape[0] > 1:
                    print("numpy array has more than 1 dimension, only the "
                          "1st dimension will be used.")
                    genelength = genelength[0]
                assert np.isreal(genelength)
                genelength = pd.DataFrame(genelength).transpose()
        return genelength

    def _count2FPKM(self, rawcount, genelength):
        """Convert gene expression data from raw count to FPKM.

        This function converts gene expression data from raw count to FPKM.

        :param rawcount: a pandas DataFrame which contains gene expression
        count data.

        :param genelength: a pandas Series, DataFrame, numpy array, or list
        which contains gene length in bp. The size of genelength should be
        equal to the number of rows in exprdata. This argument is optional.

        If using exprtype="FPKM", genelength argument is ignored. If using
        exprtype="count", the raw count will be converted to FPKM. If
        genelength is provided, the function will convert raw count to FPKM
        according to the user-supplied gene length. Otherwise, gene length
        is obtained from the internal database.

        :return: a pandas DataFrame contains FPKM.
        """

        genelength = self._processlength(rawcount, genelength)
        assert (genelength.iloc[:, 0].dropna() > 0).all(), \
            "genelength cannot contain negative value(s) or 0."

        num_nas = genelength.iloc[:, 0].isna().sum()
        if num_nas > 0:
            print("Can't find gene length for {:.2f}% genes when converting "
                  "raw count to FPKM.".format(num_nas/genelength.shape[0]*100))

        countsum = rawcount.apply(sum)
        bg = pd.DataFrame(columns=rawcount.columns, index=rawcount.index)
        bg.loc[:, :] = countsum.values

        wid = pd.DataFrame(columns=rawcount.columns, index=rawcount.index)
        wid.loc[:, :] = genelength.values

        return rawcount / (wid/1000) / (bg/1e6)

    def predict_age(self, exprdata, genelength=None, chronage=None):
        """Calculate RNA age.

        This function calculates RNA age based on pre-trained predictors.

        :param exprdata: a pandas DataFrame which contains gene expression data
        with each row represents a gene and each column represents a sample.
        Use the argument "exprtype" to specify raw count or FPKM. The index of
        "exprdata" should be gene ids and columns names of "exprdata" should
        be sample ids.

        :param genelength: a pandas Series, DataFrame, numpy array, or list
        which contains gene length in bp. The size of genelength should be
        equal to the number of rows in exprdata. This argument is optional.
        If using exprtype="FPKM", genelength argument is ignored. If using
        exprtype="count", the raw count will be converted to FPKM. If
        genelength is provided, the function will convert raw count to FPKM
        according to the user-supplied gene length. Otherwise, gene length
        is obtained from the internal database.

        :param chronage: a pandas DataFrame which contains the chronological
        age of each sample. This argument is optional. If provided, it should
        be a DataFrame with 1st column sample id and 2nd column chronological
        age. The sample order in chronage doesn't have to be in the same order
        as in exprdata. However, the samples in chronage and exprdata should
        be the same. If some samples' chronological age are not available,
        users are expected to set the chronological age in chronage to NaN.
        If chronage contains more than 2 columns, only the first 2 columns
        will be considered. If this argument is not provided, the age
        acceleration residual will not be calculated. See package tutorial
        for the definition of age acceleration residual.

        :return: a pandas DataFrame contains RNA age.

        """

        # check input:
        assert isinstance(exprdata, pd.DataFrame), \
            "exprdata should be a pandas DataFrame."
        assert exprdata.applymap(np.isreal).all().all(),\
            "Only numeric values are allowed in the exprdata DataFrame."
        assert list(exprdata.index) != list(range(exprdata.shape[0])), \
            "The index of exprdata should be gene ids."
        assert list(exprdata.columns) != list(range(exprdata.shape[1])), \
            "The column names of exprdata should be sample ids."
        assert ~np.any(exprdata.index.duplicated()), \
            "Duplicated gene names found in exprdata."
        assert (exprdata >= 0).all().all(), \
            "Gene expression data cannot contain negative value(s)."

        if chronage is not None:
            assert isinstance(chronage, pd.DataFrame), \
                "chronage should be a pandas DataFrame."
            if(chronage.shape[1] > 2):
                print("More than 2 columns are provided in chronage. "
                      "Only the first 2 columns will be used.")
            # assert ~chronage.applymap(np.isreal).all()[0], \
            #    "The 1st column in chronage should be sample ids."
            assert chronage.applymap(np.isreal).all()[1], \
                "The 2nd column in chronage should be chronological age."
            assert ~any(chronage.iloc[:, 0].duplicated()), \
                "chronage contains duplicated sample ids."
            assert set(chronage.iloc[:, 0].astype(str)) == \
                set(exprdata.columns), \
                "Samples in chronage and exprdata should be the same."
            assert ~np.any(chronage.iloc[:, 1] < 0), \
                "Chronological age contains negative value(s)."

        if self._exprtype == "count":
            exprdata = self._count2FPKM(exprdata, genelength)

        if self.idtype != "symbol":
            mg = mygene.MyGeneInfo()
            genes = list(exprdata.index)
            temp = mg.querymany(genes, scopes=self.idtype, fields='symbol',
                                species='human', returnall=True,
                                as_dataframe=True)["out"]
            temp = temp.loc[~temp["symbol"].isna(), "symbol"]
            temp = temp[~temp.index.duplicated(keep="first")]
            temp = temp.drop_duplicates(keep=False)
            genesymbol = temp[exprdata.index]
            genesymbol[genesymbol.isna()] = "unknown"
            exprdata.index = genesymbol

        location = os.path.dirname(os.path.realpath(__file__))
        if self.stype == "all":
            tempPath = os.path.join(location, "internal_data", "all",
                                    "coef_{}_{}.csv".format(self._tissue,
                                                            self._signature))
        else:
            tempPath = os.path.join(location, "internal_data", "Caucasian",
                                    "coef_{}_{}.csv".format(self._tissue,
                                                            self._signature))

        sig_internal = pd.read_csv(tempPath, index_col=0)
        genes_required = sig_internal.index[1:]
        sig_in_expr = genes_required.isin(exprdata.index)
        # full NA row
        if np.sum(~sig_in_expr) != 0:
            print("{:.2f}% genes in the gene signature are not included in "
                  "the supplied gene expression."
                  .format(np.sum(~sig_in_expr)/len(genes_required)*100))

            # impute the gene expression in the log scale
            tempmat = pd.DataFrame(columns=exprdata.columns,
                                   index=genes_required[~sig_in_expr])

            exprdata_withNA = pd.concat([exprdata, tempmat], axis=0)
            exprdata_log = np.log2(exprdata_withNA.apply(pd.to_numeric) + 1)
            ind1 = exprdata_log.isna().all(axis=1)
            ind2 = ~exprdata_log.isna().any(axis=1)
            exprdata_log.loc[(ind1 | ind2), :] = \
                exprdata_log.loc[(ind1 | ind2), :].fillna(exprdata_log.mean())
        else:
            exprdata_log = np.log2(exprdata.apply(pd.to_numeric) + 1)

        # check partial NA
        if ~exprdata_log.notnull().all().all():
            # impute the gene expression in the log scale
            imputer = KNNImputer(n_neighbors=min(10, exprdata_log.shape[1]),
                                 row_max_missing=1, col_max_missing=1)
            X_imputed = imputer.fit_transform(exprdata_log.transpose())
            exprdata_log_impute = pd.DataFrame(X_imputed).transpose()
            exprdata_log_impute.index = exprdata_log.index
            exprdata_sub = exprdata_log_impute.loc[genes_required, :]
        else:
            exprdata_sub = exprdata_log.loc[genes_required, :]

        RNAAge = exprdata_sub.apply(
            lambda x: np.sum(x.multiply(sig_internal.iloc[1:, 0])) +
            sig_internal.iloc[0, 0])
        res = pd.DataFrame(index=exprdata.columns)
        res["RNAAge"] = list(RNAAge)

        if chronage is not None:
            chronage.index = chronage.iloc[:, 0]
            res["ChronAge"] = chronage.loc[res.index].iloc[:, 1]
            # if sample size is too small, age acceleration residual
            # cannot be calculated
            if res.dropna().shape[0] > 30:
                Y = res["RNAAge"]
                X = res["ChronAge"]
                X = sm.add_constant(X)
                model = sm.OLS(Y, X).fit()
                res["AgeAccelResid"] = model.resid

        return res

    @staticmethod
    def _tissue_set():
        return set(["adipose_tissue", "adrenal_gland", "blood", "blood_vessel",
                    "brain", "breast", "colon", "esophagus", "heart", "liver",
                    "lung", "muscle", "nerve", "ovary", "pancreas",
                    "pituitary", "prostate", "salivary_gland", "skin",
                    "small_intestine", "spleen", "stomach", "testis",
                    "thyroid", "uterus", "vagina", "all_tissues"])

    @staticmethod
    def _sig_set():
        return set(["DESeq", "Pearson", "Dev", "deMagalhaes", "GenAge",
                    "GTExAge", "Peters", "all"])

    def __repr__(self):
        return """RNAAgeCalc(tissue=%r, exprtype=%r, idtype=%r, stype=%r,
signature=%r)""" % (self.tissue, self.exprtype, self.idtype,
                    self.stype, self.signature)

    def __str__(self):
        return """RNAAgeCalc
tissue: {}
exprtype: {}
idtype: {}
stype: {}
signature: {}""".format(self.tissue, self.exprtype, self.idtype,
                        self.stype, self.signature)
