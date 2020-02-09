
import pandas as pd
import os

dataPath = "/export/home/xurren/internal_data"
lengthDB = pd.read_csv(os.path.join(dataPath, "lengthDB.csv"), index_col=0)
fpkm = pd.read_csv(os.path.join(dataPath, "fpkm.csv"), index_col=0)
rawcount = pd.read_csv(os.path.join(dataPath, "rawcount.csv"), index_col=0)

tissues = ["adipose_tissue", "adrenal_gland", "blood", "blood_vessel", "brain", "breast",
           "colon", "esophagus", "heart", "liver", "lung", "muscle", "nerve", "ovary",
           "pancreas", "pituitary", "prostate", "salivary_gland", "skin", "small_intestine",
           "spleen", "stomach", "testis", "thyroid", "uterus", "vagina", "all_tissues"]
signatures = ["DESeq", "Pearson", "Dev", "deMagalhaes", "GenAge", "Kuan_Ren", "Peters", "all"]

genelist_dict = {}
for t in tissues[:-1]:
    genelist_element = {}
    for sig in signatures:
        df = pd.read_csv(os.path.join(dataPath, "coef_" + t + "_" + sig + ".csv"))
        coef_series = pd.Series(list(df.coef), index = df.geneid)
        genelist_element[sig] = coef_series
    genelist_dict[t] = genelist_element


genelist_element = {}
for sig in signatures[1:]:
    df = pd.read_csv(os.path.join(dataPath, "coef_all_tissues_" + sig + ".csv"))
    coef_series = pd.Series(list(df.coef), index = df.geneid)
    genelist_element[sig] = coef_series

genelist_dict["all_tissues"] = genelist_element


