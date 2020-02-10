import unittest
from racpy import RNAAgeCalc
from racpy import rawcount
from racpy import fpkm
from racpy.exceptions import SetterError
import pandas as pd
import numpy as np


class Test(unittest.TestCase):

    def test_assert_raises1(self):
        fpkm1 = fpkm
        age_calc = RNAAgeCalc(tissue="brain")

        with self.assertRaises(AssertionError):
            age_calc.predict_age(exprdata=1)

        with self.assertRaises(AssertionError):
            fpkm1.iloc[:, 0] = "unknown"
            age_calc.predict_age(exprdata=fpkm1)

        with self.assertRaises(AssertionError):
            # fpkm1 = fpkm
            fpkm1.index = list(range(fpkm1.shape[0]))
            age_calc.predict_age(exprdata=fpkm1)

        with self.assertRaises(AssertionError):
            # fpkm1 = fpkm
            fpkm1.columns = list(range(fpkm1.shape[1]))
            age_calc.predict_age(exprdata=fpkm1)

        with self.assertRaises(AssertionError):
            # fpkm1 = fpkm
            templist = list(fpkm1.index)
            templist[1] = templist[0]
            fpkm1.index = templist
            age_calc.predict_age(exprdata=fpkm1)

        with self.assertRaises(AssertionError):
            # fpkm1 = fpkm
            fpkm1.iloc[0, 0] = -1
            age_calc.predict_age(exprdata=fpkm1)

    def test_assert_raises2(self):
        with self.assertRaises(AssertionError):
            RNAAgeCalc(exprtype=1)

        with self.assertRaises(AssertionError):
            RNAAgeCalc(exprtype="unknown")

    def test_assert_raises3(self):
        with self.assertRaises(AssertionError):
            RNAAgeCalc(idtype=1)

        with self.assertRaises(AssertionError):
            RNAAgeCalc(idtype="unknown")

    def test_assert_raises4(self):
        chronage1 = pd.DataFrame()
        chronage1["sampleid"] = ["SRR2166176", "SRR2167642"]
        chronage1["chronage"] = [50, 60]
        age_calc = RNAAgeCalc(tissue="brain")

        with self.assertRaises(AssertionError):
            age_calc.predict_age(exprdata=fpkm, chronage=1)

        # with self.assertRaises(AssertionError):
        #    chronage1["sampleid"] = [50, 60]
        #    rnaagepy.predict_age(exprdata=fpkm, chronage=chronage1)

        with self.assertRaises(AssertionError):
            chronage1["chronage"] = ["SRR2166176", "SRR2167642"]
            age_calc.predict_age(exprdata=fpkm, chronage=chronage1)

        with self.assertRaises(AssertionError):
            chronage1["sampleid"] = ["SRR2166176", "SRR2166176"]
            age_calc.predict_age(exprdata=fpkm, chronage=chronage1)

        with self.assertRaises(AssertionError):
            chronage1["sampleid"] = ["SRR2166176", "unknown"]
            age_calc.predict_age(exprdata=fpkm, chronage=chronage1)

        with self.assertRaises(AssertionError):
            chronage1["chronage"] = [-20, 50]
            age_calc.predict_age(exprdata=fpkm, chronage=chronage1)

    def test_assert_raises5(self):
        rawcount1 = rawcount
        age_calc = RNAAgeCalc(tissue="brain", idtype="symbol")
        with self.assertRaises(AssertionError):
            age_calc._count2FPKM(rawcount=rawcount1, genelength=1)

        with self.assertRaises(AssertionError):
            genelength1 = pd.DataFrame(list(range(rawcount1.shape[0]-1)))
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = pd.DataFrame(["unknown"]*(rawcount1.shape[0]-1))
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = pd.Series(list(range(rawcount1.shape[0]-1)))
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = pd.Series(["unknown"]*(rawcount1.shape[0]-1))
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = list(range(24989))
            genelength1[0] = "unknown"
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = list(range(24988))
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = np.array(["unknown"]*rawcount1.shape[0])
            genelength1 = genelength1.reshape(1, rawcount1.shape[0])
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

        with self.assertRaises(AssertionError):
            genelength1 = list(range(24989))
            genelength1[0] = -1
            age_calc._count2FPKM(rawcount=rawcount1, genelength=genelength1)

    # def test_assert_raises6(self):
    #    res = rnaagepy.predict_age(exprdata=fpkm)
    #    with self.assertRaises(AssertionError):
    #        myplot = rnaagepy.makeplot(res)

    def test_assert_raises7(self):
        with self.assertRaises(AssertionError):
            RNAAgeCalc(stype="unknown")

    def test_getter_setter(self):
        age_calc = RNAAgeCalc(tissue="brain")
        with self.assertRaises(AssertionError):
            age_calc.exprtype = "unknown"
        age_calc.exprtype = "count"
        self.assertEqual(age_calc.exprtype, "count")
        self.assertEqual(age_calc._exprtype, "count")

        with self.assertRaises(AssertionError):
            age_calc.idtype = "unknown"
        age_calc.idtype = "refseq"
        self.assertEqual(age_calc.idtype, "refseq")
        self.assertEqual(age_calc._idtype, "refseq")

        with self.assertRaises(AssertionError):
            age_calc.stype = "unknown"
        age_calc.stype = "Caucasian"
        self.assertEqual(age_calc.stype, "Caucasian")
        self.assertEqual(age_calc._stype, "Caucasian")

        age_calc2 = RNAAgeCalc(exprtype="count", idtype="refseq")
        self.assertEqual(age_calc2.exprtype, "count")
        self.assertEqual(age_calc2._exprtype, "count")
        self.assertEqual(age_calc2.idtype, "refseq")
        self.assertEqual(age_calc2._idtype, "refseq")
        self.assertEqual(age_calc2.stype, "all")
        self.assertEqual(age_calc2._stype, "all")

    def test_getter_setter2(self):
        age_calc = RNAAgeCalc()
        self.assertEqual(age_calc._tissue, "all_tissues")
        self.assertEqual(age_calc.tissue, "all_tissues")
        self.assertEqual(age_calc._signature, "GTExAge")
        self.assertEqual(age_calc.signature, "GTExAge")
        with self.assertRaises(AssertionError):
            RNAAgeCalc(signature="unknown")
        age_calc = RNAAgeCalc(signature="DESeq")
        self.assertEqual(age_calc._signature, "Pearson")
        self.assertEqual(age_calc.signature, "Pearson")
        with self.assertRaises(AssertionError):
            RNAAgeCalc(tissue=1)
        age_calc = RNAAgeCalc(tissue="brain")
        self.assertEqual(age_calc._tissue, "brain")
        self.assertEqual(age_calc.tissue, "brain")
        self.assertEqual(age_calc._signature, "DESeq")
        self.assertEqual(age_calc.signature, "DESeq")
        age_calc = RNAAgeCalc(tissue="brain", signature="Peters")
        self.assertEqual(age_calc._signature, "Peters")
        self.assertEqual(age_calc.signature, "Peters")
        age_calc = RNAAgeCalc(tissue="notfound")
        self.assertEqual(age_calc._tissue, "all_tissues")
        self.assertEqual(age_calc.tissue, "all_tissues")
        self.assertEqual(age_calc._signature, "GTExAge")
        self.assertEqual(age_calc.signature, "GTExAge")
        with self.assertRaises(AssertionError):
            age_calc = RNAAgeCalc(tissue="notfound",
                                  signature="unknown")
        age_calc = RNAAgeCalc(tissue="notfound", signature="DESeq")
        self.assertEqual(age_calc._signature, "Pearson")
        self.assertEqual(age_calc.signature, "Pearson")
        with self.assertRaises(SetterError):
            age_calc = RNAAgeCalc(signature="DESeq")
            age_calc.signature = "Pearson"
        with self.assertRaises(SetterError):
            age_calc = RNAAgeCalc(signature="DESeq")
            age_calc.tissue = "brain"


if __name__ == '__main__':
    unittest.main()
