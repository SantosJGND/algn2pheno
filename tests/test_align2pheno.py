#!/usr/bin/python


"""
Test functions in the alignment to phenotype module.

By Carljin Boghart & Joao Santos
@INSA, @Surrey
"""

import os
from re import S

import numpy
import pandas as pd
import pytest

from modules.a2p_classes import alnmt2pheno, gpdb2lut
from modules.arg_input import get_args_a2p, get_args_scrape
from modules.scrapefunc import pokay_dl

# class test_alnmt2pheno:z

gpbd = gpdb2lut(
    "tests/data/test_align2pheno/test_align2pheno.xlsx",
    "S",
    output="test",
    odir="test",
)
a2p = alnmt2pheno(
    "tests/data/test_align2pheno/test_align2pheno.xlsx",
    "S",  # GENE PREFIX
    "HUMAN",
)


class Test_alnmt2pheno:
    def test_add_flags(self):

        dummy_data = numpy.array(
            [
                ["A12D", "High"],
                ["A12D", "Low"],
            ]
        )

        gpbd.lut = pd.DataFrame(dummy_data, columns=["Mutation", "Phenotype category"])
        gpbd.add_flags()

        assert gpbd.lut.loc[0, "Flag"] == "D"
