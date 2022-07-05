#!/usr/bin/python

"""
Parse genotype-phenotype data from a data-base, Parse sequence alignment file. Match SNVs to phenotypes, return report.

class gbdb2lut: parse phenotype effect database, 3 columns: gene, phenotype, effect. 

class alnmt2pheno: parse sequence alignment, matchh to phenotypes, return report.

By Carljin Boghart & Joao Santos
@INSA, @Surrey
"""

import logging
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import sqlalchemy
from Bio import AlignIO

pd.options.mode.chained_assignment = None  # To remove warning messages https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas

## general TODOs :


class gpdb2lut:
    def __init__(
        self,
        input,
        gene,
        sheet="S",
        gencol=1,
        phencol=7,
        table="pokay",
        odir="",
        output="",
    ):
        """
        Parse genotype-phenotype association data from input excel, tsv or sqlite database.

        :type input: str
        :type gene: str
        :type sheet: str
        :type gencol: int
        :type phencol: int
        :type table: str
        :type odir: str
        :type output: str

        :param input: input file (excel database or .tsv file or .db file)
        :param gene: gene of interest.
        :param sheet: sheet name in excel database.
        :param gencol: column number of genotype in excel database.
        :param phencol: column number of phenotype in excel database.
        :param table: table name in .db file.
        :param odir: output directory.
        :param output: output file name.

        """

        self.input = input  # "input/testdb.xlsx"
        self.gene = gene
        self.sheet = sheet  # "S"
        self.genotype = gencol  # 1
        self.phenotype = phencol  # 7
        self.dbtable = table

        if odir and odir[-1] != "/":
            odir += "/"

        self.output = odir + output  # "test2/testdb_s.csv"

    def read_in(self):
        """
        read in genotype-phenotype association data from self.input (excel database)
        :return: dataframe with 2 columns: genotype (lineage or mutation), and phenotype
        """
        input_type = os.path.splitext(self.input)[1]

        def input_db_process(input_db, gene):
            if "Gene" in input_db.columns:
                if self.gene not in input_db["Gene"].values:
                    logging.info(f"Gene {gene} not found in input database. exiting...")
                    sys.exit(1)

                input_db = input_db[input_db["Gene"] == gene].reset_index(drop=True)
            cols = ["Mutation", "Phenotype category"]
            if "Flag" in input_db.columns:
                cols += ["Flag"]

            input_db = input_db[cols]

            return input_db

        if input_type == ".xlsx":
            self.lut = pd.read_excel(
                self.input,
                sheet_name=self.sheet,
                engine="openpyxl",
                usecols=[self.genotype - 1, self.phenotype - 1],
            ).dropna(how="all")

        elif input_type in [".tsv", ".csv"]:
            input_db = (
                pd.read_csv(self.input, sep="\t")
                .dropna(how="all")
                .reset_index(drop=True)
            )

            self.lut = input_db_process(input_db, self.gene)

        elif input_type == ".db":
            if not self.dbtable:
                logging.info("No table name provided for input database. exiting...")
                sys.exit(1)
            engine = sqlalchemy.create_engine("sqlite:///{}".format(self.input))
            input_db = pd.read_sql(self.dbtable, engine)
            self.lut = input_db_process(input_db, self.gene)

        return self

    def add_flags(self):
        """
        add flags to genotype-phenotype association dataframe: D if direct association, P if partial association (i.e. as part of constellation of mutations)
        :return: look-up table with 3 columns: genotype (lineage or mutation), phenotype, flag (D or P)

        """
        if self.sheet == "Lineages":
            # Add a Flag column with default value 'D'
            self.lut.loc[:, "Flag"] = "D"
            # Deduplicate
            self.lut = self.lut.drop_duplicates()

        else:

            # self.read_ref()  # Reads in reference sequence
            self.process_mut()  # Processes mutations and adds flags accordingly

    def split_multidel(self, matchobj):
        """
        split up multi-residue deletions in the geno-pheno dataframe (incl looking up of deleted residues in reference sequence)
        :return: ", "-separated list of single residue deletions
        """
        del_positions = list(range(int(matchobj.group(1)), int(matchobj.group(2)) + 1))
        del_residues = [self.refseq[pos - 1] for pos in del_positions]
        return ", ".join(
            ["{}{}del".format(r, str(p)) for r, p in zip(del_residues, del_positions)]
        )

    def process_mut(self):
        """
        process mutation rows: split up constellations and multi-residue deletions and add flag P, or flag D for individual mutations; add protein code to mutations
        :return: look-up table with each individual mutation-phenotype combination in a separate row, with flag denoting D(irect) or P(artial) association.
        """
        # Split up multi-residue deletions
        for i, v in self.lut.iloc[:, 0].items():
            if "del" in str(v):
                self.lut.iloc[i, 0] = re.sub(
                    "del([\d]+)-([\d]+)", self.split_multidel, str(self.lut.iloc[i, 0])
                )

        # Turn mutations format from strings into lists
        self.lut.iloc[:, 0] = self.lut.iloc[:, 0].apply(
            lambda x: x.replace(" ", "").split(",")
        )

        # Add a Flag column with default value 'D'
        if "Flag" not in self.lut.columns:
            self.lut.loc[:, "Flag"] = "D"

        # Change Flag to P where part of a constellation or multi-residue deletion (length of mut list >1)
        for i, v in self.lut.iloc[:, 0].items():

            if len(v) > 1:

                if self.lut.loc[i, "Flag"] != "D":
                    self.lut.loc[i, "Flag"] = "P:" + self.lut.loc[i, "Flag"]
                else:
                    self.lut.loc[i, "Flag"] = "P"

        # Explode to turn mutations list-cols (constellations) into extra rows; deduplicate
        # print(list(self.lut.Flag))

        self.lut = self.lut.explode(self.lut.columns[0]).drop_duplicates()
        # Combine D and P, where both applicable
        self.lut = (
            self.lut.groupby([self.lut.columns[0], self.lut.columns[1]])
            .agg({"Flag": lambda x: ";".join(x)})
            .reset_index()
        )

        # Add protein code to mutation string
        self.lut.iloc[:, 0] = self.lut.iloc[:, 0].apply(lambda x: self.gene + ":" + x)

        return self

    def write_lut(self, sep="\t"):
        """
        Write look-up table to .tsv file
        :param sep: separator for .tsv file

        :return: .tsv file with geno-pheno look-up table
        """
        # Write lookup_table to csv
        if self.output:
            self.lut.to_csv(self.output + ".CLEAN.tsv", index=False, sep=sep)


class alnmt2pheno:
    def __init__(
        self,
        input,
        gene,
        reference,
        nucl=False,
        output="test",
        odir="test/",
    ):
        """
        Class for aligning and processing alignments to reference sequence, map to existing look-up table, and write output to .tsv file.

        :param input: input alignment file
        :param gene: gene name
        :param reference: reference sequence
        :param nucl: boolean, True if nucleotide alignment, False if amino acid alignment
        :param output: output file name
        :param odir: output directory

        :type input: str
        :type gene: str
        :type reference: str
        :type nucl: bool
        :type output: str
        :type odir: str


        :return: summary geno-pheno look-up table

        """

        self.input = input
        self.gene = gene
        self.reference = reference
        self.nucl = nucl
        if odir and odir[-1] != "/":
            odir += "/"
        self.odir = odir

        os.makedirs(self.odir, exist_ok=True)
        self.output = self.odir + output

    def read_input(self):
        """
        read in input (multifasta protein alignment) and geno-pheno database ("look-up table" resulting from gpdb2lut script)
        :return: self with alignment, reference sequence, ref seq index, and database
        """
        # read in alignment

        self.Alnmt = AlignIO.read(self.input, "fasta")
        ref_seq = self.reference
        # compares the full refseq identifier/description in the alnmt, to the ref_seq provided in argument.
        # This is a fix to deal with identifiers that (1) include spaces, (2) might consist of separate words of which the first is not unique.
        self.ref_seq_index = [
            i for i, s in enumerate(self.Alnmt) if s.description == ref_seq
        ]
        if len(self.ref_seq_index) == 0:
            raise ValueError(
                "Reference sequence not found in alignment. Please check alignment and reference sequence."
            )
            sys.exit(1)
        else:
            self.ref_seq_index = self.ref_seq_index[0]

        self.ref_seq = self.Alnmt[self.ref_seq_index]

    def read_db(self, database=""):
        """
        reads in database (geno-pheno look-up table)
        """
        if database:
            self.db = pd.read_csv(self.database, sep="\t")

            if self.db.shape[1] == 2:
                self.db.columns = ["Mutation", "Phenotype category"]
            if self.db.shape[1] == 3:
                self.db.columns = ["Mutation", "Phenotype category", "Flag"]

    def read_ref(self):
        """
        read in reference sequence
        :return: self with reference sequence

        """
        samples = [s.description for _, s in enumerate(self.Alnmt)]

        refseq = [
            s.seq for _, s in enumerate(self.Alnmt) if s.description == self.reference
        ]
        if len(refseq) == 0:
            logging.info("Reference sequence not found in alignment")
            self.refseq = None

        else:
            self.refseq = refseq[0]

        self.samples = samples

    def msa2snp(self, merge_dels=False):
        """
        code from https://github.com/pinbo/msa2snp/blob/master/msa2snp.py
        Copyright 2019 Junli Zhang <zhjl86@gmail.com>

        iterates over mutations only.
        generate dataframe with (1) positions that vary, (2) aa at these positions in the refseq, (3) any alternatives in the dataset (sep by /), (4-end) for each sequence the aa at these positions
        :return: self with allele overview dataframe
        """
        fasta = {s.description: s.seq for _, s in enumerate(self.Alnmt)}
        seqnames = [s.description for _, s in enumerate(self.Alnmt)]
        refid = self.reference
        refseq = str(fasta[refid])
        self.refseq = refseq

        refseq_nogap = refseq.replace("-", "")

        n = -1
        seq2 = {}  # new dictionary of rotated sequences with no-gap postions as keys

        for i in range(len(refseq)):  # i is the position of seq poition
            if refseq[i] != "-":
                n += 1
            templist = []
            for j in seqnames:  # j is seq name
                seq = fasta[j]
                templist.append(seq[i])
            if n not in seq2:
                seq2[n] = templist
            else:
                seq2[n] = [s1 + s2 for s1, s2 in zip(seq2[n], templist)]

        # output
        outlist = []

        for k in range(len(refseq_nogap)):
            seq = [
                w[0] + w[1:].replace("-", "") if len(w) > 1 else w for w in seq2[k]
            ]  # remove "-" in alleles in "A--"
            # seq= seq2
            alleles = set(seq)  # set: unique values

            if len(alleles) != 1:
                ref_allele = refseq_nogap[k]  # string
                alt_allele_set = alleles - set(ref_allele)  # set
                alt_allele = "/".join(alt_allele_set)  # string
                outlist.append(
                    [str(k + 1), ref_allele, alt_allele] + seq
                )  # K+1: starting from 1

        # merge continuous deletions: ATTT to A---
        header = ["Pos", "ref", "alt"] + seqnames

        if merge_dels:

            def allpos(alist, avalue):
                return [pos for pos, element in enumerate(alist) if element == avalue]

            outlist2 = []  # merged continuous deletions
            mm = []  # temporary merging of lists
            n = 0  # the start position of continious deletion
            for i in range(len(outlist) - 1):
                p1 = outlist[i]
                p2 = outlist[i + 1]
                if (
                    int(p1[0]) + 1 == int(p2[0])
                    and "-" in p1
                    and "-" in p2
                    and allpos(p1, "-") == allpos(p2, "-")
                ):
                    if mm:  # if mm already have merged lists
                        mm = [s1 + s2 for s1, s2 in zip(mm, p1)]
                    else:
                        mm = p1
                        n = p1[0]
                else:
                    if mm:
                        mm = [s1 + s2 for s1, s2 in zip(mm, p1)]
                        mm[0] = n
                        outlist2.append(
                            ["-" if "-" in w else w for w in mm]
                        )  # replace "-----" to "-"
                        mm = []
                    else:
                        outlist2.append(p1)
                    if i + 1 == len(outlist) - 1:
                        outlist2.append(p2)

            outlist = outlist2

        outframe = pd.DataFrame(outlist)
        outframe.columns = header
        outframe.Pos = outframe.Pos.astype(int)

        self.vcf = outframe
        self.snps = outframe.Pos
        return self

    def find_stops(self):
        """
        determine "stop" positions for each sequence, to be able to identify sequences with early stops
        ## see https://github.com/pinbo/msa2snp/blob/master/msa2snp.py for reference.
        :return: self with an array of stop positions for each sequence
        """

        seq_matrix = [list(x.seq) for x in self.Alnmt]
        seq_matrix = np.array(seq_matrix, dtype=str)
        gapm = np.array(seq_matrix == "-", dtype=int)

        gapn = np.cumsum(gapm, axis=1)
        gapn = gapn[self.ref_seq_index]

        early_array = [seq_matrix.shape[1]] * seq_matrix.shape[0]
        # early_array = np.array(early_array)

        early_stops = np.where(seq_matrix == "*")
        early_finds = []

        for x, z in enumerate(early_stops[0]):
            est = early_stops[1][x]

            early_array[early_stops[0][x]] = est + 1 - gapn[est]
            early_finds += [est + 1 - gapn[est]]
        self.early_stops = early_array

        return self

    def vcf_format(self):
        """
        construct 3 dictionaries:
        (1) pos/mut_position links each mutation (key) with the position it occurs at (value);
        (2) seq/mut_seq links each mutation (key) with list of sequences that contain it (value) - NB this includes 'alt' as sequence for pos with single alternative allele
        (3) names/rec_var is dict of dicts, links each position (outer key) with each alternative allele at this position (value/ inner key) and associated name of the mutation (value / value)
        :return: self, with a dictionary of the three above dictionaries
        """
        mut_position = {}  # dictionary w mutations as keys, and positions as values
        mut_seq = defaultdict(
            list
        )  # dictionary w mutations as keys, and list of sequences as values
        rec_var = {}

        last = 0
        for i in range(self.vcf.shape[0]):
            row = self.vcf.iloc[i]
            alts = row.alt.split("/")

            for var in alts:
                if not var:
                    continue
                if var == "-":
                    name = f"{self.gene}:{row.ref}{row.Pos}del"

                elif len(var) > len(row.ref):
                    insm = var
                    if var[0] == row.ref[0]:
                        insm = var[1:]

                    name = f"{self.gene}:ins{row.Pos}{insm}"
                else:
                    name = f"{self.gene}:{row.ref}{row.Pos}{var}"

                ## if necessary to use mutation type later separately,
                name, pos = name, row.Pos

                if pos != last:
                    rec_var[pos] = {}
                    last = pos

                rec_var[pos][var] = name
                mut_position[name] = pos
                # print(np.array(row, dtype= str) == 'alt')

                td = [(i, x) for i, x in enumerate(row) if x == "alt"]

                mut_seq[name] = [
                    self.vcf.columns[i] for i, g in enumerate(row) if g == var and i > 2
                ]

        self.mut_dict = {"pos": mut_position, "seq": mut_seq, "names": rec_var}

        return self

    def gen_mut_matrix(self, index_sort=True):
        """
        generate binary matrix: for each mutation detected in the dataset, and for each sequence in the dataset, whether the sequence contains this particular mutation.
        :return: self with this matrix
        """
        vcf2 = self.vcf.copy()
        ## deal with multiple alleles:
        vcf2.alt = vcf2.alt.apply(lambda x: x.split("/"))
        vcf2 = vcf2.explode("alt").drop_duplicates()

        vcf2 = vcf2[(vcf2.alt != "") & (vcf2.alt != "*")].reset_index(drop=True)
        boolm = np.array(vcf2.iloc[:, 3:]).T == np.array(vcf2.alt)

        vcf2.iloc[:, 3:] = np.array(boolm, dtype=int).T

        ## mut NA's with coordinates above early stops.
        logging.info("NA positions above early stops")

        for i, stop in enumerate(self.early_stops):
            rowname = vcf2.columns[3 + i]
            colh = vcf2[rowname]
            colh[(vcf2.Pos > stop)] = "NA"
            vcf2[rowname] = colh

        mut_matrix = vcf2.iloc[:, 3:].T

        mut_names = [
            self.mut_dict["names"][vcf2.Pos[x]][vcf2.alt[x]]
            for x in range(vcf2.shape[0])
        ]
        mut_matrix.columns = mut_names

        if index_sort:
            logging.info("sorting index")
            mut_matrix.rename(
                index={
                    mut_matrix.iloc[self.ref_seq_index].name: "Reference: "
                    + mut_matrix.iloc[self.ref_seq_index].name
                },
                inplace=True,
            )
            mut_matrix = mut_matrix.reindex(
                [mut_matrix.index[self.ref_seq_index]]
                + [
                    e
                    for e in list(mut_matrix.index)
                    if e != mut_matrix.index[self.ref_seq_index]
                ]
            )

        mut_matrix = self.clean_mut_matrix(
            mut_matrix, name="mut_marix", outfile="fpass.csv"
        )

        self.mut_matrix = mut_matrix

        return self

    def clean_mut_matrix(self, df, name="matrix", outfile="fpass.csv"):
        """
        remove duplicate index, duplicate columns.
        """
        logging.info(f"Cleaning {name} shape {df.shape}")
        unique_alldb = df.index.is_unique
        logging.info(f"{name} index is unique, first pass: {unique_alldb}")
        if not unique_alldb:

            dup_index = df.index.duplicated()
            logging.info(f"{sum(dup_index)} duplicated indeces found")
            dup_mut = df.columns[dup_index]
            with open(self.output + "dupindex_" + outfile, "w") as f:
                f.write("".join([str(x) for x in dup_mut]))
            df = df.loc[:, ~dup_index]

        alldb_cols_unique = df.columns.is_unique
        logging.info(f"{name} cols is unique, first pass: {alldb_cols_unique}")
        if not alldb_cols_unique:

            dup_index = df.columns.duplicated()
            logging.info(f"{sum(dup_index)} duplicated columns found")
            dup_cols = df.columns[dup_index]
            with open(self.output + "dupcols_" + outfile, "w") as f:
                f.write("".join([str(x) for x in dup_cols]))
            df = df.loc[:, ~dup_index]

        return df

    def all_mutations_matrix(self):
        """
        Create a dataframe with all mutations included in the genopheno look-up table (all db mutation matrix)
        generate binary matrix: for each mutation in the genopheno look-up table, and for each sequence in the dataset, whether the sequence contains this particular mutation.
        NB: this does not "scan for" all mutations in the genopheno look-up table; rather it copies across the mut detected in the dataset from self.mut_matrix, and automatically adds 0s for any mut not in self.mut_matrix
        :return: self with this matrix, and a list of all mutations in the genopheno look-up table
        """
        db = self.db
        all_db_mutations = sorted(
            db["Mutation"].unique()
        )  # sorting all database mutations by position
        logging.info("sorted all_db_mutations, shape {}".format(len(all_db_mutations)))
        all_db_mut_matrix = pd.DataFrame(
            index=self.mut_matrix.index, columns=all_db_mutations, data=0
        )  # default = 0
        logging.info(
            "generated binary mutations matrix, shape {}".format(
                all_db_mut_matrix.shape
            )
        )

        all_db_mut_matrix.index.name = "Sequence"

        all_db_mut_matrix.update(
            self.mut_matrix
        )  # checks for overlapping mutations and copies values from mut_matrix (detected mutations)

        self.all_db_mutations = (
            all_db_mutations  # list of all mutations in the genopheno look-up table
        )
        self.all_db_mut_matrix = all_db_mut_matrix  # binary matrix: for each mutation in the genopheno look-up table, and for each sequence in the dataset, whether the sequence contains this particular mutation.

        return self

    def flags_matrix(self):
        """
        Create a dataframe with mutations detected in the dataset, and with phenotype flags
        Generate two matrices:
        (1) self.flags: for each phenotype category, and each mutation detected in the dataset, any flags indicating nature of association (D/P)
        (2) self.allflags: for each phenotype category, and ech mutation in genopheno lookup table, any flags indicating the nature of this association (D/P)
        :return: self with flag matrices, and a list of detected mutations
        """
        db = self.db
        detected_mutations = sorted(self.mut_dict["pos"], key=self.mut_dict["pos"].get)

        flags = pd.DataFrame(
            index=db["Phenotype category"].unique(), columns=detected_mutations, data=""
        )
        flags.index.name = "Phenotype"

        mut_dbfilter = np.isin(
            np.array(detected_mutations), db.Mutation
        )  ## filter to loop only over mutations present in db.
        mint = np.array(detected_mutations)[mut_dbfilter]

        for p in db["Phenotype category"].unique():
            pheno_muts = mint[
                np.isin(mint, db[(db["Phenotype category"] == p)]["Mutation"])
            ]  # preclude if statement.
            for m in pheno_muts:
                flags.loc[[p], [m]] = (
                    db["Flag"][(db["Mutation"] == m) & (db["Phenotype category"] == p)]
                    .to_string(index=False)
                    .strip()
                )
                # Equivalent flags dataframe for all db mutations

        all_db_flags = (
            db.pivot(index="Phenotype category", columns="Mutation", values="Flag")
            .reindex(columns=self.all_db_mutations)
            .replace(float("NaN"), "")
        )

        all_db_flags.index.name = "Phenotype"
        self.allflags = all_db_flags  # matrix with phenotype categories, all mutations in the genopheno lookup table, and flags indicating nature of association (D/P/"")
        self.flags = flags  # matrix with phenotype categories, all mutations detected in the dataset, and flags indicating nature of association between them (D/P/"")
        self.detected_mutations = detected_mutations  # list of detected mutations

        return self

    def filtered_tables(self):
        """
        Generate filtered mutation and flag matrices: these only contain detected mutations with at least one phenotype flag; additionally, low coverage positions are indicated with 'nc'
        Generate copies of "full" mutation and flag matrices (for all mutations detected in dataset), but with low coverage positions indicated with 'nc'
        Set low coverage positions to 'nc' also in "all database" matrices (self.allflags and self.all_db_mut_matrix)

        :return: self with 2 dictionaries, containing filtered and full nc matrices,
        """
        flags_nc = self.flags.copy(deep=True)
        mut_matrix_nc = self.mut_matrix.copy(deep=True)

        # Change mut_matrix_nc and all_db_mut_matrix entries for low coverage positions
        regex_find = {
            False: "([A-Za-z0-9]+:[A-Z*][0-9]+)X",
            True: "([A-Za-z0-9]+:[A-Z*][0-9]+)N",
        }
        subst_nc = re.compile(
            regex_find[self.nucl]
        )  # regex pattern for apparent substitution, but actually indicator of low coverage ## how does that work?
        # --> for any other substitution at position 203: replace 0 with 'nc'

        # all_muts= np.array(list(set(self.all_db_mutations + self.detected_mutations)))
        # detected_muts= list(self.detected_mutations)
        # detected_filter= np.array(np.isin(all_muts, self.detected_mutations))
        # all_db_filter= np.isin(all_muts, self.all_db_mutations)
        # db_posx= [re.findall('[0-9]+', x)[0] for x in self.all_db_mutations]

        #
        for mutation in self.detected_mutations:
            if subst_nc.match(
                mutation
            ):  # for each 'substitution' due to low coverage in detected mutations
                pos = self.mut_dict["pos"][mutation]
                sisters = list(self.mut_dict["names"][pos].values())
                # print(list(sisters))

                pattern_same_position = re.compile(
                    subst_nc.match(mutation).group(1) + "[A-Z*]"
                )  # pattern for other substitutions at the same position

                subst_same_position = list(
                    [m for m in sisters if pattern_same_position.match(m)]
                )  # list of other substitutions at the same position, within detected mutations
                subst_same_position_all_db = list(
                    [m for m in self.all_db_mutations if pattern_same_position.match(m)]
                )  # list of other substitutions at the same position, within all database mutations

                sequences_nc = np.array(self.mut_matrix.index)[
                    self.mut_matrix[mutation] == 1
                ]

                mut_matrix_nc.loc[
                    sequences_nc, subst_same_position
                ] = "nc"  # set all combinations of sequences and substitutions to 'nc' (detected mutations)
                self.all_db_mut_matrix.loc[
                    sequences_nc, subst_same_position_all_db
                ] = "nc"  # set all combinations of sequences and substitutions to 'nc' (full db)
                # now remove columns with X in flags_nc and matrix_nc
                mut_matrix_nc = mut_matrix_nc.drop(columns=mutation)
                flags_nc = flags_nc.drop(columns=mutation)
                # print(flags_nc.columns)

        # --> for any del at 203 -> replace 0 with 'nc' or keep at 0? (i.e. how confident that D203X has existing X rather than gap?)
        # --> for any ins at 203 -> I suppose desired behaviour for D203insABC is also ‘nc’ as we can’t be sure X is 1 or more residues?

        # Filter detected mutation dataframes by flag -> detected & flagged mutations only
        flags_filtered = (
            flags_nc.replace("", float("NaN"))
            .dropna(axis=1, how="all")
            .replace(float("NaN"), "")
        )

        #
        mut_matrix_filtered = mut_matrix_nc[flags_filtered.columns]
        ##
        self.filtered = {
            "flags": flags_filtered,  # matrix with phenotype categories, detected mutations associated with them, and flags indicating nature of association (D/P/""); but with nc for low coverage positions
            "muts": mut_matrix_filtered,  # for each detected & flagged mutation, whether the sequence contains this particular mutation; but with nc for low coverage positions
        }

        self.nc = {
            "flags": flags_nc,  # matrix with phenotype categories, all mutations detected in the dataset, and flags indicating nature of association between them (D/P/""); but with nc for low coverage positions
            "muts": mut_matrix_nc,  # for each mutation detected in the dataset, and for each sequence in the dataset, whether the sequence contains this particular mutation; but with nc for low coverage positions
        }

    def prep_report(self):
        """
        prepare final report, showing for all sequences in dataset:
        (1) for each phenotype category, whether the sequence has a predicted altered phenotype (Y(es) or P(artial)) based on the genopheno associations
        (2) all mutations with phenotype flags, detected in the sequence
        (3) all mutations detected in the sequence

        # check:
        # sample Portugal/PT24068/2021 - flag "Receptor binding affinity"
        # is resulting in NA (why?)
        ### CB 20/03: is this still the case after I added "elif 'D, P' in f.values:" in the first few lines below? Or is this a different issue?

        :return: self with final report dataframe
        """
        final_report = pd.DataFrame(
            index=self.mut_matrix.index, columns=self.flags.index, data=""
        )

        mut_lists = (
            self.filtered["muts"]
            .replace(0, float("NaN"))
            .replace("nc", float("NaN"))
            .stack(dropna=True)
            .reset_index()
            .groupby("Sequence")["level_1"]
            .apply(list)
            .rename("Flagged mutations")
        )

        def summarise_flags(f):
            summary = [x.split(";") for x in f][0]

            return ";".join(set(summary))

        for seq in mut_lists.index:

            flag_summary = (
                self.filtered["flags"][mut_lists[seq]]
                .replace("", float("NaN"))
                .stack(dropna=True)
                .reset_index(name="Flag")
                .groupby("Phenotype")["Flag"]
                .apply(summarise_flags)
            )

            final_report.loc[seq][flag_summary.index] = flag_summary

        # Adding columns with Flagged mutations and All mutations
        final_report = final_report.merge(
            mut_lists.apply(";".join), left_index=True, right_index=True, how="left"
        )

        mut_lists_all = (
            self.nc["muts"]
            .replace(0, float("NaN"))
            .replace("nc", float("NaN"))
            .stack(dropna=True)
            .reset_index()
            .groupby("Sequence")["level_1"]
            .apply(list)
            .rename("All mutations")
        )

        final_report = final_report.merge(
            mut_lists_all.apply(";".join), left_index=True, right_index=True, how="left"
        )

        final_report["Nflagged"] = final_report["Flagged mutations"].apply(
            lambda x: len(str(x).split(";"))
        )

        final_report["Nmutations"] = final_report["All mutations"].apply(
            lambda x: len(str(x).split(";"))
        )

        self.final_report = final_report

        return self

    def export(self, sep="\t"):
        """
        export final report, as well as various mutation and flag tables, as .tsv files
        full_mutation table ;
        final report ;

        :param sep:

        :return:
        """
        ext_dict = {",": ".csv", "\t": ".tsv"}
        ext = ext_dict[sep]

        # print mutation matrix only (with low coverage still masquerading as substitutions)
        # Output not particularly required
        self.mut_matrix.to_csv(
            self.output + "_all_mutation_matrix" + ext, na_rep="NA", sep=sep
        )

        # full table (mutation matrix [detected mutations] + flags, with low coverage positions indicated as 'nc')
        # Keep this output for now as unclear which "full table" format (this one with "nc", or the one above with "D203X" columns) is preferred
        full_table = pd.concat([self.nc["flags"], self.nc["muts"]], axis=0)
        full_table.to_csv(
            self.output + "_all_mutation_report" + ext, na_rep="NA", sep=sep
        )

        # filtered table (mutation matrix [detected mutations that are flagged] + flags, with low coverage positions indicated as 'nc')
        # Output not particularly required but I (Carlijn) think it could be useful
        full_table = pd.concat([self.filtered["flags"], self.filtered["muts"]], axis=0)
        full_table.to_csv(
            self.output + "_flagged_mutation_report" + ext, na_rep="NA", sep=sep
        )

        # all db mut table (mutation matrix [all db mutations] + flags, with low coverage positions indicated as 'nc')
        # Output REQUIRED (request by Vitor)
        full_table = pd.concat([self.allflags, self.all_db_mut_matrix], axis=0)
        full_table.to_csv(
            self.output + "_all_db_mutations_report" + ext, na_rep="NA", sep=sep
        )

        # final report
        # Output REQUIRED (request by Vitor)
        self.final_report.to_csv(self.output + "_final_report" + ext, sep=sep)
