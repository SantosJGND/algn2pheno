#!/usr/bin/python


"""
Input arguments for the scripts in the aln2pheno package.

By Carljin Boghart & Joao Santos
@INSA, @Surrey
"""


def get_args_a2p():
    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")

        parser.add_argument(
            "--db",
            type=str,
            help="""phenotype db. if excel, provide sheet name and columns numbers \
                        for genotype and phenotype data. \
                            if not excel, provide path to db, ensure 3 columns 'Mutation', 'Phenotype category' & 'Flag'""",
        )

        parser.add_argument(
            "--sheet",
            "-s",
            type=str,
            required=False,
            default="S",
            help="Give sheet name (gene name or 'Lineages'). excel input.",
        )

        parser.add_argument(
            "--table",
            "-t",
            type=str,
            required=False,
            default="pokay",
            help="table name in db. default: pokay. sqlite input.",
        )

        parser.add_argument(
            "--gencol",
            required=False,
            type=int,
            default=1,
            help="nr of column with genotype data. excel input.",
        )
        parser.add_argument(
            "--phencol",
            required=False,
            type=int,
            default=7,
            help="nr of column with phenotype data. excel input.",
        )

        parser.add_argument(
            "-g", "--gene", required=True, help="Set gene or protein prefix"
        )

        parser.add_argument("--algn", required=True, help="Input alignment file")

        parser.add_argument(
            "-r",
            "--reference",
            required=True,
            help="Give full reference sequence as in alignment file (use quotation marks if this contains a space)",
        )
        parser.add_argument(
            "--nucl",
            action="store_true",
            default=False,
            help="provide if nucleotide alignment instead of amino acid.",
        )
        parser.add_argument(
            "--odir", "-o", default="test", type=str, help="output directory"
        )

        parser.add_argument("--output", default="test", help="Set output file prefix")
        parser.add_argument("--log", default="algn2pheno.log", help="logfile")

        parser.add_argument(
            "-f",
            "--force",
            action="store_true",
            default=False,
            help="overwrite existing files",
        )
        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)

    return args


def get_args_scrape():
    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")
        parser.add_argument(
            "--pokay", action="store_true", help="chose what db to parse and download"
        )
        parser.add_argument(
            "--medium", action="store_true", help="scrape data from medium web"
        )

        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)

    return args
