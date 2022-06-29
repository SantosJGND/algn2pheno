#!/usr/bin/python


"""
Main function for aln2pheno. 




By Carljin Boghart & Joao Santos
@INSA, @Surrey
"""

import logging
import os
import sys
import time

from aln2pheno.a2p_classes import alnmt2pheno, gpdb2lut
from aln2pheno.arg_input import get_args_a2p


def main():
    start = time.time()

    ###############################################################################
    ########  1.  SETUP  #########################################################
    args = get_args_a2p()

    if args.odir and args.odir[-1] != "/":
        args.odir += "/"

    if os.path.isdir(args.odir):
        if args.force:
            logging.info("Ouput directory exists. Overwriting.")
        else:
            logging.info("Output directory exists. Exiting. use -f to force overwrite.")
            sys.exit(1)
    else:
        os.makedirs(args.odir, exist_ok=True)

    intermediate_file_log = os.path.basename(args.db)
    intermediate_file_log = os.path.splitext(intermediate_file_log)[0]
    if os.path.splitext(args.db)[1] == ".db":
        intermediate_file_log += "." + args.table

    logging.basicConfig(
        filename=args.odir + args.log,
        filemode="w",
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )

    ###############################################################################
    ########  2.  READ ALIGNMENT  ################################################
    logging.info("Started")

    a2p = alnmt2pheno(
        args.algn,
        args.gene,
        args.reference,
        nucl=args.nucl,
        output=args.output,
        odir=args.odir,
    )

    a2p.read_input()
    a2p.read_ref()

    logging.info(f"Reference: {args.reference}")
    logging.info(f"Gene: {args.gene}")
    logging.info(f"Alignment file: {args.algn}")
    logging.info(f"Number of samples in alignment: {len(a2p.samples)}")
    logging.info(f"Alignment type: {['protein','nucleic'][int(args.nucl)]}")
    logging.info(f"Db: {args.db}")

    logging.info(f"Output: {args.odir}/{args.output}*")

    if a2p.refseq == None:
        return

    gpdb = gpdb2lut(
        args.db,
        args.gene,
        sheet=args.sheet,
        gencol=args.gencol,
        phencol=args.phencol,
        table=args.table,
        odir=args.odir,
        output=intermediate_file_log,
    )

    gpdb.refseq = a2p.refseq
    gpdb.read_in()
    gpdb.add_flags()
    gpdb.write_lut()

    a2p.read_input()

    a2p.db = gpdb.lut

    a2p.read_input()
    ###############################################################################
    ############  3.  ALIGNMENT TO PHENOTYPE  ####################################
    a2p.find_stops()
    a2p.msa2snp()
    a2p.vcf_format()
    logging.info("A total of {} variants were found.".format(len(a2p.snps)))

    a2p.gen_mut_matrix()
    try:
        a2p.all_mutations_matrix()
    except ValueError as e:
        logging.info("all_mutations_matrix failed: {}".format(e))
        elapsed = time.time() - start
        logging.info("Elapsed time: {}".format(elapsed))
        logging.info("Exiting")
        sys.exit(1)

    found_stop = (a2p.all_db_mut_matrix != 0).sum().sum()

    if found_stop == 0:
        logging.info("No mutations identified in database, please check database.")
        sys.exit(1)

    a2p.flags_matrix()
    a2p.filtered_tables()
    a2p.prep_report()
    a2p.export()

    end = time.time()

    elapsed = end - start

    logging.info(f"Finished.")
    logging.info(f"Elapsed time: {elapsed:0.4f} seconds.")


if __name__ == "__main__":

    ###########################################################################

    main()
