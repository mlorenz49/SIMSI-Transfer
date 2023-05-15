import argparse
import sys
import logging
from pathlib import Path
import os

import parsing
import pyA_test
import join_outputs
import rescoring

if __name__ == '__main__':

    # On Windows the subproccesses will import (i.e. execute) the main module at start. You need to insert an if
    # __name__ == '__main__': guard in the main module to vavoid creating subporecesses recursively.
    # https://stackoverflow.com/questions/18204782/runtimeerror-on-windows-trying-python-multiprocessing

    ####################################################################################################################
    # parsing from commandline

    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=str, help="path to raw files")
    parser.add_argument("--msms", type=str, help="path to mq file or SIMSI summary file")
    parser.add_argument("--output_pars", type=str, help="parsing results outputpath")
    parser.add_argument("--mode", type=str, help="'mq' if input is mq or 'simsi' if input is of simsi transfer")
    parser.add_argument("--output_ascore", type=str, help="ascoring results outputpath")
    parser.add_argument("--results", type=str, help="path to final results folder")

    # only used for debugging rn, normally input comes from command line as encoded above. So we are basically
    # emmulating a command line statement until we fixed the problems
    # sys.argv = ["simsi_ascore", "--raw", r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\mzML",
    #             "--msms", r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\summaries", "--output_pars",
    #             r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\psm_minimal_simsi_MD_3", "--output_ascore",
    #             r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\Acores_simsi_MD_3",
    #             "--mode", "simsi", "--results",
    #             r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results_MD_3"]

    args = parser.parse_args()

    ####################################################################################################################
    # section on logging
    # create logger
    logger = logging.getLogger()

    logger.setLevel(logging.INFO)

    # create console handler and file handler and set level to Info
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    fl = logging.FileHandler(Path(os.path.dirname(args.results)) / "simsi_ascore.log", 'w+')
    # substitute logging file path with Path(os.path.dirname(os.getcwd())) later on
    fl.setLevel(logging.INFO)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch and fl
    ch.setFormatter(formatter)
    fl.setFormatter(formatter)

    # add ch and fl to logger
    logger.addHandler(ch)
    logger.addHandler(fl)

    logger.info('entering ascoring package')

    ####################################################################################################################
    # actual code starts here

    parsing.run_msms_parser(args.raw, args.msms, args.output_pars, args.mode)

    pyA_test.run_ascore(args.raw, args.output_pars, args.output_ascore, args.mode)

    join_outputs.join_outputs(args.msms, args.output_ascore, args.results)

    rescoring.rescore_seq(args.results, args.results)

    logger.info('Exiting ascoring package')
