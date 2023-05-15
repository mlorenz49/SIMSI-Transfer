import pandas as pd
import os
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def run_msms_parser(filepath_raw, filepath_msms, filepath_output, mode):
    # dir containing the raw files
    directory_raw = filepath_raw

    # make new dir where files are saved, if dir does not exist create it
    filepath_output = Path(filepath_output)

    if not filepath_output.exists():
        os.mkdir(filepath_output)

    def raw_loop(dic, output_path=filepath_output):
        for file in os.listdir(directory_raw):
            msms_raw_idx = msms[msms[dic["raw"]] == file.replace(".mzML", "")]  # choose raw file

            # only keep relevant columns
            msms_raw_idx = msms_raw_idx[[dic["scan"], dic["rt"], dic["mo seq"], dic["charge"], dic["score"]]]

            # modify the format of modified sequence
            msms_raw_idx[dic["mo seq"]] = msms_raw_idx[dic["mo seq"]].str.replace(r"(Phospho (STY))",
                                                                                  "[79.9663]")
            msms_raw_idx[dic["mo seq"]] = msms_raw_idx[dic["mo seq"]].str.replace(r"(Oxidation (M))",
                                                                                  "[15.9949146]")
            msms_raw_idx[dic["mo seq"]] = msms_raw_idx[dic["mo seq"]].str.replace(r"(Acetyl (Protein N-term))",
                                                                                  "[42.010565]")
            msms_raw_idx[dic["mo seq"]] = msms_raw_idx[dic["mo seq"]].str.replace("_", "")

            # adding TMT as N-terminal fixed mod for all peptides, except if they are carrying a N-term. acetylation
            msms_raw_idx[dic["mo seq"]] = msms_raw_idx[dic["mo seq"]].apply(
                lambda x: "[229.162932]" + x if "[42.010565]" not in x else x)

            # change column names
            msms_raw_idx.columns = ["scan", "rt", "sequence", "charge", "percolator score"]
            msms_raw_idx.reset_index(drop=True, inplace=True)

            # save as tsv file
            msms_raw_idx.to_csv(output_path / file.replace(".mzML", ".tsv"), sep="\t")
            # return msms_raw_idx

            logger.info("succesfully saved psm_min to: " + str(output_path / file.replace(".mzML", ".tsv")))

    # mode dictionaries
    mq = {"raw": "Raw file", "scan": "Scan number", "rt": "Retention time", "mo seq": "Modified sequence",
          "charge": "Charge", "score": "Score"}

    simsi = {"raw": "Raw file", "scan": "scanID", "rt": "Retention time", "mo seq": "Modified sequence",
             "charge": "Charge", "score": "Score"}

    if mode == "mq":
        # df of maxquant msms output
        msms = pd.read_csv(filepath_msms, sep="\t", low_memory=False)
        logger.info("started parsing")
        raw_loop(mq)

    elif mode == "simsi":

        simsi_sum = Path(filepath_msms)
        logger.info("started parsing")

        for element in os.listdir(simsi_sum):

            filepath_output_string = filepath_output / element
            filepath_input_string = simsi_sum / element / str(element + "_msms.txt")

            if not filepath_output_string.exists():
                os.mkdir(filepath_output_string)

            # df of simsi transfer msms pX output file
            msms = pd.read_csv(filepath_input_string, sep="\t", low_memory=False)

            raw_loop(simsi, filepath_output_string)
    else:
        print("Not a valid psm input format")

    logger.info("finished parsing")
    logger.info("")


if __name__ == "__main__":
    run_msms_parser(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test_miniTOPAS\mzML",
                    r"C:\Users\mlorenz\Desktop\Internship\pyAscore test_miniTOPAS\summaries",
                    r"C:\Users\mlorenz\Desktop\Internship\pyAscore test_miniTOPAS\psm_minimal_simsi_CL",
                    "simsi")
