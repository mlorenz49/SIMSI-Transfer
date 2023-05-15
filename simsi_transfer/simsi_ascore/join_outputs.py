from pathlib import Path
import os
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def join_outputs(filepath_simsi_sum, filepath_ascores, filepath_output):
    logger.info("started merging ascoring results with original simsi outputfile")
    filepath_simsi_sum = Path(filepath_simsi_sum)
    filepath_ascores = Path(filepath_ascores)

    # make new dir where files are saved in case the output dir does not exist
    filepath_output = Path(filepath_output)

    if not filepath_output.exists():
        os.mkdir(filepath_output)

    for stringency in os.listdir(filepath_simsi_sum):

        # make new dir where files are saved in case the output dir does not exist
        stringency_output = filepath_output / stringency

        if not stringency_output.exists():
            os.mkdir(stringency_output)

        # loop through the dir containing all the psm files and the corresponding ascore file
        df_list = []
        for file in os.listdir(filepath_ascores / stringency):
            simsi_sum = pd.read_csv(filepath_simsi_sum / stringency / str(stringency + "_msms.txt"), sep="\t",
                                    low_memory=False, index_col=0)
            simsi_sum_raw = simsi_sum[simsi_sum["Raw file"] == file.replace(".tsv", "")]

            ascores = pd.read_csv(filepath_ascores / stringency / file, sep="\t", low_memory=False, index_col=0)
            ascores.rename(columns={"scan": "scanID"}, inplace=True)

            df_merged = pd.merge(left=simsi_sum_raw, right=ascores, on="scanID", how="outer")
            df_list.append(df_merged)

        result_stringen = pd.concat(df_list)
        result_stringen.reset_index(drop=True, inplace=True)

        # renaming the mod sequence so that its the same format as in the original simsi output file
        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].str.replace(f"[{round(79.9663)}]",
                                                                                                "(Phospho (STY))")
        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].str.replace(
            f"[{round(15.9949146)}]", "(Oxidation (M))")
        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].str.replace(
            f"n[{round(42.010565)}]", "(Acetyl (Protein N-term))")
        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].str.replace(
            f"n[{round(229.162932)}]", "")
        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].str.replace(
            f"[{round(57.021464)}]", "")
        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].str.replace(
            f"[{round(229.162932)}]", "")

        result_stringen["localized_peptide"] = result_stringen["localized_peptide"].apply(lambda x: "_" + str(x) + "_")
        result_stringen.to_csv(stringency_output / str(stringency + "_msms.txt"), sep="\t")
        logger.info("succsefully saved file to: " + str(stringency_output / str(stringency + "_msms.txt")))

    logger.info("finsihed joining outputs")
    logger.info("")


if __name__ == "__main__":
    join_outputs(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\summaries",
                 r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\Acores_simsi",
                 r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results")
