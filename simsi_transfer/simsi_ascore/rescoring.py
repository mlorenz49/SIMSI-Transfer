from pathlib import Path
import os
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def rescore_seq(filepath_results, filepath_output):
    logger.info("started rescoring modified sequence")

    filepath_results = Path(filepath_results)

    # make new dir where files are saved in case the output dir does not exist
    filepath_output = Path(filepath_output)

    if not filepath_output.exists():
        os.mkdir(filepath_output)

    for stringency in os.listdir(filepath_results):

        # make new dir where files are saved in case the output dir does not exist
        stringency_output = filepath_output / stringency

        if not stringency_output.exists():
            os.mkdir(stringency_output)

        df_raw = pd.read_csv(filepath_results / stringency / str(stringency + "_msms.txt"), sep="\t",
                             low_memory=False, index_col=0)

        # preparing the data
        df_raw.fillna("nan", inplace=True)
        df_raw["modified sequence match"] = df_raw["Modified sequence"] == df_raw["localized_peptide"]  # do sequ match
        df_raw = convert_ascore(df_raw)  # convert ascores from string to int

        # actual rescoring happens here
        rescored = []
        rescored_sequence = []
        for i in range(0, len(df_raw)):

            if contains_non_zero(df_raw.loc[df_raw.index[i], "ascores"]):  # if not only zeros in ascore continue
                if not df_raw.loc[df_raw.index[i], "modified sequence match"]:  # if seqeunces differ continue
                    rescored_sequence.append(df_raw.loc[df_raw.index[i], "localized_peptide"])
                    rescored.append(True)
                    continue
                rescored_sequence.append(df_raw.loc[df_raw.index[i], "Modified sequence"])
                rescored.append(False)
                continue

            rescored_sequence.append(df_raw.loc[df_raw.index[i], "Modified sequence"])
            rescored.append(False)

        df_raw["rescored"] = rescored  # adding a column indicating if modified sequence was rescored
        df_raw["Modified sequence rescored"] = rescored_sequence # add a column with the rescored or orig sequence
        logger.info(f"number of rescored scans/peptides: {sum(rescored)} ({round(100 / len(df_raw) * sum(rescored), 1)}"
                    f"% in stringency {str(stringency)})")

        df_raw.to_csv(stringency_output / str(stringency + "_msms.txt"), sep="\t")
        # return df_raw, rescored

    logger.info("finished rescoring modified sequence")
    logger.info("")


def convert_ascore(df_wascore):
    # first have to convert the ascore column from string into float

    ascores = df_wascore["ascores"].apply(lambda x: x.split(";"))
    ascores = ascores.apply(lambda x: [float(i) for i in x])
    df_wascore["ascores"] = ascores

    return df_wascore


def contains_non_zero(listx):
    nr_zeros = sum([i > 0.0 for i in listx])
    if nr_zeros > 0:
        # bol_zeros.append(True)
        return True
    else:
        # bol_zeros.append(False)
        return False


if __name__ == "__main__":
    rescore_seq(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results",
                r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results_rescored")
