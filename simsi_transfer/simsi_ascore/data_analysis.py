from pathlib import Path
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as sc
import numpy as np


class Datamanipulation:

    def __init__(self, result_path):
        self.result_path = Path(result_path)

        self.stringencies = []
        self.stringencies_conc = []
        for file in os.listdir(self.result_path):
            df = pd.read_csv(self.result_path / file / str(file + "_msms.txt"), sep="\t", low_memory=False,
                             index_col=0)
            self.stringencies.append(df)

            df["stringency"] = str(file)
            self.stringencies_conc.append(df)

        self.stringencies_conc = pd.concat(self.stringencies_conc)
        self.stringencies_conc.reset_index(drop=True, inplace=True)

    def histplot(self):
        # considering the peptide scores of all peptides in each stringency
        fig, axs = plt.subplots(1, len(self.stringencies))

        i = 0
        for element in self.stringencies:
            sns.histplot(data=element["pepscore"], ax=axs[i])
            i += 1

        plt.show()
        plt.close()

    def boxplot(self):
        # considering the peptide scores of all peptides in each stringency
        fig, axs = plt.subplots(1, len(self.stringencies))

        i = 0
        for element in self.stringencies:
            sns.boxplot(data=element["pepscore"], ax=axs[i])
            i += 1

        plt.show()
        plt.close()


ex1 = Datamanipulation(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results")
ex1.histplot()
ex1.boxplot()

# considering only the peptide score of peptides which were transfered #################################################

p10_trans = ex1.stringencies[0][ex1.stringencies[0]["identification"] == "t"]["pepscore"].dropna()
p15_trans = ex1.stringencies[1][ex1.stringencies[1]["identification"] == "t"]["pepscore"].dropna()
p20_trans = ex1.stringencies[2][ex1.stringencies[2]["identification"] == "t"]["pepscore"].dropna()

trans = ex1.stringencies_conc[ex1.stringencies_conc["identification"] == "t"].dropna(subset="pepscore")

# boxplots
fig, axs = plt.subplots(1, 3, figsize=(15, 10))
sns.boxplot(data=p10_trans, ax=axs[0])
axs[0].set_xlabel("p10")
sns.boxplot(data=p15_trans, ax=axs[1])
axs[1].set_xlabel("p15")
sns.boxplot(data=p20_trans, ax=axs[2])
axs[2].set_xlabel("p20")
plt.show()
plt.close()

# histplots
fig, axs = plt.subplots(1, 3, figsize=(15, 10))
sns.histplot(data=p10_trans, ax=axs[0])
axs[0].set_xlabel("p10")
sns.histplot(data=p15_trans, ax=axs[1])
axs[1].set_xlabel("p15")
sns.histplot(data=p20_trans, ax=axs[2])
axs[2].set_xlabel("p20")
plt.show()
plt.close()

fig, axs = plt.subplots(1, 2, figsize=(15, 10))
sns.histplot(data=trans, x="pepscore", hue="stringency", element="poly", ax=axs[0])
sns.histplot(data=trans, x="pepscore", hue="stringency", element="poly", stat="density", common_norm=False, ax=axs[1])
axs[0].set_title("No normalization")
axs[1].set_title("Normalization due to different sample size")
plt.show()
plt.close()

########################################################################################################################
# statistical analysis

len(p10_trans)
len(p15_trans)
len(p20_trans)

p10_trans.median()
p15_trans.median()
p20_trans.median()

# checking if data fullfiles creiteria for multiple comparison with ANOVA and TK########################################
# 1. data are independent from each other -> I htink so
# 2. data normally distributed? -> seemes to be the case

fig, axs = plt.subplots(1, 3, figsize=(20, 10))
sns.histplot(data=p10_trans, ax=axs[0])
sns.histplot(data=p15_trans, ax=axs[1])
sns.histplot(data=p20_trans, ax=axs[2])
plt.show()
plt.close()

# 3. The population standard deviations of the groups are all equal. This property is known as homoscedasticity.
np.std(p10_trans)
np.std(p15_trans)
np.std(p20_trans)
# seemes equal to me, levene test of equality of variance?
levene = sc.levene(p10_trans, p15_trans, p20_trans)

# p val is 0.146, so above 0.05 so we accept the null hypothesis that the difference in variances between the
# groups is 0

# performing the statistical testing

F, p = sc.f_oneway(p10_trans,
                   p15_trans,
                   p20_trans)

res = sc.tukey_hsd(p10_trans,
                   p15_trans,
                   p20_trans)
# print(res)
# all stringencies have singificantly different peptide scores between each other regarding transfered peptides

########################################################################################################################
# checking if loclaized sequnece strings are of the same format as in the original simsi mod seuqnce column
mq_comp = Datamanipulation(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results")
stringencies_conc = mq_comp.stringencies_conc

# fitle rout for scans/peptides which contain phospho peptides
ser = stringencies_conc["localized_peptide"] == "_nan_"
stringencies_conc_sty = stringencies_conc[~ser]

# number of identical sequences
print(sum(~(stringencies_conc_sty["Modified sequence"] == stringencies_conc_sty["localized_peptide"])))

# df containingn the scans which have different mod and loc sequences
stringencies_conc_sty_diff = stringencies_conc_sty[
    ~(stringencies_conc_sty["Modified sequence"] == stringencies_conc_sty["localized_peptide"])]

# df containing the above but only if scans have at least one ascore above 0 ###########################################

# first have to convert the ascore column from string into float
ascores = stringencies_conc_sty_diff["ascores"].apply(lambda x: x.split(";"))
ascores = ascores.apply(lambda x: [float(i) for i in x])
stringencies_conc_sty_diff["ascores"] = ascores

# test_list = stringencies_conc_sty_diff["ascores"][0]
bol_zeros = []


def contains_non_zero(listx):
    nr_zeros = sum([i > 0.0 for i in listx])
    if nr_zeros > 0:
        bol_zeros.append(True)
    else:
        bol_zeros.append(False)


stringencies_conc_sty_diff["ascores"].apply(contains_non_zero)

# the scans in this dataframe have different mod sequences according to Ascore and maxquant output
stringencies_conc_sty_diff_non0 = stringencies_conc_sty_diff[bol_zeros]

# this dataframe is the same as above but only containing transfered scans
stringencies_conc_sty_diff_non0_t = stringencies_conc_sty_diff_non0[
    stringencies_conc_sty_diff_non0["identification"] == "t"]

########################################################################################################################
# actuall analysis with above data

# comparing number of different mod sequences for STY peptides for the p10 stringency###################################
sty_p10_diff = len(stringencies_conc_sty_diff_non0[stringencies_conc_sty_diff_non0["stringency"] == "p10"])
sty_p10_tot = len(stringencies_conc_sty[stringencies_conc_sty["stringency"] == "p10"])

(100 / sty_p10_tot) * sty_p10_diff
# 15.528% of sequences have discrepancies in mod sequence in p10 for both d and t scans

# same as above but for t scans only
sty_p10_diff_t = len(stringencies_conc_sty_diff_non0_t[stringencies_conc_sty_diff_non0_t["stringency"] == "p10"])
sty_p10_tot_t = len(stringencies_conc_sty[(stringencies_conc_sty["stringency"] == "p10") & (
            stringencies_conc_sty["identification"] == "t")])

(100 / sty_p10_tot_t) * sty_p10_diff_t
# 24.24% of t scans in the p10 stringency have a discrepancy in mq and ascore modified sequences! (only considering
# sequences with at least one ascore above 0 as truely different)

# comparing number of different mod sequences for STY peptides for the p15 stringency###################################
sty_p15_diff = len(stringencies_conc_sty_diff_non0[stringencies_conc_sty_diff_non0["stringency"] == "p15"])
sty_p15_tot = len(stringencies_conc_sty[stringencies_conc_sty["stringency"] == "p15"])

(100 / sty_p15_tot) * sty_p15_diff
# 15.25% of sequences have discrepancies in mod sequence in p15!

# same as above but for t scans only
sty_p15_diff_t = len(stringencies_conc_sty_diff_non0_t[stringencies_conc_sty_diff_non0_t["stringency"] == "p15"])
sty_p15_tot_t = len(stringencies_conc_sty[(stringencies_conc_sty["stringency"] == "p15") & (
            stringencies_conc_sty["identification"] == "t")])

(100 / sty_p15_tot_t) * sty_p15_diff_t
# 24.38% of t scans in the p15 stringency have a discrepancy in mq and ascore modified sequences! (only considering
# sequences with at least one ascore above 0 as truely different)

# comparing number of different mod sequences for STY peptides for the p20 stringency###################################
sty_p20_diff = len(stringencies_conc_sty_diff_non0[stringencies_conc_sty_diff_non0["stringency"] == "p20"])
sty_p20_tot = len(stringencies_conc_sty[stringencies_conc_sty["stringency"] == "p20"])

(100 / sty_p20_tot) * sty_p20_diff
# 14.93% of sequences have discrepancies in mod sequence in p20!

# same as above but for t scans only
sty_p20_diff_t = len(stringencies_conc_sty_diff_non0_t[stringencies_conc_sty_diff_non0_t["stringency"] == "p20"])
sty_p20_tot_t = len(stringencies_conc_sty[(stringencies_conc_sty["stringency"] == "p20") & (
            stringencies_conc_sty["identification"] == "t")])

(100 / sty_p20_tot_t) * sty_p20_diff_t
# 23.32% of t scans in the p15 stringency have a discrepancy in mq and ascore modified sequences! (only considering
# sequences with at least one ascore above 0 as truely different)

# taking a look at the differnet sequences for p20
sty_p20_diff_t_df = stringencies_conc_sty_diff_non0_t[stringencies_conc_sty_diff_non0_t["stringency"] == "p20"]
sty_p20_diff_t_df.to_csv(Path(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\data_analysis\diff_sequence_p20.tsv")
                         , sep="\t")

########################################################################################################################
# taking a look into clusters and their assigned sequences from mq vs ascore comparing inter cluster t vs d scans

stringencies_conc_sty_diff_non0.to_csv(Path(r"C:\Users\mlorenz\Desktop\Internship\pyAscore "
                                            r"test\data_analysis\diff_sequence.tsv"), sep="\t")

# df containingn the scans which have different mod and loc sequences
stringencies_conc_sty["mod sequence overlapp"] = stringencies_conc_sty["Modified sequence"] == \
                                                 stringencies_conc_sty["localized_peptide"]


# df containing the above but only if scans have at least one ascore above 0 ###########################################

def ascore_non_zero(df_wascore):
    # first have to convert the ascore column from string into float
    ascores = df_wascore["ascores"].apply(lambda x: x.split(";"))
    ascores = ascores.apply(lambda x: [float(i) for i in x])
    df_wascore["ascores"] = ascores

    # test_list = df_wascore_diff["ascores"][0]
    bol_zeros = []

    def contains_non_zero(listx):
        nr_zeros = sum([i > 0.0 for i in listx])
        if nr_zeros > 0:
            bol_zeros.append(True)
        else:
            bol_zeros.append(False)

    df_wascore["ascores"].apply(contains_non_zero)
    return bol_zeros


bol_zeros = ascore_non_zero(stringencies_conc_sty)
stringencies_conc_sty_non0 = stringencies_conc_sty[bol_zeros]
