import numpy as np
import pandas as pd
import os
from pathlib import Path
import logging
import time
import warnings

import pyascore
import spec_parsers

logger = logging.getLogger(__name__)


def run_ascore(filepath_raw, filepath_psm_min, filepath_output, mode):
    # dir containing the raw files
    directory_raw = Path(filepath_raw)

    # dir containing the parsed psm files
    filepath_psm_min = Path(filepath_psm_min)

    # make new dir where files are saved in case the output dir does not exist
    filepath_output = Path(filepath_output)

    if not filepath_output.exists():
        os.mkdir(filepath_output)

    def raw_loop(raw_file, raw_file_path):
        # enter the loop extracting infromation of the individual RAW files

        # feed in spectra info
        spectra_file = raw_file_path
        spectra_parser = spec_parsers.SpectraParser(str(spectra_file), "mzML", ms_level=2)

        # To a dictionary of dictionaries:
        spectra_objects = spectra_parser.to_dict()
        # spectra_objects
        logger.info("parsed spectra from raw file for " + str(raw_file))

        for psm_stringency in os.listdir(filepath_psm_min):
            # feed in PSM info
            # static mods are carbamidomethylation of cystein aswell as TMT reporter ion being bound to lysine residues
            warnings.filterwarnings("ignore", category=UserWarning)
            psm_file = filepath_psm_min / str(psm_stringency) / raw_file.replace(".mzML", ".tsv")
            psm_parser = pyascore.IdentificationParser(str(psm_file), "percolatorTXT",
                                                       static_mods={"K": 229.162932, "C": 57.021464})

            # To a list of dictionaries:
            psm_objects = psm_parser.to_list()
            # psm_objects
            logger.info("parsed psm info for " + str(psm_stringency) + ": " + str(raw_file))

            # Ascore
            mod_mass = 79.9663304
            ascore = pyascore.PyAscore(bin_size=100.,
                                       n_top=10,
                                       mod_group="STY",
                                       mod_mass=mod_mass,
                                       mz_error=.05,
                                       fragment_types="by")

            ascore.add_neutral_loss('s', 97.9339)
            ascore.add_neutral_loss('t', 97.9339)

            pyascore_results = []

            logger.info("created pyascore object and started looping thorugh psms for " +
                        str(psm_stringency) + ": " + str(raw_file))

            # Iterate through all masses
            for psm in psm_objects:

                # Look for correct modification to be present
                mod_select = np.isclose(psm["mod_masses"], mod_mass)
                nmods = np.sum(mod_select)

                if nmods >= 1:
                    # Get the corresponsding spectrum
                    spectrum = spectra_objects[psm["scan"]]

                    # Other modifications that sould not be considered (position, masses) but needs to be included
                    # for correct assignment
                    aux_mod_pos = psm["mod_positions"][~mod_select].astype(np.uint32)
                    aux_mod_masses = psm["mod_masses"][~mod_select].astype(np.float32)

                    # Score PSMs
                    ascore.score(
                        # Array of MZ values for each peak in a spectra.
                        mz_arr=spectrum["mz_values"],
                        # Array of intensity values for each peak in a spectra.
                        int_arr=spectrum["intensity_values"],
                        # The peptide string without any modifications or n-terminal markings. (plain sequence)
                        peptide=psm["peptide"],
                        # Number of unlocalized modifications on the sequence.
                        n_of_mod=np.sum(mod_select),
                        # Maximum fragment charge to be used for score calculations.
                        max_fragment_charge=psm["charge_state"] - 1,
                        # Positions of fixed modifications. Most modification positions should start at 1 with 0 being
                        # reserved for n-terminal modifications, as seems to be the field prefered encoding.
                        aux_mod_pos=aux_mod_pos,
                        # Masses of individual fixed postion modifications.
                        aux_mod_mass=aux_mod_masses)

                    # adding alternative sites to the output not possible at the moment
                    # when joining the alternative sites for spectra with neutral loss ions the script exits
                    # alt_sites = [",".join([str(site) for site in site_list])
                    #              for site_list in ascore.alt_sites]

                    # Store scores
                    pyascore_results.append({"scan": psm["scan"],
                                             "localized_peptide": ascore.best_sequence,
                                             "pepscore": ascore.best_score,
                                             "ascores": ";".join([str(s) for s in ascore.ascores])})
                                             #"AltSites": ";".join(alt_sites)})

                    # logger.info("appended ascore into resutls list for scan: " + str(psm["scan"]))

            # save as tsv file
            filepath_output_stringency = filepath_output / str(psm_stringency)
            if not filepath_output_stringency.exists():
                os.mkdir(filepath_output_stringency)

            pyascore_results_df = pd.DataFrame(pyascore_results)
            pyascore_results_df.to_csv(filepath_output_stringency / raw_file.replace(".mzML", ".tsv"), sep="\t")
            logger.info("saved pyascore results for " + str(psm_stringency) + ": " + str(raw_file))

            # return pd.DataFrame(pyascore_results)

    start = time.time()

    if mode == "mq":
        logger.error("mq mode not avaible in the moment")
        #raw_loop()

    elif mode == "simsi":

        simsi_psm = Path(filepath_psm_min)
        logger.info("started ascoring in simsi mode")

        for raw_file in os.listdir(directory_raw):
            raw_file_path = directory_raw / str(raw_file)
            raw_loop(raw_file, raw_file_path)

    else:
        print("Not a valid psm input format")

    end = time.time()
    logger.info(f"fnished ascoring with clocked wall time of {(end - start) / 60} minutes")
    logger.info("")


if __name__ == '__main__':
    run_ascore(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\mzML",
               r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\psm_minimal_simsi_CL_full",
               r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\Acores_simsi_UT_2",
               "simsi")
