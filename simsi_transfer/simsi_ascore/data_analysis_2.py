import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import warnings
import os

from pyteomics import mzml, auxiliary, pylab_aux
import pylab
import pyascore
import simsi_ascore.spec_parsers as spec_parsers

# reading in raw file using pyteomics mzml class #######################################################################
raw = mzml.MzML(r"C:\Users\mlorenz\Desktop\Internship\pyAscore "
                r"test_miniTOPAS\mzML\03951_GA1_P043866_S00_I00_U01_R1_TMT11.mzML")

# print(next(raw))
raw.get_by_id('controllerType=0 controllerNumber=1 scan=2197')
raw.get_by_index(0)

# list(raw) # so the object from the class mzml.MzML is appearntly also an iterator (just as the object resulting from
# mzml.read(path/rawfile)

# plot spectrum
pylab.figure()
pylab_aux.plot_spectrum(raw.get_by_id('controllerType=0 controllerNumber=1 scan=2197'))
plt.show()
plt.close()

# reading in raw file using pyteomics reader() functionality ###########################################################
# format example of reading in a raw file from pyteomics:
with mzml.read(r"C:\Users\mlorenz\Desktop\Internship\pyAscore "
               r"test_miniTOPAS\mzML\03951_GA1_P043866_S00_I00_U01_R1_TMT11.mzML") as spectra:
    spectrum = next(spectra)

# plot spectrum
pylab.figure()
pylab_aux.plot_spectrum(spectrum)
plt.show()
plt.close()

# applied to my own data ###############################################################################################
# creates iterator called spectra form the mzml file using the mzml.read() function
spectra = mzml.read(r"C:\Users\mlorenz\Desktop\Internship\pyAscore "
                    r"test_miniTOPAS\mzML\03951_GA1_P043866_S00_I00_U01_R1_TMT11.mzML")

# instead of iterating over all values in the iterator using next(spectra), we can use list() on the iterrator to get
# all iterables (hence spectra) from the iterator at once
# note: doing this with list() is not really the point of using iterators, since there we dont want to generate all data
# at once

spectra_list = list(spectra)
# print(next(spectra)['id'], next(spectra)['id'])

spectra_ids = [x['id'] for x in spectra_list]
spectra_precursors = [x['precursorList'] for x in spectra_list]
spectra_precursors_counts = [x['count'] for x in spectra_precursors]
np.unique(spectra_precursors_counts, return_counts=True)

np.median(spectra_precursors_counts)
sns.histplot(spectra_precursors_counts, bins=5)
plt.show()
plt.close()

# checking the raw files from my previous test set in order to see how many precursors the scans have there

raw = mzml.MzML(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\mzML_2\03CPTAC_UCEC_P_PNNL_20170922_B1S3_f04.mzML")

spectra_list_2 = list(raw)

spectra_ids = [x['id'] for x in spectra_list_2]
spectra_precursors_2 = [x['precursorList'] for x in spectra_list_2]
spectra_precursors_2_counts = [x['count'] for x in spectra_precursors_2]
np.unique(spectra_precursors_2_counts, return_counts=True)

np.median(spectra_precursors_2_counts)
sns.histplot(spectra_precursors_2_counts, bins=5)
plt.show()
plt.close()

# in the first test data set, all scans only contained one precursor, unlike the new dataset which contains up to mul
# tiple precursors
# appearently the nr 11 here referes to the fact that we have an MS3 method going on here used for quan!: for each MS2
# scan, the top 11 peaks are tooken and fragmented in order to generate an MS3 scan used for quantification in order to
# only get the TMT tags!

# so: filter out scans which are only MS2 and look at the amount of precursors there

spectra_ms2 = [x for x in spectra_list if x["ms level"] == 2]
spectra_precursors_3 = [x['precursorList'] for x in spectra_ms2]
spectra_precursors_3_counts = [x['count'] for x in spectra_precursors_3]
np.unique(spectra_precursors_3_counts, return_counts=True)

# firas assumption was correct: we have 39972 MS2 spectra which have onl yone precursor, no scan in the dataset has more
# than one precursor

[x["ms level"] for x in spectra_list]
np.unique([x["ms level"] for x in spectra_list], return_counts=True)


# creating callable which determines if scan is MS2 or MS3 and returns True or False for _passesfilter in the specparser
def ms_lvl(spectral_object):
    if spectral_object["ms level"] == 2:
        return True
    else:
        return False


callable(ms_lvl)
# testing out the callable
ms_lvl(spectra_list[0])

ex2 = [s for s in spectra_list if ms_lvl(s)]

# my callable function is working, but the pyAscore custom filter option is trash so this isnt working
# custom fitlering only happens AFTER getting the scan, this means filtering happens only after getting the precursor

# checking raw file containing the deffective scan throwing: Process finished with exit code -1073740791 (0xC0000409)
# reading in raw file using pyteomics mzml class
spectra = mzml.MzML(
    r"C:\Users\mlorenz\Desktop\Internship\pyAscore "
    r"test_miniTOPAS\mzML_ms2_short\03951_GA4_P043866_S00_I00_U04_R1_TMT11.mzML")

spectra_list = list(spectra)
spectra_ms2 = [x for x in spectra_list if x["ms level"] == 2]
spectra_precursors = [x['precursorList'] for x in spectra_ms2]
spectra_precursors_counts = [x['count'] for x in spectra_precursors]
print(np.unique(spectra_precursors_counts, return_counts=True))

spectra_ms3 = [x for x in spectra_list if x["ms level"] == 3]
spectra_precursors = [x['precursorList'] for x in spectra_ms3]
spectra_precursors_counts = [x['count'] for x in spectra_precursors]
print(np.unique(spectra_precursors_counts, return_counts=True))

# check out deffekt scan nr. 16134
spectra.get_by_id('controllerType=0 controllerNumber=1 scan=16134')

# plot spectrum
pylab.figure()
pylab_aux.plot_spectrum(spectra.get_by_id('controllerType=0 controllerNumber=1 scan=16134'))
plt.show()
plt.close()

# suspect: only one huge peak in this MS2 scan

# checking if error in the scan is occruing in the spectrum parser from pyascore

psm_file = r"C:\Users\mlorenz\Desktop\Internship\pyAscore " \
           r"test_miniTOPAS\psm_minimal_simsi_CL\p10\03951_GA4_P043866_S00_I00_U04_R1_TMT11.tsv"

psm_parser = pyascore.IdentificationParser(str(psm_file), "percolatorTXT",
                                           static_mods={"K": 229.162932, "C": 57.021464})

warnings.filterwarnings("ignore", category=UserWarning)
psm_objects = psm_parser.to_list()
scan_16133 = [x for i, x in enumerate(psm_objects) if psm_objects[i]["scan"] == 16133]
scan_16134 = [x for i, x in enumerate(psm_objects) if psm_objects[i]["scan"] == 16134]

type(scan_16133[0]['peptide'])
type(scan_16134[0]['peptide'])
# scans seemes to be fine

# checking if error in the scan is occruing in the raw file parser from pyascore

spectra_file = r"C:\Users\mlorenz\Desktop\Internship\pyAscore " \
               r"test_miniTOPAS\mzML_ms2_short\03951_GA4_P043866_S00_I00_U04_R1_TMT11.mzML"

spectra_parser = pyascore.SpectraParser(str(spectra_file), "mzML")

# To a list of dictionaries:
spectra_objects = spectra_parser.to_list()
spectra_objects[0]

# spectrum from the deffective scan also seemes to be fine

# testing how to filter out MS3 scans using the pyascore spectraparser #################################################

spectra_file = r"C:\Users\mlorenz\Desktop\Internship\pyAscore " \
               r"test_miniTOPAS\mzML_ms2\03951_GA1_P043866_S00_I00_U01_R1_TMT11.mzML"

spectra_parser = pyascore.SpectraParser(str(spectra_file), "mzML", ms_level=1)

# To a list of dictionaries:
spectra_objects = spectra_parser.to_list()
spectra_objects[0]
spectra_objects_ms1 = [x for i, x in enumerate(spectra_objects) if spectra_objects[i]['ms_level'] == 1]

# playing around with pyascore spectrum parser im changing: ############################################################

spectra_file = r"C:\Users\mlorenz\Desktop\Internship\pyAscore " \
               r"test_miniTOPAS\mzML\03951_GA1_P043866_S00_I00_U01_R1_TMT11.mzML"
reader = mzml.MzML(spectra_file)

extractor = spec_parsers.MzMLExtractor()

reader.map(extractor.extract)

extractor.extract(next(reader))
# I changed the extractor.extract function so that it checks the ms level after reading it form the scan and fills in
# NaNs for the rest of the scan info before getting to the _get_precursor function which always threw the error
