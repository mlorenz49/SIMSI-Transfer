import logging

from . import parsing
from . import pyA_test
from . import join_outputs
from . import rescoring

logger = logging.getLogger(__name__)


def main(raw_folders, simsi_msms_folders, output_pars,  mode, output_ascore, results):

    logger.info('entering ascoring package')

    ####################################################################################################################
    # actual code starts here

    parsing.run_msms_parser(raw_folders, simsi_msms_folders, output_pars, mode)

    pyA_test.run_ascore(raw_folders, output_pars, output_ascore, mode)

    join_outputs.join_outputs(simsi_msms_folders, output_ascore, results)

    rescoring.rescore_seq(results, results)

    logger.info('Exiting ascoring package')


# if __name__ == "__main__":
#     main(r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\mzML",
#          r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\summaries",
#          r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\psm_minimal_simsi_MD_3",
#          r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\Acores_simsi_MD_3",
#          "simsi",
#          r"C:\Users\mlorenz\Desktop\Internship\pyAscore test\results_MD_3")
