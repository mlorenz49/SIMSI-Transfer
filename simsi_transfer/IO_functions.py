import os
import logging
import glob
import multiprocessing
from pathlib import Path
import argparse

import pandas as pd

from .processing_functions import purge_mrc_files


logger = logging.getLogger(__name__)


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def export_summary_file(summary_file, mainpath, pval, state):
    """
    Generate and return summary dataframes, a merge of msms.txt identifications and MaRaCluster clustering results with
    one line for each ms2 scan generated during acquisition
    :param summary_file: summary of msms.txt, msmsscans.txt, and MaRaCluster cluster.csv files
    :param mainpath: Path to main processing folder
    :param pval: MaRaCluster stringency threshold; a higher value (e.g. p30) is more stringent than a lower one (e.g.
    p15), but results in less transfers
    :param state: set to 'merged', 'transferred', or 'filtered' for merged files, msmsScans, and msms files respectively
    """
    sumpath = mainpath / Path('summaries')
    if not os.path.exists(sumpath):
        os.makedirs(sumpath)
    pval_path = sumpath / Path(pval)
    if not os.path.exists(pval_path):
        os.makedirs(pval_path)
    if state == 'merged':
        summary_file.to_csv(pval_path / Path(f'{pval}_summary.csv'), sep='\t', index=False, na_rep='NaN')
    elif state == 'transferred':
        summary_file = summary_file.sort_values(by=['Raw file', 'scanID'])
        summary_file.to_csv(pval_path / Path(f'{pval}_msmsScans.csv'), sep='\t', index=False, na_rep='NaN')
    elif state == 'filtered':
        summary_file = summary_file.sort_values(by=['Sequence', 'Modified sequence'])
        summary_file.to_csv(pval_path / Path(f'{pval}_msms.csv'), sep='\t', index=False, na_rep='NaN')
    else:
        raise ValueError("Incorrect state given, please set state parameter to 'merged', 'transferred', or 'filtered'")


def export_simsi_evidence_file(evidence_file, mainpath, pval):
    sumpath = mainpath / Path('summaries')
    pval_path = sumpath / Path(pval)
    evidence_file.to_csv(pval_path / Path(f'{pval}_evidence.csv'), sep='\t', index=False, na_rep='NaN')


def open_msmsscans_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated msmsscans.txt dataframe
    """
    cols = ['Raw file', 'Scan number', 'm/z', 'Mass', 'Retention time', 'Precursor full scan number',
            'Reporter intensity 1', 'Reporter intensity 2',
            'Reporter intensity 3', 'Reporter intensity 4', 'Reporter intensity 5', 'Reporter intensity 6',
            'Reporter intensity 7', 'Reporter intensity 8', 'Reporter intensity 9', 'Reporter intensity 10',
            'Reporter intensity 11', 'Reporter intensity corrected 1', 'Reporter intensity corrected 2',
            'Reporter intensity corrected 3', 'Reporter intensity corrected 4', 'Reporter intensity corrected 5',
            'Reporter intensity corrected 6', 'Reporter intensity corrected 7', 'Reporter intensity corrected 8',
            'Reporter intensity corrected 9', 'Reporter intensity corrected 10', 'Reporter intensity corrected 11']
    try:
        msmsscans = pd.read_csv(mainpath / Path('msmsScans.txt'), sep='\t', usecols=cols).rename(
            columns={'Scan number': 'scanID'})
        tmt = 11
    except ValueError:
        cols.remove('Reporter intensity 11')
        cols.remove('Reporter intensity corrected 11')
        msmsscans = pd.read_csv(mainpath / Path('msmsScans.txt'), sep='\t', usecols=cols).rename(
            columns={'Scan number': 'scanID'})
        tmt = 10
    return msmsscans, tmt


def open_msms_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated msms.txt dataframe
    """
    msmstxt = pd.read_csv(mainpath / Path('msms.txt'), sep='\t',
                          usecols=['Raw file', 'Scan number', 'Sequence', 'Modified sequence', 'Length',
                                   'Modifications', 'Missed cleavages', 'Proteins', 'Gene Names', 'Protein Names',
                                   'Charge', 'Mass error [ppm]', 'PIF', 'Precursor Intensity', 'PEP', 'Score',
                                   'Delta score', 'Reverse']).rename(
        columns={'Scan number': 'scanID'})
    return msmstxt


def open_evidence_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    evidence = pd.read_csv(mainpath / Path('evidence.txt'), sep='\t',
                           usecols=['Sequence', 'Modified sequence', 'Leading proteins', 'Raw file', 'Charge',
                                    'Calibrated retention time', 'Retention time', 'Retention length',
                                    'Calibrated retention time start', 'Calibrated retention time finish',
                                    'Retention time calibration', 'Type', 'Intensity', 'Reverse'])
    return evidence


def open_summary_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated msmsscans.txt dataframe
    """
    try:
        return pd.read_csv(mainpath / Path('summary.txt'), sep='\t',
                           usecols=['Raw file', 'Experiment', 'Fraction'])
    except ValueError:
        temp = pd.read_csv(mainpath / Path('summary.txt'), sep='\t',
                           usecols=['Raw file', 'Experiment'])
        temp['Fraction'] = 1
        return temp


def open_maracluster_clusters(mainpath, pval):
    maracluster_df = pd.read_csv(mainpath / Path(f'MaRaCluster.clusters_{pval}.tsv'),
                                 sep='\t', names=['Raw file', 'scanID', 'clusterID'])
    maracluster_df['Raw file'] = maracluster_df['Raw file'].apply(purge_mrc_files)
    return maracluster_df


def parse_args(argv):
    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_txt_folder', default=None, metavar="DIR", required=True,
                       help='''Path to MaxQuant combined/txt output folder.''')

    apars.add_argument('--raw_folder', default=None, metavar="DIR", required=True,
                       help='''Full path to folder containing .raw or .mzML files.''')

    apars.add_argument('--stringencies', default="20,15,10", metavar="S",
                       help='''Clustering thresholds at which to produce cluster files, listed as comma separated list. 
                               The higher the stringency value, the more strict the clustering is.
                          ''')

    apars.add_argument('--output_folder', default="./simsi_output", metavar="DIR",
                       help='''Full path to desired SIMSI output folder.''')

    apars.add_argument('--num_threads', type=int, default=multiprocessing.cpu_count(), metavar='N',
                       help='''Number of threads, by default this is equal to the number of CPU cores available on the device.
                          ''')

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return Path(args.mq_txt_folder), Path(args.raw_folder), parse_stringencies(args.stringencies), Path(
        args.output_folder), args.num_threads


def parse_stringencies(stringencies):
    if stringencies == '':
        stringencies = [20, 15, 10]
    else:
        try:
            stringencies = [int(i) for i in stringencies.split(',')]
        except ValueError:
            logger.error('This is not a stringency list. Please input the list in the following format: 30,25,20,15,10,5')

    return stringencies


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
