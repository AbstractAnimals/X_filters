# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 15:08:00 2017

@author: manager
"""

#!/usr/bin/python
import csv
import argparse # package to help with argument parsing
import sys
import os
import logging # module to enable logging
import numpy as np
import time

# import custom functions
from X_filtering_functions import *


# only run following code if was called from command line
if __name__ == '__main__':
    # setup argument parser
    parser = argparse.ArgumentParser(description="Script to filter vcf files.")
    parser.add_argument('-m', '--meta-input', type=str, default='', help='Name of input meta file.')
    parser.add_argument('-i', '--input', type=str, default='', help='Name of input subset vcf file.')
    parser.add_argument('-o', '--output', type=str, default='', help='Name of output image.')
    parser.add_argument('--log-file', type=str, default='', help='File to write log information to (uses stdout if none specified).')
  
    # parse command line arguments
    opts = parser.parse_args(sys.argv[1:])
    
    
     # setup logging, this will log anything info level or above
    logging.basicConfig(filename=(opts.log_file if opts.log_file != '' else None), filemode='a', level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    logging.info('Start of plotting.')
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

    except ImportError as ex:
        logging.error('Problem importing plotting library.' + str(ex))
        sys.exit(-1)
    
    # getting SNP ids of interest from vcf file
    SNP_IDs = get_SNP_IDs_from_VCF(opts.input)
    # get the (normalised) coverage of SNPs of interest
    male_coverages, female_coverages = get_coverages_from_meta(opts.meta_input, SNP_IDs)
    
    # plot the coverages found 
    plt.figure()
    plt.hist(male_coverages, bins=100, label='male', alpha=0.5)
    plt.hist(female_coverages, bins=100, label='female', alpha=0.5)
    plt.legend(loc='best')
    plt.xlabel('coverage')
    plt.title('Dist of coverage for %s'%key)
    plt.savefig('%s_coverage_dist.png' % (opts.output))
    
    