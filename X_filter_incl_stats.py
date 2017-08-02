#!/usr/bin/python
import os
import sys
import logging
from optparse import OptionParser
from collections import defaultdict
#from utils import utils_logging
from scipy import stats
import numpy as np
import csv
import argparse # package to help with argument parsing
import time

# setup argument parser
parser = argparse.ArgumentParser(description="Script to filter vcf files.")
parser.add_argument('-i', '--input', type=str, default='', help='Name of input vcf file.')
parser.add_argument('-o', '--output', type=str, default='', help='Name of output vcf file.')
parser.add_argument('-m', '--meta-output', type=str, default='', help='Name of meta data output file.')
parser.add_argument('-gq', '--gq-threshold', type=int, default=20, help='GQ threshold to use (default=20).')
parser.add_argument('--fold-change-margin', type=float, default=0.2, help='Margin around 2.0 to use for fold change check (default=0.2).')
parser.add_argument('--log-file', type=str, default='', help='File to write log information to (uses stdout if none specified).')
parser.add_argument('-r', '--reverse', action='store_true', help='Interchange the labels on the males and females.')

# parse command line arguments
opts = parser.parse_args(sys.argv[1:])

# config values
individual_start_col = 9


# initialise counters
removed = 0
total = 0

# setup logging, this will log anything info level or above
logging.basicConfig(filename=(opts.log_file if opts.log_file != '' else None), filemode='a', level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

logging.info('Start of filtering.')


total_sample_read_depth = total_read_dp_per_individual(opts.input, individual_start_col, opts.gq_threshold)
#print(total_sample_read_depth)


# open files for reading
try:
    f = open(opts.input, "r")
    fw = open(opts.output, "w")
    fm = open(opts.meta_output, "w")
    logging.info('Opened input file %s' % opts.input)
    logging.info('Opened output file %s' % opts.output)
    logging.info('Opened output meta file %s' % opts.meta_output)
    csv_reader = csv.reader(f, delimiter="\t")
    csv_writer = csv.writer(fw, delimiter="\t")
    csv_meta_writer = csv.writer(fm, delimiter="\t")
    s = time.time()
    for row in csv_reader:
        if row[0].startswith("##"):
            continue
        if row[0].startswith("#"): # header column
            headers = row
            individuals = headers[individual_start_col:]
            male_cols, female_cols = find_genders(individuals, offset=individual_start_col, reverse=opts.reverse)
            csv_writer.writerow(row)
            csv_meta_writer.writerow(["locus", "position", "is_male_heterozygote", "n_male_homozygote", "n_male_heterozygote",
                                      "n_female_homozygote", "n_female_heterozygote", "n_gq_filtered", "male_mean_coverage",
                                      "female_mean_coverage", "fold_change", "fold_change_in_range", "t_stat_eq", "pvalue_eq"])
        else:
            total += 1
            gq_filtered = filter_by_gq(row, opts.gq_threshold, offset=individual_start_col)  # filter individuals where gq is less than given threshold                            
            n_hm_male, n_ht_male, n_hm_female, n_ht_female = count_zygote_gt_type(row, female_cols, male_cols)
            is_male_heterozygote = at_least_one_heterozygote(row, male_cols)
            male_mean_coverage, female_mean_coverage, fold_change = calc_coverage_and_fold_change(row, male_cols, female_cols, total_sample_read_depth, normalise=True)
            t_stat_eq, pvalue_eq = stats.ttest_ind(male_dps,female_dps)    
            if fold_change is None:
                fold_change_in_range = None
            else:
                if (2.0 - opts.fold_change_margin) < fold_change < (2.0 + opts.fold_change_margin):
                    fold_change_in_range = True
                else:
                    fold_change_in_range = False
                    
               
            if not is_male_heterozygote and pvalue_eq<.05 and fold_change_in_range:
                csv_writer.writerow(row)
            else:
                removed += 1
            csv_meta_writer.writerow([row[0], row[1], is_male_heterozygote, n_hm_male, n_ht_male, n_hm_female, n_ht_female, gq_filtered,
                                      male_mean_coverage, female_mean_coverage, fold_change, fold_change_in_range, t_stat_eq, pvalue_eq])
    e = time.time()
    logging.info('Filtered %d/%d records leaving %d in %.2f seconds.' % (removed, total, total - removed, e-s))
    f.close()
    fw.close()
    fm.close()
except IOError as ioerror:
    logging.error('Problem opening files: ' + str(ioerror))
