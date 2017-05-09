#!/usr/bin/python
import csv
import argparse # package to help with argument parsing
import sys
import os
import logging # module to enable logging
import numpy as np

# setup argument parser
parser = argparse.ArgumentParser(description="Script to filter vcf files.")
parser.add_argument('-i', '--input', type=str, default='', help='Name of input vcf file.')
parser.add_argument('-o', '--output', type=str, default='', help='Name of output vcf file.')
parser.add_argument('-m', '--meta-output', type=str, default='', help='Name of meta data output file.')
parser.add_argument('-gq', '--gq-threshold', type=int, default=20, help='GQ threshold to use (default=20).')
parser.add_argument('--fold-change-margin', type=float, default=0.2, help='Margin around 2.0 to use for fold change check (default=0.2).')

# parse command line arguments
opts = parser.parse_args(sys.argv[1:])

# config values
individual_start_col = 9

# define useful functions
def find_genders(x, offset):
    males = []
    females = []
    for i, ind in enumerate(x):
        if ind[2].lower() == 'm':
            males.append(i + offset)
        elif ind[2].lower() == 'f':
            females.append(i + offset)
    return males, females

def is_heterozygote(snp_info):
    first_part = snp_info[:snp_info.find(':')]
    if first_part == '0/1':
        return True
    elif first_part == '0/0' or first_part == '1/1':
        return False
    else:
        return None

def at_least_one_heterozygote(row, cols):
    for col in cols:
        if is_heterozygote(row[col]):
            return True
    return False

def is_gq_greater_than(snp_info, gq_threshold):
    parts = snp_info.split(':')
    try:
        gq_value = int(parts[-1]) # we expect GQ to be the last field
    except Exception as exception:
        logging.debug('GQ value not integer valued: %s of %s.' % (parts[-1], snp_info))
        return False
    return gq_value >= gq_threshold

def count_zygote_gt_type(row, male_cols, female_cols):
    """
    Count the number of male, female, homozygote and heterozygote genotypes.
    """
    n_hm_male, n_ht_male, n_hm_female, n_ht_female = 0, 0, 0, 0
    for i in male_cols:
        is_ht = is_heterozygote(row[i])
        if is_ht:
            n_ht_male += 1
        elif is_ht == False:
            n_hm_male += 1

    for i in female_cols:
        is_ht = is_heterozygote(row[i])
        if is_ht:
            n_ht_female += 1
        elif is_ht == False:
            n_hm_female += 1

    return n_hm_male, n_ht_male, n_hm_female, n_ht_female

def filter_by_gq(row, gq_threshold, offset, empty_str='.'):
    n = 0
    for i in range(offset, len(row)):
        if not is_gq_greater_than(row[i], gq_threshold):
            n += 1
            row[i] = empty_str
    return n

def dp_values(row, cols, dp_idx=2):
    # get a list of the dp values for the given rows
    dp_values = []
    for col in cols:
        parts = row[col].split(':')
        if len(parts) == 8:
            try:
                dp_value = int(parts[dp_idx])
                dp_values.append(dp_value)
            except Exception as ex:
                logging.warning('DP value not an integer %s.' % parts[dp_idx])
    return dp_values

def calc_coverage_and_fold_change(row, male_cols, female_cols, normalise=True):
    male_dps = np.array(dp_values(row, male_cols), np.float)
    female_dps = np.array(dp_values(row, female_cols), np.float)
    if (len(male_dps) == 0) or (len(female_dps) == 0):
        return None, None, None
    total_dp = np.sum(male_dps) + np.sum(female_dps)
    if total_dp == 0:
        logging.error('Total coverage depth is 0.')
        return None, None, None
    # normalise the depts
    if normalise:
        male_dps /= total_dp
        female_dps /= total_dp
    male_mean_coverage = np.mean(male_dps)
    female_mean_coverage = np.mean(female_dps)
    fold_change = female_mean_coverage/male_mean_coverage
    return male_mean_coverage, female_mean_coverage, fold_change


# initialise counters
removed = 0
total = 0

# setup logging, this will log anything info level or above
logging.basicConfig(filename='X_filtering.log', filemode='a', level=logging.INFO, format='%(asctime)s %(name)s: %(message)s')

logging.info('Starting of filtering.')

# open files for reading
with open(opts.input, "r") as f, open(opts.output, "w") as fw, open(opts.meta_output, "w") as fm:
    logging.info('Opened input file %s' % opts.input)
    logging.info('Opened output file %s' % opts.output)
    logging.info('Opened output meta file %s' % opts.meta_output)
    csv_reader = csv.reader(f, delimiter="\t")
    csv_writer = csv.writer(fw, delimiter="\t")
    csv_meta_writer = csv.writer(fm, delimiter="\t")
    for row in csv_reader:
        if row[0].startswith("##"):
            continue
        if row[0].startswith("#"): # header column
            headers = row
            individuals = headers[individual_start_col:]
            male_cols, female_cols = find_genders(individuals, offset=individual_start_col)
            csv_writer.writerow(row)
            csv_meta_writer.writerow(["locus", "position", "is_male_heterozygote", "n_male_homozygote", "n_male_heterozygote",
                                      "n_female_homozygote", "n_female_heterozygote", "n_gq_filtered", "male_mean_coverage",
                                      "female_mean_coverage", "fold_change", "fold_change_in_range"])
        else:
            total += 1
            gq_filtered = filter_by_gq(row, opts.gq_threshold, offset=individual_start_col)  # filter individuals where gq is less than given threshold
            n_hm_male, n_ht_male, n_hm_female, n_ht_female = count_zygote_gt_type(row, male_cols, female_cols)
            is_male_heterozygote = at_least_one_heterozygote(row, male_cols)
            male_mean_coverage, female_mean_coverage, fold_change = calc_coverage_and_fold_change(row, male_cols, female_cols, normalise=True)
            if fold_change is None:
                fold_change_in_range = None
            else:
                if (2.0 - opts.fold_change_margin) < fold_change < (2.0 + opts.fold_change_margin):
                    fold_change_in_range = True
                else:
                    fold_change_in_range = False
            if not is_male_heterozygote and fold_change_in_range:
                csv_writer.writerow(row)
            else:
                removed += 1
            csv_meta_writer.writerow([row[0], row[1], is_male_heterozygote, n_hm_male, n_ht_male, n_hm_female, n_ht_female, gq_filtered,
                                      male_mean_coverage, female_mean_coverage, fold_change, fold_change_in_range])
