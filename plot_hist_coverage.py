#!/usr/bin/python
import csv
import argparse # package to help with argument parsing
import sys
import os
import logging # module to enable logging
import numpy as np
import time

# setup argument parser
parser = argparse.ArgumentParser(description="Script to filter vcf files.")
parser.add_argument('-i', '--input', type=str, default='', help='Name of input vcf file.')
parser.add_argument('-o', '--output-id', type=str, default='', help='Identifying string to put in output filenames.')
parser.add_argument('-gq', '--gq-threshold', type=int, default=20, help='GQ threshold to use (default=20).')
parser.add_argument('--fold-change-margin', type=float, default=0.2, help='Margin around 2.0 to use for fold change check (default=0.2).')
parser.add_argument('--log-file', type=str, default='', help='File to write log information to (uses stdout if none specified).')
parser.add_argument('-r', '--reverse', action='store_true', help='Interchange the labels on the males and females.')

# parse command line arguments
opts = parser.parse_args(sys.argv[1:])

# config values
individual_start_col = 9

# define useful functions
def find_genders(x, offset, reverse=False):
    males = []
    females = []
    for i, ind in enumerate(x):
        if ind[2].lower() == 'm':
            males.append(i + offset)
        elif ind[2].lower() == 'f':
            females.append(i + offset)
    if reverse:
        return females, males
    else:
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
logging.basicConfig(filename=(opts.log_file if opts.log_file != '' else None), filemode='a', level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

logging.info('Start of plotting.')

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

except ImportError as ex:
    logging.error('Problem importing plotting library.' + str(ex))
    sys.exit(-1)

# keep fold change for each group
fold_changes = {}
fold_changes['all'] = []
fold_changes['filtered'] = []
fold_changes['excluded'] = []

coverages = {}
coverages['all'] = {'male': [], 'female': []}
coverages['filtered'] = {'male': [], 'female': []}
coverages['excluded'] = {'male': [], 'female': []}

fold_change_issue = 0

# open files for reading
try:
    f = open(opts.input, "r")
    logging.info('Opened input file %s' % opts.input)
    csv_reader = csv.reader(f, delimiter="\t")
    s = time.time()
    for row in csv_reader:
        if row[0].startswith("##"):
            continue
        if row[0].startswith("#"): # header column
            headers = row
            individuals = headers[individual_start_col:]
            male_cols, female_cols = find_genders(individuals, offset=individual_start_col, reverse=opts.reverse)
        else:
            total += 1
            gq_filtered = filter_by_gq(row, opts.gq_threshold, offset=individual_start_col)  # filter individuals where gq is less than given threshold
            n_hm_male, n_ht_male, n_hm_female, n_ht_female = count_zygote_gt_type(row, male_cols, female_cols)
            is_male_heterozygote = at_least_one_heterozygote(row, male_cols)
            male_mean_coverage, female_mean_coverage, fold_change = calc_coverage_and_fold_change(row, male_cols, female_cols, normalise=True)
            if fold_change is None:
                fold_change_in_range = None
                fold_change_issue += 1 
            else:
                fold_changes['all'].append(fold_change)
                coverages['all']['male'].append(male_mean_coverage)
                coverages['all']['female'].append(female_mean_coverage)
                if (2.0 - opts.fold_change_margin) < fold_change < (2.0 + opts.fold_change_margin):
                    fold_change_in_range = True
                    fold_changes['filtered'].append(fold_change)
                    coverages['filtered']['male'].append(male_mean_coverage)
                    coverages['filtered']['female'].append(female_mean_coverage)
                else:
                    fold_change_in_range = False
                    fold_changes['excluded'].append(fold_change)
                    coverages['excluded']['male'].append(male_mean_coverage)
                    coverages['excluded']['female'].append(female_mean_coverage)

            if not is_male_heterozygote and fold_change_in_range:
                pass
            else:
                removed += 1
    e = time.time()
    logging.info('Filtered %d/%d records leaving %d in %.2f seconds.' % (removed, total, total - removed, e-s))
    f.close()
    
    # now plot the distributions
    for i, key in enumerate(fold_changes.keys()):
        plt.subplot(3,1,i+1)
        plt.hist(fold_changes[key], range=(0.0, 5.0), bins=20)
        plt.xlabel('fold change')
        plt.title('Dist of fold change for %s'%key)
    output_id = '_id_%s_' % (opts.output_id)
    plt.savefig('%s%sfold_change_dist.png' % (os.path.basename(opts.input), output_id))

    plt.figure()
    for i, key in enumerate(coverages.keys()):
        plt.subplot(3,1,i+1)
        plt.hist(coverages[key]['male'], bins=20, label='male', alpha=0.8)
        plt.hist(coverages[key]['female'], bins=20, label='female', alpha=0.8)
        plt.legend(loc='best')
        plt.xlabel('coverage')
        plt.title('Dist of fold change for %s'%key)
    plt.savefig('%s%scoverage_dist.png' % (os.path.basename(opts.input), output_id))


except IOError as ioerror:
    logging.error('Problem opening files: ' + str(ioerror))
