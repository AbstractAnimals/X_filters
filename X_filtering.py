#!/usr/bin/python
import csv
import argparse # package to help with argument parsing
import sys
import os
import logging # module to enable logging
import numpy as np
import time

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
        if len(parts) > dp_idx:
            try:
                dp_value = int(parts[dp_idx])
                dp_values.append(dp_value)
            except Exception as ex:
                logging.warning('DP value not an integer %s.' % parts[dp_idx])
        else: dp_values.append(0)
    return dp_values

def total_read_dp_per_individual(input_file, individual_start_col, gq_threshold):
    f = open(input_file, "r")
    logging.info('Opened input file %s' % input_file)
    csv_reader = csv.reader(f, delimiter="\t")
    s = time.time()
    all_samples_total_coverages={}
    for row in csv_reader:
        if row[0].startswith("##"):
            continue
        if row[0].startswith("#"): # header column
            headers = row
            individuals = headers[individual_start_col:]
            for i, individual in enumerate(individuals):
                all_samples_total_coverages[i + individual_start_col]=0
        else:
            gq_filtered = filter_by_gq(row, gq_threshold, offset=individual_start_col)
            for col_idx in range(individual_start_col, len(row)):
                individual_dp = dp_values(row, [col_idx])[0]
                all_samples_total_coverages[col_idx] += individual_dp
    return all_samples_total_coverages

def calc_coverage_and_fold_change(row, male_cols, female_cols, total_sample_read_depth, normalise=True):
    male_dps = np.array(dp_values(row, male_cols), np.float)
    female_dps = np.array(dp_values(row, female_cols), np.float)
    male_dp_total = np.array(list(map(lambda x: total_sample_read_depth[x], male_cols)))
    female_dp_total = np.array(list(map(lambda x: total_sample_read_depth[x], female_cols)))
    if np.any(male_dp_total == 0) or np.any(female_dp_total == 0):
        logging.error('Total coverage depth is 0.')
        return None, None, None

    # normalise the depths
    if normalise:
        male_dps = male_dps/male_dp_total * 1000000
        female_dps = female_dps/female_dp_total * 1000000

    male_mean_coverage = np.mean(male_dps)
    female_mean_coverage = np.mean(female_dps)
    fold_change = female_mean_coverage/male_mean_coverage
    return male_mean_coverage, female_mean_coverage, fold_change

# only run following code if was called from command line
if __name__ == '__main__':
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
                                        "female_mean_coverage", "fold_change", "fold_change_in_range"])
            else:
                total += 1
                gq_filtered = filter_by_gq(row, opts.gq_threshold, offset=individual_start_col)  # filter individuals where gq is less than given threshold
                n_hm_male, n_ht_male, n_hm_female, n_ht_female = count_zygote_gt_type(row, female_cols, male_cols)
                is_male_heterozygote = at_least_one_heterozygote(row, male_cols)
                male_mean_coverage, female_mean_coverage, fold_change = calc_coverage_and_fold_change(row, male_cols, female_cols, total_sample_read_depth, normalise=True)
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
        e = time.time()
        logging.info('Filtered %d/%d records leaving %d in %.2f seconds.' % (removed, total, total - removed, e-s))
        f.close()
        fw.close()
        fm.close()
    except IOError as ioerror:
        logging.error('Problem opening files: ' + str(ioerror))
