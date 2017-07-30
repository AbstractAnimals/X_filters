# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 14:43:31 2017

@author: manager
"""
import logging
import numpy as np
import time
import csv
import sys
import os


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