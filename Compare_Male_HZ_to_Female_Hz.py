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
parser.add_argument('-o', '--output', type=str, default='', help='Output prefix to use.')
parser.add_argument('-gq', '--gq-threshold', type=int, default=20, help='GQ threshold to use (default=20).')
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

def is_gq_greater_than(snp_info, gq_threshold):
    parts = snp_info.split(':')
    try:
        gq_value = int(parts[-1]) # we expect GQ to be the last field
    except Exception as exception:
        logging.debug('GQ value not integer valued: %s of %s.' % (parts[-1], snp_info))
        return False
    return gq_value >= gq_threshold

def count_zygote_gt_by_gender(row, male_cols, female_cols):
    """
    Count the number of male, female, homozygote and heterozygote genotypes.
    """
    n_male, n_hz_male, n_female, n_hz_female = 0, 0, 0, 0
    for i in male_cols:
        is_hz = is_heterozygote(row[i])
        if is_hz is not None:
            n_male += 1
            if is_hz:
                n_hz_male += 1

    for i in female_cols:
        is_hz = is_heterozygote(row[i])
        if is_hz is not None:
            n_female += 1
            if is_hz:
                n_hz_female += 1

    if n_male == 0 or n_female == 0:
        return None, None
    else:
        return  n_hz_male/float(n_male), n_hz_female/float(n_female)

def filter_by_gq(row, gq_threshold, offset, empty_str='.'):
    n = 0
    for i in range(offset, len(row)):
        if not is_gq_greater_than(row[i], gq_threshold):
            n += 1
            row[i] = empty_str
    return n

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


# open files for reading
try:
    f = open(opts.input, "r")
    logging.info('Opened input file %s' % opts.input)
    csv_reader = csv.reader(f, delimiter="\t")
    male_hzs, female_hzs = [], []
    consensus, position = [], []
    s = time.time()
    for row in csv_reader:
        if row[0].startswith("##"):
            continue
        if row[0].startswith("#"): # header column
            headers = row
            individuals = headers[individual_start_col:]
            male_cols, female_cols = find_genders(individuals, offset=individual_start_col, reverse=opts.reverse)
        else:
            gq_filtered = filter_by_gq(row, opts.gq_threshold, offset=individual_start_col)  # filter individuals where gq is less than given threshold
            male_hz, female_hz = count_zygote_gt_by_gender(row, male_cols, female_cols)
            if male_hz is not None and female_hz is not None:
                male_hzs.append(male_hz)
                female_hzs.append(female_hz)
                consensus.append(row[0])
                position.append(row[1])

    e = time.time()
    logging.info('Found values for %d snps.' % (len(male_hzs)))
    f.close()

    plt.plot(female_hzs, male_hzs, marker='o', ls='none')
    plt.xlabel('female')
    plt.ylabel('male')
    plt.title('Heterzygote proportion')
    plt.savefig('%s_heterozygote_scatter.png' % (opts.output))

    with open('%s_heterozygote_values.csv' % (opts.output), 'w') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(['consensus', 'position', 'male_hz', 'female_hz'])
        for i in range(len(male_hzs)):
            csv_writer.writerow([consensus[i], position[i], male_hzs[i], female_hzs[i]])

except IOError as ioerror:
    logging.error('Problem opening files: ' + str(ioerror))
