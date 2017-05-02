#!/usr/bin/python
import csv

#definte filenames
#filename = "/storage/home/users/pm74/SexLinkedLoci/Manuscript/VCF_files/T.oc_assmb_T.oc_georder_0.8_allo_sym_no_T.marini_rm_PLf1_PLm12_PLf9_DVf9.vcf.recode.vcf"
folder = '/home/manager/pop_genomics/2nd_dataset_newer_samtools/Sexlinked_loci_Judith_recall'
filename = folder + "/T.OC_assmb_Sex_linked_candidates_MAF_0.05_maxmiss_0.5_rm_indv_greater_than0.2subset_T.OC_allo_sym.vcf.recode.vcf"
output_vcf_file=folder + "/T.oc_filtered_X.vcf"
output_meta_file=folder + "/T.oc_meta_filter_file.txt"

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
    
# initialise counters
removed = 0
total = 0

# open files for reading
with open(filename, "r") as f, open(output_vcf_file, "w") as fw, open(output_meta_file, "w") as fm:
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
            csv_meta_writer.writerow(["locus", "position", "is_male_heterozygote"])
        else:
            total += 1
            is_male_heterozygote = at_least_one_heterozygote(row, male_cols)
            if not is_male_heterozygote:
                csv_writer.writerow(row)
            else:
                removed += 1
            csv_meta_writer.writerow([row[0], row[1], is_male_heterozygote])
