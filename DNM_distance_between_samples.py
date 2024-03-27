
# takes as input a multisample vcf file of de novo mutations and identifies variants at close proximity between samples
# usage: python DNM_distance_between_samples.py -v <vcf_file> -o <output_folder> 

import argparse
import os

parser = argparse.ArgumentParser(description='DNM clustering')
parser.add_argument('-v', '--vcf', help='vcf file', required=True)
parser.add_argument('-o', '--output', help='output folder', required=True)
args = parser.parse_args()

vcf = args.vcf
output = args.output

def get_samples(lines):
    samples = []
    for line in lines:
        if line.startswith('#CHROM'):
            samples = line.strip().split('\t')[9:]
    return samples

def get_variants_by_sample(lines, samples):
    # create a dictionary with samples as keys and a list of mutations as values, in the format of a vcf line
    mutations = {}
    for sample in samples:
        mutations[sample] = []
    for line in lines:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            for i, sample in enumerate(samples):
                genotype = line[9 + i].split(':')[0]
                if genotype not in ['0/0', '0|0', './.', '.|.']:
                    # on append chr, pos, ref, alt, ID, variant type sous la forme d'un dictionnaire
                    dict = {'chr': line[0], 'pos': line[1], 'ref': line[3], 'alt': line[4], 'ID': line[2], 'type': "SNV" if len(line[3]) == 1 and len(line[4]) == 1 else "INDEL"}
                    mutations[sample].append(dict)
    return mutations

def get_chromosomes(lines):
    chromosomes = []
    for line in lines:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            if line[0] not in chromosomes:
                chromosomes.append(line[0])
    return chromosomes

# if arguments not provided, print help
if not vcf or not output:
    parser.print_help()
    exit()

# read vcf file
with open(vcf, 'r') as f:
    lines = f.readlines()

# get samples
samples = get_samples(lines)

# get chromosomes
chromosomes = get_chromosomes(lines)

# create output folder
if not os.path.exists(output):
    os.makedirs(output)

# create a dictionary with samples as keys and a list of mutations as values
variants_by_sample = get_variants_by_sample(lines, samples)

# append dist_to_closest_variant_in_samples to each variant
for sample in samples:
    for variant in variants_by_sample[sample]:
        variant['dist_to_closest_variant_in_samples'] = -1
        variant['oter_sample'] = ''
        variant['other_variant'] = ''
        for other_sample in samples:
            if other_sample != sample:
                for other_variant in variants_by_sample[other_sample]:
                    if other_variant['chr'] == variant['chr']:
                        dist = abs(int(variant['pos']) - int(other_variant['pos']))
                        if variant['dist_to_closest_variant_in_samples'] == -1 or dist < variant['dist_to_closest_variant_in_samples']:
                            variant['dist_to_closest_variant_in_samples'] = dist
                            variant['other_sample'] = other_sample
                            variant['other_variant'] = other_variant['ID']


# save results in output folder
output_file = output + "/" + vcf.split('/')[-1].replace('.vcf', '_DNM_distance_between_samples.txt')
with open(output_file, 'w') as f:
    # header
    f.write('sample\tchr\tpos\tref\talt\tID\ttype\tdist_to_closest_variant_in_samples\tother_sample\tother_variant\n')
    for sample in samples:
        for variant in variants_by_sample[sample]:
            f.write(sample + '\t' + variant['chr'] + '\t' + variant['pos'] + '\t' + variant['ref'] + '\t' + variant['alt'] + '\t' + variant['ID'] + '\t' + variant['type'] + '\t' + str(variant['dist_to_closest_variant_in_samples']) + '\t' + variant['other_sample'] + '\t' + variant['other_variant'] + '\n')
print('Results saved in ' + output_file)
