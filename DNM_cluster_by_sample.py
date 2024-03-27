
# takes as input a multisample vcf file of de novo mutations and identifies clusters
# usage: python DNM_cluster_by_sample.py -v <vcf_file> -o <output_folder> -c <max_distance_to_call_a_cluster_in_a_sample>

import argparse
import os

parser = argparse.ArgumentParser(description='DNM clustering')
parser.add_argument('-v', '--vcf', help='vcf file', required=True)
parser.add_argument('-o', '--output', help='output folder', required=True)
parser.add_argument('-c', '--cluster_max_size_between_mut', help='max distance to call a cluster in a sample', required=True)
args = parser.parse_args()

vcf = args.vcf
output = args.output
cluster_max_distance = int(args.cluster_max_size_between_mut)

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
if not vcf or not output or not cluster_max_distance:
    parser.print_help()
    exit()

# read vcf file
with open(vcf, 'r') as f:
    lines = f.readlines()

# get samples
samples = get_samples(lines)

# create a dictionary with samples as keys and a list of mutations as values
variants_by_sample = get_variants_by_sample(lines, samples)

# get chromosomes
chromosomes = get_chromosomes(lines)

class Cluster:
    def __init__(self, sample, chr, start, end, variants):
        self.sample = sample
        self.chr = chr
        self.start = start
        self.end = end
        self.variants = variants

    def get_var_count(self):
        return len(self.variants)
    
    def get_size(self):
        return self.end - self.start


# create a dictionary with samples as keys and a list of clusters as values
clusters = []

# print clusters by sample
for sample in samples:
    for chr in chromosomes:
        variants_to_analyse = [variant for variant in variants_by_sample[sample] if variant['chr'] == chr]
        # reorder variants by position
        variants_to_analyse = sorted(variants_to_analyse, key=lambda x: int(x['pos']))
        # on ajoute none en distance to previous et none en distance to next pour chaque variant
        for i, variant in enumerate(variants_to_analyse):
            variant['distance_to_previous'] = None
            variant['distance_to_next'] = None

        for i, variant in enumerate(variants_to_analyse):
            if i > 0:
                variant['distance_to_previous'] = int(variant['pos']) - int(variants_to_analyse[i - 1]['pos'])
            if i < len(variants_to_analyse) - 1:
                variant['distance_to_next'] = int(variants_to_analyse[i + 1]['pos']) - int(variant['pos'])

        # on affiche les variants pour lesquels la distance est inférieure à la distance maximale
        cluster_ongoing = False
        for i, variant in enumerate(variants_to_analyse):
            if variant['distance_to_previous'] and variant['distance_to_previous'] <= cluster_max_distance:
                # deux possiblités : soit on commence un cluster, soit on ajoute le variant au cluster en cours
                if cluster_ongoing:
                    # on ajoute le variant au cluster
                    end = int(variant['pos'])
                    clusters[-1].end = end
                    clusters[-1].variants.append(variant)
                else:
                    cluster_ongoing = True
                    # on génère un cluster avec comme start la position du variant précédent
                    start = int(variants_to_analyse[i - 1]['pos'])
                    end = int(variant['pos'])
                    variants = [variants_to_analyse[i - 1], variant]
                    clusters.append(Cluster(sample, chr, start, end, variants))
            else:
                cluster_ongoing = False


# create output folders for invividual cluster files
if not os.path.exists(output):
    os.makedirs(output)

if not os.path.exists(output + '/clusters_by_sample'):
    os.makedirs(output + '/clusters_by_sample')

# generate indivdual_files_by_clusters
for cluster in clusters:
    cluster_filename = output + '/clusters_by_sample/' + cluster.sample + '_' + cluster.chr + '_' + str(cluster.start) + '_' + str(cluster.end) + "_" + str(cluster.get_var_count()) + "variants_" + str(cluster.get_size()) + 'bp.tsv'
    with open(cluster_filename, 'w') as f:
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + cluster.sample + '\n')
        for variant in cluster.variants:
            f.write(variant['chr'] + '\t' + variant['pos'] + '\t' + variant['ID'] + '\t' + variant['ref'] + '\t' + variant['alt'] + '\t.\t.\t.\t.\t.\n')

# create a summary file
summary_filename = output + '/clusters_summary.txt'
with open(summary_filename, 'w') as f:
    f.write('sample\tchr\tstart\tend\tsize\tvar_count\n')
    for cluster in clusters:
        f.write(cluster.sample + '\t' + cluster.chr + '\t' + str(cluster.start) + '\t' + str(cluster.end) + '\t' + str(cluster.get_size()) + '\t' + str(cluster.get_var_count()) + '\n')

print('Done!')
