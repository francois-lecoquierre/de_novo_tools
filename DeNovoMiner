
import os
import sys
import argparse
import yaml
import subprocess
import pandas as pd


# args management
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="vcf file", required=True)
parser.add_argument("-p", "--ped", help="ped file", required=True)
parser.add_argument("-c", "--config", help="config file", required=True)
parser.add_argument("-o", "--output_dir", help="output directory", required=True)
args = parser.parse_args()

if not os.path.isfile(args.vcf):
    print("vcf file does not exist")
    parser.print_help()
    sys.exit(1)
if not os.path.isfile(args.ped):
    print("ped file does not exist")
    parser.print_help()
    sys.exit(1)
if not os.path.isfile(args.config):
    print("config file does not exist")
    parser.print_help()
    sys.exit(1)
if not args.output_dir:
    print("output dir is missing")
    parser.print_help()
    sys.exit(1)


class Trio:
    def __init__(self, family_id, child_id, father_id, mother_id, multi_vcf, list_of_individuals_from_same_family_in_ped_file, tmp_dir, output_dir, params_yaml):
        self.family_id = family_id
        self.child_id = child_id
        self.father_id = father_id
        self.mother_id = mother_id
        self.child_index_in_vcf = None
        self.father_index_in_vcf = None
        self.mother_index_in_vcf = None
        self.trio_is_ok_for_analysis = False
        self.trio_problems = []
        self.multi_vcf = multi_vcf
        self.list_of_individuals_from_same_family_in_ped_file = list_of_individuals_from_same_family_in_ped_file
        self.params_yaml = params_yaml
        self.params_dict = self.yaml_to_dict()
        self.tmp_dir = tmp_dir
        self.output_dir = output_dir
        self.dnv_count = None
        self.snv_count = None
        self.indel_count = None
        self.check_trio_is_ok_for_analysis()
        self.get_trio_index_in_multi_vcf()
        self.validate_parameters()
        self.call_de_novo()
        # self.print_trio_parameters()

    def get_indiv_from_vcf(self):
        # returns the list of individuals in the multi vcf
        result = subprocess.run(['bcftools', 'query', '-l', self.multi_vcf], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Erreur lors de l'exécution de bcftools: {result.stderr}")
            return None
        samples = result.stdout.split('\n')
        samples = [sample for sample in samples if sample]
        return samples

    def get_trio_index_in_multi_vcf(self):
        # returns the index of each individual to be used by bcftools
        if self.trio_is_ok_for_analysis:
            indiv = self.get_indiv_from_vcf()
            self.child_index_in_vcf = indiv.index(self.child_id)
            self.father_index_in_vcf = indiv.index(self.father_id)
            self.mother_index_in_vcf = indiv.index(self.mother_id)
        else:
            return None
        

    def check_trio_is_ok_for_analysis(self):
        # verify that the trio is defined in the ped file
        if not self.father_id:
            self.trio_is_ok_for_analysis = False
            self.trio_problems.append("father not defined in ped file")
        if not self.mother_id:
            self.trio_is_ok_for_analysis = False
            self.trio_problems.append("mother not defined in ped file")
        if not self.child_id:
            self.trio_is_ok_for_analysis = False
            self.trio_problems.append("child not defined in ped file")
        
        # verify that the trio is in the multi vcf
        indiv = self.get_indiv_from_vcf()
        if self.child_id not in indiv:
            self.trio_is_ok_for_analysis = False
            self.trio_problems.append("child is not in vcf")
        if self.father_id not in indiv:
            self.trio_is_ok_for_analysis = False
            self.trio_problems.append("father is not in vcf")
        if self.mother_id not in indiv:
            self.trio_is_ok_for_analysis = False
            self.trio_problems.append("mother is not in vcf")
        
        # return True if no problems
        if not self.trio_problems:
            self.trio_is_ok_for_analysis = True
            return True, []
        
    
    def print_trio_parameters(self):
        print("*********************************************")
        print("New trio")
        print("*********************************************")
        print("family_id : ", self.family_id)
        print("child_id : ", self.child_id)
        print("father_id : ", self.father_id)
        print("mother_id : ", self.mother_id)
        print("trio_is_ok_for_analysis : ", self.trio_is_ok_for_analysis)
        print("trio_problems : ", self.trio_problems)
        print("multi_vcf : ", self.multi_vcf)
        print("list_of_individuals_from_same_family_in_ped_file : ", self.list_of_individuals_from_same_family_in_ped_file)
        print("child_index_in_vcf : ", self.child_index_in_vcf)
        print("father_index_in_vcf : ", self.father_index_in_vcf)
        print("mother_index_in_vcf : ", self.mother_index_in_vcf)
    

    def print_filtration_parameters(self):
        print("*********************")
        print("Filtration parameters")
        print("*********************")
        for key, value in self.params_dict.items():
            print(key, " : ", value)
        print("family_id : ", self.family_id)
        print("child_id : ", self.child_id)
        print("father_id : ", self.father_id)
        print("mother_id : ", self.mother_id)
        

    def print_variant_counts(self):
        print("*********************************************")
        print("Variant counts for trio " + self.family_id + " (child : " + self.child_id + ")" )
        print("dnv_count : ", self.dnv_count)
        print("snv_count : ", self.snv_count)
        print("indel_count : ", self.indel_count)

    def yaml_to_dict(self):
        with open(self.params_yaml, 'r') as f:
            params_dict = yaml.safe_load(f)
        return params_dict

    def validate_parameters(self):
        # check that the dict has the right keys
        keys = ["dp_trio_min", "gq_trio_min", "VAF_prob_min", "delta_VAF_vs_parents_min", "ad1_plus_ad2_divided_by_dp_min"]
        for key in keys:
            if key not in self.params_dict:
                print("key " + key + " is missing in the config file")
                sys.exit(1)
        
        # check that the values are integers or floats
        for key, value in self.params_dict.items():
            if not isinstance(value, int) and not isinstance(value, float):
                print("value for key " + key + " is not an integer or a float")
                sys.exit(1)

    def call_de_novo(self):
        # this method uses bcftools to filter the multi vcf and get the de novo variants
        
        def get_variant_count(vcf_file):
            # uses grep to get the number of variants in a vcf file
            # command = ["grep", "-v", "^#", vcf_file, "|", "wc", "-l"]
            command = "grep -v '^#' " + vcf_file + " | wc -l"
            print("Counting variants in " + vcf_file)
            os.system(command)
        
        def get_variant_count_with_python(vcf_file):
            print("Counting variants in " + vcf_file)
            with open(vcf_file, "r") as f:
                lines = f.readlines()
                count = 0
                for line in lines:
                    if not line.startswith("#"):
                        count += 1
                return count
        
        def double_quote(string):
            return '"' + string + '"'

        
        def filter_variants(filtration_synopsis, filtration_ID, sub_command, input):
            # sub_command is the -i option of bcftools view
            print("***")
            print("Filtering process : " + filtration_synopsis)
            # if the child_id is not in the input file, we add it to the output file name
            if not self.child_id in input:
                output_file = self.tmp_dir + "/" + os.path.splitext(os.path.basename(input))[0] + "." + self.child_id + "." + filtration_ID + ".vcf"
            else:
                output_file = self.tmp_dir + "/" + os.path.splitext(os.path.basename(input))[0] + "." + filtration_ID + ".vcf"
            # command=["bcftools", "view", "-i",  sub_command, input, ">", output_file]
            command = "bcftools view -i " + sub_command + " " + input + " > " + output_file
            print("Running command : \n" + command)
            # output,error  = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            os.system(command)
            print("Filtration complete")
            get_variant_count(output_file)
            return output_file  
        
        def split_snvs_and_indels(input_file, output_dir):
            output_file_snv = output_dir + "/" + os.path.splitext(os.path.basename(input_file))[0] + ".snv.vcf"
            output_file_indel = output_dir + "/" + os.path.splitext(os.path.basename(input_file))[0] + ".indel.vcf"
            command_snv = ["bcftools", "view", "-v", "snps", "-o", output_file_snv, input_file]
            command_indel = ["bcftools", "view", "-v", "indels", "-o", output_file_indel, input_file]
            print("Running command : \n" + " ".join(command_snv))
            output,error  = subprocess.Popen(command_snv, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            print("Running command : \n" + " ".join(command_indel))
            output,error  = subprocess.Popen(command_indel, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            get_variant_count(output_file_snv)
            get_variant_count(output_file_indel)
            return output_file_snv, output_file_indel
        
        # parameters from the config file
        dp_trio_min = self.params_dict["dp_trio_min"]
        gq_trio_min = self.params_dict["gq_trio_min"]
        VAF_prob_min = self.params_dict["VAF_prob_min"]
        delta_VAF_vs_parents_min = self.params_dict["delta_VAF_vs_parents_min"]
        ad1_plus_ad2_divided_by_dp_min = self.params_dict["ad1_plus_ad2_divided_by_dp_min"]

        self.print_filtration_parameters()

        # Step 1 : filter on genotype : keep variants that are alt in the child and ref in the parents
        sub_command = 'GT[' + str(self.child_index_in_vcf) + "]='alt' && GT[" + str(self.father_index_in_vcf) + "]='ref' && GT[" + str(self.mother_index_in_vcf) + "]='ref'"
        command_synopsis = "Filtration on genotype : alt in prob, ref in parents"
        variants = filter_variants(command_synopsis, "GT", double_quote(sub_command), self.multi_vcf)

        # Step 2 : filter on depth : keep variants that have a depth greater than dp_trio_min in prob and parents
        sub_command = "FORMAT/DP[" + str(self.father_index_in_vcf) + "] > " + str(dp_trio_min) + " && FORMAT/DP[" + str(self.mother_index_in_vcf) + "] > " + str(dp_trio_min) + " && FORMAT/DP[" + str(self.child_index_in_vcf) + "] > " + str(dp_trio_min)
        command_synopsis = "Filtration on depth : dp > " + str(dp_trio_min) + " in prob and parents"
        variants = filter_variants(command_synopsis, "DP", double_quote(sub_command), variants)

        # Step 3 : filter on GQ : keep variants that have a GQ greater than gq_trio_min in prob and parents
        sub_command = "FORMAT/GQ[" + str(self.father_index_in_vcf) + "] > " + str(gq_trio_min) + " && FORMAT/GQ[" + str(self.mother_index_in_vcf) + "] > " + str(gq_trio_min) + " && FORMAT/GQ[" + str(self.child_index_in_vcf) + "] > " + str(gq_trio_min)
        command_synopsis = "Filtration on GQ : GQ > " + str(gq_trio_min) + " in prob and parents"
        variants = filter_variants(command_synopsis, "GQ", double_quote(sub_command), variants)

        # Step 4 : filter on VAF : keep variants that have an VAF greater than vaf_prob_min in the child
        sub_command = "(AD[" + str(self.child_index_in_vcf) + ":1] / FORMAT/DP[" + str(self.child_index_in_vcf) + "]) > " + str(VAF_prob_min)
        command_synopsis = "Filtration on VAF : VAF > " + str(VAF_prob_min) + " in prob"
        variants = filter_variants(command_synopsis, "VAF", double_quote(sub_command), variants)
        
        # Step 5 : filter on VAF prob / VAF parents > delta_VAF_vs_parents_min
        # sub_command = "(AD[" + str(self.child_index_in_vcf) + ":1] / FORMAT/DP[" + str(self.child_index_in_vcf) + "]) / (AD[" + str(self.father_index_in_vcf) + ":1] / FORMAT/DP[" + str(self.father_index_in_vcf) + "]) > " + str(delta_VAF_vs_parents_min) + " & (AD[" + str(self.child_index_in_vcf) + ":1] / FORMAT/DP[" + str(self.child_index_in_vcf) + "]) / (AD[" + str(self.mother_index_in_vcf) + ":1] / FORMAT/DP[" + str(self.mother_index_in_vcf) + "]) > " + str(delta_VAF_vs_parents_min)
        # sub_command = "(sum(AD[" + str(self.child_index_in_vcf) + ":1]) / sum(FORMAT/DP[" + str(self.child_index_in_vcf) + "])) / (sum(AD[" + str(self.father_index_in_vcf) + ":1]) / sum(FORMAT/DP[" + str(self.father_index_in_vcf) + "])) > " + str(delta_VAF_vs_parents_min) + " & (sum(AD[" + str(self.child_index_in_vcf) + ":1]) / sum(FORMAT/DP[" + str(self.child_index_in_vcf) + "])) / (sum(AD[" + str(self.mother_index_in_vcf) + ":1]) / sum(FORMAT/DP[" + str(self.mother_index_in_vcf) + "])) > " + str(delta_VAF_vs_parents_min)
        VAF_child = "(AD[" + str(self.child_index_in_vcf) + ":1] / FORMAT/DP[" + str(self.child_index_in_vcf) + "])"
        VAF_father = "(AD[" + str(self.father_index_in_vcf) + ":1] / sum(FORMAT/DP[" + str(self.father_index_in_vcf) + "]))"
        VAF_mother = "(AD[" + str(self.mother_index_in_vcf) + ":1] / sum(FORMAT/DP[" + str(self.mother_index_in_vcf) + "]))"
        sub_command = VAF_child + " / " + VAF_father + " > " + str(delta_VAF_vs_parents_min) + " && " + VAF_child + " / " + VAF_mother + " > " + str(delta_VAF_vs_parents_min)
        command_synopsis = "Filtration on delta_VAF_vs_parents : VAF prob / VAF parents > " + str(delta_VAF_vs_parents_min)
        variants = filter_variants(command_synopsis, "VAF_vs_parents", double_quote(sub_command), variants)
    
        # Step 6 : filter on ((AD[id:0] + AD[id:1]) / DP) > ad1_plus_ad2_divided_by_dp_min : keep variants that have an ad1_plus_ad2_divided_by_dp_min greater than ad1_plus_ad2_divided_by_dp_min in the child and parents
        sub_command = "(AD[" + str(self.child_index_in_vcf) + ":0] + AD["  + str(self.child_index_in_vcf) + ":1]) / FORMAT/DP[" + str(self.child_index_in_vcf) + "] > " + str(ad1_plus_ad2_divided_by_dp_min) + " && (AD[" + str(self.father_index_in_vcf) + ":0] + AD["  + str(self.father_index_in_vcf) + ":1]) / FORMAT/DP[" + str(self.father_index_in_vcf) + "] > " + str(ad1_plus_ad2_divided_by_dp_min) + " && (AD[" + str(self.mother_index_in_vcf) + ":0] + AD["  + str(self.mother_index_in_vcf) + ":1]) / FORMAT/DP[" + str(self.mother_index_in_vcf) + "] > " + str(ad1_plus_ad2_divided_by_dp_min)
        command_synopsis = "Filtration on ad1_plus_ad2_divided_by_dp_min : (AD[id:0] + AD[id:1]) / DP > " + str(ad1_plus_ad2_divided_by_dp_min)
        variants = filter_variants(command_synopsis, "dp_sum", double_quote(sub_command), variants)
        self.de_novo_variant_vcf = variants

        # Variant counts
        self.dnv_count = get_variant_count_with_python(variants)

        # Split snvs and indels
        snvs, indels = split_snvs_and_indels(variants, self.output_dir)

        # Variant counts
        self.snv_count = get_variant_count_with_python(snvs)
        self.indel_count = get_variant_count_with_python(indels)

    def compare_with_truth(self, truth_vcf_with_path, truth_bed_with_path, reference_with_path):
        truth_vcf_folder = os.path.dirname(truth_vcf_with_path)
        truth_vcf = os.path.basename(truth_vcf_with_path)
        truth_bed_folder = os.path.dirname(truth_bed_with_path)
        truth_bed = os.path.basename(truth_bed_with_path)
        reference_folder = os.path.dirname(reference_with_path)
        reference = os.path.basename(reference_with_path)
        vcf_1_folder = os.path.dirname(self.de_novo_variant_vcf)
        vcf_1 = os.path.basename(self.de_novo_variant_vcf)
        command = "docker run -it -v `pwd`:/data -v " + truth_vcf_folder + ":/truth_vcf -v " + truth_bed_folder + ":/truth_bed -v " + reference_folder + ":/reference -v " + vcf_1_folder + ":/vcf_1 pkrusche/hap.py /opt/hap.py/bin/som.py /data/" + truth_vcf + " /data/" + truth_bed + " /data/" + reference + " /data/" + vcf_1 + " -o /data/output"
        print("Running command : \n" + command)
        exit()


def create_trio_objects_from_ped_file(ped_file, multi_vcf, output_folder, params_yaml):
    # returns a list of trio objects
    
    def validate_ped_format(ped_file):
        # ped must have at least 4 columns : family_id, child_id, father_id, mother_id
        with open(ped_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 4:
                        print("ped file must have at least 4 columns, found " + str(len(parts)))
                        sys.exit(1)

    def get_all_individuals_from_families(ped_file):
        # this function returns a dictionnary with family_id as key and a list of individuals as value
        # example : {'fam1': ['ind1', 'ind2', 'ind3'], 'fam2': ['ind4', 'ind5', 'ind6']}
        individuals_from_families = {}
        with open(ped_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    # field 1 is family id
                    # if present, fields 2, 3 and 4 are individuals from the family
                    parts = line.rstrip().split("\t")
                    family_id = parts[0]
                    if family_id not in individuals_from_families:
                        individuals_from_families[family_id] = []
                    if len(parts) >= 2:
                        individuals_from_families[family_id].append(parts[1])
                    if len(parts) >= 3:
                        individuals_from_families[family_id].append(parts[2])
                    if len(parts) >= 4:
                        individuals_from_families[family_id].append(parts[3])
        return individuals_from_families
    
    validate_ped_format(ped_file)

    # Tmp folder
    if output_folder.endswith("/"):
        output_folder = output_folder[:-1]
    tmp_dir = output_folder + "/tmp"
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    # List of trios
    trio_list = []
    with open(ped_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                family_id, proband_id, father_id, mother_id = '', '', '', ''
                parts = line.rstrip().split("\t")
                if len(parts) >= 1:
                    family_id = parts[0]
                if len(parts) >= 2:
                    proband_id = parts[1]
                if len(parts) >= 3:
                    father_id = parts[2]
                if len(parts) >= 4:
                    mother_id = parts[3]
                list_of_individuals_from_same_family_in_ped_file = get_all_individuals_from_families(ped_file)[family_id]
                trio = Trio(family_id, proband_id, father_id, mother_id, multi_vcf, list_of_individuals_from_same_family_in_ped_file, tmp_dir, output_folder, params_yaml)
                trio_list.append(trio)
    return trio_list


def write_summary_statistics(list_of_trios):
    # this function creates a dataframe from the list of trios
    data = []
    for trio in list_of_trios:
        data.append([trio.family_id, trio.child_id, trio.father_id, trio.mother_id, trio.dnv_count, trio.snv_count, trio.indel_count])
    df = pd.DataFrame(data, columns=["family_id", "child_id", "father_id", "mother_id", "#all_candidate", "#snv_candidates", "#indel_candidates"])
    output_file = args.output_dir + "/DeNovoMiner_counts_" + os.path.splitext(os.path.basename(args.vcf))[0] + ".csv"
    df.to_csv(output_file, index=False)
    print("Results written to " + output_file)
    print(df)



##################################################### MAIN #####################################################

print("Starting DeNovoMiner.py")

list_of_trios=create_trio_objects_from_ped_file(args.ped, args.vcf, args.output_dir, args.config)
write_summary_statistics(list_of_trios)






