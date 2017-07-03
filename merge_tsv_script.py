

import os
import sys
import re
#python.exe script_name.py "path_to_gtf" "path_to_tsv_folder"
# globals
path_to_gtf = sys.argv[1]
path_to_tsv_folder = sys.argv[2]

# functions


def create_convert_dict(path_to_gtf=path_to_gtf):
    convert_dict = {}
    with open(path_to_gtf, 'r') as annotation:
        for line in annotation:
            gtf_file_string = line.split()
            gene_id = gtf_file_string[9][1:-2]
            transcript_id = gtf_file_string[11][1:-2]
            if gene_id in convert_dict:
                convert_dict[transcript_id].append(gene_id)
            else:
                convert_dict[transcript_id] = [gene_id]
    return convert_dict


def write_convert_dict(convert_dict, output_file='convert_dict.txt'):
    with open(output_file, 'w') as out_file:
        for keys in convert_dict.keys():
            r = convert_dict.get(keys)
            t = str(keys) + str('\t') + str(r) + str('\n')
            out_file.write(t)


def get_tsv_files_list(path_to_tsv_folder=path_to_tsv_folder):
    tsv_files_list = []
    for d, dirs, files in os.walk(path_to_tsv_folder):
        for f in files:
            if f.endswith('.tsv'):
                tsv_file = os.path.join(d, f)
                tsv_files_list.append(tsv_file)
    return tsv_files_list


def print_count_gene(count_gene, tsv_file_path):
    p = re.compile(r'SRR[0-9]*')
    output_file = p.search(tsv_file_path).group() + '.result'
    with open(output_file, 'w') as out:
        out.write('target_id\ttpm\n')
        for gene_id, counts in count_gene.items():
            out.write('{0}\t{1}\n'.format(gene_id, counts))


def count_gene_tpm(tsv_file_path):
    """ Returns tuple"""
    count_gene = {}
    with open(tsv_file_path, 'r') as table:
        for line in table:
            if not line.startswith("target_id"):
                tsv_file_line = line.split()
                transcript_id_from_tsv = tsv_file_line[0]
                gene_id = convert_dict[transcript_id_from_tsv][0]
                tpm = tsv_file_line[3]
                if gene_id in count_gene:
                    count_gene[gene_id] += float(tpm)
                else:
                    count_gene[gene_id] = float(tpm)

    return count_gene


# begin
convert_dict = create_convert_dict(path_to_gtf)
write_convert_dict(convert_dict)
tsv_files_list = get_tsv_files_list()

for tsv_file in tsv_files_list:
    print_count_gene(count_gene_tpm(tsv_file), tsv_file)

# end
