import argparse
import os
import itertools
import sys
from pathlib import Path

def read_mapping_file(mapping_name):
    batch_mappings_id={}
    #mappingname =  "mapping"+str(batch_id)+".txt"
    incoming_ctr=0
    #print("MAPPING",batch_id,mappingname)
    with open(mapping_name) as g:
        for id, reads in itertools.zip_longest(*[g] * 2):
            #print(id, reads)
            inter_id = id.replace('\n', '')
            #id = int(inter_id.replace('consensus', ''))
            # print("ID",id)
            #print(reads)
            reads = reads.replace('\n', '')
            reads = reads.replace('[', '')
            reads = reads.replace(']', '')
            # reads=reads.replace('\'','')
            reads = reads.replace("'", "")
            readslist = reads.split(", ")
            #print(readslist)
            # print(len(readslist))
            batch_mappings_id[id] = readslist

            incoming_ctr+=len(readslist)
            #print(batch_mappings_id)
    #if batch_id == 5:
        #print(batch_mappings_id)
    print("INCOMING",incoming_ctr, ", ",mapping_name)
    return incoming_ctr


def main(args):
    print("Hello World")

    directory = args.fastq_folder #os.fsencode(args.fastq_folder)
    pat = Path(directory)
    file_list = list(pat.rglob('*_mapping.txt'))
    incoming_cter_glob=0
    for file in file_list:
        if not file=="/home/alexanderpetri/isONform_analysis/Para_out_500_September/transcriptome_mapping.txt":
            incoming_cter_glob+=read_mapping_file(file)
    print("High Abundance",incoming_cter_glob)
    other_file_list = list(pat.rglob('*_mapping_low_abundance.txt'))
    incoming_cter_glob_low = 0
    for o_file in other_file_list:
        incoming_cter_glob_low+=read_mapping_file(o_file)
    print("LOW ABUNDANCE",incoming_cter_glob_low)
    sum=incoming_cter_glob+incoming_cter_glob_low
    print("SUM",sum)

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description="De novo clustering of long-read transcriptome reads",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('--fastq_folder', type=str, default=False,
                            help='Path to input fastq folder with reads in clusters')
        args = parser.parse_args()
        main(args)
        """if len(sys.argv) == 1:
            parser.print_help()
            sys.exit()"""