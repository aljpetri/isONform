import itertools
import pickle
import os

import shutil
from modules import IsoformGeneration
from modules import Parallelization_side_functions

class Read:
    def __init__(self, sequence, reads, merged):
        self.sequence = sequence
        self.reads = reads
        self.merged = merged

def generate_consensus_path(work_dir, mappings1, mappings2, all_sequences, spoa_count):
    """ This method is used to generate the consensus file needed for the consensus generation
    INPUT:  work_dir  : The working directory in which to store the file
    OUTPUT: spoa_ref:   The consensus
    """
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    seq_count = 0
    for id in mappings1:
        if seq_count < spoa_count:
            sequence = all_sequences[id]
            reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
            seq_count += 1
    for id in mappings2:
        if seq_count < spoa_count:
            if not (id in all_sequences):
                sequence = id
                reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
                seq_count += 1
    reads_path.close()
    spoa_ref = IsoformGeneration.run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"))
    return spoa_ref


""" This function reads the consensus file and saves the infos in batch_reads_id
INPUT:  work_dir  : The working directory in which to store the file
OUTPUT: spoa_ref:   The consensus
"""


def read_spoa_file(batch_id, cl_dir):
    filename = "spoa" + str(batch_id) + "merged.fasta"
    batch_reads_id = {}
    with open(os.path.join(cl_dir, filename)) as f:
        for id, sequence in itertools.zip_longest(*[f] * 2):
            inter_id = id.replace('\n', '')
            new_id = int(inter_id.replace('>consensus', ''))
            sequence = sequence.replace('\n', '')
            batch_reads_id[new_id] = sequence
    return batch_reads_id


def read_mapping_file(batch_id, cl_dir):
    batch_mappings_id = {}
    mappingname = "mapping" + str(batch_id) + ".txt"
    incoming_ctr = 0
    with open(os.path.join(cl_dir, mappingname)) as g:
        for id, reads in itertools.zip_longest(*[g] * 2):
            inter_id = id.replace('\n', '')
            id = int(inter_id.replace('consensus', ''))
            reads = reads.replace('\n', '')
            reads = reads.replace('[', '')
            reads = reads.replace(']', '')
            reads = reads.replace("'", "")
            readslist = reads.split(", ")
            batch_mappings_id[id] = readslist
            incoming_ctr += len(readslist)
    return batch_mappings_id


def read_batch_file(batch_id, all_infos_dict, all_reads_dict, cl_dir):
    batchfilename = str(batch_id) + "_batchfile.fa"
    with open(os.path.join(cl_dir, batchfilename)) as h:
        for id, seq in itertools.zip_longest(*[h] * 2):
            id = id.replace('\n', '')
            id = id.replace('>', '')
            seq = seq.replace('\n', '')
            all_reads_dict[id] = seq
            all_infos_dict[batch_id] = {}


def write_final_output(all_infos_dict, outfolder, iso_abundance, cl_dir, folder, write_fastq):
    write_low_abundance = False
    support_name = "support_" + str(folder) + ".txt"
    other_support_name = "support_" + str(folder) + "low_abundance.txt"
    if write_fastq:
        consensus_name = "cluster" + str(folder) + "_merged.fq"
        other_consensus_name = "cluster" + str(folder) + "_merged_low_abundance.fq"
    else:
        consensus_name = "cluster" + str(folder) + "_merged.fa"
        other_consensus_name = "cluster" + str(folder) + "_merged_low_abundance.fa"
    mapping_name = "cluster" + str(folder) + "_mapping.txt"
    other_mapping_name = "cluster" + str(folder) + "_mapping_low_abundance.txt"
    support_file = open(os.path.join(outfolder, support_name), "w")
    other_support_file = open(os.path.join(outfolder, other_support_name), "w")
    consensus_file = open(os.path.join(outfolder, consensus_name), "w")
    other_consensus = open(os.path.join(outfolder, other_consensus_name), 'w')
    mapping_file = open(os.path.join(outfolder, mapping_name), 'w')
    other_mapping = open(os.path.join(outfolder, other_mapping_name), 'w')
    skipped_reads = {}
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            new_id = str(folder) + "_" + str(batchid) + "_" + str(id)
            if not infos.merged:
                if len(infos.reads) >= iso_abundance or iso_abundance == 1:
                    mapping_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].reads))
                    if write_fastq:
                        consensus_file.write("@{0}\n{1}\n+\n{2}\n".format(new_id, all_infos_dict[batchid][id].sequence,
                                                                      "+" * len(all_infos_dict[batchid][id].sequence)))
                    else:
                        consensus_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
                    support_file.write("{0}: {1}\n".format(new_id, len(all_infos_dict[batchid][id].reads)))
                else:
                    if write_low_abundance:
                        if write_fastq:
                            other_consensus.write("@{0}\n{1}\n+\n{2}\n".format(new_id, all_infos_dict[batchid][id].sequence,
                                                                           "+" * len(
                                                                               all_infos_dict[batchid][id].sequence)))
                        else:
                            other_consensus.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
                        if new_id in all_infos_dict:
                            other_support_file.write("{0}: {1}\n".format(new_id, len(all_infos_dict[new_id].reads)))
                        else:
                            other_support_file.write("{0}: {1}\n".format(new_id, 1))
                        other_mapping.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].reads))
    if write_low_abundance:
        for skipfile in os.listdir(cl_dir):
            if skipfile.startswith("skip"):
                with open(os.path.join(cl_dir, skipfile)) as h:
                    for id, seq in itertools.zip_longest(*[h] * 2):
                        id = id.replace('\n', '')
                        id = id.replace('>', '')
                        seq = seq.replace('\n', '')
                        skipped_reads[id] = seq

                for acc, seq in skipped_reads.items():
                    other_consensus.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, "+" * len(seq)))
                    other_support_file.write("{0}: {1}\n".format(acc, len(seq)))
                    other_mapping.write(">{0}\n{1}\n".format(acc, acc))

    consensus_file.close()
    mapping_file.close()
    other_consensus.close()
    other_mapping.close()


# TODO: add the rest of variables for this method and move filewriting to wrappermethod


def actual_merging_process(all_infos_dict, delta, delta_len,
                           delta_iso_len_3, delta_iso_len_5, max_seqs_to_spoa, all_batch_sequences, work_dir):
    all_infos_list = sorted(all_infos_dict.items(), key=lambda x: x[0], reverse=True)
    for b_i, (batchid, id_dict) in enumerate(all_infos_list[:len(all_infos_list) - 1]):
        batch_id_list = sorted(id_dict.items(), key=lambda x: len(x[1].sequence))
        for b_j, (batchid2, id_dict2) in enumerate(all_infos_list[b_i + 1:]):
            batch_id_list2 = sorted(id_dict2.items(), key=lambda x: len(x[1].sequence))
            if not batchid2 <= batchid:
                for t_id1, (id, infos) in enumerate(batch_id_list):
                    if not infos.merged:
                        for t_id2, (id2, infos2) in enumerate(batch_id_list2):
                            if len(infos2.sequence) >= len(infos.sequence):
                                batch_id_long = batchid2
                                batch_id_short = batchid
                                id_long = id2
                                infos_long = infos2
                                id_short = id
                                infos_short = infos
                            else:
                                batch_id_long = batchid
                                batch_id_short = batchid2
                                id_long = id
                                infos_long = infos
                                id_short = id2
                                infos_short = infos2

                            if infos_long.merged:
                                continue
                            if DEBUG:
                                print("bid", batchid, ": ", id, "bid2", batchid2, ": ", id2)
                            if not infos_short.merged:
                                consensus1 = infos_long.sequence
                                consensus2 = infos_short.sequence
                                good_to_pop = IsoformGeneration.align_to_merge(consensus1, consensus2, delta, delta_len,
                                                                                delta_iso_len_3,
                                                                               delta_iso_len_5)
                                if DEBUG:
                                    print(good_to_pop)
                                if good_to_pop:
                                    # if the first combi has more than 50 reads support
                                    if len(infos_long.reads) > 50:
                                        all_infos_dict[batch_id_long][
                                            id_long].reads = infos_long.reads + infos_short.reads
                                        infos_short.merged = True
                                    else:
                                        mappings1 = infos_long.reads
                                        mappings2 = infos_short.reads
                                        new_consensus = generate_consensus_path(work_dir, mappings1, mappings2,
                                                                                all_batch_sequences, max_seqs_to_spoa)
                                        all_infos_dict[batch_id_long][id_long].sequence = new_consensus
                                        all_infos_dict[batch_id_long][
                                            id_long].reads = infos_long.reads + infos_short.reads
                                        all_infos_dict[batch_id_short][id_short].merged = True


def join_back_via_batch_merging(outdir, delta, delta_len, delta_iso_len_3,
                                delta_iso_len_5, max_seqs_to_spoa, iso_abundance,write_fastq):
    print("Batch Merging")
    unique_cl_ids = set()
    subfolders = [f.path for f in os.scandir(outdir) if f.is_dir()]
    # iterate over all folders in the out directory
    for folder in subfolders:
        tmp_fname = folder.split('/')
        fname = tmp_fname[-1]
        cl_id = fname
        unique_cl_ids.add(cl_id)
    # enter the folder containing all output for each cluster
    for cl_id in unique_cl_ids:
        #print(cl_id)
        all_infos_dict = {}
        batch_reads = {}
        batch_mappings = {}
        all_reads_dict = {}
        cl_dir = os.path.join(outdir, cl_id)
        # iterate over all batchfiles that were generated into the cluster's folder
        for batchfile in os.listdir(cl_dir):
            if batchfile.endswith("_batch"):
                tmp_bname = batchfile.split('/')
                tmp_bname2 = tmp_bname[-1].split('_')
                batch_id = int(tmp_bname2[0])
                #print(batch_id)
                batchfilename = str(batch_id) + "_batch"
                batch_file = open(os.path.join(cl_dir, batchfilename), 'rb')
                all_reads_dict[batch_id] = pickle.load(batch_file)
                for key in all_reads_dict.keys():
                    all_infos_dict[key] = {}
                spoa_name = "spoa" + str(batch_id)
                #print(os.path.join(cl_dir, spoa_name))
                spoa_file = open(os.path.join(cl_dir, spoa_name), 'rb')
                batch_reads[batch_id] = pickle.load(spoa_file)
                map_name = "mapping" + str(batch_id)
                map_file = open(os.path.join(cl_dir, map_name), 'rb')
                batch_mappings[batch_id] = pickle.load(map_file)
        # collect the information so that we have one data structure containing all infos
        for b_id, value in batch_reads.items():
            all_infos_batch = {}
            for cons_id, seq in value.items():
                bm_this_batch = batch_mappings[b_id]
                reads = bm_this_batch[cons_id]
                read_mapping = Read(seq, reads, False)
                all_infos_batch[cons_id] = read_mapping
            all_infos_dict[b_id] = all_infos_batch

        nr_reads = 0
        for b_id, b_infos in all_infos_dict.items():
            for c_id, c_infos in b_infos.items():
                if not c_infos.merged:
                    nr_reads += len(c_infos.reads)
        # perform the merging step during which all consensuses are compared and if possible merged
        actual_merging_process(all_infos_dict, delta, delta_len,
                               delta_iso_len_3, delta_iso_len_5, max_seqs_to_spoa, all_reads_dict, outdir)
        nr_reads = 0
        for b_id, b_infos in all_infos_dict.items():
            for c_id, c_infos in b_infos.items():
                if not c_infos.merged:
                    nr_reads += len(c_infos.reads)
            write_final_output(all_infos_dict, outdir, iso_abundance, cl_dir, cl_id, write_fastq)
        #shutil.rmtree(os.path.join(outdir,cl_id))

DEBUG = False


def main():
    # outfolder= "/home/alexanderpetri/isONform_analysis/Paraout_500_cl0"
    # outfolder = "/home/alexanderpetri/isONform_analysis/Para_out_500batchadd"
    # outfolder="/home/alexanderpetri/isONform_analysis/Batch_parallel_testing"
    # outfolder = "/home/alexanderpetri/isONform_analysis/Para_out_500batch"
    outfolder = "/home/alexanderpetri/Desktop/IsONform_test_results/SIRV100k/552997d_t2"
    delta = 0.15
    delta_len = 5
    merge_sub_isoforms_3 = True
    merge_sub_isoforms_5 = True
    delta_iso_len_3 = 30
    delta_iso_len_5 = 50
    max_seqs_to_spoa = 200
    iso_abundance = 4
    join_back_via_batch_merging(outfolder, delta, delta_len, delta_iso_len_3, delta_iso_len_5,
                                max_seqs_to_spoa, iso_abundance)

    Parallelization_side_functions.generate_full_output(outfolder)
    # merge_batches(max_batchid, tmp_work_dir, outfolder, all_reads,merge_sub_isoforms_3,merge_sub_isoforms_5, delta, delta_len,max_seqs_to_spoa, delta_iso_len_3, delta_iso_len_5,iso_abundance)
    print("removing temporary workdir")
    # shutil.rmtree(work_dir)


if __name__ == "__main__":
    main()
