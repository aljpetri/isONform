import itertools
from consensus import *
#from recordclass import recordclass
from IsoformGeneration import align_to_merge
import tempfile
import pickle
from recordclass import *
from Parallelization_side_functions import *
from modules import  help_functions
""" This method is used to generate the consensus file needed for the consensus generation
INPUT:  work_dir  : The working directory in which to store the file
OUTPUT: spoa_ref:   The consensus
"""
def generate_consensus_path(work_dir,mappings1,mappings2, all_sequences,spoa_count):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    seq_count=0
    for id in mappings1:
        if seq_count<spoa_count:
                sequence=all_sequences[id]
                reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
                seq_count+=1
    for id in mappings2:
        if seq_count<spoa_count:
            if id in all_sequences:
                sequence=all_sequences[id]
            else:
                sequence=id
                reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
                seq_count += 1
    reads_path.close()
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
    return spoa_ref
def read_spoa_file(batch_id,cl_dir):
    filename ="spoa"+str(batch_id)+"merged.fa"
    batch_reads_id={}
    with open( os.path.join(cl_dir, filename)) as f:
        for id, sequence in itertools.zip_longest(*[f] * 2):
            #print(id, sequence)
            inter_id = id.replace('\n', '')
            new_id = int(inter_id.replace('>consensus', ''))
            sequence = sequence.replace('\n', '')
            # print("Seq_len",len(sequence))
            batch_reads_id[new_id] = sequence
    return batch_reads_id
def read_mapping_file(batch_id,cl_dir):
    batch_mappings_id={}
    mappingname =  "mapping"+str(batch_id)+".txt"
    incoming_ctr=0
    #print("MAPPING",batch_id,mappingname)
    with open(os.path.join(cl_dir,mappingname)) as g:
        for id, reads in itertools.zip_longest(*[g] * 2):
            # print(id, reads)
            inter_id = id.replace('\n', '')
            id = int(inter_id.replace('consensus', ''))
            # print("ID",id)
            reads = reads.replace('\n', '')
            reads = reads.replace('[', '')
            reads = reads.replace(']', '')
            # reads=reads.replace('\'','')
            reads = reads.replace("'", "")
            readslist = reads.split(", ")
            # print(readslist)
            # print(len(readslist))
            batch_mappings_id[id] = readslist

            incoming_ctr+=len(readslist)
            #print(batch_mappings_id)
    #if batch_id == 5:
        #print(batch_mappings_id)
    print("INCOMING",incoming_ctr)
    return batch_mappings_id
def read_batch_file(batch_id,all_infos_dict,all_reads_dict,cl_dir):
    batchfilename = str(batch_id)+"_batchfile.fa"
    with open(os.path.join(cl_dir,batchfilename)) as h:
        for id, seq in itertools.zip_longest(*[h] * 2):
            #print(id, seq)
            id = id.replace('\n', '')
            id = id.replace('>', '')
            # id=int(inter_id.replace('',''))
            seq = seq.replace('\n', '')
            #print(seq)
            all_reads_dict[id] = seq
            all_infos_dict[batch_id] = {}

#TODO:write skipped files into overall output
def write_final_output(all_infos_dict,outfolder,iso_abundance,cl_dir,folder):
    nr_reads = 0
    for b_id, b_infos in all_infos_dict.items():
        for c_id, c_infos in b_infos.items():
            if not c_infos.merged:
                nr_reads += len(c_infos.reads)
    print("Nr reads in all_infos_dict just before writing",nr_reads)
    print("Writing file")
    mapping_out_cter=0
    #print(folder)
    consensus_name = "cluster"+str(folder)+"_merged.fq"
    other_consensus_name = "cluster"+str(folder)+"_merged_low_abundance.fq"
    mapping_name = "cluster"+str(folder)+"_mapping.txt"
    other_mapping_name = "cluster"+str(folder)+"_mapping_low_abundance.txt"
    #print(consensus_name)
    consensus_file = open(os.path.join(outfolder, consensus_name), "w")
    other_consensus = open(os.path.join(outfolder, other_consensus_name), 'w')
    mapping_file = open(os.path.join(outfolder, mapping_name), 'w')
    other_mapping = open(os.path.join(outfolder, other_mapping_name), 'w')
    skipped_reads = {}
    read_cter_mapping=0
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            # print(id, " ", all_infos_dict[batchid][id].merged)
            if not infos.merged:
                new_id = str(batchid) + "_" + str(id)
                mapping_out_cter += len(all_infos_dict[batchid][id].reads)
                if len(infos.reads) >= iso_abundance or iso_abundance==1:
                    #print(len(all_infos_dict[batchid][id].reads))

                    #print(all_infos_dict[batchid][id].reads)
                    read_cter_mapping += len(all_infos_dict[batchid][id].reads)
                    mapping_file.write(">{0}\n{1}\n".format(new_id,  all_infos_dict[batchid][id].reads))
                    consensus_file.write("@{0}, support: {1}\n{2}\n+\n{2}\n".format(new_id, len(all_infos_dict[batchid][id].reads), all_infos_dict[batchid][id].sequence,
                                                                      "+" * len(all_infos_dict[batchid][id].sequence)))
                else:
                        other_consensus.write("@{0}\n{1}\n+\n{2}\n".format(new_id, all_infos_dict[batchid][id].sequence,
                                                                           "+" * len(
                                                                               all_infos_dict[batchid][id].sequence)))
                        other_mapping.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].reads))
    nr_skipped=0
    for skipfile in os.listdir(cl_dir):
        if skipfile.startswith("skip"):
                this_skip = 0
                with open(os.path.join(cl_dir,skipfile)) as h:
                    for id, seq in itertools.zip_longest(*[h] * 2):

                        id = id.replace('\n', '')
                        id = id.replace('>', '')
                        # id=int(inter_id.replace('',''))
                        seq = seq.replace('\n', '')
                        # print(seq)
                        skipped_reads[id] = seq
                        this_skip+=1
                #print("All reads",all_reads)
                #print(this_skip)
                nr_skipped+=this_skip
                for acc,seq in skipped_reads.items():
                    other_consensus.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq,"+"*len(seq)))
                    other_mapping.write(">{0}\n{1}\n".format(acc, acc))
    print("NR out mappings",str(mapping_out_cter))
    print("Skipped:",nr_skipped)
    #print("RCM",read_cter_mapping)
    consensus_file.close()
    mapping_file.close()
    other_consensus.close()
    other_mapping.close()
#TODO: add the rest of variables for this method and move filewriting to wrappermethod
def actual_merging_process(all_infos_dict,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,all_batch_sequences,work_dir):
    cter = 0
    seq_count_bef=0
    #print("NEWAID",all_infos_dict[5])
    for batchid,id_dict in all_infos_dict.items():
        for batchid2, id_dict2 in all_infos_dict.items():
            if not batchid2 <= batchid:# and not batchid2==batchid:
                #print("bid",batchid,"bid2",batchid2)
                for id,infos in id_dict.items():
                    id_merged=False
                    #print("IM",infos.merged)
                    #print("BID", str(batchid))
                    #print("BID2", str(batchid2))
                    #print("IDDICT", id_dict)
                    if not infos.merged:
                        #print("IDDICT", id_dict)
                        #print("IM", infos.merged)
                        for id2, infos2 in id_dict2.items():
                            if infos.merged:
                                id_merged=True
                                continue
                            #print("bid", batchid,": ",id, "bid2", batchid2,": ",id2)
                            if not infos2.merged:
                                cter+=1
                                consensus1=infos.sequence
                                consensus2=infos2.sequence
                                #print(consensus1, "\n")
                                #print(consensus2, "\n")
                                good_to_pop=align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
                                if good_to_pop:
                                    #print("Before")
                                    #print(len(all_infos_dict[batchid][id].reads))
                                    #print(len(all_infos_dict[batchid2][id2].reads))
                                    #print("Merging", batchid, "_", id, " and ", batchid2, "_", id2)
                                    #if the first combi has more than 50 reads support
                                    if len(infos.reads) > 50:
                                        #if the length of the first combis sequence is greater or equal to the length of the second combis sequence
                                        if len(infos.sequence) >= len(infos2.sequence):
                                            #print(len(all_infos_dict[batchid][id].reads))
                                            #print(len(all_infos_dict[batchid2][id2].reads))
                                            all_infos_dict[batchid][id].reads = infos2.reads + infos.reads
                                            #print(len(all_infos_dict[batchid][id].reads))
                                            #print(len(all_infos_dict[batchid2][id2].reads))
                                            infos2.merged = True
                                            print("FIRST")
                                            #all_infos_dict[batchid][id].reads=infos.reads+infos2.reads
                                        else: #the second sequence is longer than the first
                                            #add the reads of the first consensus to the second and mark the first combi as marked
                                            all_infos_dict[batchid2][id2].reads = infos2.reads+infos.reads
                                            infos.merged=True
                                            print("SECOND")

                                    else:
                                        #print("Else")
                                        # TODO generate consensus and add all infos to longer read id
                                        if len(infos.sequence) >= len(
                                                infos2.sequence):

                                            mappings1 = infos.reads
                                            mappings2 = infos2.reads

                                            new_consensus=generate_consensus_path(work_dir,mappings1,mappings2, all_batch_sequences,max_seqs_to_spoa)
                                            #print("CONS1",all_infos_dict[batchid][id])
                                            all_infos_dict[batchid][id].sequence = new_consensus
                                            all_infos_dict[batchid][id].reads = infos.reads + infos2.reads
                                            all_infos_dict[batchid2][id2].merged = True
                                            print("Third")
                                        else:
                                            #print("BM",batch_mappings)
                                            mappings1 = infos2.reads
                                            mappings2 = infos.reads
                                            #print("R1",reads1)
                                            #print("M1",mappings1)
                                            #print("R2", reads2)
                                            new_consensus=generate_consensus_path(work_dir,mappings1,mappings2, all_batch_sequences,max_seqs_to_spoa)
                                            #print("CONS2",all_infos_dict[batchid2][id2])
                                            all_infos_dict[batchid2][id2].sequence = new_consensus
                                            all_infos_dict[batchid2][id2].reads = infos2.reads + infos.reads
                                            all_infos_dict[batchid][id].merged = True
                                            print("fourth")
                            #print("After", batchid,id)
                            #print(len(all_infos_dict[batchid][id].reads))
                            #print(len(all_infos_dict[batchid2][id2].reads))
                            #print(all_infos_dict[batchid2][id2].merged)
                            #print(all_infos_dict[batchid][id].merged)
                            nr_reads = 0
                            for b_id, b_infos in all_infos_dict.items():
                                for c_id, c_infos in b_infos.items():
                                    if not c_infos.merged:
                                        nr_reads += len(c_infos.reads)
                            if not nr_reads==24135:
                                print("NR_reads in all_infos_dict", nr_reads)
    #print("AT pos",len(all_infos_dict[2][67].reads))
    print("Combi count ", cter)
    print("SCB",str(seq_count_bef))
def join_back_via_batch_merging(outdir,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,iso_abundance):
    print("Batch Merging")
    ingoing_mapping_read_cter=0
    #print(outdir)
    Read = recordclass('Read', "sequence reads merged")
    #print(outdir, tmp_work_dir)
    unique_cl_ids = set()
    full_mapping_sum=0
    subfolders = [f.path for f in os.scandir(outdir) if f.is_dir()]
    #print(subfolders)
    #iterate over all folders in the out directory
    for folder in subfolders:
        # print(file)
        tmp_fname = folder.split('/')
        fname=tmp_fname[-1]
        #print("FNAME",fname)

        cl_id = fname # file.split('_')
        unique_cl_ids.add(cl_id)
        #print(unique_cl_ids)
        #enter the folder containing all output for each cluster
    for cl_id in unique_cl_ids:
            #print("CLID",str(cl_id))
            all_infos_dict = {}
            batch_reads = {}
            batch_mappings = {}
            all_reads_dict = {}
            #out_pattern=os.path.join(outdir,"isoform")
            cl_dir = os.path.join(outdir, cl_id)
            #print("CLDIR",cl_dir)
            #iterate over all batchfiles that were generated into the cluster's folder
            for batchfile in os.listdir(cl_dir):
                #print(batchfile)
                if batchfile.endswith("_batchfile.fa"):
                    #print(batchfile)
                    tmp_bname = batchfile.split('/')
                    #print(tmp_bname)
                    tmp_bname2=tmp_bname[-1].split('_')

                    batch_id=int(tmp_bname2[0])

                    #mkdir_p(out_pattern)
                    #read all files needed to perform the batch merging and store the respective infos into all_infos_dict as well as all_reads_dict
                    read_batch_file(batch_id, all_infos_dict, all_reads_dict, cl_dir)
                    batch_reads[batch_id]=read_spoa_file(batch_id,cl_dir)
                    #print("BatchReads",batch_reads)
                    batch_mappings[batch_id]=read_mapping_file(batch_id,cl_dir)
                    ingoing_mapping_read_cter += len(batch_mappings[batch_id])
                    #print("ALL INFOS",all_infos_dict)
            #collect the information so that we have one data structure containing all infos

                #print("B_reads",batch_reads)
                    #if batch_id==5:
                        #print("B_Map",batch_mappings)

            for b_id, value in batch_reads.items():
                all_infos_batch={}
                for cons_id, seq in value.items():
                        #print("K",b_id,"V",value)
                        bm_this_batch=batch_mappings[b_id]
                        reads = bm_this_batch[cons_id]
                        full_mapping_sum+=len(reads)
                        # print("READS",reads)
                        read_mapping = Read(seq, reads, False)
                        #print("Infos",b_id,cons_id)
                        #if b_id==5:
                            #print("Mappiong",cons_id,"::",read_mapping)
                        all_infos_batch[cons_id] = read_mapping
                all_infos_dict[b_id]=all_infos_batch
            #print("AID",all_infos_dict[5])
            nr_reads=0
            for b_id, b_infos in all_infos_dict.items():
                for c_id, c_infos in b_infos.items():
                    if not c_infos.merged:
                        nr_reads += len(c_infos.reads)
            #perform the merging step during which all consensuses are compared and if possible merged

            actual_merging_process(all_infos_dict,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,all_reads_dict,outdir)
            nr_reads=0
            for b_id, b_infos in all_infos_dict.items():
                for c_id, c_infos in b_infos.items():
                    if not c_infos.merged:
                        nr_reads += len(c_infos.reads)
            print("NR_reads in all_infos_dict", nr_reads)
            #write the final output into files
            #print("AID2",all_infos_dict)
            write_final_output(all_infos_dict,outdir,iso_abundance,cl_dir,cl_id)
    print("FUll mapping sum",full_mapping_sum)
    #rint("#mappings incoming",ingoing_mapping_read_cter)
        #shutil.rmtree(batch_id)
def main():
    #outfolder= "/home/alexanderpetri/isONform_analysis/Paraout_500_cl0"
    #outfolder = "/home/alexanderpetri/isONform_analysis/Para_out_500batchadd"
    #outfolder="/home/alexanderpetri/isONform_analysis/Batch_parallel_testing"
    outfolder = "/home/alexanderpetri/isONform_analysis/Para_out_500batch"
    outfolder="/home/alexanderpetri/isONform_analysis/Para_out_500_September"
    delta=0.10
    delta_len=5
    merge_sub_isoforms_3=True
    merge_sub_isoforms_5 = True
    delta_iso_len_3= 30
    delta_iso_len_5 = 50
    max_seqs_to_spoa=200
    iso_abundance=1
    join_back_via_batch_merging(outfolder,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,iso_abundance)
    generate_full_output(outfolder)
    #merge_batches(max_batchid, tmp_work_dir, outfolder, all_reads,merge_sub_isoforms_3,merge_sub_isoforms_5, delta, delta_len,max_seqs_to_spoa, delta_iso_len_3, delta_iso_len_5,iso_abundance)
    print("removing temporary workdir")
    #shutil.rmtree(work_dir)

if __name__ == "__main__":
    main()