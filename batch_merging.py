import itertools
from consensus import *
from recordclass import recordclass
from IsoformGeneration import align_to_merge,parse_cigar_diversity_isoform_level
import tempfile
import pickle

""" This method is used to generate the consensus file needed for the consensus generation

INPUT:  work_dir  : The working directory in which to store the file


OUTPUT: spoa_ref:   The consensus
"""
def generate_consensus_path(work_dir,mappings1,mappings2, reads1,reads2,spoa_count):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    seq_count=0
    for id in mappings1:
        if seq_count<spoa_count:
            sequence=reads1[id]
            reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
    for id in mappings2:
        if seq_count<spoa_count:
            sequence=reads2[id]
            reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
    reads_path.close()
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
    return spoa_ref
"""The last step of our algorithm: We take all isoforms that were generated for each batch and align them to merge same isoforms
"""
def merge_batches(max_batchid,work_dir, outfolder,all_reads,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa, delta_iso_len_3, delta_iso_len_5,iso_abundance):
    all_infos_dict={}
    #Todo: add working directory
    Read = recordclass('Read',"sequence reads merged")
    for batchid in range(0,max_batchid+1):
        batch_reads={}
        batch_mappings={}
        batch_sequences={}
        filename="spoa" + str(batchid) + "merged.fa"
        batchfilename=str(batchid)+"_batchfile.fa"
        mappingname= "mapping" + str(batchid) + ".txt"
        print("File: ",filename)
        all_infos_dict[batchid]={}
        with open(os.path.join(outfolder,filename)) as f:
            for id, sequence in itertools.zip_longest(*[f] * 2):
                #print(id, sequence)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('>consensus',''))
                sequence = sequence.replace('\n', '')
                #print("Seq_len",len(sequence))
                batch_reads[id]=sequence
        #print(batch_reads)
        with open(os.path.join(outfolder,mappingname)) as g:
            for id, reads in itertools.zip_longest(*[g] * 2):
                #print(id, reads)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('consensus',''))
                reads = reads.replace('\n', '')
                batch_mappings[id]=reads
                #print("reads_len",len(reads))
        with open(os.path.join(outfolder,batchfilename))as h:
            for id, seq in itertools.zip_longest(*[h] * 2):
                #print(id, reads)
                id=id.replace('\n','')
                #id=int(inter_id.replace('',''))
                seq = seq.replace('\n', '')
                batch_sequences[id]=seq
        for key, value in batch_reads.items():
            reads=batch_mappings[key]
            read_mapping=Read(value,reads,False)
            all_infos_dict[batchid][key]=read_mapping
    print("Len of batchid 0:",len(all_infos_dict[0]))
    print("Len of batchid 1:", len(all_infos_dict[1]))
    cter=0
    for batchid,id_dict in all_infos_dict.items():
        for batchid2, id_dict2 in all_infos_dict.items():
            if not batchid2 <= batchid:# and not batchid2==batchid:
                for id,infos in id_dict.items():
                    if not infos.merged:
                        for id2, infos2 in id_dict2.items():
                            if not infos2.merged:
                                cter+=1
                                consensus1=infos.sequence
                                consensus2=infos2.sequence
                                #print(consensus1, "\n")
                                #print(consensus2, "\n")
                                good_to_pop=align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
                                if good_to_pop:
                                    if len(infos.reads)>50:
                                        if len(all_infos_dict[batchid][id].sequence)>=len(all_infos_dict[batchid2][id2].sequence):
                                            all_infos_dict[batchid][id].reads=infos.reads+infos2.reads
                                            all_infos_dict[batchid2][id2].merged=True
                                        else:
                                            all_infos_dict[batchid2][id2].reads=infos2.reads +infos.reads
                                            all_infos_dict[batchid][id].merged=True
                                        print("Merging",batchid,"_",id," and ",batchid2,"_",id2)
                                    else:
                                        print("Else")
                                        # TODO generate consensus and add all infos to longer read id
                                        if len(all_infos_dict[batchid][id].sequence) >= len(
                                                all_infos_dict[batchid2][id2].sequence):
                                            reads1=all_reads[batchid]
                                            mappings1=batch_mappings[batchid][id]
                                            mappings2 = batch_mappings[batchid2][id2]
                                            reads2=all_reads[batchid2]
                                            new_consensus=generate_consensus_path(work_dir,mappings1,mappings2, reads1,reads2,max_seqs_to_spoa)
                                            print("CONS1",all_infos_dict[batchid][id])
                                            all_infos_dict[batchid][id].sequence =new_consensus
                                            all_infos_dict[batchid][id].reads = infos.reads + infos2.reads
                                            all_infos_dict[batchid2][id2].merged = True
                                        else:
                                            reads1 = all_reads[batchid2]
                                            mappings1 = batch_mappings[batchid2][id2]
                                            mappings2 = batch_mappings[batchid][id]
                                            reads2 = all_reads[batchid]
                                            new_consensus=generate_consensus_path(work_dir,mappings1,mappings2, reads1,reads2,max_seqs_to_spoa)
                                            print("CONS2",all_infos_dict[batchid2][id2])
                                            all_infos_dict[batchid2][id2].sequence = new_consensus
                                            all_infos_dict[batchid2][id2].reads = infos2.reads + infos.reads
                                            all_infos_dict[batchid][id].merged = True

    #for m_batch in merged_batches
    #print("MB",merged_batches," ",len(merged_batches))
    #print(all_infos_dict)
    print("Combi count ",cter)
    print("Writing file")
    consensus_name = "cluster_merged.fa"
    mapping_name="cluster_mapping.txt"

    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    mapping_file = open(os.path.join(outfolder, mapping_name), 'w')
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            print(id," ",all_infos_dict[batchid][id].merged)
            if not all_infos_dict[batchid][id].merged:
                if len(all_infos_dict[batchid][id].reads)>=iso_abundance:
                    new_id=str(batchid)+"_"+str(id)
                    mapping_file.write(">{0}\n{1}\n".format(new_id,all_infos_dict[batchid][id].reads))
                    consensus_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
    consensus_file.close()
    mapping_file.close()

def main():
    outfolder = ""
    delta=0.10
    delta_len=20
    merge_sub_isoforms_3=True
    merge_sub_isoforms_5 = True
    delta_iso_len_3= 30
    delta_iso_len_5 = 50
    max_seqs_to_spoa=200
    iso_abundance=1
    max_batchid=24
    work_dir = tempfile.mkdtemp()
    print("Temporary workdirektory:", work_dir)
    with open(os.path.join(outfolder, "all_batches_reads.txt"), 'rb') as file:
        all_reads = pickle.load(file)

    merge_batches(max_batchid, work_dir, outfolder, all_reads,merge_sub_isoforms_3,merge_sub_isoforms_5, delta, delta_len,max_seqs_to_spoa, delta_iso_len_3, delta_iso_len_5,iso_abundance)
    print("removing temporary workdir")
    shutil.rmtree(work_dir)

if __name__ == "__main__":
    main()