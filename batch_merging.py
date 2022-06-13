from collections import namedtuple
import itertools
from consensus import *
from recordclass import recordclass
from IsoformGeneration import align_to_merge,parse_cigar_diversity_isoform_level
"""The last step of our algorithm: We take all isoforms that were generated for each batch and align them to merge same isoforms
"""
def merge_batches(max_batchid, outfolder,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa,iso_abundance, delta_iso_len_3, delta_iso_len_5,all_batch_reads):
    all_infos_dict={}
    merged_id=0
    merged_batches={}
    Read = recordclass('Read',"sequence reads merged")
    for batchid in range(0,max_batchid+1):
        batch_reads={}
        batch_mappings={}
        filename="spoa" + str(batchid) + "merged.fa"
        mappingname= "mapping" + str(batchid) + ".txt"
        print("File: ",filename)
        all_infos_dict[batchid]={}
        with open(filename) as f:
            for id, sequence in itertools.zip_longest(*[f] * 2):
                #print(id, sequence)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('>consensus',''))
                sequence = sequence.replace('\n', '')
                #print("Seq_len",len(sequence))
                batch_reads[id]=sequence
        #print(batch_reads)
        with open(mappingname) as g:
            for id, reads in itertools.zip_longest(*[g] * 2):
                #print(id, reads)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('consensus',''))
                reads = reads.replace('\n', '')
                batch_mappings[id]=reads
                #print("reads_len",len(reads))
        for key, value in batch_reads.items():
            reads=batch_mappings[key]
            read_mapping=Read(value,reads,False)
            all_infos_dict[batchid][key]=read_mapping
    used_combis={}
    print("Len of batchid 0:",len(all_infos_dict[0]))
    print("Len of batchid 1:", len(all_infos_dict[1]))
    cter=0
    for batchid,id_dict in all_infos_dict.items():

        for id,infos in id_dict.items():
            #if batchid ==0 and id==2:
            #print("BATCH",all_infos_dict[batchid])
            if not infos.merged:
                for batchid2,id_dict2 in all_infos_dict.items():
                    if not batchid2<=batchid:
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
                                        print(True)
                                        supporting_reads=infos.reads+infos2.reads
                                        entries=Read(consensus1,supporting_reads,False)
                                        merged_batches[merged_id]=entries
                                        merged_id+=1
                                        if len(all_infos_dict[batchid][id].sequence)>=len(all_infos_dict[batchid2][id2].sequence):
                                            all_infos_dict[batchid][id].reads=reads+infos2.reads
                                        #all_infos_dict[batchid][id]._replace(merged = True)
                                            all_infos_dict[batchid2][id2].merged=True
                                        else:
                                            #sequence=all_infos_dict[batchid2][id2].sequence
                                            #reads=reads+infos.reads
                                            #merged=False
                                            #new_Tuple=Read(sequence,reads,merged)
                                            #all_infos_dict[batchid2][id2]=new_Tuple
                                            all_infos_dict[batchid2][id2].reads=reads + infos.reads
                                            # all_infos_dict[batchid][id]._replace(merged = True)
                                            #other_seq=all_infos_dict[batchid][id].sequence
                                            #other_new_Tuple = Read(other_seq, [], True)
                                            #all_infos_dict[batchid][id]=other_new_Tuple
                                            all_infos_dict[batchid][id].merged=True
                                        print("Merging",batchid,"_",id," and ",batchid2,"_",id2)
                                    else:
                                        print("CONS1",all_infos_dict[batchid][id])
                                        print("CONS2",all_infos_dict[batchid2][id2])
    #for m_batch in merged_batches
    print("MB",merged_batches," ",len(merged_batches))
    #print(all_infos_dict)
    print("Combi count ",cter)
    print("Writing file")
    consensus_name = "cluster_merged.fa"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            print(id," ",all_infos_dict[batchid][id].merged)
            if not all_infos_dict[batchid][id].merged:
                new_id=str(batchid)+"_"+str(id)
                consensus_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
    consensus_file.close()
def main():
    outfolder = "out_local"
    delta=0.10
    k_size=20
    delta_len=k_size
    merge_sub_isoforms_3=True
    merge_sub_isoforms_5 = True
    delta_iso_len_3= 30
    delta_iso_len_5 = 50
    max_seqs_to_spoa=200
    iso_abundance=1
    max_batchid=1
    merge_batches(max_batchid, outfolder,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa,iso_abundance, delta_iso_len_3, delta_iso_len_5)
if __name__ == "__main__":
    main()