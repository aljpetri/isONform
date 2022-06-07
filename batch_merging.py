from collections import namedtuple
import itertools
from consensus import *
from IsoformGeneration import align_to_merge,parse_cigar_diversity_isoform_level
"""The last step of our algorithm: We take all isoforms that were generated for each batch and align them to merge same isoforms
"""
def merge_batches(max_batchid, outfolder,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa,iso_abundance, delta_iso_len_3, delta_iso_len_5):
    all_infos_dict={}
    Read = namedtuple('Read',"sequence reads merged")
    for batchid in range(0,max_batchid+1):
        batch_reads={}
        batch_mappings={}
        filename="spoa" + str(batchid) + "merged.fa"
        mappingname= "mapping" + str(batchid) + ".txt"
        print("File: ",filename)
        all_infos_dict[batchid]={}
        with open(filename) as f:
            for id, sequence in itertools.zip_longest(*[f] * 2):
                print(id, sequence)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('>consensus',''))
                sequence = sequence.replace('\n', '')
                print("Seq_len",len(sequence))
                batch_reads[id]=sequence
        print(batch_reads)
        with open(mappingname) as g:
            for id, reads in itertools.zip_longest(*[g] * 2):
                print(id, reads)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('consensus',''))
                reads = reads.replace('\n', '')
                batch_mappings[id]=reads
                print("reads_len",len(reads))
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
            if not infos.merged:
                for batchid2,id_dict2 in all_infos_dict.items():
                    if not batchid2<=batchid:
                        for id2, infos2 in id_dict2.items():
                            if not infos2.merged:
                                cter+=1
                                consensus1=infos.sequence
                                consensus2=infos2.sequence
                                print(consensus1, "\n")
                                print(consensus2, "\n")
                                good_to_pop=align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
                                if good_to_pop:
                                    if len(infos.reads)>50:
                                        all_infos_dict[batchid][id]._replace(merged = True)
                                        all_infos_dict[batchid2][id2]._replace(merged = True)
                                    print("True")
    #print(all_infos_dict)
    print("Combi count ",cter)
    print("DONE")
def main():
    outfolder = "out_local"
    delta=0.10
    k_size=20
    delta_len=k_size
    batch_id=0
    merge_sub_isoforms_3=True
    merge_sub_isoforms_5 = True
    delta_iso_len_3=30
    delta_iso_len_5 = 50
    max_seqs_to_spoa=200
    iso_abundance=1
    max_batchid=1
    merge_batches(max_batchid, outfolder,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa,iso_abundance, delta_iso_len_3, delta_iso_len_5)
if __name__ == "__main__":
    main()