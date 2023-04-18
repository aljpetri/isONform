import itertools
from modules.consensus import *
from recordclass import recordclass
from IsoformGeneration import align_to_merge
import tempfile
import pickle
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


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
                sequence=all_sequences[id]
                reads_path.write(">{0}\n{1}\n".format(str(id), sequence))
                seq_count += 1
    reads_path.close()
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
    return spoa_ref
"""The last step of our algorithm: We take all isoforms that were generated for each batch and align them to merge same isoforms
"""
def merge_batches(max_batchid,work_dir, outfolder,all_reads,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa, delta_iso_len_3, delta_iso_len_5,iso_abundance,rc_threshold):
    all_infos_dict={}
    #if max_batchid==0:
    #    print("returning from batch-merging")
    #    return
    #Todo: add working directory
    Read = recordclass('Read',"sequence reads merged")
    seq_count=0
    all_batch_sequences={}
    #print(list(all_reads.keys()))
    for batchid in range(0,max_batchid+1):
        batch_reads={}
        batch_mappings={}

        filename="spoa" + str(batchid) + "merged.fasta"
        batchfilename=str(batchid)+"_batchfile.fa"
        mappingname= "mapping" + str(batchid) + ".txt"
        #print("File: ",filename)
        all_infos_dict[batchid]={}
        fq_res=os.path.join(outfolder,filename)
        #fastq = { acc : (seq,qual) for acc, (seq,qual) in readfq(open(fq_res, 'r'))}
        #print(fastq)
        with open(os.path.join(outfolder,filename)) as f:
            for id, sequence in itertools.zip_longest(*[f] * 2):
                #print(id, sequence)
                inter_id=id.replace('\n','')
                #print(inter_id)
                id=int(inter_id.replace('>consensus',''))
                sequence = sequence.replace('\n', '')
                #print("Seq_len",len(sequence))
                batch_reads[id]=sequence
                seq_count+=1
        #print(batch_reads)
        with open(os.path.join(outfolder,mappingname)) as g:
            for id, reads in itertools.zip_longest(*[g] * 2):
                #print(id, reads)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('consensus',''))
                #print("ID",id)
                reads = reads.replace('\n', '')
                reads=reads.replace('[','')
                reads = reads.replace(']', '')
                #reads=reads.replace('\'','')
                reads = reads.replace("'", "")
                readslist=reads.split(", ")
                #print(readslist)
                #print(len(readslist))
                batch_mappings[id]=readslist
        #print("BATCHLEn",len(batch_mappings))
        #print("NR sequences",seq_count)
                #print("reads_len",len(reads))
        with open(os.path.join(outfolder,batchfilename))as h:
            for id, seq in itertools.zip_longest(*[h] * 2):
                #print(id, reads)
                id=id.replace('\n','')
                id=id.replace('>','')
                #id=int(inter_id.replace('',''))
                seq = seq.replace('\n', '')
                all_batch_sequences[id]=seq
        #TODO: got a key error here. Check this script and correct what is wrong
        #print("Keys")
        #print(batch_mappings.keys())
        #print(batch_reads.keys())
        #print(len(batch_mappings),len(batch_reads))
        for key, value in batch_reads.items():
            if key in batch_mappings:
                reads=batch_mappings[key]
                #print("READS",reads)
                read_mapping=Read(value,reads,False)
                all_infos_dict[batchid][key]=read_mapping
    #print("Len of batchid 0:",len(all_infos_dict[0]))
    #print("Len of batchid 1:", len(all_infos_dict[1]))
    cter=0
    #print("count of input sequences:",str(seq_count))
    #print(read_mapping)
    for batchid,id_dict in all_infos_dict.items():
        for batchid2, id_dict2 in all_infos_dict.items():
            if not batchid2 <= batchid:# and not batchid2==batchid:
                #print("bid",batchid,"bid2",batchid2)
                for id,infos in id_dict.items():
                    if not infos.merged:
                        for id2, infos2 in id_dict2.items():
                            if not infos2.merged:
                                cter+=1
                                consensus1=infos.sequence
                                consensus2=infos2.sequence
                                #print("ID",id)
                                #print(consensus1, "\n")
                                #print("ID2",id2)
                                #print(consensus2, "\n")
                                good_to_pop=align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
                                #print(good_to_pop)
                                if good_to_pop:
                                    if len(infos.reads) > 50:
                                        if len(infos.sequence)>=len(infos2.sequence):
                                            infos.reads=infos.reads+infos2.reads
                                            infos2.merged=True
                                        else:
                                            infos2.reads=infos2.reads +infos.reads
                                            infos.merged=True
                                        #print("Merging",batchid,"_",id," and ",batchid2,"_",id2)
                                    else:
                                        #print("Else")
                                        # TODO generate consensus and add all infos to longer read id
                                        if len(infos.sequence) >= len(
                                                infos2.sequence):

                                            mappings1 = infos.reads
                                            mappings2 = infos2.reads

                                            new_consensus=generate_consensus_path(work_dir,mappings1,mappings2, all_batch_sequences,max_seqs_to_spoa)
                                            #print("CONS1",all_infos_dict[batchid][id])
                                            infos.sequence = new_consensus
                                            infos.reads = infos.reads + infos2.reads
                                            infos2.merged = True
                                        else:
                                            #print("BM",batch_mappings)
                                            mappings1 = infos2.reads
                                            mappings2 = infos.reads
                                            #print("R1",reads1)
                                            #print("M1",mappings1)
                                            #print("R2", reads2)
                                            new_consensus=generate_consensus_path(work_dir,mappings1,mappings2, all_batch_sequences,max_seqs_to_spoa)
                                            #print("CONS2",all_infos_dict[batchid2][id2])
                                            infos2.sequence = new_consensus
                                            infos2.reads = infos2.reads + infos.reads
                                            infos.merged = True

    #for m_batch in merged_batches
    #print("MB",merged_batches," ",len(merged_batches))
    #print(all_infos_dict)
    #print("Combi count ",cter)
    print("Writing file")
    consensus_name = "cluster_merged.fastq"
    other_consensus_name="cluster_merged_low_abundance.fastq"
    mapping_name="cluster_mapping.txt"
    other_mapping_name="cluster_mapping_low_abundance.txt"
    support_name = "support.txt"
    other_support_name = "support_low_abundance.txt"

    support_file = open(os.path.join(outfolder, support_name), "w")
    other_support_file = open(os.path.join(outfolder, other_support_name), "w")
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    other_consensus =open(os.path.join(outfolder, other_consensus_name), 'w')
    mapping_file = open(os.path.join(outfolder, mapping_name), 'w')
    other_mapping=open(os.path.join(outfolder, other_mapping_name), 'w')
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            #print(id," ",all_infos_dict[batchid][id].merged)
            if not all_infos_dict[batchid][id].merged:
                new_id = str(batchid) + "_" + str(id)
                if len(all_infos_dict[batchid][id].reads) >= iso_abundance:
                    mapping_file.write(">{0}\n{1}\n".format(new_id,all_infos_dict[batchid][id].reads))
                    consensus_file.write("@{0}\n{1}\n+\n{2}\n".format(new_id, all_infos_dict[batchid][id].sequence,"+" * len(all_infos_dict[batchid][id].sequence)))
                    support_file.write("{0}: {1}\n".format(new_id, len(all_infos_dict[batchid][id].reads)))
                else:
                    other_consensus.write("@{0}\n{1}\n+\n{2}\n".format(new_id, all_infos_dict[batchid][id].sequence,"+" * len(all_infos_dict[batchid][id].sequence)))
                    other_mapping.write(">{0}\n{1}\n".format(new_id,all_infos_dict[batchid][id].reads))
                    other_support_file.write("{0}: {1}\n".format(new_id, len(all_infos_dict[batchid][id].reads)))
    consensus_file.close()
    mapping_file.close()
    other_consensus.close()
    other_mapping.close()
def main():
    outfolder = "100kSIRV/test8"
    delta=0.10
    delta_len=5
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