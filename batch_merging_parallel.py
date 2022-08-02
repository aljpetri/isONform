import itertools
from consensus import *
from recordclass import recordclass
from IsoformGeneration import align_to_merge
import tempfile
import pickle
from Parallelization_side_functions import *
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
def read_spoa_file(batch_id,cl_dir):
    filename ="spoa"+str(batch_id)+"merged.fa"
    batch_reads={}
    with open( os.path.join(cl_dir, filename)) as f:
        for id, sequence in itertools.zip_longest(*[f] * 2):
            # print(id, sequence)
            inter_id = id.replace('\n', '')
            id = int(inter_id.replace('>consensus', ''))
            sequence = sequence.replace('\n', '')
            # print("Seq_len",len(sequence))
            batch_reads[id] = sequence
    return batch_reads
def read_mapping_file(batch_id,cl_dir):
    batch_mappings = {}

    mappingname =  "mapping"+str(batch_id)+".txt"
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
            batch_mappings[id] = readslist
    return batch_mappings
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
"""The last step of our algorithm: We take all isoforms that were generated for each batch and align them to merge same isoforms
"""
def merge_batches(max_batchid,work_dir, outfolder,all_reads,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa, delta_iso_len_3, delta_iso_len_5,iso_abundance):
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

        filename="spoa" + str(batchid) + "merged.fa"
        batchfilename="0_batchfile.fa"
        mappingname= "mapping0.txt"
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

        #print(batch_mappings)
        for key, value in batch_reads.items():
            reads=batch_mappings[key]
            #print("READS",reads)
            read_mapping=Read(value,reads,False)
            all_infos_dict[batchid][key]=read_mapping
    #print("Len of batchid 0:",len(all_infos_dict[0]))
    #print("Len of batchid 1:", len(all_infos_dict[1]))
    cter=0
    print("count of input sequences:",str(seq_count))
    #print(read_mapping)
    for batchid,id_dict in all_infos_dict.items():
        for batchid2, id_dict2 in all_infos_dict.items():
            if not batchid2 <= batchid:# and not batchid2==batchid:
                print("bid",batchid,"bid2",batchid2)
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
                                    if len(infos.reads) > 50:
                                        if len(infos.sequence)>=len(infos2.sequence):
                                            infos.reads=infos.reads+infos2.reads
                                            infos2.merged=True
                                        else:
                                            infos2.reads=infos2.reads +infos.reads
                                            infos.merged=True
                                        print("Merging",batchid,"_",id," and ",batchid2,"_",id2)
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
    print("Combi count ",cter)
    print("Writing file")
    consensus_name = "cluster_merged.fa"
    other_consensus_name="cluster_merged_low_abundance.fa"
    mapping_name="cluster_mapping.txt"
    other_mapping_name="cluster_mapping_low_abundance.txt"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    other_consensus =open(os.path.join(outfolder, other_consensus_name), 'w')
    mapping_file = open(os.path.join(outfolder, mapping_name), 'w')
    other_mapping=open(os.path.join(outfolder, other_mapping_name), 'w')
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            print(id," ",all_infos_dict[batchid][id].merged)
            if not all_infos_dict[batchid][id].merged:
                new_id = str(batchid) + "_" + str(id)
                if len(all_infos_dict[batchid][id].reads) >= iso_abundance:
                    mapping_file.write(">{0}\n{1}\n".format(new_id,all_infos_dict[batchid][id].reads))
                    consensus_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
                else:
                    other_consensus.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
                    other_mapping.write(">{0}\n{1}\n".format(new_id,all_infos_dict[batchid][id].reads))
    consensus_file.close()
    mapping_file.close()
    other_consensus.close()
    other_mapping.close()
#TODO:write skipped files into overall output
def write_final_output(all_infos_dict,outfolder,iso_abundance):
    print("Writing file")
    consensus_name = "cluster_merged.fa"
    other_consensus_name = "cluster_merged_low_abundance.fa"
    mapping_name = "cluster_mapping.txt"
    other_mapping_name = "cluster_mapping_low_abundance.txt"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    other_consensus = open(os.path.join(outfolder, other_consensus_name), 'w')
    mapping_file = open(os.path.join(outfolder, mapping_name), 'w')
    other_mapping = open(os.path.join(outfolder, other_mapping_name), 'w')
    for batchid, id_dict in all_infos_dict.items():
        for id, infos in id_dict.items():
            #print(id, " ", all_infos_dict[batchid][id].merged)
            if not all_infos_dict[batchid][id].merged:
                new_id = str(batchid) + "_" + str(id)
                if len(all_infos_dict[batchid][id].reads) >= iso_abundance:
                    mapping_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].reads))
                    consensus_file.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
                else:
                    other_consensus.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].sequence))
                    other_mapping.write(">{0}\n{1}\n".format(new_id, all_infos_dict[batchid][id].reads))
    consensus_file.close()
    mapping_file.close()
    other_consensus.close()
    other_mapping.close()
#TODO: add the rest of variables for this method and move filewriting to wrappermethod
def actual_merging_process(all_infos_dict,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,all_batch_sequences,work_dir):
    cter = 0
    for batchid,id_dict in all_infos_dict.items():
        for batchid2, id_dict2 in all_infos_dict.items():
            if not batchid2 <= batchid:# and not batchid2==batchid:
                print("bid",batchid,"bid2",batchid2)
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
                                    if len(infos.reads) > 50:
                                        if len(infos.sequence)>=len(infos2.sequence):
                                            infos.reads=infos.reads+infos2.reads
                                            infos2.merged=True
                                        else:
                                            infos2.reads=infos2.reads +infos.reads
                                            infos.merged=True
                                        print("Merging",batchid,"_",id," and ",batchid2,"_",id2)
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
    print("Combi count ", cter)
#TODO: Finish implementation of this method and replace call of other merge_function by call to this
def join_back_via_batch_merging(tmp_work_dir, outdir, split_mod, residual,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,iso_abundance):
    print("Batch Merging")
    Read = recordclass('Read', "sequence reads merged")
    print(outdir, tmp_work_dir)
    unique_cl_ids = set()
    all_infos_dict={}
    subfolders = [f.path for f in os.scandir(outdir) if f.is_dir()]
    for folder in subfolders:
        print("Folder",folder)
        file = os.fsdecode(folder)
        # print(file)
        tmp_fname = folder.split('/')
        fname=tmp_fname[-1]
        print("FNAME",fname)

        cl_id = fname # file.split('_')
        unique_cl_ids.add(cl_id)
        print(unique_cl_ids)
        for cl_id in unique_cl_ids:
            batch_reads = {}
            batch_mappings = {}
            all_reads_dict = {}
            out_pattern=os.path.join(outdir,"isoform")
            cl_dir = os.path.join(outdir, cl_id)
            print("CLDIR",cl_dir)
            for batchfile in os.listdir(cl_dir):
                print(batchfile)
                if batchfile.endswith("_batchfile.fa"):
                    print(batchfile)
                    tmp_bname = batchfile.split('/')
                    #print(tmp_bname)
                    tmp_bname2=tmp_bname[-1].split('_')

                    batch_id=int(tmp_bname2[0])
                    # print(type(tmp_work_dir), type(cl_id))
                    #batches_pattern = os.path.join(os.fsdecode(outdir), cl_id + '_*')
                    # print("joining all", out_pattern, "from", batches_pattern)
                    mkdir_p(out_pattern)
                    #print("BPattern",batches_pattern)
                    error_file = open(os.path.join(out_pattern, 'cat.stderr'), 'w')
                    outfilename = os.path.join(out_pattern, 'isoforms.fastq')

                    #with open(outfilename, 'wb') as outfile:

                    #TODO:open needed files, perform the batch_merging as before
                        #tmp_file_name=file_id.split("/")
                        #tmp_imp_infos=tmp_file_name[-1].split("_")
                        #clu_id=tmp_imp_infos[0]
                        #batch_id=tmp_imp_infos[1]
                        #print(file_id)
                        #print(batch_id)
                        #TODO:write all reads' infos into one big dictionary
                    #print("Reading", batch_id)
                    read_batch_file(batch_id, all_infos_dict, all_reads_dict, cl_dir)
                    batch_reads=read_spoa_file(batch_id,cl_dir)
                    #print("BatchReads",batch_reads)
                    batch_mappings=read_mapping_file(batch_id,cl_dir)
                    #print("ALL INFOS",all_infos_dict)
            for key, value in batch_reads.items():
                            #print(key,value)
                            reads = batch_mappings[key]
                            # print("READS",reads)
                            read_mapping = Read(value, reads, False)
                            #print("Infos",batch_id,key)
                            all_infos_dict[batch_id][key] = read_mapping
        #filename = os.path.join(str(cl_id), 'corrected_reads.fastq')
        #if filename == outfilename:
                            # don't want to copy the output into the output
        #    continue
        #with open(filename, 'rb') as readfile:
        #                print("Hello World")
        actual_merging_process(all_infos_dict,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,max_seqs_to_spoa,all_reads_dict,outdir)
        write_final_output(all_infos_dict,outdir,iso_abundance)
                        #shutil.rmtree(batch_id)
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