def merge_batches(max_batchid, outfolder,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa,iso_abundance, delta_iso_len_3, delta_iso_len_5):
    all_reads_dict={}#dictionary holding all reads from all the batches
    all_mappings_dict={}
    for batchid in range(0,max_batchid+1):
        batch_reads={}
        batch_mappings={}
        filename="spoa" + str(batchid) + ".fa"
        mappingname= "mapping" + str(batchid) + ".txt"
        print("File: ",filename)
        #file = open("spoa" + str(batchid) + ".fa", 'r')
        import itertools
        all_reads_dict[batchid]={}
        all_mappings_dict[batchid]={}
        with open(filename) as f:
            for id, sequence in itertools.zip_longest(*[f] * 2):
                print(id, sequence)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('>consensus',''))

                sequence = sequence.replace('\n', '')
                print("Seq_len",len(sequence))
                batch_reads[id]=sequence
        with open(mappingname) as g:
            for id, reads in itertools.zip_longest(*[g] * 2):
                print(id, reads)
                inter_id=id.replace('\n','')
                id=int(inter_id.replace('>consensus',''))

                reads = reads.replace('\n', '')
                batch_mappings[id]=reads
                print("reads_len",len(reads))

    print(all_reads_dict)
def main():
    outfolder = "out_local"
    delta=0.05
    k_size=20
    delta_len=2*k_size
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