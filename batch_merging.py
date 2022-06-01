def merge_batches(max_batchid, outfolder,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,max_seqs_to_spoa,iso_abundance, delta_iso_len_3, delta_iso_len_5):
    all_reads_dict={}#dictionary holding all reads from all the batches
    for batchid in range(0,max_batchid+1):
        filename="spoa" + str(batchid) + ".fa"
        print("File: ",filename)
        #file = open("spoa" + str(batchid) + ".fa", 'r')
        import itertools
        with open(filename) as f:
            for line1, line2 in itertools.zip_longest(*[f] * 2):
                print(line1, line2)
                print("Seq_len",len(line2))


def main():
    outfolder = "out_local"
    delta=0.10
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