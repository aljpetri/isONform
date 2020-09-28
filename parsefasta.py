
import fastaparser

class Read(object):
    id = ""
    description = " "
    sequence = ""

    # The class "constructor" - It's actually an initializer
    def __init__(self, id, description, sequence):
        self.id = id
        self.description = description
        self.sequence = sequence

def make_read(id, description, sequence):
    read = Read(id, description, sequence)
    return read
def read_fasta_and_remove_deletions():
    with open("/home/alexanderpetri/Desktop/RAWDATA_PhD1/RBMY.fa") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        newparser=[]
        #for seq in parser:
            # seq is a FastaSequence object
            #print('ID:', seq.id)
            #print('Description:', seq.description)
            #print('Sequence:', seq.sequence_as_string())
            #print()
        for seq in parser:
            prevseq=seq.sequence_as_string()
            newstr = prevseq.replace("-", "")

            read=make_read(seq.id,seq.description,newstr)
            newparser.append(read)
    print("Deletions from reads were removed")
    return newparser


