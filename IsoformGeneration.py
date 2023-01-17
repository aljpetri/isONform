import _pickle as pickle
import collections
from sys import stdout
import networkx as nx
import subprocess
import os, sys
from EqualityObject import *
import itertools
import tempfile
from consensus import *
import matplotlib.pyplot as plt
from collections import namedtuple
from GraphGeneration import *
def draw_Graph(DG):
    # defines the graph layout to use spectral_layout. Pos gives the position for each node
    pos = nx.spectral_layout(DG)
    # draws the network in the given layout. Only nodes are labeled, the edges are directed but not labeled
    nx.draw_networkx(DG, pos, font_weight='bold')
    # add labels for the edges and draw them into the graph
    # labels = nx.get_edge_attributes(DG, 'weight')
    # nx.draw_networkx_edge_labels(DG,pos, edge_labels=labels)
    plt.show()
"""Method to delete read information from nodes
INPUT:      DG                  Directed Graph
            node                node for which read information is deleted
            supported_reads     list of reads which we want to delete from the nodes 
OUPUT:      reads:              dictionary, which does not contain the read information anymore
"""
def remove_reads_from_node(DG,node,supported_reads):
    #print("Node",node)
    reads=DG._node[node]['reads']
    for read in supported_reads:
        if read in reads.keys():
            del reads[read]
    return reads
"""Method to delete read information from edges
INPUT:      DG                  Directed Graph
            edge                edge for which read information is deleted
            supported_reads     list of reads which we want to delete from the nodes 
OUPUT:      reads:              dictionary, which does not contain the read information anymore
"""
def remove_reads_from_edge(DG,edge,supported_reads):

    if DEBUG:
        print("DEBUG")
        print("Edge", edge)
        #print(DG.edges())
        #print(list(DG.nodes.data()))
        #if 120 in supported_reads:
        #    print("Removing 120")
    reads = DG[edge[0]][edge[1]]['edge_supp']
    for read in supported_reads:
        if read in reads:
            reads.remove(read)
        else:
            if DEBUG:
                print("Strange")
        return reads

"""
Method to make sure that an isoform only contains reads which do actually end with this node

INPUT   edgelist            list of edges starting at current node
        DG                  Networkx DigraphObject
        supported_reads     List of reads which support the path up to this point

OUTPUT: supported_reads:    List of reads which support the path and do not have any further nodes before t
"""
def subtract_wrong_reads(edgelist,supported_reads,DG):
    #supported_reads_t=supported_reads.copy()
    #print("Subtracting wrong reads from ")
    #print(supported_reads)
    reads_to_remove=[]
    for edge in edgelist:
        other_node = edge[1]
        if not other_node == 't':
            #print(other_node)
            other_node_reads = DG._node[other_node]['reads']
            #print(other_node_reads)
            for read in supported_reads:
                #print(read)
                # if a read is in both the current node and the subsequent node we are currently looking at
                if read in other_node_reads.keys():
                    #print("Deleting read ")
                    #print(read)
                    reads_to_remove.append(read)
    #print("Reads to remove:")
    #print(reads_to_remove)
    for remread in reads_to_remove:
        if remread in supported_reads:
            supported_reads.remove(remread)
    #print("Supported reads after")
    #print(supported_reads)
    return supported_reads

"""Method to generate the final isoforms by iterating through the graph structure
INPUT:      DG          Directed Graph
            support       list of reads 
OUPUT:      filename    file which contains all the final isoforms
"""
def compute_equal_reads2(DG,support):
    #path_and_support will hold the infos concerning the found paths
    node_support_left=set(support)
    visited_nodes_isoforms={}
    isoforms={}
    all_supp=set(support)
    #we iterate as long as still not all support was allocated to a path
    while node_support_left:
        #print("Node_support_left",node_support_left)
        node = "s"
        #current_node_support = node_support_left
        read=node_support_left.pop()
        current_node_support = node_support_left
        current_node_support.add(read)
        visited_nodes = []
        #As long as we have not visited the bubble end node we continue walking through our graph
        while node!="t":
            visited_nodes.append(node)
            out_edges=DG.out_edges(node)
            next_found=False
            for edge in out_edges:
                #if DEBUG:
                    #print("edge",edge)
                edge_supp = DG[edge[0]][edge[1]]['edge_supp']
                if read in edge_supp:
                    node=edge[1]
                    current_node_support=current_node_support.intersection(edge_supp)
                    next_found=True
                    break
            if not next_found:
                break

        if current_node_support:
                id = list(current_node_support)[0]
                isoforms[id]=list(current_node_support)
                node_support_left-=current_node_support
                visited_nodes_isoforms[id]=visited_nodes
        else:
            print("no current_node_support")
    print("Found ",len(isoforms)," isoforms")
    return isoforms,visited_nodes_isoforms



    """calls spoa and returns the consensus sequence for the given reads"""

def run_spoa(reads, spoa_out_file, spoa_path):
    with open(spoa_out_file, "w") as output_file:
        #print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call(["spoa", reads, "-l", "0", "-r", "0", "-g", "-2"],
                                  stdout=output_file, stderr=null)
        #print('Done.')
        stdout.flush()
    #output_file.close()
    l = open(spoa_out_file, "r").readlines()
    output_file.close()
    consensus = l[1].strip()
    del l

    null.close()
    return consensus

    # retreives the sequences which are to be aligned using spoa and writes the consensus into a file
    """
        curr_best_seqs is an array with q_id, pos1, pos2
        the current read is on index 0 in curr_best_seqs array
    """
def generate_isoform_using_spoa(curr_best_seqs,reads, work_dir,outfolder,batch_id, max_seqs_to_spoa=200,iso_abundance=1):
    print("Generating the Isoforms")
    mapping = {}
    consensus_name="spoa"+str(batch_id)+"merged.fasta"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    seq_counter=0
    mapping_cter = 0
    for key,value in curr_best_seqs.items():
       # print("CBS",curr_best_seqs)
    #for equalreads in curr_best_seqs:
        name = 'consensus' + str(value[0])
        #name = 'consensus' + str(equalreads[0])
        mapping[name] = []
        reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
        #if len(equalreads) == 1:
        #seq_counter+=len(value)

        if len(value)>=iso_abundance:
            if len(value) == 1:
                seq_counter+=1
            #rid = equalreads[0]
                rid = key
                singleread = reads[rid]
            #print(singleread)
                seq = singleread[1]
            # name='consensus'+str(rid)
                mapping[name].append(singleread[0])
                consensus_file.write(">{0}\n{1}\n".format(name, seq))
                reads_path.close()
            else:
            #print("Equalreads has different size")
            #for i, q_id in enumerate(equalreads):
                for i, q_id in enumerate(value):
                    seq_counter+=1
                    singleread = reads[q_id]
                    seq = singleread[1]
                #print(seq)
                    mapping[name].append(singleread[0])
                    #if i > max_seqs_to_spoa:
                    #    break
                #print("read ", q_id,": ",seq)
                    if i<max_seqs_to_spoa:
                        reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
                reads_path.close()
                spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
            #print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
                consensus_file.write(">{0}\n{1}\n".format(name, spoa_ref))
        #else:
            #print("NOADD",key,value)
    #print(mapping)

    #print("Mapping has length "+str(len(mapping)))

        #print("Mapping",mapping)
    mapping_name="mapping"+str(batch_id)+".txt"
    mappingfile = open(os.path.join(outfolder, mapping_name), "w")
        #print(mapping)
    for id, readlist in mapping.items():

            mapping_cter+=len(readlist)
            if len(readlist)>=iso_abundance:
                mappingfile.write("{0}\n{1}\n".format(id, readlist))

    mappingfile.close()
    # consensus_file.write(">{0}\n{1}\n".format('consensus', spoa_ref))

    consensus_file.close()
    #print(seq_counter," sequences, ",mapping_cter," mappings")


def generate_consensuses(curr_best_seqs,reads,id,id2,work_dir,max_seqs_to_spoa,called_consensuses):

    consensus_infos = {}
    consensuses={}
    #print("CBS",curr_best_seqs)
    if not (id in called_consensuses):
        consensus_infos[id]=curr_best_seqs[id]
        #print("ConInfos id", consensus_infos[id])
    else:
        consensuses[id]=called_consensuses[id]

    if not (id2 in called_consensuses):
        consensus_infos[id2]=curr_best_seqs[id2]
    else:
        consensuses[id2] = called_consensuses[id2]

    #print("CInfos",consensus_infos)
    if bool(consensus_infos):
        for key,value in consensus_infos.items():
            reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
            #print("k",key,"vv",value)
            if len(value) == 1:
                rid = key
                singleread = reads[rid]
                seq = singleread[1]
                consensuses[key]=seq
                #consensus_file.write(">{0}\n{1}\n".format(name, seq))
                reads_path.close()
            else:
                for i, q_id in enumerate(value):
                    singleread = reads[q_id]
                    seq = singleread[1]
                    if i < max_seqs_to_spoa:
                    #print("read ", q_id,": ",seq)
                        reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))

                reads_path.close()
                spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
                    #print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
                consensuses[key] = spoa_ref
    return consensuses
    #print(mapping)
def parse_cigar_diversity_isoform_level(cigar_tuples, delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,overall_len):
    #print("Overall_length",overall_len)
    miss_match_length = 0
    alignment_len = 0
    # print("Now we are parsing....")
    #print(cigar_tuples)
    too_long_indel = False
    three_prime=True
    five_prime=True
    this_start_pos=0
    max_pos=0
    for i, elem in enumerate(cigar_tuples):
        this_start_pos=alignment_len
        max_pos=i
        cig_len = elem[0]
        #print("cigar_len", cig_len)
        cig_type = elem[1]
        alignment_len += cig_len
        if (cig_type != '=') and (cig_type != 'M'):
            # we want to add up all missmatches to compare to sequence length
            miss_match_length += cig_len
            if cig_len <= delta_len:
                continue
            # we also want to make sure that we do not have too large internal sequence differences

            #if the user wants us to merge the 5'end
            #make sure that the startposition of the cigar tuple is in delta_iso_len
            if merge_sub_isoforms_5 and this_start_pos <delta_iso_len_5:
                #startpos=this_start_pos+cig_len
                #print("startpos",startpos)
                #we only want to merge if the cigar tuple does not span into the read by more than delta_len + delta_iso_len_5 base pairs
                if this_start_pos+cig_len<delta_iso_len_5+delta_len:
                    continue
            #we now want that the cigar string starts at most delta_iso_len+delta_len base pairs before the end of the read
            if merge_sub_isoforms_3 and this_start_pos > ((overall_len- delta_iso_len_3)-delta_len):
                continue
            #the user wanted us to neither merge at 3' nor at 5', but we found a cigar tuple entry longer than delta_len indicating a missmatch
            return False

    #We calculate the diversity of our alignment
    diversity = (miss_match_length / alignment_len)
    max_bp_diff = max(delta * alignment_len, delta_len)
    mod_div_rate = max_bp_diff / alignment_len
    #print("3'",three_prime," 5'",five_prime)

    #print("diversity", diversity, "mod_div", mod_div_rate)
    #we additionally make sure that the two consensuses are not too diverse
    diversity_bool=diversity <= mod_div_rate
    #if any of the three parameters we look at tells us not to merge we do not merge
    if diversity_bool:  # delta_perc:
        #print("Div:",diversity_bool,"3' ",three_prime," 5' ",five_prime)
        return True
    else:
        #print("no pop due to diversity")
        return False
def get_overall_alignment_len(cigar_tuples):
    overall_len=0
    for i, elem in enumerate(cigar_tuples):
        overall_len += elem[0]
    return overall_len
def align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5):
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(consensus1, consensus2,
                                                                                       match_score=2,
                                                                                       mismatch_penalty=-2,
                                                                                       opening_penalty=12, gap_ext=1)
    overall_len=get_overall_alignment_len(cigar_tuples)
    #print(cigar_string)
    #print(s1_alignment)
    #print(s2_alignment)
    #print(cigar_tuples)
    #print(overall_len)
    good_to_pop = parse_cigar_diversity_isoform_level(cigar_tuples, delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,overall_len)
    #print(good_to_pop)
    return good_to_pop



def generate_new_full_consensus(id,id2,reads,curr_best_seqs,work_dir,max_seqs_to_spoa):
    #print("again best seqs")
    #print(curr_best_seqs)
    consensus_infos1 = curr_best_seqs[id]
    consensus_infos2 = curr_best_seqs[id2]
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    all_consensus_infos=consensus_infos1+consensus_infos2
    if len(all_consensus_infos) > max_seqs_to_spoa:
        all_consensus_infos=all_consensus_infos[0:max_seqs_to_spoa]
    #print("ACI",all_consensus_infos)
    for q_id in all_consensus_infos:
        singleread = reads[q_id]
        seq = singleread[1]

        # print("read ", q_id,": ",seq)
        reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
    reads_path.close()
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
    return spoa_ref
def generate_all_consensuses(all_consensuses,alternative_consensuses,curr_best_seqs,reads,work_dir,max_seqs_to_spoa):
    print("Generating")
    for id,seqs in curr_best_seqs.items():
        reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
        # print("k",key,"vv",value)
        if len(seqs) == 1:
            #print(id)
            rid = id
            singleread = reads[rid]
            seq = singleread[1]
            #all_consensuses[id] = seq
            consensus_tuple = (id, seq)
            alternative_consensuses.append(consensus_tuple)
            # consensus_file.write(">{0}\n{1}\n".format(name, seq))
            reads_path.close()
        else:
            for i, q_id in enumerate(seqs):
                singleread = reads[q_id]
                seq = singleread[1]
                if i < max_seqs_to_spoa:
                    # print("read ", q_id,": ",seq)
                    reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
            reads_path.close()
            spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
            # print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
            all_consensuses[id] = spoa_ref
            consensus_tuple=(id,spoa_ref)
            alternative_consensuses.append(consensus_tuple)
def add_merged_reads(curr_best_seqs, id2,id1):
    consensus_to_update=curr_best_seqs[id2]
    consensus_updates=curr_best_seqs[id1]
    #print("TOUPDATE:",consensus_to_update)
    #print("UPDATES:",consensus_updates)
    new_consensus=consensus_to_update+consensus_updates
    curr_best_seqs[id2]=new_consensus
    curr_best_seqs.pop(id1)
"""
This method is used to find out how similar the consensuses are (by figuring out how much of their intervals are shared.
INPUT: isoform_paths: List object of all nodes visited for each isoform
        outfolder:         The folder in which we want to write our infos
        batch_id:           The id of the current batch to be able to tell apart the infos of different batches
    The method produces two files:  A similarity file giving the similarity values for each pair of consuensuses.
                                    A path file giving the paths of all final isoforms
"""
def merge_consensuses(curr_best_seqs,work_dir,isoform_paths,outfolder,delta,delta_len,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,reads,max_seqs_to_spoa,delta_iso_len_3,delta_iso_len_5):
    #print("merging consensuses")
    new_best_seqs={}
    #print("best_seqs")
    #print(curr_best_seqs)
    #for key,value in reads.items():
    #    print(key,value)
    #merge_set=set()
    new_consensuses={}
    all_consensuses={}
    alternative_consensuses=[]
    generate_all_consensuses(all_consensuses,alternative_consensuses,curr_best_seqs,reads,work_dir,max_seqs_to_spoa)
    alternative_consensuses.sort(key=lambda x: len(x[1]))
    #for consensus in alternative_consensuses:
        #print(consensus[0],len(consensus[1]))
    for i, consensus in enumerate(alternative_consensuses[:len(alternative_consensuses)-1]):
        #print("con",consensus)
        id1=consensus[0]
        #if id1 in merge_set:
        #    continue
        if not id1 in new_consensuses:
            seq1=consensus[1]
        else:
            seq1=new_consensuses[id1]
        for consensus2 in alternative_consensuses[i+1:]:
            id2=consensus2[0]
            #if id2 in merge_set:
            #    continue
            #if not id1==id2: #todo get rid of this line
            #print(id1,id2)
            seq2=consensus2[1]
            #as soon as we have a length difference larger than delta_iso_len_3+delta_iso_len_5 we break out of the inner loop
            if len(seq2)-len(seq1)>delta_iso_len_3+delta_iso_len_5:
                break
            consensus1 = seq1
            consensus2 = seq2
            merge_consensuses = align_to_merge(consensus1, consensus2, delta, delta_len, merge_sub_isoforms_3,
                                                   merge_sub_isoforms_5, delta_iso_len_3, delta_iso_len_5)
            if merge_consensuses:
                #print("WEMERGE")
                #merge_set.add(id2)
                first_consensus=[item for item in alternative_consensuses if item[0] == id2][0]
                #second_consensus=[item for item in alternative_consensuses if item[0] == id2][0]
                if len(curr_best_seqs[id2]) > 50:
                    print("merge id1 into id2",id," ", id2)
                    add_merged_reads(curr_best_seqs, id2,id1)
                    new_consensuses[id2]=first_consensus[1]
                else:
                    new_consensuses[id2]=generate_new_full_consensus(id1,id2,reads,curr_best_seqs,work_dir,max_seqs_to_spoa)
                    add_merged_reads(curr_best_seqs, id2, id1)
                break

    for id,support in curr_best_seqs.items():
        if not id in new_consensuses:
            seq=reads[id][1]
            new_consensuses[id]=seq
    return new_consensuses



def generate_isoforms_new(equal_reads, reads, work_dir, outfolder, batch_id, max_seqs_to_spoa,iso_abundance,new_consensuses):
    print("Generating the Isoforms-merged")
    mapping = {}
    consensus_name = "spoa" + str(batch_id) + "merged.fasta"
    support_name = "support" + str(batch_id) + ".txt"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    support_file = open(os.path.join(outfolder, support_name), 'w')
    seq_counter = 0
    other_mapping_cter = 0
    # we iterate over all items in curr_best_seqs
    for key, value in equal_reads.items():
        # print(key,value)
        # print(len(value))
        seq_counter += len(value)
        name = 'consensus' + str(key)
        # We generate a list as value for mapping[key]
        if name not in mapping:
                    mapping[name] = []
                    # iterate over all reads in value and add them to mapping
        for i, q_id in enumerate(value):
                    # print("Q_ID",q_id)
                    singleread = reads[q_id]
                    mapping[name].append(singleread[0])
                # we do not have to calculate consensus as already done
        sequence = new_consensuses[key]
        consensus_file.write(">{0}\n{1}\n".format(name, sequence))
        #consensus_file.write("@{0}\n{1}\n+\n{2}\n".format(name, sequence, "+" * len(sequence)))

        other_mapping_cter = 0
        for key, value in mapping.items():
            other_mapping_cter += len(value)

        # we have to treat the mapping file generation a bit different as we have to collect also the supporting reads of the subconsensuses
    mapping_name = "mapping" + str(batch_id) + ".txt"
    mappingfile = open(os.path.join(outfolder, mapping_name), "w")
    mapping_cter = 0
    for id, mapped_read_list in mapping.items():
            mapping_cter += len(mapped_read_list)
            mappingfile.write("{0}\n{1}\n".format(id, mapped_read_list))
            support_file.write("{0}: {1}\n".format(id, len(mapped_read_list)))

    consensus_file.close()



"""
Wrapper method used for the isoform generation
"""
def generate_isoforms(DG,all_reads,reads,work_dir,outfolder,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,delta_iso_len_3,delta_iso_len_5,iso_abundance,max_seqs_to_spoa=200):
    #print("s", DG.nodes["s"]['reads'])
    #print("t", DG.nodes["t"]['reads'])
    #print("Importante")
    #print(merge_sub_isoforms_3)
    #print(merge_sub_isoforms_5)
    equal_reads,isoform_paths=compute_equal_reads2(DG,reads)
    #equal_reads_name='equal_reads_'+str(batch_id)+'.txt'
    #print("EQUALS",equal_reads)
    #print("BID",batch_id)
    #print("s",DG.nodes["s"]['reads'])
    #print("t",DG.nodes["t"]['reads'])
    #with open(os.path.join(outfolder, equal_reads_name), 'wb') as file:
    #    file.write(pickle.dumps(equal_reads))
    #isoform_paths_name='isoform_paths_'+str(batch_id)+'.txt'
    #with open(os.path.join(outfolder, isoform_paths_name), 'wb') as file:
    #    file.write(pickle.dumps(isoform_paths))
    #with open('all_reads.txt', 'wb') as file:
    #    file.write(pickle.dumps(all_reads))
        #old_outfolder=os.path.join(outfolder,"out")
    #calculate_isoform_similarity(equal_reads,work_dir,isoform_paths,outfolder,delta,delta_len,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,all_reads,max_seqs_to_spoa,delta_iso_len_3=0,delta_iso_len_5=0)
    #generate_isoform_using_spoa(equal_reads,all_reads, work_dir,outfolder,batch_id, max_seqs_to_spoa,iso_abundance)
    #print(merge_sub_isoforms_3)
    #merge_sub_isos_5=(merge_sub_isoforms_5=='True')
    #merge_sub_isos_3 = (merge_sub_isoforms_3 == 'True')
    #print("DO WE MERGE?",merge_sub_isoforms_5,merge_sub_isoforms_3)
    if merge_sub_isoforms_5 or merge_sub_isoforms_3:
        #print("MergingTrue")
        #print("EQUALS",equal_reads)
        #for c_id, supp_reads in equal_reads.items():
            #for supp_read in supp_reads:
                #print("{0}_{1}: {2} ".format(c_id,supp_read, all_reads[supp_read]))

        new_consensuses=merge_consensuses(equal_reads, work_dir, isoform_paths, outfolder, delta, delta_len, batch_id,
                                 merge_sub_isoforms_3, merge_sub_isoforms_5, all_reads, max_seqs_to_spoa,
                                 delta_iso_len_3, delta_iso_len_5)

        generate_isoforms_new(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa,
                                           iso_abundance,new_consensuses)
    else:
        generate_isoform_using_spoa(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance)


DEBUG=True

def main():
    outfolder="100kSIRV/20722_abundance2_par/22"
    batch_id = 0
    print("Hello World")
    file = open(os.path.join(outfolder,'equal_reads_'+str(batch_id)+'.txt'), 'rb')
    equal_reads = pickle.load(file)
    file.close()
    file = open(os.path.join(outfolder,'isoform_paths_'+str(batch_id)+'.txt'), 'rb')
    isoform_paths = pickle.load(file)
    file.close()

    file2 = open(os.path.join(outfolder, 'all_reads_' + str(batch_id)+".txt" ), 'rb')
    #file2 = open(os.path.join(outfolder,'all_reads_'+str(batch_id)+'.txt'), 'rb')
    all_reads = pickle.load(file2)
    file2.close()
    work_dir = tempfile.mkdtemp()
    outfolder = "out_local"
    delta=0.10
    delta_len=5

    merge_sub_isoforms_3=True
    merge_sub_isoforms_5 = True
    delta_iso_len_3=30
    delta_iso_len_5 = 50
    max_seqs_to_spoa=200
    iso_abundance=1

    #if merge_sub_isoforms_5 or merge_sub_isoforms_3:
    #    merged_dict,merged_consensuses,called_consensuses,consensus_map = calculate_isoform_similarity(equal_reads, work_dir, isoform_paths, outfolder, delta, delta_len, batch_id,
    #                             merge_sub_isoforms_3, merge_sub_isoforms_5, all_reads, max_seqs_to_spoa,
    #                             delta_iso_len_3, delta_iso_len_5)
    #    generate_isoforms_new()
    #    generate_isoform_using_spoa_merged(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance,merged_dict,merged_consensuses,called_consensuses,consensus_map)
    #else:
    #    generate_isoform_using_spoa(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance)

if __name__ == "__main__":
    main()
    """
    main concept: 
    while reads_for_isoforms:
    start at s and iterate through the nodes up to t.
    Greedy approach:
    get nextnode by getting all out_edges.
    follow the edge which has
    max(current_node['reads']and nextnode['reads'])
    meaning which has the maximum amount of elements which are also in current_node
    all reads that are not in nextnode are deleted from supported_reads
    One isoform:
    Supported_reads as soon as edge[1]==t
    Solve isoform by using spoa.
    """



