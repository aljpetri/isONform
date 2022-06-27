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
    #print("Edge",edge)
    reads = DG[edge[0]][edge[1]]['edge_supp']
    #print(reads)
    #draw_Graph(DG)
    for read in supported_reads:
        if read in reads:
            reads.remove(read)
    return reads
"""Method to delete nodes and edges which do are not supported by any reads anymore
INPUT:      DG                  Directed Graph
            visitee_nodes       Nodes which make up an isoform and from which we delete the reads
            visited_edges       Edges which make up an isoform and from which we delete the reads
            supported_reads     list of reads which we want to delete from the nodes 
OUPUT:      reads:              dictionary, which does not contain the read information anymore
"""
def clean_graph(DG,visited_nodes,visited_edges,supported_reads):
    update_dict={}
    #edge_update_dict={}
    for edge in visited_edges:
        #print("visited_edge:",edge)
        new_reads = remove_reads_from_edge(DG, edge, supported_reads)
        #print("New_reads",new_reads)
        if new_reads:
            edge_tuple=(edge[0],edge[1])
            update_dict[edge_tuple]=new_reads
        else:
            #print("removing edge",edge)
            DG.remove_edge(edge[0],edge[1])
    nx.set_node_attributes(DG, update_dict, 'reads')
    for node in visited_nodes:
        if DG.degree(node)==0:
            DG.remove_node(node)
            #print("removing node", node)
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
"""
Method to find the edge, which is supported by the maximum amount of nodes, used to tell which node we look into next

INPUT   DG                  Networkx DigraphObject
        current_node        The current start node
        supported_reads     List of reads which support the path up to this point
        edge_attr           Dict holding all edges as key and their respective support as values

OUTPUT: next_node:          the node which has the maximum support
        support_list:       List of supporting reads
"""
def get_best_supported_edge_node(DG,current_node,supported_reads,edge_attr):
    edgelist = list(DG.out_edges(current_node))
    #print("now at",current_node)
    #print("Edgelist",edgelist)
    #print("initial supported reads",supported_reads)
    final_support=[]
    similarity_val=0
    next_node=""
    #iterate over all possible next nodes (other_node)
    for edge in edgelist:
        #print("current node")
        #print(current_node)
        supp_reads=supported_reads
        #print("Initial supp")
        #print(supp_reads)
        edge_reads=edge_attr[edge]
        #print("edge_reads",edge_reads)
        shared_reads=list(set(supp_reads).intersection(edge_reads))
        #print("Shared REads",shared_reads)
        if len(shared_reads)>similarity_val:
                #print("SIM")
                #print(similarity_val)
                similarity_val=len(shared_reads)
                #print(similarity_val)
                final_support=shared_reads
                next_node=edge[1]
    return (next_node,final_support)
"""Method to generate the final isoforms by iterating through the graph structure
INPUT:      DG          Directed Graph
            reads       list of reads 
OUPUT:      filename    file which contains all the final isoforms
"""
def compute_equal_reads(DG,reads):
    startnode = 's'
    startreads=DG._node['s']['reads']
    #print("Startreads",startreads)
    #endnode='t'
    visited_nodes_for_isoforms= {}
    supported_reads=[]
    reads_for_isoforms=reads
    isoforms= {}
    edge_attr=nx.get_edge_attributes(DG,"edge_supp")
    #print("EdgeAttri")
    #print(edge_attr)
    #while still reads have to be assigned to an isoform
    while(reads_for_isoforms):
        #print("RFI",reads_for_isoforms)
        current_node=startnode

        supported_reads=reads_for_isoforms
        reached_t=False
        visited_nodes = []
        visited_edges = []
        #While the end of the read was not reached, iterate over subsequent nodes
        while(not reached_t):
            #add current node to the list of visited_nodes
            visited_nodes.append(current_node)
            prev_node=current_node
            #print("CurrnodebefMethod",current_node)
            #print()
            current_node,supported_reads=get_best_supported_edge_node(DG,current_node,supported_reads,edge_attr)
            #print("current node returned by get best supported edge node", current_node)
            edge_tup=(prev_node,current_node)
            #print("edge_tup",edge_tup)
            visited_edges.append(edge_tup)
            if not supported_reads:
                break
            #print("Supported:")

            if current_node=="t":
                reached_t=True
            #print("Still supported")
            #print(*supported_reads)
            #print("after")
            #print(current_node)
            #if(support_list):
            #supported_reads=list(support_list)
            #print("Current Node: "+current_node)
        #print("Cleaning graph")
        #print("visited_edges:",visited_edges)
        #print(DG.edges(data=True))
        #print("Supported", supported_reads)
        clean_graph(DG,visited_nodes,visited_edges,supported_reads)
        #print(DG.edges(data=True))
        id=supported_reads[0]
        isoforms[id]=supported_reads
        #reads_for_isoforms[supported_reads[0]]=visited_nodes
        visited_nodes_for_isoforms[id]=visited_nodes
        #print(reads_for_isoforms)
        for sup_read in supported_reads:
            #print(sup_read)
            reads_for_isoforms.remove(sup_read)
        #print("Isoforms")
        #print(isoforms)
        if DEBUG==True:
            print("VisitedNodes" , id,"has ", len(visited_nodes)," elements: ",visited_nodes)
        #print(visited_nodes)
    return isoforms,visited_nodes_for_isoforms

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
    consensus_name="spoa"+str(batch_id)+".fa"
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
    print(seq_counter," sequences, ",mapping_cter," mappings")

def generate_isoform_using_spoa_merged(curr_best_seqs, reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance,merged_dict,merged_consensuses,called_consensuses,consensus_map):
    print("Generating the Isoforms-merged")
    mapping = {}
    consensus_name = "spoa" + str(batch_id) + "merged.fa"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')
    seq_counter=0
    other_mapping_cter=0
    #we iterate over all items in curr_best_seqs
    for key, value in curr_best_seqs.items():
        #print(key,value)
        #print(len(value))
        seq_counter += len(value)
        name = 'consensus' + str(key)
        print("Sequence count:",seq_counter)
        print("Mapping count:",other_mapping_cter)
        #print(key,value)
        #print(key,len(value))
        #print(merged_dict)
        #if our key is in merged_consensuses, we know that we merged this consensus during get_isoform_similarity
        if key in merged_consensuses:
            #If the key is in merged_dict i.e.  key is still a consensus id
            if key in merged_dict:
                #We generate a list as value for mapping[key]
                if name not in mapping:
                    mapping[name] = []
                    #iterate over all reads in value and add them to mapping
                for i, q_id in enumerate(value):
                    #print("Q_ID",q_id)
                    singleread = reads[q_id]
                    mapping[name].append(singleread[0])
                #we do not have to calculate consensus as already done
                sequence=merged_dict[key].id_seq
                consensus_file.write(">{0}\n{1}\n".format(name, sequence))
            else:
                name = 'consensus' + str(consensus_map[key])
                for i, q_id in enumerate(value):
                    #print("Q_ID", q_id)
                    singleread = reads[q_id]
                    if name not in mapping:
                        mapping[name] = []

                    mapping[name].append(singleread[0])
        #there is no else case because we skip keys that are no merged_dict keys
        else:
            # This is exactly what we do in the normal generate_isoform_using_spoa

            # name = 'consensus' + str(equalreads[0])
            reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
            if name not in mapping:
                mapping[name]=[]
            if len(value) >= iso_abundance:
                if len(value) == 1:
                    # rid = equalreads[0]
                    rid = key
                    #print("R_ID", rid)
                    singleread = reads[rid]
                    # print(singleread)
                    seq = singleread[1]
                    mapping[name].append(singleread[0])
                    consensus_file.write(">{0}\n{1}\n".format(name, seq))
                    reads_path.close()
                else:
                    # print("Equalreads has different size")
                    # for i, q_id in enumerate(equalreads):
                    for i, q_id in enumerate(value):
                        #print("Q_IDs", q_id)
                        singleread = reads[q_id]
                        seq = singleread[1]
                        mapping[name].append(singleread[0])
                        if i < max_seqs_to_spoa:
                            reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
                    reads_path.close()
                    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
                    # print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
                    consensus_file.write(">{0}\n{1}\n".format(name, spoa_ref))
            #else:
                #print("HW")
            # print(mapping)

            #print("Mapping has length "+str(len(mapping)))
        other_mapping_cter=0
        for key, value in mapping.items():
            other_mapping_cter+=len(value)
        # we have to treat the mapping file generation a bit different as we have to collect also the supporting reads of the subconsensuses
    mapping_name = "mapping" + str(batch_id) + ".txt"
    mappingfile = open(os.path.join(outfolder, mapping_name), "w")
    mapping_cter=0
    for id, mapped_read_list in mapping.items():
            if not id in merged_consensuses:
                if len(mapped_read_list) >= iso_abundance:
                    mapping_cter+=len(mapped_read_list)
                    mappingfile.write("{0}\n{1}\n".format(id, mapped_read_list))
            else:
                    if id in merged_dict:
                        full_read_list=[]
                        full_read_list.extend(mapping[id])
                        for cons_id in merged_dict[id].otherIDs:
                            #print("CID",mapping[cons_id])
                            full_read_list.extend(mapping[cons_id])
                            #mapping_cter+=1
                            #single_read=reads[cons_id]
                            #full_read_list.append(single_read[0])
                        mappingfile.write("{0}\n{1}\n".format(id, full_read_list))
    print(seq_counter, " sequences, ", mapping_cter, " mappings")
    mappingfile.close()
    # consensus_file.write(">{0}\n{1}\n".format('consensus', spoa_ref))

    consensus_file.close()

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
    miss_match_length = 0
    alignment_len = 0
    # print("Now we are parsing....")
    # print(cigar_tuples)
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
            # we also want to make sure that we do not have too large internal sequence differences
            if cig_len > delta_len:
                if i == 0:
                    if merge_sub_isoforms_5:
                        if cig_len<delta_iso_len_5:
                            five_prime=True
                        else:
                            five_prime=False
                elif i == (len(cigar_tuples)-1):
                    if merge_sub_isoforms_3:
                        if this_start_pos > (overall_len- delta_iso_len_3):
                            three_prime=True
                        else:
                            three_prime=False
                elif merge_sub_isoforms_3:
                    if this_start_pos > (overall_len - delta_iso_len_3):
                        three_prime = True
                    else:
                        three_prime = False
                else:
                    return False
    print("Tuples",len(cigar_tuples),"max_pos",max_pos)
    diversity = (miss_match_length / alignment_len)
    max_bp_diff = max(delta * alignment_len, delta_len)
    mod_div_rate = max_bp_diff / alignment_len
    print("3'",three_prime," 5'",five_prime)

    print("diversity", diversity, "mod_div", mod_div_rate)
    if diversity <= mod_div_rate and three_prime and five_prime:  # delta_perc:
        return True
    else:
        # print("no pop due to diversity")
        return False
def get_overall_alignment_len(cigar_tuples):
    overall_len=0
    for i, elem in enumerate(cigar_tuples):
        overall_len += elem[0]
    return overall_len
def align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5):
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(consensus1, consensus2,
                                                                                       match_score=2,
                                                                                       mismatch_penalty=-8,
                                                                                       opening_penalty=12, gap_ext=1)
    overall_len=get_overall_alignment_len(cigar_tuples)
    print(cigar_string)
    print(cigar_tuples)
    print(overall_len)
    good_to_pop = parse_cigar_diversity_isoform_level(cigar_tuples, delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5,overall_len)
    print(good_to_pop)
    return good_to_pop



def generate_new_full_consensus(id,id2,reads,curr_best_seqs,work_dir,max_seqs_to_spoa):
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
"""
This method is used to find out how similar the consensuses are (by figuring out how much of their intervals are shared.
INPUT: isoform_paths: List object of all nodes visited for each isoform
        outfolder:         The folder in which we want to write our infos
        batch_id:           The id of the current batch to be able to tell apart the infos of different batches
    The method produces two files:  A similarity file giving the similarity values for each pair of consuensuses.
                                    A path file giving the paths of all final isoforms
"""

def calculate_isoform_similarity(curr_best_seqs,work_dir,isoform_paths,outfolder,delta,delta_len,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,reads,max_seqs_to_spoa,delta_iso_len_3,delta_iso_len_5):
    print("calculating similarity")
    merged_dict={}
    Merged_consensus = namedtuple("Merged_consensus", "id_seq otherIDs")
    #consensus_name = "inter_spoa" + str(batch_id) + ".fa"
    called_consensuses={}
    #consensus_file_inter = open(os.path.join(outfolder, consensus_name), 'w')
    eq_file_name="similarity_batch_"+str(batch_id)+".txt"
    path_file_name="path_"+str(batch_id)+".txt"
    similarity_file = open(os.path.join(outfolder, eq_file_name), 'w')
    path_file=open(os.path.join(outfolder,path_file_name),"w")
    consensus_map={}
    #merged consensuses holds all ids for which a merge was done
    merged_consensuses=set()
    for id, path in isoform_paths.items():
        path_set=set(path)
        path_file.write("{0}: {1}\n".format(id, path))
        l1 = len(path)
        if DEBUG:
            print(id, ": ", path)
        for id2,path2 in isoform_paths.items():
            if not id==id2:
                path2_set=set(path2)
                equal_elements=path_set.intersection(path2_set)
                equal_elements_nr=len(equal_elements)
                l2=len(path2)
                #get equality value for path 1: The equality value is a number between 0 and 1 denoting the portion of shared nodes for each consensus pair
                equality_1=equal_elements_nr/l1
                equality_2=equal_elements_nr/l2
                if equality_2>=equality_1:
                    if DEBUG:
                        print("combi ",id,", ",id2,", eq2: ",equality_2,"\n")
                        # if equality_2 > 0.8:
                    if True:

                        similarity_file.write("{0} subisoform of {1} indicated by {2}\n".format(id2, id, equality_2))
                        if id2 not in merged_consensuses:
                            consensuses=generate_consensuses(curr_best_seqs,reads,id,id2,work_dir,max_seqs_to_spoa,called_consensuses)
                        #print(consensuses)
                            consensus1=consensuses[id]
                            consensus2=consensuses[id2]
                            merge_consensuses=align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
                        #if we do not merge the consensuses, we add them to called_consensuses to be able to retreive them if needed
                            if not merge_consensuses:
                                called_consensuses[id]=consensus1
                                called_consensuses[id2]=consensus2

                            else:
                                # We have to figure out whether id is in merged_dict
                                if not(id in merged_dict):
                                    if len(reads[id]) > 50:
                                        this_merged=Merged_consensus(consensus1,[id2])

                                    else:
                                        full_consensus=generate_new_full_consensus(id,id2,reads,curr_best_seqs,work_dir,max_seqs_to_spoa)
                                        this_merged = Merged_consensus(full_consensus, [id2])
                                    merged_dict[id]=this_merged

                                else:
                                    #print(merged_dict[id].other_IDs)
                                    merged_dict[id].otherIDs.append(id2)
                                    #print(merged_dict[id].other_IDs)
                                merged_consensuses.add(id)
                                merged_consensuses.add(id2)
                                consensus_map[id2]=id
                        elif id2 in merged_dict:
                            id_list_2 = []
                            id_list_2.append(id2)
                            for idee in merged_dict[id2].otherIDs:
                                id_list_2.append(idee)
                            all_mergeable=True
                            for id2_entry in id_list_2:
                                consensuses = generate_consensuses(curr_best_seqs, reads, id, id2_entry, work_dir,
                                                                   max_seqs_to_spoa, called_consensuses)
                                consensus1 = consensuses[id]
                                consensus2 = consensuses[id2_entry]
                                merge_consensuses = align_to_merge(consensus1, consensus2, delta, delta_len,
                                                                   merge_sub_isoforms_3, merge_sub_isoforms_5,
                                                                   delta_iso_len_3, delta_iso_len_5)
                                if not merge_consensuses:
                                    all_mergeable=False
                                    break
                            if all_mergeable:
                                if id in merged_dict:
                                    print("IDLIST2",id_list_2)
                                    merged_dict[id].otherIDs.extend(id_list_2)
                                else:
                                    this_merged = Merged_consensus(consensus1, id_list_2)
                                    merged_dict[id] = this_merged
                                    merged_consensuses.add(id)
                                    merged_dict.pop(id2)
                                for id2_entry in id_list_2:
                                    consensus_map[id2_entry] = id
                if DEBUG:
                    print("Equality of ",id," vs ",id2,": ",equality_1)
                similarity_file.write("{0},{1}: {2} \n".format(id, id2,equality_1))
    path_file.close()
    return merged_dict, merged_consensuses, called_consensuses, consensus_map
DEBUG=False
"""
Wrapper method used for the isoform generation
"""
def generate_isoforms(DG,all_reads,reads,work_dir,outfolder,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,delta,delta_len,delta_iso_len_3,delta_iso_len_5,iso_abundance,max_seqs_to_spoa=200):
    equal_reads,isoform_paths=compute_equal_reads(DG,reads)
    equal_reads_name='equal_reads_'+str(batch_id)+'.txt'
    with open(os.path.join(outfolder, equal_reads_name), 'wb') as file:
        file.write(pickle.dumps(equal_reads))
    isoform_paths_name='isoform_paths_'+str(batch_id)+'.txt'
    with open(os.path.join(outfolder, isoform_paths_name), 'wb') as file:
        file.write(pickle.dumps(isoform_paths))
    with open('all_reads.txt', 'wb') as file:
        file.write(pickle.dumps(all_reads))
    #calculate_isoform_similarity(equal_reads,work_dir,isoform_paths,outfolder,delta,delta_len,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,all_reads,max_seqs_to_spoa,delta_iso_len_3=0,delta_iso_len_5=0)
    #generate_isoform_using_spoa(equal_reads,all_reads, work_dir,outfolder,batch_id, max_seqs_to_spoa,iso_abundance)
    if merge_sub_isoforms_5 or merge_sub_isoforms_3:
        merged_dict,merged_consensuses,called_consensuses,consensus_map = calculate_isoform_similarity(equal_reads, work_dir, isoform_paths, outfolder, delta, delta_len, batch_id,
                                 merge_sub_isoforms_3, merge_sub_isoforms_5, all_reads, max_seqs_to_spoa,
                                 delta_iso_len_3, delta_iso_len_5)
        generate_isoform_using_spoa_merged(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance,merged_dict,merged_consensuses,called_consensuses,consensus_map)
    else:
        generate_isoform_using_spoa(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance)




def main():
    outfolder="100kSIRV/test8"
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
    if merge_sub_isoforms_5 or merge_sub_isoforms_3:
        merged_dict,merged_consensuses,called_consensuses,consensus_map = calculate_isoform_similarity(equal_reads, work_dir, isoform_paths, outfolder, delta, delta_len, batch_id,
                                 merge_sub_isoforms_3, merge_sub_isoforms_5, all_reads, max_seqs_to_spoa,
                                 delta_iso_len_3, delta_iso_len_5)
        generate_isoform_using_spoa_merged(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance,merged_dict,merged_consensuses,called_consensuses,consensus_map)
    else:
        generate_isoform_using_spoa(equal_reads, all_reads, work_dir, outfolder, batch_id, max_seqs_to_spoa, iso_abundance)

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



