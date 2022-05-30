import _pickle as pickle
import collections
from sys import stdout
import networkx as nx
import subprocess
import os, sys
from EqualityObject import *
import itertools
import matplotlib.pyplot as plt
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
    #print("cleaning the graph")
    #print("Visited Edges",visited_edges)
    #print("Visited nodes",visited_nodes)

    #print("supported_reads",supported_reads)
    #all_edges_dict={}
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
        isoforms[supported_reads[0]]=supported_reads
        #reads_for_isoforms[supported_reads[0]]=visited_nodes
        visited_nodes_for_isoforms[supported_reads[0]]=visited_nodes
        #print(reads_for_isoforms)
        for sup_read in supported_reads:
            #print(sup_read)
            reads_for_isoforms.remove(sup_read)
        #print("Isoforms")
        #print(isoforms)
        if DEBUG==True:
            print("VisitedNodes" , supported_reads[0],"has ", len(visited_nodes)," elements: ",visited_nodes)
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

    mapping = {}
    consensus_name="spoa"+str(batch_id)+".fa"
    consensus_file = open(os.path.join(outfolder, consensus_name), 'w')

    for key,value in curr_best_seqs.items():
    #for equalreads in curr_best_seqs:
        name = 'consensus' + str(value[0])
        #name = 'consensus' + str(equalreads[0])
        mapping[name] = []
        reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
        #if len(equalreads) == 1:
    #TODO: add information of how many reads support this isoform by consensusx_support
        if len(value)>=iso_abundance:
            if len(value) == 1:
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
                    singleread = reads[q_id]
                    seq = singleread[1]
                #print(seq)
                    mapping[name].append(singleread[0])
                    if i > max_seqs_to_spoa:
                        break
                #print("read ", q_id,": ",seq)
                    reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
                reads_path.close()
                spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
            #print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
                consensus_file.write(">{0}\n{1}\n".format(name, spoa_ref))
    #print(mapping)

    #print("Mapping has length "+str(len(mapping)))

        #print("Mapping",mapping)
        mapping_name="mapping"+str(batch_id)+".txt"
        mappingfile = open(os.path.join(outfolder, mapping_name), "w")
        for id, seq in mapping.items():
            if len(seq)>=iso_abundance:
                mappingfile.write("{0}\n{1}\n".format(id, seq))
    mappingfile.close()
    # consensus_file.write(">{0}\n{1}\n".format('consensus', spoa_ref))

    consensus_file.close()
    # for i, (q_id, pos1, pos2) in  enumerate(grouper(curr_best_seqs, 3)):
    #    seq = reads[q_id][1][pos1: pos2 + k_size]
    #    if i > max_seqs_to_spoa:
    #        break
    #    reads_path.write(">{0}\n{1}\n".format(str(q_id), seq))
    # reads_path.close()
    #print("Isoforms generated")
def generate_consensuses(curr_best_seqs,reads,id,id2,work_dir,max_seqs_to_spoa):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    consensus_infos = {}
    consensuses={}
    consensus_infos[id]=curr_best_seqs[id]
    consensus_infos[id2]=curr_best_seqs[id2]
    for key,value in consensus_infos.items():
            if len(value) == 1:
                rid = key
                singleread = reads[rid]
                seq = singleread[1]
                consensuses[key]=seq
                #consensus_file.write(">{0}\n{1}\n".format(name, seq))
                #reads_path.close()
            else:
                for i, q_id in enumerate(value):
                    singleread = reads[q_id]
                    seq = singleread[1]

                    if i > max_seqs_to_spoa:
                        break
                #print("read ", q_id,": ",seq)
                    reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
                reads_path.close()
                spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
            #print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
                consensuses[key]=spoa_ref
    return consensuses
    #print(mapping)
def parse_cigar_diversity_isoform_level(cigar_tuples, delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5):
    miss_match_length = 0
    alignment_len = 0
    # print("Now we are parsing....")
    # print(cigar_tuples)
    too_long_indel = False
    three_prime=False
    five_prime=False
    for i, elem in enumerate(cigar_tuples):

        cig_len = elem[0]
        cig_type = elem[1]
        alignment_len += cig_len
        if (cig_type != '=') and (cig_type != 'M'):
            # we want to add up all missmatches to compare to sequence length
            miss_match_length += cig_len
            # we also want to make sure that we do not have too large internal sequence differences
            if cig_len > 2 * delta_len:
                if i == 1:
                    if merge_sub_isoforms_3:
                        if cig_len<delta_iso_len_3:
                            three_prime=True
                else:
                    return False
                if i == (len(cigar_tuples) - 1):
                    if merge_sub_isoforms_5:
                        if cig_len < delta_iso_len_5:
                            five_prime=True
                else:
                    # print("ELE",elem)
                    # print("No pop due to delta_len")
                    return False
    diversity = (miss_match_length / alignment_len)

    max_bp_diff = max(delta * alignment_len, delta_len)
    mod_div_rate = max_bp_diff / alignment_len
    if DEBUG:
        print("diversity", diversity, "mod_div", mod_div_rate)
    if diversity <= mod_div_rate and three_prime and five_prime:  # delta_perc:
        return True
    else:
        # print("no pop due to diversity")
        return False
def align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5):
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(consensus1, consensus2,
                                                                                       match_score=2,
                                                                                       mismatch_penalty=-8,
                                                                                       opening_penalty=12, gap_ext=1)
    good_to_pop = parse_cigar_diversity_isoform_level(cigar_tuples, delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
    if good_to_pop:
        print("YESS")
"""
This method is used to find out how similar the consensuses are (by figuring out how much of their intervals are shared.
INPUT: isoform_paths: List object of all nodes visited for each isoform
        outfolder:         The folder in which we want to write our infos
        batch_id:           The id of the current batch to be able to tell apart the infos of different batches
    The method produces two files:  A similarity file giving the similarity values for each pair of consuensuses.
                                    A path file giving the paths of all final isoforms
"""
def calculate_isoform_similarity(curr_best_seqs,work_dir,isoform_paths,outfolder,delta,delta_len,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,reads,max_seqs_to_spoa,delta_iso_len_3=0,delta_iso_len_5=0):
    consensus_name = "inter_spoa" + str(batch_id) + ".fa"
    called_consensuses={}
    consensus_file_inter = open(os.path.join(outfolder, consensus_name), 'w')
    eq_file_name="similarity_batch_"+str(batch_id)+".txt"
    path_file_name="path_"+str(batch_id)+".txt"
    similarity_file = open(os.path.join(outfolder, eq_file_name), 'w')
    path_file=open(os.path.join(outfolder,path_file_name),"w")
    mapping={}
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
                if equality_2>equality_1:
                    if equality_2>0.7:
                        similarity_file.write("{0} subisoform of {1} indicated by {3}\n".format(id2, id,equality_2))
                        res=False
                        #the following code figures out whether path2 is a subisoform of path by finding out whether path2 is a sublist of path
                        for idx in range(len(path) - len(path2) + 1):
                            if path[idx: idx + len(path2)] == path2:
                                res = True
                                similarity_file.write("Sublist TRUE\n")
                                break
                        if res:

                            consensuses=generate_consensuses(curr_best_seqs,reads,id,id2,work_dir,max_seqs_to_spoa)
                            consensus1=consensuses[id]
                            consensus2=consensuses[id2]
                            merge_consensuses=align_to_merge(consensus1,consensus2,delta,delta_len,merge_sub_isoforms_3,merge_sub_isoforms_5,delta_iso_len_3,delta_iso_len_5)
                            if not merge_consensuses:
                                called_consensuses[id]=consensus1
                                called_consensuses[id2]=consensus2
                            else:
                                if len(reads[id])>50:


                #    if equality_1==1.0:
                #        similarity_file.write("{0} subisoform of {1}\n".format(id, id2))
                #    min_equality=equality_2
                #    max_equality=equality_1
                if DEBUG:
                    print("Equality of ",id," vs ",id2,": ",equality_1)
                similarity_file.write("{0},{1}: {2} \n".format(id, id2,equality_1))

DEBUG=False
"""
Wrapper method used for the isoform generation
"""
def generate_isoforms(DG,all_reads,reads,work_dir,outfolder,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,max_seqs_to_spoa=200,iso_abundance=1,delta_iso_len_3=0,delta_iso_len_5=0):
    equal_reads,isoform_paths=compute_equal_reads(DG,reads)
    calculate_isoform_similarity(equal_reads,work_dir,isoform_paths,outfolder,delta,delta_len,batch_id,merge_sub_isoforms_3,merge_sub_isoforms_5,all_reads,max_seqs_to_spoa,delta_iso_len_3=0,delta_iso_len_5=0)
    generate_isoform_using_spoa(equal_reads,all_reads, work_dir,outfolder,batch_id, max_seqs_to_spoa,iso_abundance)








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



