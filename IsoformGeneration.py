import _pickle as pickle
import collections
from sys import stdout
import networkx as nx
import subprocess
import os, sys
from EqualityObject import *
import itertools
import matplotlib.pyplot as plt
"""Method to delete read information from nodes
INPUT:      DG                  Directed Graph
            node                node for which read information is deleted
            supported_reads     list of reads which we want to delete from the nodes 
OUPUT:      reads:              dictionary, which does not contain the read information anymore
"""
def remove_reads_from_node(DG,node,supported_reads):
    reads=DG.nodes[node]['reads']
    print("Before")
    print(reads)
    for read in supported_reads:
        if read in reads.keys():
            del reads[read]
        else:
            print("ERROR: "+str(read)+" is not in "+str(reads))
    return reads
def remove_edges_for_node(DG,node):
    print("Removing edges from graph")
def clean_graph(DG,visited_nodes,supported_reads):
    print("cleaning the graph")
    print("Visited nodes")
    print(visited_nodes)
    print("supported_reads")
    print(supported_reads)
    update_dict={}
    for node in visited_nodes:
        new_reads=remove_reads_from_node(DG, node, supported_reads)
        if new_reads:
            update_dict[node]=new_reads
            nx.set_node_attributes(DG,update_dict,'reads')
            new_reads=DG.nodes[node]['reads']
            print("after")
            print(new_reads)
        else:
            print("Removing node "+str(node))
            DG.remove_node(node)
"""
Method to make sure that an isoform only contains reads which do actually end with this node

INPUT   edgelist            list of edges starting at current node
        DG                  Networkx DigraphObject
        supported_reads     List of reads which support the path up to this point

OUTPUT: supported_reads:    List of reads which support the path and do not have any further nodes before t
"""
def subtract_wrong_reads(edgelist,supported_reads,DG):
    #supported_reads_t=supported_reads.copy()
    print("Subtracting wrong reads from ")
    print(supported_reads)
    reads_to_remove=[]
    for edge in edgelist:
        other_node = edge[1]
        if not other_node == 't':
            print(other_node)
            other_node_reads = DG._node[other_node]['reads']
            print(other_node_reads)
            for read in supported_reads:
                print(read)
                # if a read is in both the current node and the subsequent node we are currently looking at
                if read in other_node_reads.keys():
                    print("Deleting read ")
                    print(read)
                    reads_to_remove.append(read)
    print("Reads to remove:")
    print(reads_to_remove)
    for remread in reads_to_remove:
        if remread in supported_reads:
            supported_reads.remove(remread)
    print("Supported reads after")
    print(supported_reads)
    return supported_reads
"""
Method to find the edge, which is supported by the maximum amount of nodes, used to tell which node we look into next

INPUT   DG                  Networkx DigraphObject
        current_node        The current start node
        supported_reads     List of reads which support the path up to this point
        delta_len           Parameter to distinguish isoforms by length differences

OUTPUT: next_node:          the node which has the maximum support
        support_list:       List of supporting reads
"""
def get_best_supported_edge_node(DG,current_node,supported_reads,delta_len):
    edgelist = list(DG.out_edges(current_node))

    #print(current_node)
    print("Edgelist")
    print(edgelist)
    #print("initial supported reads")
    #print(supported_reads)
    next_node=None
    curr_node_reads = DG._node[current_node]['reads']

    #print("CurrnodeREads")
    #print(curr_node_reads)
    next_node_eq = []
    #iterate over all possible next nodes (other_node)
    for edge in edgelist:
        other_node = edge[1]
        one_next_node=[]
        #print("OTHERNode")
        #print(other_node)
        if not other_node == 't':
            #print(other_node)

            other_node_reads = DG._node[other_node]['reads']

            #print(other_node_reads)
            #iterate over all reads which make up the path up to current_node
            for read in supported_reads:
                #print(read)
                #if a read is in both the current node and the subsequent node we are currently looking at
                if read in other_node_reads.keys():
                    #print("Current node"+current_node)
                    start_tup = curr_node_reads[read]
                    end_tup = other_node_reads[read]
                    len=end_tup[0]-start_tup[1]
                    len_found=False
                    if one_next_node:
                        for eq_obj in one_next_node:
                            cur_len=eq_obj.get_length()
                            if abs(len-cur_len)< delta_len:
                                eq_obj.add_supported_read(read)
                                eq_obj.increment_eq()
                                #print(eq_obj)
                                #eq_obj.print_readlist()
                                len_found=True
                        if not len_found:
                            eq_obj = EqualityObject(other_node,len, 0, read)
                            one_next_node.append(eq_obj)
                    else:
                        eq_obj = EqualityObject(other_node,len,0,read)
                        one_next_node.append(eq_obj)
        else:
            supported_reads_t=supported_reads.copy()
            supported_reads_t=subtract_wrong_reads(edgelist,supported_reads_t,DG )
            print(supported_reads_t)
            if supported_reads_t:
                print("Next node: t")
                next_node = "t"
                return(next_node,supported_reads_t)
        for node in one_next_node:
            next_node_eq.append(node)
    max_support=-1
    #print("Next node eq")
    for eq_obj in next_node_eq:
        #print(eq_obj)
        #print("max_support "+str(max_support))
        cur_equality=eq_obj.get_equality()
        #print("cur_eq"+str(cur_equality))
        if cur_equality>max_support:
            max_support=cur_equality
            supported_reads=eq_obj.get_support()
            next_node=eq_obj.get_node()
    #print("Sup")
    #print(supported_reads)
    return (next_node,supported_reads)
"""Method to generate the final isoforms by iterating through the graph structure
INPUT:      DG          Directed Graph
            reads       list of reads 
OUPUT:      filename    file which contains all the final isoforms
"""
def compute_equal_reads(DG,reads,delta_len):
    startnode = 's'
    startreads=DG._node['s']['reads']
    #print("Startreads")
    #print(startreads)
    #endnode='t'
    supported_reads=[]
    reads_for_isoforms=reads
    isoforms= {}
    #while still reads have to be assigned to an isoform
    while(reads_for_isoforms):
        current_node=startnode
        supported_reads=reads_for_isoforms
        reached_t=False
        visited_nodes = []
        #While the end of the read was not reached, iterate over subsequent nodes
        while(not reached_t):
            #add current node to the list of visited_nodes
            visited_nodes.append(current_node)
            #print("supporting_edge")
            #print(supporting_edge)
            #print("CurrnodebefMethod")
            #print(current_node)
            current_node,supported_reads=get_best_supported_edge_node(DG,current_node,supported_reads,delta_len)
            if current_node=="t":
                reached_t=True
            #print(support_list)
            #print("after")
            #print(current_node)
            #if(support_list):
            #supported_reads=list(support_list)
            #print("Current Node: "+current_node)
        print("Cleaning graph")
        clean_graph(DG,visited_nodes,supported_reads)
        isoforms[supported_reads[0]]=supported_reads
        print(reads_for_isoforms)
        for sup_read in supported_reads:
            print(sup_read)
            reads_for_isoforms.remove(sup_read)
        print("Isoforms")
        print(isoforms)
        print("VisitedNodes")
        print(visited_nodes)
    return isoforms

    """calls spoa and returns the consensus sequence for the given reads"""

def run_spoa(reads, spoa_out_file, spoa_path):
    with open(spoa_out_file, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call(["/home/alexanderpetri/spoa/build/bin/spoa", reads, "-l", "0", "-r", "0", "-g", "-2"],
                                  stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    del l
    output_file.close()
    null.close()
    return consensus

    # generates the sequences which are to be aligned using spoa and writes the consensus into a file
    """
        curr_best_seqs is an array with q_id, pos1, pos2
        the current read is on index 0 in curr_best_seqs array
    """
def generate_isoform_using_spoa(curr_best_seqs,reads, work_dir,outfolder, max_seqs_to_spoa=200):


    mapping = {}
    consensus_file = open(os.path.join(outfolder, "spoa.fa"), 'w')

    #curr_best_seqs = curr_best_seqs[0:3]
    for key,value in curr_best_seqs.items():
    #for equalreads in curr_best_seqs:
        name = 'consensus' + str(value[0])
        #name = 'consensus' + str(equalreads[0])
        mapping[name] = []
        reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
        #if len(equalreads) == 1:
        if len(value) == 1:
            #rid = equalreads[0]
            rid = key
            singleread = reads[rid]
            # print(singleread)
            seq = singleread[1]
            # name='consensus'+str(rid)
            mapping[name].append(singleread[0])
            consensus_file.write(">{0}\n{1}\n".format(name, seq))
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
                reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
            reads_path.close()
            spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
            #print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
            consensus_file.write(">{0}\n{1}\n".format(name, spoa_ref))
    # print(mapping)

    # print("Mapping has length "+str(len(mapping)))
    mappingfile = open(os.path.join(outfolder, "mapping.txt"), "w")
    for id, seq in mapping.items():
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
    print("Isoforms generated")


"""
Wrapper method used for the isoform generation
"""
def generate_isoforms(DG,all_reads,reads,work_dir,outfolder,delta_len,max_seqs_to_spoa=200):
    equal_reads=compute_equal_reads(DG,reads,delta_len)
    generate_isoform_using_spoa(equal_reads,all_reads, work_dir,outfolder, max_seqs_to_spoa=200)








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



