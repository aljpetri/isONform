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
    print("Node",node)
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
    print("Edge",edge)
    reads = DG[edge[0]][edge[1]]['edge_supp']
    print(reads)
    for read in supported_reads:
        if read in reads:
            reads.remove(read)
    return reads
"""Method to delete nodes and edges which do not support any reads any more
INPUT:      DG                  Directed Graph
            visitee_nodes       Nodes which make up an isoform and from which we delete the reads
            visited_edges       Edges which make up an isoform and from which we delete the reads
            supported_reads     list of reads which we want to delete from the nodes 
OUPUT:      reads:              dictionary, which does not contain the read information anymore
"""
def clean_graph(DG,visited_nodes,visited_edges,supported_reads):
    #print("cleaning the graph")
    print("Visited Edges",visited_edges)
    print("Visited nodes")
    print(visited_nodes)
    #print("supported_reads")
    #print(supported_reads)
    all_edges_dict={}
    update_dict={}
    edge_update_dict={}
    for edge in visited_edges:
        print("visited_edge:",edge)
        new_reads = remove_reads_from_edge(DG, edge, supported_reads)
    #for node in visited_nodes:
     #   print(node)
      #  new_reads=remove_reads_from_node(DG, node, supported_reads)

        print("New_reads",new_reads)
        if new_reads:
            print("true")
            edge_tuple=(edge[0],edge[1])
            update_dict[edge_tuple]=new_reads
            nx.set_node_attributes(DG,update_dict,'reads')
        else:
            print("removing edge",edge)
            DG.remove_edge(edge[0],edge[1])
    for node in visited_nodes:
        if DG.degree(node)==0:
            DG.remove_node(node)
            print("removing node", node)
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
    #print("now at")
    #print(current_node)
    print("Edgelist")
    print(edgelist)
    #print("initial supported reads")
    #print(supported_reads)
    final_support=[]
    #print("CurrnodeREads")
    #print(curr_node_reads)
    similarity_val=0
    next_node=""
    #iterate over all possible next nodes (other_node)
    for edge in edgelist:
        #print("current node")
        #print(current_node)
        supp_reads=supported_reads
        print("Initial supp")
        print(supp_reads)
        edge_reads=edge_attr[edge]
        print("edge_reads",edge_reads)
        shared_reads=list(set(supp_reads).intersection(edge_reads))
        print("Shared REads")
        print(shared_reads)
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
    #print("Startreads")
    #print(startreads)
    #endnode='t'
    supported_reads=[]
    reads_for_isoforms=reads
    isoforms= {}
    edge_attr=nx.get_edge_attributes(DG,"edge_supp")
    #print("EdgeAttri")
    #print(edge_attr)
    #while still reads have to be assigned to an isoform
    while(reads_for_isoforms):
        print("RFI",reads_for_isoforms)
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
            #print("supporting_edge")
            #print(supporting_edge)
            #print("CurrnodebefMethod")
            #print(current_node)
            current_node,supported_reads=get_best_supported_edge_node(DG,current_node,supported_reads,edge_attr)
            print("current node returned by get best supported edge node", current_node)
            edge_tup=(prev_node,current_node)
            visited_edges.append(edge_tup)
            if not supported_reads:
                break
            #print("Supported:")
            #print(supported_reads)
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
        print("visited_edges:",visited_edges)
        print(DG.edges(data=True))
        clean_graph(DG,visited_nodes,visited_edges,supported_reads)
        print(DG.edges(data=True))
        isoforms[supported_reads[0]]=supported_reads
        #print(reads_for_isoforms)
        for sup_read in supported_reads:
            #print(sup_read)
            reads_for_isoforms.remove(sup_read)
        #print("Isoforms")
        #print(isoforms)
        #print("VisitedNodes")
        #print(visited_nodes)
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
def generate_isoform_using_spoa(curr_best_seqs,reads, work_dir,outfolder, max_seqs_to_spoa=200):
    #print("reads")
    #print(reads)
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
    #print("Isoforms generated")


"""
Wrapper method used for the isoform generation
"""
def generate_isoforms(DG,all_reads,reads,work_dir,outfolder,max_seqs_to_spoa=200):
    equal_reads=compute_equal_reads(DG,reads)
    generate_isoform_using_spoa(equal_reads,all_reads, work_dir,outfolder, max_seqs_to_spoa)








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



