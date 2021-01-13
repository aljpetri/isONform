import _pickle as pickle
import collections
from sys import stdout
import networkx as nx
import subprocess
import os, sys
import itertools
import matplotlib.pyplot as plt
"""Method to generate the final isoforms by iterating through the graph structure
INPUT:      DG          Directed Graph
            reads       list of reads 
OUPUT:      filename    file which contains all the final isoforms
"""
def compute_equal_reads(DG,reads):
    startnode = 's'
    #startreads=DG._node['s']['reads']
    #print(startreads)
    #endnode='t'
    supported_reads=[]
    reads_for_isoforms=reads
    visited_nodes=[]
    isoforms= {}
    #while still reads have to be assigned to an isoform
    while(reads_for_isoforms):
        current_node=startnode
        supported_reads=reads_for_isoforms
        reached_t=False
        #While the end of the read was not reached iterate over subsequent nodes
        while(not reached_t):
            edgelist=list(DG.out_edges(current_node))
            #print(edgelist)
            maximum_support = 0
            #supporting_edge=None
            #current_node_reads=DG.nodes[current_node]['reads']
            visited_nodes.append(current_node)
            for edge in edgelist:
                #print(edge)
                equality_counter=0
                other_node=edge[1]
                if not other_node=='t':
                    other_node_reads=DG.nodes[other_node]['reads']
                    #print("other node reads")
                    #print(other_node_reads)
                    support_list=[]
                    #print("Supported Reads")
                    #print(supported_reads)
                    for read in supported_reads:
                        read_id=read
                        #print("ReadID:"+str(read_id))
                        read_id_both_reads=[item for item in other_node_reads if item[0] == read_id]
                        if read_id_both_reads:
                            equality_counter=equality_counter+1
                            #print(equality_counter)
                            support_list.append(read_id)

                    if equality_counter>maximum_support:
                        maximum_support=equality_counter
                        supporting_edge=edge
                        #print("supporting_edge")
                        #print(supporting_edge)
                        supported_reads = support_list
                else:
                    current_node="t"
                    reached_t=True
                    break
            #print("supporting_edge")
            #print(supporting_edge)
            current_node=supporting_edge[1]
            #print("Current Node: "+current_node)
        isoforms[supported_reads[0]]=supported_reads
        print(reads_for_isoforms)
        for sup_read in supported_reads:
            print(sup_read)
            reads_for_isoforms.remove(sup_read)
            #startreads.remove(sup_read)
        #for vis_node in visited_nodes:
        #    node_att_prev=DG._node[vis_node]["reads"]
        #    print("Node att prev")
        #    print(node_att_prev)
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


"""Wrapper method used for the isoform generation
"""
def generate_isoforms(DG,all_reads,reads,work_dir,outfolder,max_seqs_to_spoa=200):
    equal_reads=compute_equal_reads(DG,reads)
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



