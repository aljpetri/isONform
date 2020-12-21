import _pickle as pickle
import collections
from sys import stdout
import networkx as nx
import itertools
import matplotlib.pyplot as plt
"""Method to generate the final isoforms by iterating through the graph structure
INPUT:      DG          Directed Graph
            reads       list of reads 
OUPUT:      filename    file which contains all the final isoforms
"""
def generate_isoforms(DG,reads):
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
            print("Current Node: "+current_node)
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
    print("Final Isoforms:")
    print(isoforms)








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

"""function to merge consecutive nodes, if they contain the same reads to simplify the graph
    INPUT: DG Directed Graph
    OUTPUT: DG: Directed Graph with merged nodes
"""
def merge_nodes(DG):
    #iterate over the edges to find all pairs of nodes
    edgesView=DG.edges.data()
    for ed_ge in edgesView:
        startnode = ed_ge[0]
        endnode = ed_ge[1]
        #we only need to know the out degree of the start node and the end degree of the end node
        start_out_degree = DG.out_degree(startnode)
        end_in_degree = DG.in_degree(endnode)
        #if the degrees are both equal to 1 and if none of the nodes is s or t
        if(start_out_degree==end_in_degree==1and startnode!="s"and endnode!="t"):
            #print("Merging nodes "+startnode+" and "+endnode)
            #use the builtin function to merge nodes, prohibiting self_loops decreases the amount of final edges
            DG=nx.contracted_nodes(DG, startnode,endnode,self_loops=False)

