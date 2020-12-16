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
    startnode = "s"
    endnode="t"
    supported_reads=[]
    reads_for_isoforms=reads
    visited_nodes=[]
    while(reads_for_isoforms):
        current_node=startnode
        reached_t=False
        while(not reached_t):
            edgelist=list(DG.out_edges(current_node))
            print(edgelist)
            maximum_support=0
            supporting_edge=None
            current_node_reads=DG.nodes[current_node]['reads']
            for edge in edgelist:
                other_node=ed
                other_node_reads=DG.nodes[current_node]['reads']


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

