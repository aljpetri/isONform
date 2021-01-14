import networkx as nx
from collections import Counter
from main import draw_Graph
#TODO find out how to find bubbles->implement pop_bubbles
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
    return(DG)

"""
function to find cycles in the graph (, which denote repetitive regions)
INPUT: DG Directed Graph
OUTPUT: List_of_cycles: A list holding all cycles present in DG
"""
def find_repetative_regions(DG):
    # collect the cycles in the graph (denoting repetative regions) using the builtin networkx function
    altcyc = list(nx.simple_cycles(DG))
    print("Alternative cycles:")
    print(altcyc)
    #data structure which holds all the cycles in the graph
    list_of_cycles = []
    #iterate over the cycles to retrive all nodes which are part of the cycles
    for comp in altcyc:
        #if len(comp) > 1:
        intermediate_cycle = []
        #iterate over the nodes in each cycle to get a better output format (start, nodename,end)
        for node_i in comp:
            intermediate = tuple(map(int, node_i.split(', ')))
            print(intermediate)
            # print(type(intermediate))
            intermediate_cycle.append((node_i,intermediate))
            # real_cycles.append(intermediate_sorted)
            # print(intermediate_sorted)
            #if not node_i in cycle_nodes:
            #    cycle_nodes.append(node_i)
            #sort the nodes in a cycle by start coordinates to simplify the resolving of cycles later on
        cycle_sorted = sorted(intermediate_cycle, key=lambda x: x[0])
        #print("Cycle_sorted type")
        #print(str(type(cycle_sorted)))
        #print(cycle_sorted)
        list_of_cycles.append(cycle_sorted)
    if list_of_cycles:
        print("Found repetative region in reads")
        for cyc in list_of_cycles:
            print(cyc)
    else:
        print("No cycles found in the graph")
    return (list_of_cycles)

def find_min_max(contains_s,contains_t,listofnodes):
    if contains_t:
        max_element = "t"
        intermediate_sorted = sorted(listofnodes, key=lambda x: x[0])
        min_element = intermediate_sorted[0]
        print(min_element)
    if contains_s:
        intermediate_sorted = sorted(listofnodes, key=lambda x: x[0])
        max_element = intermediate_sorted[-1]
        min_element = "s"
    if not (contains_t or contains_s):
        intermediate_sorted = sorted(listofnodes, key=lambda x: x[0])
        min_element = intermediate_sorted[0]
        max_element = intermediate_sorted[-1]
    print("intermediate_sorted")
    print(intermediate_sorted)
    return(min_element,max_element)


"""Helper method for find_bubbles: This method figures out which are the different pathsin between min_node and max_node
INPUT:      cycle:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element: The node deemed to be the starting node of the bubble (given as tuple)
            max_element: The node deemed to be the end node of the bubble (given as tuple)
            DG:         the directed graph we want to pop the bubble in
            contains_s: A boolean value denoting whether the cycle contains node "s"
            contains_t: A boolean value denoting whether the cycle contains node "t"

"""
#TODO:correct this method
def get_path_nodes(cycle,min_element,max_element,DG,contains_s,contains_t):
    if contains_s:
        min_node="s"
        min_readlist = DG.nodes[min_node]['reads']
    else:
        min_node = str(min_element[0]) + ", " + str(min_element[1]) + ", " + str(min_element[2])
        inter_min_readlist = DG.nodes[min_node]['reads']
        min_readlist=[i[0] for i in inter_min_readlist]
    if contains_t:
        max_node="t"
        max_readlist = DG.nodes[max_node]['reads']
    else:
        max_node = str(max_element[0]) + ", " + str(max_element[1]) + ", " + str(max_element[2])
        inter_max_readlist = DG.nodes[max_node]['reads']
        max_readlist = [i[0] for i in inter_max_readlist]
    shared_reads=[]
    for read in min_readlist:
        if read in max_readlist:
            shared_reads.append(read)
    min_edges=DG.out_edges(min_node)
    max_edges=DG.in_edges(max_node)
    min_node_out=[]
    max_node_in=[]
    for edge in min_edges:
        if edge[1] in cycle:
            min_node_out.append(edge[1])
    for edge in max_edges:
        if edge[0] in cycle:
            max_node_in.append(edge[0])
    if( len(min_node_out) >2 or len(max_node_in)>2):
        print("Think again!")
        print(min_node_out)
        print(max_node_in)
    else:
        print("Seems correct!")
        print(min_node_out)
        print(max_node_in)
    out_path_reads=[]
    in_path_reads=[]
    #TODO: something in here does not work out correctly yet
    for out_node in min_node_out:
        inter_out_readlist=list(DG.nodes[out_node]['reads'])
        out_readlist = [i[0] for i in inter_out_readlist]
        for read in out_readlist:
            if read in shared_reads:
                out_path_reads.append(out_readlist)
    for in_node in max_node_in:
        inter_in_readlist = list(DG.nodes[in_node]['reads'])
        in_readlist = [i[0] for i in inter_in_readlist]
        for read in in_readlist:
            if read in shared_reads:
                in_path_reads.append(in_readlist)
    print("out_path_reads")
    print(out_path_reads)
    print("in_path_reads")
    print(in_path_reads)
    return min_node_out,max_node_in

#def measure_path_length_differences(min_node_out,max_node_in,delta_len):
    #TODO:until here the whole algo seems to work correctly. Now however, we have to think about missing exons
#    min_path_1_reads=DG.nodes[min_node_out[0]]['reads']
#    print("min_path_1_reads")
#    print(min_path_1_reads)
#    min_path_2_reads=DG.nodes[min_node_out[1]]['reads']
#    print("MinReadlist")
#    print(min_readlist)
#    print("MaxReadlist")
#    print(max_readlist)
""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
            max_bubblesize  parameter giving the maximum size of a bubble to still pop it (TODO: properly define max_bubblesize
            delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """
def find_bubbles(DG,max_bubblesize,delta_len):
    #get undirected version of the graph to find bubbles
    UG=DG.to_undirected()

    #bubble_limit=2*max_bubblesize+2
    #print("Bubble limit: "+str(bubble_limit))
    #collect the bubbles in the graph (Bubbles denote possible mutations in the minimizers)
    list_of_bubbles = []
    #find all cycles in the undirected graph->bubbles
    listofcycles = nx.cycle_basis(UG)
    print(listofcycles)
    #iterate over the different bubbles
    for listint,cycle in enumerate(listofcycles):
        listofnodes=[]
        contains_s=False
        contains_t=False
        #iterate over the nodes
        for node in cycle:
            #treat nodes different if they are s or t (just because of how they were generated)
            if "s" in node:
                s_tuple=("s",0,0)
                contains_s=True

            elif "t" in node:
                contains_t=True
                t_tuple = ("t", 0, 0)
            #for all other nodes the attributes (startnode,endnode,read_id) are recovered from the nodeids
            else:
                mytuple = tuple(map(int, node.split(', ')))
                listofnodes.append(mytuple)
        #find the minimum and maximum for each bubble
        min_element,max_element=find_min_max(contains_s,contains_t,listofnodes)
        #find the nodes
        get_path_nodes(cycle,min_element,max_element,DG,contains_s,contains_t)
        print("Cycle")
        print(cycle)
        print("MinElement")
        print(min_element)
        print("MaxElement")
        print(max_element)
                #print(listofnodes)
            #print (min(listofnodes, key = lambda x: x[0]))
            #print(max(listofnodes, key=lambda x: x[1]))
    #        if not(len(listofnodes)> bubble_limit):
    #            possible_bubble_sorted = sorted(listofnodes, key=lambda x: x[0])
    #            list_of_bubbles.append(possible_bubble_sorted)
    #    else:
    #        listofcycles.pop(listint)

    #for cycle in list_of_cycles:
    #    list_of_bubbles= [x for x in list_of_bubbles if x!= cycle]
    #for bubble in list_of_bubbles:
    #    endtuple=bubble[-1]
    #    starttuple=bubble[0]
    #    endname=str(endtuple[0])+','+str(endtuple[1])
    #    startname = str(starttuple[0]) + ',' + str(starttuple[1])
    #    list_of_bubbles_sliced=bubble[1:-1]
        #for bub_node in list_of_bubbles_sliced:
        #readlistandpos=DG.nodes[startname]['reads']
        #print(readlistandpos)
   # print(list_of_bubbles)
    #print("Long cycles done")
    #print("List of cycles:")
    #print(listofcycles)

    #print(mytuple)

    #simpleDG = nx.contracted_nodes(DG, 1, 3)
"""Overall method used to simplify the graph
During this method: - Cycles (denoting repetative regions) in the graph are identified  and resolved
                    - Nodes are added (using delta_len, to make sure, no INDELS were overseen)
                    - Bubbles are identified and if possible popped             
INPUT:  DG: Directed Graph before simplifications,
        max_bubblesize: maximum number of elements making up a bubble
        delta_len: Maximum length differences in between two reads
OUTPUT: DG: Graph after simplification took place    """
def simplifyGraph(DG,max_bubblesize,delta_len):
    #remove edges which yield self loops, not sure yet whether it makes sense to remove or if needed
    print("self loops")
    print(list(nx.selfloop_edges(DG)))
    list_of_cycles=find_repetative_regions(DG)
    find_bubbles(DG, max_bubblesize, delta_len)
    DG=merge_nodes(DG)
    return DG
