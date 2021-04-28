import networkx as nx
from collections import Counter
from main import draw_Graph

# TODO find out how to find bubbles->implement pop_bubbles
"""function to merge consecutive nodes, if they contain the same reads to simplify the graph
    INPUT: DG Directed Graph
    OUTPUT: DG: Directed Graph with merged nodes
"""


def merge_nodes(DG):
    # iterate over the edges to find all pairs of nodes
    edgesView = DG.edges.data()
    for ed_ge in edgesView:
        startnode = ed_ge[0]
        endnode = ed_ge[1]
        # we only need to know the out degree of the start node and the end degree of the end node
        start_out_degree = DG.out_degree(startnode)
        end_in_degree = DG.in_degree(endnode)
        # if the degrees are both equal to 1 and if none of the nodes is s or t
        if (start_out_degree == end_in_degree == 1 and startnode != "s" and endnode != "t"):
            # print("Merging nodes "+startnode+" and "+endnode)
            # use the builtin function to merge nodes, prohibiting self_loops decreases the amount of final edges
            DG = nx.contracted_nodes(DG, startnode, endnode, self_loops=False)
    return (DG)


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
    # data structure which holds all the cycles in the graph
    list_of_cycles = []
    # iterate over the cycles to retrive all nodes which are part of the cycles
    for comp in altcyc:
        # if len(comp) > 1:
        intermediate_cycle = []
        # iterate over the nodes in each cycle to get a better output format (start, nodename,end)
        for node_i in comp:
            intermediate = tuple(map(int, node_i.split(', ')))
            print(intermediate)
            # print(type(intermediate))
            intermediate_cycle.append((node_i, intermediate))
            # real_cycles.append(intermediate_sorted)
            # print(intermediate_sorted)
            # if not node_i in cycle_nodes:
            #    cycle_nodes.append(node_i)
            # sort the nodes in a cycle by start coordinates to simplify the resolving of cycles later on
        cycle_sorted = sorted(intermediate_cycle, key=lambda x: x[0])
        # print("Cycle_sorted type")
        # print(str(type(cycle_sorted)))
        # print(cycle_sorted)
        list_of_cycles.append(cycle_sorted)
    if list_of_cycles:
        print("Found repetative region in reads")
        for cyc in list_of_cycles:
            print(cyc)
    else:
        print("No cycles found in the graph")
    return (list_of_cycles)


"""Helper method for find_bubbles: This method finds the minimum and maximum nodes in the bubble
INPUT:      cycle:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element: The node deemed to be the starting node of the bubble (given as tuple)
            max_element: The node deemed to be the end node of the bubble (given as tuple)
            DG:         the directed graph we want to pop the bubble in
            contains_s: A boolean value denoting whether the cycle contains node "s"
            contains_t: A boolean value denoting whether the cycle contains node "t"


def find_min_max(contains_s,contains_t,listofnodes,DG):
    #TODO: this method is not able to sort nodes which have the same startpos. This is wrong Idea: if samestart: access nodes and find out which is first
    #This could also result in an error for max_element
    if contains_t:
        max_element = "t"
        #get_min_element(listofnodes,DG)
    if contains_s:

        min_element = "s"
    if not (contains_t or contains_s):
        intermediate_sorted = sorted(listofnodes, key=lambda x: x[0])
        min_element = intermediate_sorted[0]
        max_element = intermediate_sorted[-1]
    print("intermediate_sorted")
    print(intermediate_sorted)
    return(min_element,max_element)
"""

# def measure_path_length_differences(min_node_out,max_node_in,delta_len):
# TODO:until here the whole algo seems to work correctly. Now however, we have to think about missing exons
#    min_path_1_reads=DG.nodes[min_node_out[0]]['reads']
#    print("min_path_1_reads")
#    print(min_path_1_reads)
#    min_path_2_reads=DG.nodes[min_node_out[1]]['reads']
#    print("MinReadlist")
#    print(min_readlist)
#    print("MaxReadlist")
#    print(max_readlist)
"""Helper function for find_bubbles, called by get_min_max. 
Returns the in degree of a node with respect to edges which are inside of a bubble
INPUT:      DG:     The graph for which the element is to be found
            listofnodes:A list containing all nodes which are part of the bubble
            node:   The current node
OUTPUT:     in_degree: An integer, indicating the number of incoming edges, which are part of the bubble.
"""


def get_in_degree(DG, listofnodes, node):
    # retreive all in edges of node
    in_edges_all = DG.in_edges(node)
    in_degree = 0
    # if in_edges_all contains any elements
    if in_edges_all:
        # for each edge in in_edges_all
        for edge in in_edges_all:
            # if the other node is in listofnodes
            if edge[0] in listofnodes:
                # increase in_degree by 1
                in_degree = in_degree + 1
    return in_degree


"""Helper function for find_bubbles, called by get_min_max. 
Returns the out degree of a node with respect to edges which are inside of a bubble
INPUT:      DG:     The graph for which the element is to be found
            listofnodes:A list containing all nodes which are part of the bubble
            node:   The current node
OUTPUT:     out_degree: An integer, indicating the number of outgoing edges, which are part of the bubble.
"""


def get_out_degree(DG, listofnodes, node):
    # retreive all out edges of node
    out_edges_all = DG.out_edges(node)
    out_degree = 0
    # if out_edges_all contains any elements
    if out_edges_all:
        for edge in out_edges_all:
            # if the other node is in listofnodes
            if edge[1] in listofnodes:
                # increase out_degree by 1
                out_degree = out_degree + 1
    return out_degree


"""Helper method for find_bubbles
Returns the minimum element (=the source) of the bubble as well as the maximum element(=the sink) of the bubble
INPUT:      contains_s: Boolean value indication whether the bubble contains the source node s
            contains_t: Boolean value indication whether the bubble contains the sink node t
            listofnodes:A list containing all nodes which are part of the bubble
            DG:     The graph for which the element is to be found
OUTPUT:     min_element: the source node of the bubble(local source)
            max_element: the sink node of the bubble(local sink)
"""


def get_min_max(contains_s, contains_t, listofnodes, DG):
    print("getminmax")
    if contains_t:
        max_element = "t"
    else:
        for node in listofnodes:
            in_deg = get_in_degree(DG, listofnodes, node)
            out_deg = get_out_degree(DG, listofnodes, node)
            if in_deg == 2 and out_deg == 0:
                max_element = node
            elif in_deg == 2 or out_deg == 0:
                print("ERROR: in_deg= " + str(in_deg) + ", out_deg=" + str(out_deg))
    if contains_s:
        min_element = "s"
    else:
        for node in listofnodes:
            in_deg = get_in_degree(DG, listofnodes, node)
            out_deg = get_out_degree(DG, listofnodes, node)
            if in_deg == 0 and out_deg == 2:
                min_element = node
            elif in_deg == 0 or out_deg == 2:
                print("ERROR: in_deg= " + str(in_deg) + ", out_deg=" + str(out_deg))
    return min_element, max_element

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

"""Helper method for find_bubbles: This method figures out which reads belong to which part from bubble source to bubble sink
INPUT:      cycle:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element: The node deemed to be the starting node of the bubble (given as tuple)
            max_element: The node deemed to be the end node of the bubble (given as tuple)
            DG:         the directed graph we want to pop the bubble in
            contains_s: A boolean value denoting whether the cycle contains node "s"
            contains_t: A boolean value denoting whether the cycle contains node "t"

"""


# TODO:correct this method
def get_path_reads(cycle, min_element, max_element, DG, contains_s, contains_t):
    #find min_node and max_node as well as their read lists: We need if else as "s" and "t" have slightly different read data structures
    if contains_s:
        min_node = "s"
        min_readlist = DG.nodes[min_node]['reads']
    else:
        inter_min_readlist = DG.nodes[min_element]['reads']
        min_readlist = [i[0] for i in inter_min_readlist]
    if contains_t:
        max_node = "t"
        max_readlist = DG.nodes[max_node]['reads']
    else:
        inter_max_readlist = DG.nodes[max_element]['reads']
        max_readlist = [i[0] for i in inter_max_readlist]
    print("MinreadList")
    print(min_readlist)
    print("MaxreadList")
    print(max_readlist)
    #find all reads which are present in min_node as well as in max_node as those are the reads we need to figure out whether their lengths are equal
    shared_reads = []
    for read in min_readlist:
        if read in max_readlist:
            shared_reads.append(read)
    print("SHARED")
    print(shared_reads)
    #find the nodes which are directly connected to min_node (min_node_out) and to max_node(max_node_in) this enables the finding of which reads we have to compare
    min_edges = DG.out_edges(min_element)
    max_edges = DG.in_edges(max_element)
    min_node_out = []
    max_node_in = []
    for edge in min_edges:
        if edge[1] in cycle:
            min_node_out.append(edge[1])
    for edge in max_edges:
        if edge[0] in cycle:
            max_node_in.append(edge[0])
    if (len(min_node_out) > 2 or len(max_node_in) > 2):
        print("Think again!")
        print(min_node_out)
        print(max_node_in)
    else:
        print("Seems correct!")
        print(min_node_out)
        print(max_node_in)
    out_path_reads = []
    in_path_reads = []
    path_starts={}
    #print("MinNodeOut")
    #print(min_node_out)
    #print("MaxNodeIn")
    #print(max_node_in)
    # TODO: something in here does not work out correctly yet
    for out_node in min_node_out:
        inter_out_readlist = list(DG.nodes[out_node]['reads'])
        #print("interoutreadlist")
        #print(inter_out_readlist)
        out_readlist = [i[0] for i in inter_out_readlist]
        #print("Outreadlist")
        #print(out_readlist)
        out_path_reads_list=[]
        for read in out_readlist:
            if read in shared_reads:
                out_path_reads_list.append(read)
        path_starts[out_node]=out_path_reads_list
    print("in_path_reads")
    print(in_path_reads)
    print("Pathstarts")
    print(path_starts)
    #iterate over all shared reads and get the pathlength for each.
    path_lengths={}
    for key,value in path_starts.items():
        path_length_list=[]
        for r_id in value:
            #TODO: write something meaningful here to get the length of the path for each read
            read_length=9
            path_length_list.append(read_length)
        path_lengths[key]=path_length_list
    return path_starts


""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
            max_bubblesize  parameter giving the maximum size of a bubble to still pop it (TODO: properly define max_bubblesize
            delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """
def merge_bubble_nodes(DG,path_reads,min_element,max_element):
    print("Hello World")
def find_and_pop_bubbles(DG, max_bubblesize, delta_len):
    # get undirected version of the graph to find bubbles
    UG = DG.to_undirected()

    # bubble_limit=2*max_bubblesize+2
    # print("Bubble limit: "+str(bubble_limit))
    # collect the bubbles in the graph (Bubbles denote possible mutations in the minimizers)
    list_of_bubbles = []
    # find all cycles in the undirected graph->bubbles
    listofcycles = nx.cycle_basis(UG)
    print(listofcycles)
    # iterate over the different bubbles
    for listint, cycle in enumerate(listofcycles):
        listofnodes = []
        listofnormalnodes = []
        contains_s = False
        contains_t = False
        # iterate over the nodes
        for node in cycle:
            # treat nodes different if they are s or t (just because of how they were generated)
            if "s" in node:
                s_tuple = ("s", 0, 0)
                contains_s = True

            elif "t" in node:
                contains_t = True
                t_tuple = ("t", 0, 0)
            # for all other nodes the attributes (startnode,endnode,read_id) are recovered from the nodeids
            else:
                mytuple = tuple(map(int, node.split(', ')))
                listofnodes.append(mytuple)
                listofnormalnodes.append(node)
        print("ListofNodes")
        print(listofnodes)
        # find the minimum and maximum for each bubble
        min_element, max_element = get_min_max(contains_s, contains_t, listofnormalnodes, DG)
        # find the nodes which are on either path to be able to tell apart the paths on the bubble
        (path_reads,min_readlist,max_readlist)=get_path_reads(listofnormalnodes, min_element, max_element, DG, contains_s, contains_t)
        # TODO compare the reads of both paths: If they differ by less than delta_len: Pop the bubble
        merge_bubble_nodes(DG,path_reads,min_element,max_element)
        print("Cycle")
        print(cycle)
        print("MinElement")
        print(min_element)
        print("MaxElement")
        print(max_element)
        # print(listofnodes)
        # print (min(listofnodes, key = lambda x: x[0]))
        # print(max(listofnodes, key=lambda x: x[1]))
    #        if not(len(listofnodes)> bubble_limit):
    #            possible_bubble_sorted = sorted(listofnodes, key=lambda x: x[0])
    #            list_of_bubbles.append(possible_bubble_sorted)
    #    else:
    #        listofcycles.pop(listint)

    # for cycle in list_of_cycles:
    #    list_of_bubbles= [x for x in list_of_bubbles if x!= cycle]
    # for bubble in list_of_bubbles:
    #    endtuple=bubble[-1]
    #    starttuple=bubble[0]
    #    endname=str(endtuple[0])+','+str(endtuple[1])
    #    startname = str(starttuple[0]) + ',' + str(starttuple[1])
    #    list_of_bubbles_sliced=bubble[1:-1]
    # for bub_node in list_of_bubbles_sliced:
    # readlistandpos=DG.nodes[startname]['reads']
    # print(readlistandpos)


# print(list_of_bubbles)
# print("Long cycles done")
# print("List of cycles:")
# print(listofcycles)

# print(mytuple)

# simpleDG = nx.contracted_nodes(DG, 1, 3)
"""Overall method used to simplify the graph
During this method: - Bubbles are identified and if possible popped     
                    - Nodes are merged        
INPUT:  DG: Directed Graph before simplifications,
        max_bubblesize: maximum number of elements making up a bubble
        delta_len: Maximum length differences in between two reads
OUTPUT: DG: Graph after simplification took place    """


def simplifyGraph(DG, max_bubblesize, delta_len):
    print("Simplifying the Graph (Merging nodes, popping bubbles)")
    # remove edges which yield self loops, not sure yet whether it makes sense to remove or if needed
    print("self loops")
    print(list(nx.selfloop_edges(DG)))
    list_of_cycles = find_repetative_regions(DG)
    print(list_of_cycles)
    print(DG.edges())
    #find_and_pop_bubbles(DG, max_bubblesize, delta_len)
    DG = merge_nodes(DG)
    return DG
