import networkx as nx
from collections import Counter
from main import draw_Graph
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
    DG=merge_nodes(DG)
    return DG

""" Method to simplify the graph. In the first step the cycles in the graph are looked up and the cyles are told apart from the bubbles
    INPUT:  DG  Directed Graph
            max_bubblesize  parameter giving the maximum size of a bubble to still pop it (TODO: properly define max_bubblesize
            delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """
def find_bubbles(DG,max_bubblesize,delta_len):

    UG=DG.to_undirected()
    #otherfun=list(nx.simple_cycles(DG))
    #print("other functions' cycles:\n")
    #print(otherfun)
    bubble_limit=2*max_bubblesize+2
    #print("Bubble limit: "+str(bubble_limit))
    list_of_cycles=find_repetative_regions(DG)
    #collect the bubbles in the graph and tell them apart from the cycles (Bubbles denote possible mutations in the minimizers)
    print("List of possible bubbles:")
    list_of_bubbles = []
    listofcycles = nx.cycle_basis(UG)
    for listint,cycle in enumerate(listofcycles):
        listofnodes=[]
        if len(cycle)>1:
            #print(cycle)
            for node in cycle:
                mytuple = tuple(map(int, node.split(', ')))
                #intermediate_sorted = sorted(mytuple, key=lambda x: x[0])
                listofnodes.append(mytuple)
                #print(listofnodes)
            #print (min(listofnodes, key = lambda x: x[0]))
            #print(max(listofnodes, key=lambda x: x[1]))
            if not(len(listofnodes)> bubble_limit):
                possible_bubble_sorted = sorted(listofnodes, key=lambda x: x[0])
                list_of_bubbles.append(possible_bubble_sorted)
        else:
            listofcycles.pop(listint)

    for cycle in list_of_cycles:
        list_of_bubbles= [x for x in list_of_bubbles if x!= cycle]
    for bubble in list_of_bubbles:
        endtuple=bubble[-1]
        starttuple=bubble[0]
        endname=str(endtuple[0])+','+str(endtuple[1])
        startname = str(starttuple[0]) + ',' + str(starttuple[1])
        list_of_bubbles_sliced=bubble[1:-1]
        #for bub_node in list_of_bubbles_sliced:
        #readlistandpos=DG.nodes[startname]['reads']
        #print(readlistandpos)
    print(list_of_bubbles)
    print("Long cycles done")
    #print("List of cycles:")
    #print(listofcycles)

    #print(mytuple)

    #simpleDG = nx.contracted_nodes(DG, 1, 3)