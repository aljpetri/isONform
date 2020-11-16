import networkx as nx
from collections import Counter
def iterate_edges_to_add_nodes(DG,delta_len,list_of_cycles):
    print(delta_len)
    edgesView=DG.edges.data()
    print("EdgesView:")
    print(edgesView)
    nodelist=list(DG.nodes(data=True))
    #print("Nodelist:")
    #print(nodelist)
    #print("Nodelist done")
    #for edge in list(DG.edges(data=True)):
    #    startnode=edge[0]
    #    endnode=edge[1]
    #    print("Startlist")
    #    #print(startlist)
    print("In between methods")
    for ed_ge in edgesView:
        startnode=ed_ge[0]
        endnode=ed_ge[1]
        startlist=DG.nodes[startnode]['reads']
        endlist=DG.nodes[endnode]['reads']
        #print(startnode)
        #print("Startlist")
        #print(startlist)
        add_nodes(startnode,startlist,endnode,endlist,delta_len,list_of_cycles)
def add_nodes(startnode,startlist,endnode,endlist,delta_len,cycle_nodes):
    #TODO: Find out whether tuple is in a cycle. If yes take the positions with the minimum distance to measure if over delta_len
    startnode_in_cyc=False
    endnode_in_cyc = False
    print(cycle_nodes)
    print(startnode)
    for cyc_node in cycle_nodes:
        if cyc_node == startnode:
            startnode_in_cyc = True
        if cyc_node== endnode:
            endnode_in_cyc=True
    if startnode_in_cyc:
        if endnode_in_cyc:
            print("startnode and endnode true")
            endcount_list = Counter(elem[0] for elem in endlist)
            endcounts = endcount_list.items()
            print("endcount_list:")
            print(endcounts)
            endreads=[]
            for r_id, count in endcounts:
                if count > 1:
                   endreads.append(r_id)
            print("endreads:")
            print(endreads)
        print("startnode true")
        startcount_list=Counter(elem[0] for elem in startlist)
        startcounts=startcount_list.items()
        print("startcount_list:")
        print(startcounts)
        startreads = []
        for r_id, count in startcounts:
            if count > 1:
                startreads.append(r_id)
        print("startreads:")
        print(startreads)
    elif endnode_in_cyc:
        endcount_list = Counter(elem[0] for elem in endlist)
        endcounts = endcount_list.items()
        print("endcount_list:")
        print(endcounts)
        endreads = []
        for r_id, count in endcounts:
            if count > 1:
                endreads.append(r_id)
        print("endreads:")
        print(endreads)
    else:
        print("none in cycle")
    #Output = [item for item in list_of_cycles if item[0] == 3]

    #TODO: 4 cases: 1st: startnode in list_of_cycles 2nd:endnode in list of cycles 3rd:both in listofcycles, 4th:none in listofcycles
    #Output = [item for item in list_of_cycles if item[0] == 3]
    for starttuple in startlist:
        start_r_id = starttuple[0]
        endtuple_list = [tup for tup in endlist if tup[0] == start_r_id]

"""
function to find cycles in the graph (, which denote repetitive regions)
INPUT: DG Directed Graph
OUTPUT: List_of_cycles: A list holding all cycles present in DG
"""
def find_repetative_regions(DG):
    altcyc = nx.simple_cycles(DG)
    print("Alternative cycles:")
    list_of_cycles = []
    cycle_nodes=[]
    # collect the cycles in the graph (denoting repetative regions)
    for comp in altcyc:
        if len(comp) > 1:
            intermediate_cycle = []
            for node_i in comp:
                intermediate = tuple(map(int, node_i.split(', ')))
                print(intermediate)
                # print(type(intermediate))
                intermediate_cycle.append((node_i,intermediate))
                # real_cycles.append(intermediate_sorted)
                # print(intermediate_sorted)
                if not node_i in cycle_nodes:
                    cycle_nodes.append(node_i)
            cycle_sorted = sorted(intermediate_cycle, key=lambda x: x[0])
            print(cycle_sorted)
            list_of_cycles.append(cycle_sorted)
    if list_of_cycles:
        print("Found repetative region in reads")
        for cyc in list_of_cycles:
            print(cyc)
    return (list_of_cycles,cycle_nodes)
def resolve_cycles(cycle_nodes,list_of_cycles,DG):
    for cycle in list_of_cycles:
        reads_in_cycle=[]
        for node in cycle_nodes:
            readlist = (node,DG.nodes[node]['reads'])
            reads_in_cycle.append(readlist)
    for i,readlist in enumerate(reads_in_cycle):
        node=readlist[0]
        list_of_reads=readlist[1]
        read_count_list = Counter(elem[0] for elem in list_of_reads)
        read_counts = read_count_list.items()
        print("readcount_list:")
        print(read_counts)
        multiple_reads = []
        max_count=0
        for r_id, count in read_counts:
            if count > 1:
                if count>max_count:
                    max_count=count
                multiple_reads.append(r_id)
        if max_count>2:
            raise NameError('Cycle too complex!!!')
        print("Maximum count: " + str(max_count))
        print("Multiple_reads: ")
        print(multiple_reads)
        print(node)
        print(list_of_reads)
        if multiple_reads:
            first_mult_read=[item for item in list_of_reads if item[0] == multiple_reads[0]]
            print("Firstmultread")
            print(first_mult_read)
        DG.nodes[node]['reads']

"""Overall method used to simplify the graph
During this method: - Cycles (denoting repetative regions) in the graph are identified  and resolved
                    - Nodes are added (using delta_len, to make sure, no INDELS were overseen
                    - Bubbles are identified and if possible popped             
INPUT:  DG: Directed Graph before simplifications,
        max_bubblesize: maximum number of elements making up a bubble
        delta_len: Maximum length differences in between two reads
OUTPUT: DG: Graph after simplification took place    """
def simplifyGraph(DG,max_bubblesize,delta_len):
    #remove edges which yield self loops, not sure yet whether it makes sense to remove or if needed
    DG.remove_edges_from(nx.selfloop_edges(DG))
    list_of_cycles,cycle_nodes=find_repetative_regions(DG)
    resolve_cycles(cycle_nodes,list_of_cycles,DG)
    #iterate_edges_to_add_nodes(DG,delta_len,cycle_nodes)


def find_bubbles(DG,max_bubblesize,delta_len):
    """
            Method to simplify the graph. In the first step the cycles in the graph are looked up and the cyles are told apart from the bubbles

        """
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