import networkx as nx
from collections import Counter
from main import draw_Graph
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


    #Output = [item for item in list_of_cycles if item[0] == 3]
    for starttuple in startlist:
        start_r_id = starttuple[0]
        endtuple_list = [tup for tup in endlist if tup[0] == start_r_id]

"""Sub method to bend the edges to the correct nodes after adding nodes
INPUT:  in_edges:   list of edges going into the node before
        out_edges:  list of edges going out of the node before
        DG:         graph with new nodes not connected to other nodes yet and the previous state concerning the looked at node
        new_nodes:  list of nodes which were added to the graph and need edges
        node:       node which was split
OUTPUT: DG          graph with the correct edges assigned to each node
"""
def generate_correct_edges(in_edges,out_edges,DG,new_nodes,node):
    print("Generate_correct_edges")
    #TODO: iterate over all nodes and get needed info to test which edge is connected to which node
    original_node_positions=DG.nodes[node]['reads']
    #print("OriginalNodeInfos")
    #print(original_node_positions)
    new_node_positions = {}
    edges_to_remove=[]
    #iterate through the new nodes and retreive the information of which reads have which positions in the node
    for nnode in new_nodes:
        #print("Nnode")
        #print(nnode)
        infos_nnode=DG.nodes[nnode[0]]['reads']
        new_node_positions[nnode[0]]=infos_nnode
    #print("Othernode_positions:")
    #print(new_node_positions)
    pos_distances_in = {}
    #Implementation copied from out_edges TODO: Fix everything to apply to in_edges

    for edge in in_edges:
        othernode = edge[0]
        print(othernode)
        other_infos = DG.nodes[othernode]['reads']
        pos_distances_in[othernode] = []
        # generate the distances for the initial node to have a reference for the new nodes
        for tuple in other_infos:
            r_id = tuple[0]
            startpos = tuple[1]
            known_tuples = [item for item in original_node_positions if item[0] == r_id]
            origin_tuple = known_tuples[0]
            endpos = origin_tuple[2]
            distance = abs(startpos - endpos)
            distance_tuple = (node, r_id, distance, startpos)
            pos_distances_in[othernode].append(distance_tuple)
        # iterate through the entries in new_node_positions (all entries appointed to new nodes)
        for node_id, read_infos in new_node_positions.items():
            # iterate through the entries for one node
            for r_info_tuple in read_infos:
                dist_tupl = [item for item in pos_distances_in[othernode] if item[1] == r_info_tuple[0]]
                for dist_tup in dist_tupl:
                    distance = abs(dist_tup[3] - r_info_tuple[2])
                    if distance < dist_tup[2]:
                        pos_distances_in[othernode].remove(dist_tup)
                        new_dist_tup = (node_id, r_info_tuple[0], distance, dist_tup[3])
                        pos_distances_in[othernode].append(new_dist_tup)
    #print("New Output")
    for node_dist_id, node_distance_infos in pos_distances_in.items():
        print("Node ID")
        print(node_dist_id)
        print("infos")
        print(node_distance_infos[0][0])
        startcount_list = list(Counter(elem[0] for elem in node_distance_infos))
        print(startcount_list)
        if not (node in startcount_list):
            print("Removing edge " + node_dist_id + ";" + node)
            DG.remove_edge(node_dist_id,node)
        for id in startcount_list:
            if not id == node:
                print("Adding edge " + node_dist_id+ ";" + id)
                DG.add_edge(node_dist_id,id)
        for distance_tuple in node_distance_infos:
            print(distance_tuple)
    #First and hopefully correct implementation
    pos_distances_out = {}
    for edge in out_edges:
        othernode = edge[1]
        print(othernode)
        other_infos = DG.nodes[othernode]['reads']
        pos_distances_out[othernode]=[]
        #generate the distances for the initial node to have a reference for the new nodes
        for tuple in other_infos:
            r_id=tuple[0]
            startpos=tuple[1]
            known_tuples=[item for item in original_node_positions if item[0] == r_id]
            origin_tuple=known_tuples[0]
            endpos=origin_tuple[2]
            distance=abs(startpos-endpos)
            distance_tuple=(node,r_id,distance,startpos,)
            pos_distances_out[othernode].append(distance_tuple)
        print("Pos_distances_out")
        print(pos_distances_out)
        #iterate through the entries in new_node_positions (all entries appointed to new nodes)
        for node_id,read_infos in new_node_positions.items():
            #iterate through the entries for one node
            for r_info_tuple in read_infos:
                dist_tupl=[item for item in pos_distances_out[othernode] if item[1] == r_info_tuple[0]]
                print(dist_tupl)
                for dist_tup in dist_tupl:
                #dist_tup=dist_tupl[0]
                    print(r_info_tuple)
                    distance=abs(dist_tup[3]-r_info_tuple[2])
                    if distance<dist_tup[2]:
                        pos_distances_out[othernode].remove(dist_tup)
                        new_dist_tup=(node_id,r_info_tuple[0],distance,dist_tup[3])
                        pos_distances_out[othernode].append(new_dist_tup)
    print("New Input")
    for node_dist_id,node_distance_infos in pos_distances_out.items():
        print("Node ID")
        print(node_dist_id)
        print("infos")
        print(node_distance_infos[0][0])
        endcount_list = list(Counter(elem[0] for elem in node_distance_infos))
        print(endcount_list)
        if not(node in endcount_list):
            print("Removing edge "+node+";"+node_dist_id)
            DG.remove_edge(node, node_dist_id)
        for id in endcount_list:
            if not id==node:
                print("Adding edge " + id + ";" + node_dist_id)
                DG.add_edge(id, node_dist_id)
        for distance_tuple in node_distance_infos:
            print(distance_tuple)


    print("Pos_distances_out")
    print(pos_distances_out)
"""
Method used to break the cycles which are present in the graph
INPUT: DG Directed Graph
OUTPUT: DG Directed Graph altered to not contain cycles anymore(this means that nodes which occur more then once are split )
"""
#TODO: Add resolving method for self cycles (directly repeating intervals)
def resolve_cycles(list_of_cycles,DG):
    self_cycles=list(nx.selfloop_edges(DG))
    #TODO:write method to get rid of self cycles
    for self_cyc in self_cycles:
        print("Self cycles found in the graph:")
        print(self_cycles)
    #retreive all nodes which are present in cycles
    for cycle in list_of_cycles:
        print("Cycle:")
        print(cycle)
        reads_in_cycle=[]
        cycle_nodes=[]
        #generate a set reads_in_cycle which contains sets of node-names and the reads which the node contains
        for node in cycle:
            readlist = (node[0],DG.nodes[node[0]]['reads'])
            reads_in_cycle.append(readlist)
            if not node[0] in cycle_nodes:
                cycle_nodes.append(node[0])
    print("Nodes in cycle:")
    print(reads_in_cycle)
    print("Cycle_Nodes:")
    print(cycle_nodes)
    #iterate through all nodes of all cycles in the graph
    for i,readlist in enumerate(reads_in_cycle):
        node=readlist[0]
        list_of_reads=readlist[1]
        read_count_list = Counter(elem[0] for elem in list_of_reads)
        read_counts = read_count_list.items()
        #print("readcount_list:")
        #print(read_counts)
        multiple_reads = []
        max_count=0
        #scan through read_counts, (which denotes the counts of how often an interval is present in a read)
        #find the maximum number of ocurances of a certain interval in one read
        for r_id, count in read_counts:
            if count > 1:
                if count>max_count:
                    max_count=count
                multiple_reads.append(r_id)
        #if max_count>2:
        #    raise NameError('Cycle too complex!!!')
        #print("Maximum count: " + str(max_count))
        #prints the list of reads which have multiple occurances of a certain interval
        #print("Multiple_reads: ")
        #print(multiple_reads)
        #print(node)
        #print(list_of_reads)
        multi_info={}
        new_nodes=[]
        if multiple_reads:
            for one_multi_read in multiple_reads:
                mult_read_info=[item for item in list_of_reads if item[0] == one_multi_read]
                mult_read_pos=[]
                for info in mult_read_info:
                    mult_read_pos.append((info[1],info[2]))
                multi_info[one_multi_read]=mult_read_pos
            for i in range(0,max_count):
                if i>0:
                    node_info=[]
                    naming =None
                    for key,value in multi_info.items():
                        tuple_to_del=(key,value[i][0],value[i][1])
                        node_info.append(tuple_to_del)
                        list_of_reads.remove(tuple_to_del)
                        #print("Tuple to del:")
                        #print(tuple_to_del)
                    naming = node_info[0]
                    name=str(naming[1]) + ", " + str(naming[2]) + ", " + str(naming[0])
                    newinfos=(naming[1],naming[2],name)
                    #TODO add node again
                    DG.add_node(name, reads=node_info)
                    new_nodes.append((name,naming))
                    print("Adding node "+name)
                    print("List of reads")
                    print(list_of_reads)
                    print("Node Info")
                    print(node_info)
            print("New nodes:")
            print(new_nodes)
            print("Outgoing Edges")
            #find the edges which point FROM the repeating node
            out_edges = list(DG.out_edges(node))
            in_edges = list(DG.in_edges(node))
            generate_correct_edges(in_edges,out_edges,DG,new_nodes,node)



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
    #DG.remove_edges_from(nx.selfloop_edges(DG))
    list_of_cycles=find_repetative_regions(DG)
    #resolve_cycles(list_of_cycles,DG)
    #iterate_edges_to_add_nodes(DG,delta_len,cycle_nodes)
    #paths = list(nx.all_simple_paths(DG, source="s", target="t"))
    #print(list(paths))
    #print("Found "+str(len(paths))+ "simple paths")
    #print("merging of nodes")
    #merge_nodes(DG)
    #draw_Graph(DG)
    #print("Graph traits after merge_nodes:")
    #print("#Nodes for DG: " + str(DG.number_of_nodes()) + " , #Edges for DG: " + str(DG.number_of_edges()))
    #nodelist = list(DG.nodes)
    #print(nodelist)
    #edgelist=list(DG.edges.data())
    #print(edgelist)


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