import networkx as nx
from collections import Counter
from consensus import *
import matplotlib.pyplot as plt
# TODO find out how to find bubbles->implement pop_bubbles
"""function to merge consecutive nodes, if they contain the same reads to simplify the graph
    INPUT: DG Directed Graph
    OUTPUT: DG: Directed Graph with merged nodes
"""
def draw_Graph(DG):
    # defines the graph layout to use spectral_layout. Pos gives the position for each node
    pos = nx.spectral_layout(DG)
    # draws the network in the given layout. Only nodes are labeled, the edges are directed but not labeled
    nx.draw_networkx(DG, pos, font_weight='bold')
    # add labels for the edges and draw them into the graph
    # labels = nx.get_edge_attributes(DG, 'weight')
    # nx.draw_networkx_edge_labels(DG,pos, edge_labels=labels)
    plt.show()
def generate_subgraphs(DG,bubbles):
    for bubble in bubbles:
        SG=DG.subgraph(bubble)
        draw_Graph(SG)
def generate_subgraph(DG,bubbles):
    #for bubble in bubbles:
        SG=DG.subgraph(bubbles)
        draw_Graph(SG)
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


"""Helper method for get_bubble_start_end: This method finds the minimum and maximum nodes in the bubble
The method iterates through the nodes in a bubble and collects all their out_nodes. For The minimum node there will not be an out_node
INPUT:      listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            DG:         the directed graph we want to pop the bubble in
OUTPUT:     startnode: The node deemed to be the starting node of the bubble (given as tuple)
            """
def find_bubble_start(DG,listofnodes):
    real_bubble=True
    out= set()
    outedges=[]
    for node in listofnodes:
        out_edges_all=DG.out_edges(node)
        for out_edge in out_edges_all:
            if (out_edge[1] in listofnodes):
                out.add(out_edge[1])
                outedges.append(out_edge)
    set_a = set(listofnodes)
    set_b = out
    subtraction=set_a-set_b
    list_of_strings = [str(s) for s in subtraction]
    if len(list_of_strings)>1:
        real_bubble=False
    bubble_start_reads=set()
    for out_edge in outedges:
        edge_infos=DG[out_edge[0]][out_edge[1]]["edge_supp"]
        for entry in edge_infos:
            bubble_start_reads.add(entry)
    startnode = " ".join(list_of_strings)
    print("start_reads",bubble_start_reads)
    return  real_bubble, startnode,bubble_start_reads
"""Helper method for get_bubble_start_end: This method finds the maximum node in the bubble
The method iterates through the nodes in a bubble and collects all their in_nodes. For The maximum node there will not be an in_node
INPUT:      listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            DG:         the directed graph we want to pop the bubble in
OUTPUT:     endnode: The node deemed to be the end node of the bubble (given as tuple)
            """
def find_bubble_end(DG,listofnodes):
    real_bubble=True
    in_nodes=set()
    inedges=[]
    bubble_end_reads=set()
    for node in listofnodes:
        in_edges_all=DG.in_edges(node)
        for in_edge in in_edges_all:
            if in_edge[0] in listofnodes:
                in_nodes.add(in_edge[0])
                inedges.append(in_edge)
    set_a = set(listofnodes)
    set_b = in_nodes
    subtraction = set_a - set_b
    list_of_strings = [str(s) for s in subtraction]
    if len(list_of_strings)>1:
        real_bubble=False
    for in_edge in inedges:
        if(in_edge[1]==list_of_strings[0]):
            edge_infos=DG[in_edge[0]][in_edge[1]]["edge_supp"]
            for entry in edge_infos:
                bubble_end_reads.add(entry)
    endnode= " ".join(list_of_strings)
    return real_bubble,endnode,bubble_end_reads


"""Helper method for find_bubbles
Returns the minimum element (=the source) of the bubble as well as the maximum element(=the sink) of the bubble
INPUT:      listofnodes:A list containing all nodes which are part of the bubble
            DG:     The graph for which the element is to be found
OUTPUT:     min_element: the source node of the bubble(local source)
            max_element: the sink node of the bubble(local sink)
            contains_s: Boolean value indication whether the bubble contains the source node s
            contains_t: Boolean value indication whether the bubble contains the sink node t
"""
def get_bubble_start_end(DG, listofnodes):
    real_bubble_s,startnode,bubble_start_reads=find_bubble_start(DG,listofnodes)
    real_bubble_e,endnode,bubble_end_reads=find_bubble_end(DG,listofnodes)
    real_bubble=real_bubble_e and real_bubble_s
    shared_reads=list(bubble_start_reads.intersection(bubble_end_reads))
    print("shared_reads:", shared_reads)
    return startnode,endnode,real_bubble,shared_reads


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

"""def linearize_bubble(DG,path_reads,min_element,max_element ):
    node_distances={}
    edge_support
    for node in bubble_nodes:
        node_distances[node]=get_avg_distance_to_start()
        out_edgelist = list(DG.out_edges(node))
        in_edgelist = list(DG.in_edges(node))
        for edge in out_edgelist:
            if edge[1] in bubble_nodes:
                #delete edge
        for edge in in_edgelist:
            if edge[0] in bubble_nodes:
                #delete edge
    curr_node=min_element
    while not curr_node= max_element:
            DG.add_edge(curr_node, next_node, length=length)
            edge_support[curr_node, next_node] = path_reads
    nx.set_edge_attributes(DG, edge_support, name='edge_supp')
    print("Bubble popped by linearization")
    return(DG)"""

        #sort nodes by node_distances
def generate_consensus_via_spoa(read_ids, startpositions,endpositions,curr_best_seqs, reads,work_dir, outfolder, max_seqs_to_spoa=200):
    def generate_isoform_using_spoa(curr_best_seqs, reads, work_dir, outfolder, max_seqs_to_spoa=200):
        # print("reads")
        # print(reads)
        mapping = {}
        consensus_file = open(os.path.join(outfolder, "spoa.fa"), 'w')

        # curr_best_seqs = curr_best_seqs[0:3]
        for key, value in curr_best_seqs.items():
            # for equalreads in curr_best_seqs:
            name = 'consensus' + str(value[0])
            # name = 'consensus' + str(equalreads[0])
            mapping[name] = []
            reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
            # if len(equalreads) == 1:
            if len(value) == 1:
                # rid = equalreads[0]
                rid = key
                singleread = reads[rid]
                # print(singleread)
                seq = singleread[1]
                # name='consensus'+str(rid)
                mapping[name].append(singleread[0])
                consensus_file.write(">{0}\n{1}\n".format(name, seq))
                reads_path.close()
            else:
                # print("Equalreads has different size")
                # for i, q_id in enumerate(equalreads):
                for i, q_id in enumerate(value):
                    singleread = reads[q_id]
                    seq = singleread[1]
                    # print(seq)
                    mapping[name].append(singleread[0])
                    if i > max_seqs_to_spoa:
                        break
                    reads_path.write(">{0}\n{1}\n".format(singleread[0], seq))
                reads_path.close()
                spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
                # print("spoa_ref for " + name + " has the following form:" + spoa_ref[0:25])
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
        # print("Isoforms generated")

def linearize_bubble(DG,path_reads,min_element,max_element ):
    print("Hello World")

""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
           delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """
#TODO: We have to generate a consensus of the reads in one bubble path to be able to have a pairwise sequence alignment between the bubble paths
def align_and_linearize_bubble_nodes(DG,path_reads,min_element,max_element,delta_len,all_reads):
    #TODO: generate consensus sequences via spoa
    consensus1=""
    consensus2=""
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score=parasail_alignment(consensus1,consensus2,match_score=2,missmatch_penalty=-2,opening_penalty=3,gap_ext=1)
    #TODO: save consensus into each readfile (could make sense at least)
    #TODO implement linearize_bubble to alter DG in order to pop the bubble
    DG=linearize_bubble(DG,path_reads,min_element,max_element )
    return DG


def find_bubbles(DG):
    #draw_Graph(DG)
    # get undirected version of the graph to find bubbles
    UG = DG.to_undirected()
    # collect the bubbles in the graph (Bubbles denote possible mutations in the minimizers)
    # find all cycles in the undirected graph->bubbles
    list_of_bubbles =nx.cycle_basis(UG)
    return list_of_bubbles

"""function to find shared reads between the bubble_start and bubble_end. Those reads are the ones we will loook at to figure out how similar the paths are
    INPUT:      DG:              The graph for which the element is to be found
                bubble_start:    the source node of the bubble(local source)
                bubble_end:      the sink node of the bubble(local sink)
    OUTPUT:     shared_reads:    a list of reads which are in bubble_start as well as in bubble_end

"""
"""def find_shared_reads(DG,bubble_start, bubble_end):
    all_nodes_reads=DG.nodes(data='reads')
    
    start_reads=all_nodes_reads[bubble_start]
    print("startreads",start_reads)
    end_reads=all_nodes_reads[bubble_end]
    print("end_reads", end_reads)
    start_read_list=start_reads.keys()
    end_read_list=end_reads.keys()
    shared_reads=list(set(start_read_list).intersection(end_read_list))
    print("Shared REads ",shared_reads)
    return shared_reads"""
def get_path_nodes(cycle, min_element, max_element, DG):
    # find the nodes which are directly connected to min_node (min_node_out) and to max_node(max_node_in) this enables the finding of which reads we have to compare
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
    return min_node_out, max_node_in
def get_path_reads(DG,min_node_out,shared_reads):
    path_starts = {}
    for out_node in min_node_out:
        inter_out_readlist = DG.nodes[out_node]['reads']
        print(inter_out_readlist)
        out_readlist = [i for i in inter_out_readlist]
        out_path_reads_list = []
        for read in out_readlist:
            if read in shared_reads:
                out_path_reads_list.append(read)
        path_starts[out_node] = out_path_reads_list
    print("Pathstarts")
    print(path_starts)
    return path_starts
"""Helper method which tests if a found path is viable (meaning if the structure denoted as bubble is an actual bubble or not)
This is done by following the edges and trying to verify that a path is actually supported by at least one read from bubble_start to bubble_end
INPUT:          DG: The directed graph
                path_start: the start node of the path
                initial_support: a list of reads supporting the edge bubble_start,path_start
                cyclce: List of nodes which make up the "bubble"
                bubble_end: The node deemed to be the sink of the "bubble"
OUTPUT:         A boolean value telling us whether there is a set of reads that consistently supports the path
"""
def test_path_viability(DG,path_start,initial_support,cycle,bubble_end):
    curr_support=initial_support
    curr_node=path_start
    #we only iterate until we have reached bubble_end
    while curr_node != bubble_end:
        further_node=False
        curr_out=DG.out_edges(curr_node)
        #we have to test for all possible out edges to be sure to find the path
        for out_edge in curr_out:
            next_node=out_edge[1]
            #we only continue investigating the next_node if it is part of the bubble
            if next_node in cycle:
                #get the edges which support the edge to next_node
                next_support=DG[curr_node][next_node]["edge_supp"]
                #figure out whether there is any overlap between the current support and next_support
                intersect_supp=list(set(curr_support).intersection(next_support))
                #if we have an overlap, this means that at least one read from our initial_support als supports this edge
                if intersect_supp:
                    #we have found a node to continue
                    further_node=True
                    #set the node to be our current_node
                    curr_node=next_node
                    #as we only expect one node to be viable: break the for loop to continue in the next node
                    break
        #we did not find any node we could continue with->our path is not the path of a real bubble
        if not further_node:
            return False
    #we were able to reach bubble_end->we just confirmed that we have a bubble path
    return True
"""Helper method for find_bubbles: This method figures out which reads belong to which part from bubble source to bubble sink
INPUT:      cycle:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element: The node deemed to be the starting node of the bubble (given as tuple)
            max_element: The node deemed to be the end node of the bubble (given as tuple)
            DG:         the directed graph we want to pop the bubble in
            contains_s: A boolean value denoting whether the cycle contains node "s"
            contains_t: A boolean value denoting whether the cycle contains node "t"

"""

def get_path_reads_length(cycle, min_element, max_element, DG,shared_reads):
    #We want to find the nodes, which denote the start points for each path(as we have to find out which reads are in which path)
    min_node_out, max_node_in=get_path_nodes(cycle, min_element, max_element, DG)
    #Now we want to get the actual reads for each path
    path_starts=get_path_reads(DG,min_node_out,shared_reads)
    #iterate over all shared reads and get the pathlength for each
    print("Path_starts ",path_starts)
    path_lengths={}
    for key,value in path_starts.items():
        is_viable_path=test_path_viability(DG,key,value,cycle,max_element)
        if not(is_viable_path):
            generate_subgraph(DG,cycle)
            continue
        id_counter=0
        read_length=0
        for r_id in value:
            #get the path lengths (done by average length)
                #in order to calculate the length of the path, we have to subtract the end coordinste from the min node from the start coordinate of the max node
                max_node_infos=DG.nodes[max_element]['reads']
                min_node_infos=DG.nodes[min_element]['reads']
                positions=max_node_infos[r_id]
                otherpos=min_node_infos[r_id]
                end_pos=otherpos[1]
                start_pos=positions[0]
                read_length+=(start_pos-end_pos)
                id_counter+=1
        path_length=read_length/id_counter
        path_lengths[key]=path_length
    return path_lengths

"""
function which finds the bubbles and if feasible pops them in the graph

INPUT:
    DG:         The directed graph in which we pop the bubbles
    delta_len   Maximum length difference for which the paths are merged ie the bubble is popped
    all_reads   dictionary containing all the reads(string sequence) and their ids
    
    
OUTPUT:
    DG:         Directed Graph containing the popped bubbles"""

def find_and_pop_bubbles(DG, delta_len,all_reads,reads_for_isoforms):
    bubbles=find_bubbles(DG)
    print(type(bubbles))
    nr_bubbles=len(bubbles)
    print("Found "+str(nr_bubbles)+" bubbles in our graph")
    print("Bubbles: ",bubbles)
    #just for debugging and curiosity reasons: We introduce an integer counting the number of pseudobubbles (not true bubbles)
    filter_count=0
    #generate_subgraphs(DG, bubbles)
    # iterate over the different bubbles
    for listint, bubble_nodes in enumerate(bubbles):
        print("bubble:",bubble_nodes)
        # find the minimum and maximum for each bubble
        bubble_start, bubble_end, real_bubble,shared_reads= get_bubble_start_end(DG, bubble_nodes)
        # we have to figure out whether we are having a real bubble, this is indicated by the boolean value real_bubble
        #If the bubble is not a true bubble: Skip this "bubble"
        if not real_bubble:
            filter_count+=1
            print("Filtered ",filter_count," bubbles out")
            continue
        #get the reads which are in bubble_start and bubble_end
        #shared_reads=find_shared_reads(DG, bubble_start, bubble_end)
        # find the nodes which are on either path to be able to tell apart the paths on the bubble
        readlen_dict=get_path_reads_length(bubble_nodes, bubble_start, bubble_end, DG,shared_reads)
        #print("Listof normal",listofnormalnodes)
        print("min",bubble_start)
        print("max",bubble_end)
        print("readlendict",readlen_dict)
        # TODO compare the reads of both paths: If they differ by less than delta_len: Pop the bubble
        lengthdiff=None
        for node,length in readlen_dict.items():
            if not lengthdiff:
                lengthdiff = length
            else:
                lengthdiff = abs(lengthdiff-length)
        if lengthdiff < delta_len:
            DG=align_and_linearize_bubble_nodes(DG, readlen_dict, bubble_start, bubble_end,delta_len,all_reads)

    return DG
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


def simplifyGraph(DG, delta_len,all_reads,reads_for_isoforms):
    print("Simplifying the Graph (Merging nodes, popping bubbles)")
    # remove edges which yield self loops, not sure yet whether it makes sense to remove or if needed
    #print("self loops")
    #print(list(nx.selfloop_edges(DG)))
    list_of_cycles = find_repetative_regions(DG)
    print(list_of_cycles)
    #print(DG.edges())
    s_reads = DG.nodes["s"]['reads']
    print(type(s_reads))
    DG = merge_nodes(DG)
    find_and_pop_bubbles(DG, delta_len, all_reads,reads_for_isoforms)

    return DG
