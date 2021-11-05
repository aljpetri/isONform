import networkx as nx
from collections import Counter, namedtuple
from consensus import *
import matplotlib.pyplot as plt
from IsoformGeneration import *
from functools import cmp_to_key


#TODO: possibily better to use hybrid approach to finding node positions: calculations as well as finding the index of the minimizer in the original read
"""Helper function used to plot the graph. Taken from GraphGeneration.
    INPUT: DG   Directed Graph to plot
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


"""Helper method used to generate subgraphs for a list of subnodes
INPUT:      DG          Directed Graph to plot
            bubbles:    list of bubblenodelists to be plotted
"""


def generate_subgraphs(DG, bubbles):
    for bubble in bubbles:
        SG = DG.subgraph(bubble)
        draw_Graph(SG)


"""Helper method used to generate a subgraph for a list of nodes
INPUT:              DG          Directed Graph to plot
                    bubbles:    list of bubblenodelists to be plotted
"""


def generate_subgraph(DG, bubble):
    # for bubble in bubbles:
    SG = DG.subgraph(bubble)
    draw_Graph(SG)


"""Method to reduce the number of nodes in our graph
INPUT:      DG          Directed Graph 
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


"""
function to detect cycles in the graph (which denote repetitive regions)
INPUT:      DG             Directed Graph
OUTPUT: List_of_cycles     A list holding all cycles present in DG
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


"""
This function was taken from https://networkx.org/documentation/stable/_modules/networkx/algorithms/cycles.html#cycle_basis
As the original implementation relies on sets and therefore yields nondeterministic results, the code is altered to yield a deterministic behavior.
"""


def cycle_basis(G, root=None):
    """ Returns a list of cycles which form a basis for cycles of G.
    A basis for cycles of a network is a minimal collection of
    cycles such that any cycle in the network can be written
    as a sum of cycles in the basis.  Here summation of cycles
    is defined as "exclusive or" of the edges. Cycle bases are
    useful, e.g. when deriving equations for electric circuits
    using Kirchhoff's Laws.
    Parameters
    ----------
    G : NetworkX Graph
    root : node, optional
       Specify starting node for basis.
    Returns
    -------
    A list of cycle lists.  Each cycle list is a list of nodes
    which forms a cycle (loop) in G.
    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_cycle(G, [0, 1, 2, 3])
    >>> nx.add_cycle(G, [0, 3, 4, 5])
    >>> print(nx.cycle_basis(G, 0))
    [[3, 4, 5, 0], [1, 2, 3, 0]]
    Notes
    -----
    This is adapted from algorithm CACM 491 [1]_.
    References
    ----------
    .. [1] Paton, K. An algorithm for finding a fundamental set of
       cycles of a graph. Comm. ACM 12, 9 (Sept 1969), 514-518.
    See Also
    --------
    simple_cycles
    """
    # replaced gnodes to be implemented by a sorted list rather than a set
    gnodes = sorted(list(G.nodes()))
    cycles = []
    while gnodes:  # loop over connected components
        if root is None:
            root = gnodes.pop()
        stack = [root]
        pred = {root: root}
        used = {root: list()}  # we use a list instead of a set
        while stack:  # walk the spanning tree finding cycles
            z = stack.pop()  # use last-in so cycles easier to find
            zused = used[z]
            for nbr in G[z]:
                if nbr not in used:  # new node
                    pred[nbr] = z
                    stack.append(nbr)
                    used[nbr] = {z}
                elif nbr == z:  # self loops
                    cycles.append([z])
                elif nbr not in zused:  # found a cycle
                    pn = used[nbr]
                    cycle = [nbr, z]
                    p = pred[z]
                    while p not in pn:
                        cycle.append(p)
                        p = pred[p]
                    cycle.append(p)
                    cycles.append(cycle)
                    used[nbr].add(z)
        # altered to work with a list: we only remove the elements which were part of gnodes
        for i in pred:
            if i in gnodes:
                gnodes.remove(i)
        # gnodes -= sorted(list(pred)) #original code
        root = None
    return cycles

"""Helper method that is used to detect bubbles in our graph by using the cycle_basis method
INPUT: DG   our directed graph
OUTPUT: list_of_bubbles   a list containing all bubbles found (represented as lists themselves)
"""
def find_bubbles(DG):
    # draw_Graph(DG)
    # get undirected version of the graph to find bubbles
    UG = DG.to_undirected()
    # collect the bubbles in the graph (Bubbles denote possible mutations in the minimizers)
    # find all cycles in the undirected graph->bubbles
    # list_of_bubbles =nx.cycle_basis(UG)
    list_of_bubbles = cycle_basis(UG)
    return list_of_bubbles


"""
Helper method used to find all the reads, which are in startnode and may be part of a bubble(at least the next node is part of the bubble)
INPUT:          listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
                DG:         the directed graph we want to pop the bubble in
                startnode: The node deemed to be the starting node of the bubble (given as tuple)  
OUTPUT:         start_reads: set holding all reads which are supporting startnode and an out_edge which lies in the bubble
"""


def get_start_reads(DG, startnode, listofnodes):
    out_edges = DG.out_edges(startnode)
    start_reads = set()
    # print("getstartreads")
    for o_edge in out_edges:
        # print(o_edge)
        nextnode = o_edge[1]
        if nextnode in listofnodes:
            edge_infos = DG[o_edge[0]][nextnode]['edge_supp']
            # print("edge_infos",edge_infos)
            start_reads.update(edge_infos)
    # print(start_reads)
    return start_reads


"""
Helper method used to find all the reads which are in endnode and may be part of a bubble(at least the previous node is part of the bubble)
INPUT:          listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
                DG:         the directed graph we want to pop the bubble in
                endnode: The node deemed to be the ending node of the bubble (given as tuple)  
OUTPUT:         end_reads: set holding all reads which are supporting endnode and an in_edge which lies in the bubble
"""


def get_end_reads(DG, endnode, listofnodes):
    in_edges = DG.in_edges(endnode)
    end_reads = set()
    # print("getstartreads")
    for i_edge in in_edges:
        if i_edge[1] in listofnodes:
            end_reads.update(DG[i_edge[0]][i_edge[1]]['edge_supp'])
    # print(end_reads)
    return end_reads


"""Helper method for get_bubble_start_end: This method finds the minimum and maximum nodes in the bubble
The method iterates through the nodes in a bubble and collects all their out_nodes. For The minimum node there will not be an out_node
INPUT:      listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            DG:         the directed graph we want to pop the bubble in
OUTPUT:     startnode: The node deemed to be the starting node of the bubble (given as tuple)
            still_viable: Boolean value denoting whether the cycle is still a contender for being a bubble
            start_reads: All reads that are supporting the bubble_start_node as well as one edge towards a bubble node
            """
"""
#get all incoming edges for node x. If start of incoming edge in bubble-> this is not the source node of the bubble
def find_bubble_start(DG,listofnodes):
    #startnodes: a dictionary holding all nodes that
    startnodes= set()
    viable=False
    #final_startnodes=[]
    #startnode=None
    bubble_nodes=set(listofnodes)
    #we iterate over all nodes in our bubble
    for n in listofnodes:
        print("n",n)
        if len(startnodes)>1:
            print(startnodes)
            return startnodes
        in_edges=DG.in_edges(n)
        for edge in in_edges:
            #if the source of the edge is in the bubble
            if (edge[0] in bubble_nodes):
                print(edge[0])
                #set the indicator value to be TRUE
                #TODO:THis simply does not work->no double breaking out of loops
                break
        startnodes.add(n)
        print(startnodes)
    return startnodes"""
"""
Helper method for get_bubble_start_end: This method finds the maximum node in the bubble
The method iterates through the nodes in a bubble and collects all their in_nodes. For The maximum node there will not be an in_node
INPUT:      listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            DG:         the directed graph we want to pop the bubble in
OUTPUT:     endnode: The node deemed to be the end node of the bubble (given as tuple)
            still_viable: Boolean value denoting whether the cycle is still a contender for being a bubble
            end_reads: All reads that are supporting the bubble_end_node as well as one edge from a bubble node
            """
#get all incoming edges for node x. If start of incoming edge in bubble-> this is not the source node of the bubble
def find_bubble_start(DG,listofnodes):
    #startnodes: a dictionary holding all nodes that
    startnodes= {}
    viable=False
    final_startnodes=[]
    startnode=None
    #we iterate over all nodes in our bubble
    for n in listofnodes:
        #we use a boolean value to indicate whether n has an in_edge from any other bubble_node
        gets_list_node = False
        in_edges=DG.in_edges(n)
        for edge in in_edges:
            #if the source of the edge is in the bubble
            if (edge[0] in listofnodes):
                #set the indicator value to be TRUE
               gets_list_node=True
        startnodes[n]=gets_list_node
    print(startnodes)
    for key,value in startnodes.items():
        if not value:
            final_startnodes.append(key)
    print(final_startnodes)
    return final_startnodes
"""
#get all incoming edges for node x. If start of incoming edge in bubble-> this is not the source node of the bubble
def find_bubble_end(DG,listofnodes):
    #endnodes: a dictionary holding all nodes that
    endnodes= set()

    #final_startnodes=[]
    #startnode=None
    bubble_nodes=set(listofnodes)
    #we iterate over all nodes in our bubble
    for n in listofnodes:
        if len(endnodes)>1:
            return endnodes
        #we use a boolean value to indicate whether n has an in_edge from any other bubble_node
        out_edges=DG.out_edges(n)
        for edge in out_edges:
            #if the source of the edge is in the bubble
            if (edge[1] in bubble_nodes):
                #this node is clearly not an endnode. We therefore skip the rest of the edge iterations
                break
        endnodes.add(n)
    return endnodes"""
def find_bubble_end(DG,listofnodes):

    endnodes= {}
    viable=False
    final_endnodes=[]
    endnode=None
    for n in listofnodes:
        gets_list_node = False
        out_edges=DG.out_edges(n)
        for edge in out_edges:
            if (edge[1] in listofnodes):
               gets_list_node=True
        endnodes[n]=gets_list_node
    print(endnodes)
    for key,value in endnodes.items():
        if not value:
            final_endnodes.append(key)
    print(final_endnodes)

    return final_endnodes

"""

Returns the reads that are supporting startnode as well as endnode of the bubble 
INPUT:      listofnodes:A list containing all nodes which are part of the bubble
            DG:     The graph for which the element is to be found
OUTPUT:     min_element: the source node of the bubble(local source)
            max_element: the sink node of the bubble(local sink)
            contains_s: Boolean value indication whether the bubble contains the source node s
            contains_t: Boolean value indication whether the bubble contains the sink node t
"""

def get_shared_reads(DG,startnode,endnode):
    bubble_start_reads=set(DG.nodes[startnode]['reads'])
    bubble_end_reads = set(DG.nodes[endnode]['reads'])
    #real_bubble = real_bubble_e and real_bubble_s
    print("bubble_start_and_end")
    print(bubble_start_reads)
    print(bubble_end_reads)
    shared_reads = list(bubble_start_reads.intersection(bubble_end_reads))
    print("shared_reads:", shared_reads)
    return shared_reads



"""Helper method used to generate the consensus sequence for each bubble path
INPUT:      work_dir:               The work directory we are in
            consensus_attributes:   list of tuples containing the read_id and the positions that make up the limits of our subsequences
            reads:                  dictionary containing the actual reads (the sequences as well as their original ids)
            k_size:                 the length of parameter k needed to extract the correct subsequence length
OUTPUT:     spoa_ref:               the consensus sequence generated from the reads

"""
#TODO: for one read the positions do not fit(it has a shorter sequence for our consensus->somewhere we have the wrong position)
def generate_consensus_path(work_dir, consensus_attributes, reads, k_size):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    seq_infos={}
    endseqlist = []
    for i, (q_id, pos1, pos2) in enumerate(consensus_attributes, 3):
        print("read ",q_id)
        print("Printing full seq:")
        print( reads[q_id][1])
        if pos2 == 0:
            pos2 = len(reads[q_id][1]) - k_size
        seq = reads[q_id][1][pos1: pos2 + k_size]
        seq_infos[q_id]=(pos1,pos2+k_size,seq)
        # startseq=reads[q_id][1][pos1:pos1+k_size]
        endseq = reads[q_id][1][pos2:pos2 + k_size]
        # startseqlist.append(startseq)
        endseqlist.append(endseq)
        print(q_id, "from ", pos1, "to", pos2 + k_size, ": ", seq)
        reads_path.write(">{0}\n{1}\n".format(str(q_id) + str(pos1) + str(pos2), seq))
    # print("start",startseqlist)
    print("end", endseqlist)
    reads_path.close()
    # print(reads_path.name)
    # sys.exit()
    #print("seq_infos",seq_infos)
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
    print("spoa_ref", spoa_ref)
    return spoa_ref,seq_infos

"""Helper method used to test whether two sequences are close enough to pop their underlying bubble
INPUT:      cigar_string    a string denoting the cigar output of the alignment
            delta_len       a threshold which we use to distinguish between mutations and errors
OUTPUT:     good_to_pop     boolean value indicating our finding
"""
def parse_cigar_differences(cigar_string, delta_len):
    good_to_pop = True
    print("Now we are parsing....")
    print(cigar_string)
    for i,elem in enumerate(cigar_string):
        cig_len = elem[0]
        cig_type = elem[1]
        # all matches are absolutely fine
        if (cig_type != '=') and (cig_type != 'M'):
            print(cig_type)
            print(cig_len)
            # we have a non match, now we have to figure out whether this is due to an exon or indel (everything with len>delta_len is defined as exon)
            if cig_len > delta_len:
                    # we found an exon being the difference between the paths->no bubble popping feasible
                    good_to_pop = False
                    print(cigar_string)
    print(good_to_pop)
    return good_to_pop
#TODO fully comment this header



"""helper method used to add the end minimizer positions to a given dictionary
INPUT:      positions:              dictionary containing the read id as key and the positions(start,end) as tuple
            node_end_minimizers:    The end minimizer dictionary we would like to extend
            node:                   the node for which we would like to add the end_minimizer information 
"""
def add_end_minimizer_info(positions,node_end_minimizers,node):
    curr_dict={}
    for key,value in positions.items():
        curr_dict[key]=value[1]
    node_end_minimizers[node]=curr_dict


def get_dist_to_prev(DG,prev_node,curr_node):
    print("getting dist to prev")
    print("prev_reads",prev_node)
    print("this_reads",curr_node)
    prev_supp=DG.nodes[prev_node]['reads']
    curr_supp=DG.nodes[curr_node]['reads']
    print(prev_supp)
    print(curr_supp)
    intersect_supp = list(set(prev_supp).intersection(curr_supp))
    sum=0
    print(intersect_supp)
    #we did not yet come up with anything better than calculating the average distance for all reads supporting both nodes
    for i,read in enumerate(intersect_supp):
        print(i+1)
        prev_pos=prev_supp[read].end_mini_start
        curr_pos=curr_supp[read].end_mini_start
        sum+=(curr_pos-prev_pos)
    avg_dist=sum/(i+1)
    return avg_dist


def new_distance_to_start(DG,pnl_start,curr_node,consensus_log):
    print("consensus_log",consensus_log)
    print("currnode",curr_node)
    print(str(pnl_start))
    print("newdistancetostart")
    seq=consensus_log[pnl_start]
    node_seq=DG.nodes[curr_node]['end_mini_seq']
    print("Is ",node_seq ," in ",seq," ?")
    possible_pos=seq.find(node_seq)
    print(possible_pos)
    return possible_pos



"""
Helper method utilized by linearize bubbles which removes the edges we have to get rid of(all edges which connect 2 nodes which are part of the current bubble
    INPUT:      DG      Our directed Graph
                path_reads:         The reads which we are looking at to get the distances
                bubble_start        The start node of our bubble
                bubble_end:         The end node of our bubble
                path_nodes          A list of nodes which make up a path in our bubble
    OUTPUT:     edges_to_delete         A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
"""
def remove_edges(DG, path_reads, bubble_start, bubble_end, path_nodes,support_dict,consensus_log):
    #print("Path nodes", path_nodes)
    #print("PAth_reads", path_reads)
    #print("Bubble_start", bubble_start)
    #print("support_dict",support_dict)
    node_distances = {}
    # we store all the infos for the deleted edges into edges_to_delete
    edges_to_delete = {}

    # we have to delete the edges which connect the nodes that are inside the bubble
    for startnode_id, path_node_list in path_nodes.items():
        prev_node = bubble_start
        curr_node = bubble_start
        print("PathNodeList", path_node_list)
        entry = DG.get_edge_data(prev_node, startnode_id)
        edges_to_delete[prev_node,startnode_id]=entry
        if not(path_node_list):
            pnl_start=bubble_end
        else:
            pnl_start=path_node_list[0]
        for index, path_node in enumerate(path_node_list):
            #print(index, ", ", path_noded)
            #if pnl_start!=bubble_end:

            inter_dist = new_distance_to_start(DG, pnl_start, path_node, consensus_log)
            if inter_dist == -1:
                dist_to_prev = get_dist_to_prev(DG, prev_node, curr_node)
                #we found the distance to the previous node, however we are still missing the distance of the previous node to s
                if prev_node == bubble_start:
                    prev_to_start_dist = 0
                else:
                    prev_to_start_dist = node_distances[prev_node]
                dist = prev_to_start_dist+dist_to_prev
            else:
                dist = inter_dist
            node_distances[path_node] = dist
            if path_node != path_node_list[-1]:
                edges_to_delete[path_node, path_node_list[index + 1]] = DG[path_node][path_node_list[index + 1]]
            else:
                entry = DG.get_edge_data(path_node, bubble_end)
                edges_to_delete[path_node, bubble_end] = entry
            curr_node=path_node
    #print("Edges To Delete", edges_to_delete)
    for edge, edge_infos in edges_to_delete.items():
        #print(edge)
        DG.remove_edge(edge[0], edge[1])
    #print("node_distances ", node_distances)
    return edges_to_delete, node_distances



def get_avg_interval_length(DG,node):
    print("calculating average length for node ", node)
    node_support=DG.nodes[node]['reads']
    sum = 0
    i = 0
    for r_id, positions in node_support.items():
        print(positions)
        i += 1
        sum += (positions.end_mini_start-positions.start_mini_end)
        print(sum)
    finalresult=int(sum/(i+1))
    #finalresult=100000
    print("finalresult",finalresult)
    return finalresult



"""Helper method: Adds additional support values for the bubble_nodes(needed to get a consistent graph)
INPUT:      DG:
            nextnode:       the node of which we need the support
OUTPUT:     additional_support:  dictionary of additional support values
"""

def additional_node_support(DG, new_support,other_support, node_dist_dict, s_infos, node, prevnode_this_path,other_prevnode,seq_infos,bubble_start):
    Read_infos = namedtuple('Read_Infos',
                            'start_mini_end end_mini_start original_support')
    if other_prevnode==bubble_start:
        other_dist=0
    else:
        other_dist=node_dist_dict[other_prevnode]
    this_dist=node_dist_dict[node]
    print("other:dist",other_dist)
    #this will contain all reads that are being to be added to the node with their respective (virtual) positions
    additional_support = {}
    print("Node_dist", this_dist)
    print("node",node)
    print("s_infos",s_infos)
    this_reads=DG.nodes[node]['reads']
    avg_len=get_avg_interval_length(DG,node)
    #a lot of print statements used to debug
    print("avg_len",avg_len)
    print("new_support ",new_support)
    print("prevnode_this_path", prevnode_this_path)
    print("prevnode_other_path",other_prevnode)

    for r_id in new_support:
        print("r_id",r_id)
        #figure out whether the r_id is part of this bubble_path
        if r_id not in this_reads:
            # what we do if the read meets the bubble and is not in this_reads
            # figure out if r_id is also represented by bubble_start
            if r_id not in s_infos.keys():
                # r_id is not in s_infos (meaning the read is not supporting bubble_start), if it would be in s_infos, we do not have to recompute anything

                #get the position of the read in previous node
                #calculate relative distance to start (relative to previous node)
                #start=position+relative dist
                #end=start+avg_length
                print("r_id in other path")
                print("other prevnode",other_prevnode)
                previous_other_path_reads=DG.nodes[other_prevnode]['reads']
                if r_id in previous_other_path_reads:
                    pos_info_tup=previous_other_path_reads[r_id]
                    prev_end=pos_info_tup.end_mini_start
                    relative_dist= int(this_dist-other_dist)
                    print("relative_dist",relative_dist)
                    newend=prev_end+relative_dist
                    newstart=int(newend-avg_len)
                    additional_support[r_id] = Read_infos(newstart, newend, False)
            #the read is in this path so we do not have to do a lot here
            else:

                # print("r_id",r_id)
                print("S_Infos", s_infos[r_id])
                #start_pos=s_infos[r_id].start_mini_end
                end_pos=s_infos[r_id].end_mini_start
                #start_pos, end_pos = s_infos[r_id]
                newend=int(this_dist) + end_pos
                #print("newstart",newstart)
                newstart=newend-avg_len
                additional_support[r_id] = Read_infos(newstart, newend,False)

    print("Additional node_supp after", additional_support)
    return additional_support


"""Helper method used to merge two dicts with the values being lists 
INPUT: dict1: the first dict
       dict2: the second dict
OUTPUT: merged_dict: dictionary containing all keys in dict1 and dict2
"""
def merge_two_dicts(dict1, dict2):
    merged_dict = {}
    #print("dict1", dict1)
    #print("dict2", dict2)
    for key, value in dict1.items():
        merged_dict[key] = value
    for key2, val2 in dict2.items():
        if key2 not in merged_dict:
            merged_dict[key2] = val2
    return merged_dict

#TODO:this may not be sufficient
def sort_by_node_label(seq):
    seq_duplicates = sorted(seq, key=lambda x: x[1])
    return seq_duplicates
def find_next_node_of_this_path(actual_node_distances,path1):
    print("path1",path1)
    print("actual_node_distances(AND)", actual_node_distances)
    for node in actual_node_distances:

        print("Node_which_we_found",node)
        if node[0] in path1:
            print("found ",node, "in both ")
            return node[0]
"""adds new edges to the graph 
INPUT:          DG      our graph object
                edges_to_delete     list of edges that were deleted from the graph
                bubble_start        start of the bubble in question
                bubble_end          end of the bubble in question
                all_shared          all reads that are shared between start and end node
                path_nodes          all nodes that are part of the bubble which are neither start nor endnode
                actual_node_distances   dictionary having the path nodes and their respective distance to bubble_start


the function does not output anything but updated the graph object DG 
"""

    # print("FinalEdges",DG.edges(data=True))
def find_real_nextnode(nextnode1,nextnode2,node_dist,bubble_end):
    if nextnode1 == bubble_end:
        return nextnode2
    elif nextnode2 == bubble_end:
        return nextnode1
    elem1 = node_dist[nextnode1]
    elem2 = node_dist[nextnode2]
    if elem1 < elem2:
        return nextnode1
    elif elem2 < elem1:
        return nextnode2
    else:
        print("We have to come up with something!")
        print("is ",nextnode1," < ",nextnode2,"? As we got ",elem1,"== ",elem2)
        return nextnode1
    #TODO: we can just take the 0th element of the path bc this is what we did before. Nextnode1 is located at pos 1
def get_next_node(path, overall_nextnode,bubble_end):
    print("get_ndex_node")
    print(str(type(path)))
    print(path)
    thisidx=path.index(overall_nextnode)
    print("thisidx",thisidx)
    if len(path)<2:
        nextnode=bubble_end
    else:
        nextnode=path[thisidx+1]
    return nextnode
def prepare_adding_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes,node_dist,seq_infos):  # ,node_dist):
    counter = 0
    path1 = []
    path2 = []
    edge_params={}
    #actual_node_distances=[(k, v) for k, v in node_dist.items()]
    #print("actual_node_distances", actual_node_distances)
    #actual_node_distances.sort(key=lambda tup: tup[1])
    print("Bubble_end_node", bubble_end)
    print("EdgestoDelete", edges_to_delete)
    #print("Allnodeslen", len(actual_node_distances))
    #TODO: assert that allnodes and pathnodes have the same order->less errors
    #print("allnodes", actual_node_distances)
    print("pathnodes", path_nodes)
    #sort_by_dist_to_s(actual_node_distances,DG)

    #actual_node_distances=sort_by_node_label(actual_node_distances)

    #print("actual_node_distance s1", actual_node_distances)
    #sorted(actual_node_distances, key=cmp_to_key(compare))
    #print("actual_node_distance s2", actual_node_distances)

    # we assign both paths to variables to make them easier accessible.
    #TODO: change to use queue instead of list from collection import queue ->popleft
    for id, path in path_nodes.items():
        print(id, path)
        if counter == 0:
            for p in path:
                path1.append(p)
        else:
            for p in path:
                path2.append(p)
        counter = counter + 1
    #this dictionary contains the reads with their respected positions that are added to the nodes
    new_node_supp_dict = {}
    print("path1",path1)
    print("path2",path2)

    s_infos = DG.nodes[bubble_start]['reads']
    print("s INFOS", s_infos, "of node ", bubble_start)
    prevnode1 = bubble_start
    prevnode2 = bubble_start

    if path1:
        nextnode1 = path1[0]
        print("nextnode1",nextnode1)

    else:
        nextnode1=bubble_end
        #print("Creepy ands", actual_node_distances)
    if path2:
        nextnode2 = path2[0]
        print("nextnode2", nextnode2)
    else:
        nextnode2 = bubble_end
        #print("Creepy ands", actual_node_distances)
    prevnode = bubble_start
    while path1 or path2:
        print("P1",path1)
        print("P2", path2)
        new_edge_supp1 = edges_to_delete[prevnode1, nextnode1]['edge_supp']
        print("NES1 from ",prevnode1," to ",nextnode1,": ",new_edge_supp1)
        new_edge_supp2 = edges_to_delete[prevnode2, nextnode2]['edge_supp']
        print("NES1 from ", prevnode2, " to ", nextnode2, ": ", new_edge_supp2)
        full_edge_supp = new_edge_supp1 + new_edge_supp2
        print("FES: ", full_edge_supp)
        overall_nextnode=find_real_nextnode(nextnode1,nextnode2,node_dist,bubble_end)
        #TODO:seperate the finding of next node from the actual adding of edges to make it easier to hunt down bugs
        if overall_nextnode in path1:
            nextnode1=get_next_node(path1, overall_nextnode,bubble_end)
            print("nextnode1",nextnode1)
            additional_supp = additional_node_support(DG, new_edge_supp2, new_edge_supp1, node_dist, s_infos, overall_nextnode,
                                                      prevnode1, prevnode2, seq_infos, bubble_start)
            print("additionalNodeSupport1 for",overall_nextnode)
            print("path 1 before remove", path1)
            path1.remove(overall_nextnode)
            print("path 1 after remove",path1)
            prevnode1=overall_nextnode
            curr_node=prevnode1
        elif overall_nextnode in path2:
            nextnode2=get_next_node(path2,overall_nextnode,bubble_end)
            additional_supp = additional_node_support(DG, new_edge_supp1, new_edge_supp2, node_dist, s_infos, overall_nextnode,
                                                      prevnode2, prevnode1, seq_infos, bubble_start)
            print("additionalNodeSupport2 for", overall_nextnode)
            print("path 2 before remove", path2)
            path2.remove(overall_nextnode)
            print("path 2 after remove", path2)
            prevnode2=overall_nextnode
            curr_node = prevnode2
        else:
            print("ERRORRRRRR", overall_nextnode," neither in ",nextnode1," nor in ",nextnode2)
        old_node_supp = DG.nodes[overall_nextnode]['reads']
        new_node_supp_dict[overall_nextnode] = merge_two_dicts(additional_supp, old_node_supp)
        nx.set_node_attributes(DG, new_node_supp_dict, "reads")
        edge_params[prevnode,overall_nextnode]=full_edge_supp
        DG.add_edge(prevnode, overall_nextnode, edge_supp=full_edge_supp)

        print("Adding edge from ", prevnode, "to ", overall_nextnode)
        prevnode=curr_node

    new_edge_supp1 = edges_to_delete[prevnode1, bubble_end]['edge_supp']
    new_edge_supp2 = edges_to_delete[prevnode2, bubble_end]['edge_supp']
    full_edge_supp = new_edge_supp1 + new_edge_supp2
    DG.add_edge(prevnode, bubble_end, edge_supp=full_edge_supp)
    print("Adding edge fa from ", prevnode, "to ", bubble_end)
    # this is the main part of the linearization. We iterate over all_nodes and try to find out which path the nodes belong to.
    # This info is needed as we need the current state ob both paths to add the correct edge_support and node_support to the graph

""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
           delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """


def align_bubble_nodes(delta_len, all_reads, consensus_infos, work_dir, k_size):
    print("aligning")
    consensus_list = []
    consensus_log= {}
    seq_infos={}
    for path_node, consensus_attributes in consensus_infos.items():
        print("consensus", consensus_attributes)
        print("path_node",path_node)
        if len(consensus_attributes) > 1:
            con,seq_infos_from_fun = generate_consensus_path(work_dir, consensus_attributes, all_reads, k_size)

            consensus_log[path_node]=con
            seq_infos.update(seq_infos_from_fun)
            consensus_list.append(con)
        else:
            (q_id, pos1, pos2) = consensus_attributes[0]
            print("consensus_attributes", q_id, ", ", pos1, ", ", pos2)
            con = all_reads[q_id][1][pos1: pos2 + k_size]
            seq_infos[q_id]=(pos1,pos2+k_size,con)
            consensus_log[path_node]=con
        consensus_list.append(con)
    print(consensus_list)
    print("consensus_log",consensus_log)
    consensus1 = consensus_list[0]
    consensus2 = consensus_list[1]
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(consensus1, consensus2,
                                                                                       match_score=2,
                                                                                       missmatch_penalty=-2,
                                                                                       opening_penalty=3, gap_ext=1)
    print(s1_alignment)
    print(s2_alignment)
    print(cigar_string)
    print(cigar_tuples)
    print("cigar done")
    good_to_pop = parse_cigar_differences(cigar_tuples, delta_len)
    cigar_alignment=(s1_alignment,s2_alignment)
    consensuses=(consensus1,consensus2)
    return good_to_pop,cigar_alignment,seq_infos,consensus_log


def get_path_nodes(cycle, bubble_start, DG):
    # find the nodes which are directly connected to min_node (min_node_out) and to max_node(max_node_in) this enables the finding of which reads we have to compare
    min_edges = DG.out_edges(bubble_start)
    min_node_out = []
    for edge in min_edges:
        if edge[1] in cycle:
            min_node_out.append(edge[1])
    return min_node_out


"""Helper method for find_bubbles: This method figures out which reads belong to which part from bubble source to bubble sink
INPUT:      cycle:          A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element:    The node deemed to be the starting node of the bubble (given as tuple)
            max_element:    The node deemed to be the end node of the bubble (given as tuple)
            DG:             the directed graph we want to pop the bubble in
            contains_s:     A boolean value denoting whether the cycle contains node "s"
            contains_t:     A boolean value denoting whether the cycle contains node "t"
OUTPUT: path_starts:        Dictionary holding the information about which reads support the source node (later used as initial_supp in test_path_viability)
"""


def get_path_reads(DG, min_node_out, shared_reads):
    path_starts = {}

    for out_node in min_node_out:
        print(out_node)
        inter_out_readlist = DG.nodes[out_node]['reads']
        print(inter_out_readlist)
        #out_readlist = [ for i in inter_out_readlist]
        out_path_reads_list = []
        for read in inter_out_readlist:
            if read in shared_reads:
                out_path_reads_list.append(read)
        path_starts[out_node] = out_path_reads_list
    print("Pathstarts")
    print(path_starts)

    return path_starts


"""Helper method which tests if a found path is viable (meaning if the structure denoted as bubble is an actual bubble or not)
This is done by following the edges and trying to verify that a path is actually supported by at least one read from bubble_start to bubble_end
INPUT:          DG:                 The directed graph
                path_start:         the start node of the path
                initial_support:    a list of reads supporting the edge bubble_start,path_start
                cyclce:             List of nodes which make up the "bubble"
                bubble_end:         The node deemed to be the sink of the "bubble"
OUTPUT:         is_viable:          A boolean value telling us whether there is a set of reads that consistently supports the path
                visited_nodes:      a list of nodes in the order they were visited
                intersect supp:    the set of reads that are supporting the full path 

"""


def test_path_viability(DG, path_start, initial_support, cycle,bubble_start, bubble_end):
    curr_support = initial_support
    curr_node = path_start
    visited_nodes = []
    intersect_supp = []
    #this is a somewhat iffy situation: bubble_start is directly connected to bubble_end for this path
    if path_start==bubble_end:
        print("This is something!")
        print("visited_nodes", visited_nodes)
        #we have to access the edge connecting the nodes to find the correct support of this path
        curr_support=DG[bubble_start][curr_node]["edge_supp"]
        return (True, visited_nodes,curr_support)
    # we only iterate until we have reached bubble_end
    while curr_node != bubble_end:
        visited_nodes.append(curr_node)
        further_node = False
        curr_out = DG.out_edges(curr_node)
        # we have to test for all possible out edges to be sure to find the path
        for out_edge in curr_out:
            next_node = out_edge[1]

            # we only continue investigating the next_node if it is part of the bubble
            if next_node in cycle:
                # get the edges which support the edge to next_node
                next_support = DG[curr_node][next_node]["edge_supp"]
                # figure out whether there is any overlap between the current support and next_support
                intersect_supp = list(set(curr_support).intersection(next_support))
                #print("intersect supp",intersect_supp)
                # if we have an overlap, this means that at least one read from our initial_support also supports this edge
                if intersect_supp:
                    # we have found a node to continue
                    further_node = True
                    # set the node to be our current_node
                    curr_node = next_node
                    # update the set of still supporting reads
                    curr_support = intersect_supp
                    #print("curr_support")
                    # as we only expect one node to be viable: break the for loop to continue in the next node
                    break
        # we did not find any node we could continue with->our path is not the path of a real bubble
        if not further_node:
            is_viable = False
            return is_viable, visited_nodes, intersect_supp
    # we were able to reach bubble_end->we just confirmed that we have a bubble path
    is_viable = True
    print("Visited_nodes ", visited_nodes)
    #print("intersect_supp", intersect_supp)
    return (is_viable, visited_nodes, intersect_supp)


"""Helper method used for the bubble detection. This function finds the nodes which mark the start points of the bubble paths
INPUT:      Cycle:          List of nodes which make up the bubble
            bubble_start:   Source node of the bubble
            bubble_end:     Sink node of the bubble
            DG:             The directed graph
            shared_reads:   list of reads that are supporting the source node as well as the sink node of the bubble
OUTPUT:     path_starts:    

"""


def get_path_starts(cycle, bubble_start, DG, shared_reads):
    # We want to find the nodes, which denote the start points for each path(as we have to find out which reads are in which path)
    min_node_out = get_path_nodes(cycle, bubble_start, DG)
    print("minnodeout",min_node_out)
    # Now we want to get the actual reads for each path
    path_starts = get_path_reads(DG, min_node_out, shared_reads)
    # iterate over all shared reads and get the pathlength for each
    print("Path_starts ", path_starts)
    return path_starts


"""Helper method for find_bubbles: This method figures out which reads belong to which part from bubble source to bubble sink
INPUT:      cycle:          A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element:    The node deemed to be the starting node of the bubble (given as tuple)
            max_element:    The node deemed to be the end node of the bubble (given as tuple)
            DG:             the directed graph we want to pop the bubble in
            contains_s:     A boolean value denoting whether the cycle contains node "s"
            contains_t:     A boolean value denoting whether the cycle contains node "t"
"""


def get_consensus_positions(r_ids, bubble_start, bubble_end, DG, shared_reads):
    read_list = []
    max_node_infos = DG.nodes[bubble_end]['reads']
    min_node_infos = DG.nodes[bubble_start]['reads']
    for r_id in r_ids:
        if r_id in shared_reads:

            bubble_end_pos = max_node_infos[r_id]
            bubble_start_pos = min_node_infos[r_id]
            start_of_bubble = bubble_start_pos.end_mini_start
            end_of_bubble = bubble_end_pos.end_mini_start
            entry = (r_id, start_of_bubble, end_of_bubble)
            read_list.append(entry)
    return read_list


"""
function which finds the bubbles and if feasible pops them in the graph
INPUT:
    DG:         Our directed graph 
    delta_len   Maximum length difference for which the paths are merged ie the bubble is popped
    all_reads   dictionary containing all the reads(string sequence) and their ids


OUTPUT:
    bubble_state, popable_bubbles,no_pop_tup_list"""

#TODO:split up into two functions this would require a lot of work and shared parameters. MAybe better to leave as is
def find_poppable_bubbles(DG, delta_len, all_reads, work_dir, k_size, bubbles):
    popable_bubbles = []
    bubble_state = []
    no_pop_tup_list=[]
    # just for debugging and curiosity reasons: We introduce an integer counting the number of pseudobubbles (not true bubbles)
    filter_count = 0
    Supported = namedtuple('Supported', 'bubble_nodes, bubble_supported')
    # iterate over the different bubbles
    for listint, bubble_nodes in enumerate(bubbles):
        viable_bubble = False
        #print("bubble:", bubble_nodes)
        # find the local source and local sink for each bubble, if we find more than one source or sink, the bubble is not a supported bubble. However supported bubbles additionally have to be supported by at least two reads
        bubble_start_list =find_bubble_start(DG,bubble_nodes)
        bubble_end_list =find_bubble_end(DG,bubble_nodes)
        bubble_start=None
        bubble_end=None
        #if the bubbles have more start and/or endnodes than one it is not viable->Don't continue investigating
        if not (len(bubble_start_list) ==1 and len(bubble_end_list)==1):
            filter_count += 1
            supported = Supported(bubble_nodes, False)
            bubble_state.append(supported)
            print("Filtered ", filter_count, " bubbles out")
            continue
        else:
            bubble_start=bubble_start_list[0]
            bubble_end=bubble_end_list[0]
        print("bubble_start", bubble_start)
        print("bubble_end", bubble_end)
        #get the reads which are supporting start and endnode of the bubble. Those are the reads we are most interested in (all or a subset of those reads actually support the bubble)
        shared_reads=get_shared_reads(DG,bubble_start,bubble_end)

        # we have to figure out whether we are having a sufficiently supported_bubble, meaning at least two reads are supporting bubble_start and bubble_end
        # If the bubble is not a supported bubble: Skip the further investigation
        if len(shared_reads)<2:
            filter_count += 1
            supported = Supported(bubble_nodes, False)
            bubble_state.append(supported)
            print("Filtered ", filter_count, " bubbles out")
            continue

        path_starts = get_path_starts(bubble_nodes, bubble_start, DG, shared_reads)#TODO: this may not
        pre_consensus_infos = {}
        path_nodes_dict = {}
        no_viab_bubble = False
        support_dict = {}
        #print("Pathstarts", path_starts)
        #key: path_start: the id of the node which denotes the start of our path, value: reads supporting the node
        for key, value in path_starts.items():
            (is_viable_path, path_nodes, support) = test_path_viability(DG, key, value, bubble_nodes, bubble_start,bubble_end)
            print("path nodes for ",key,": ",path_nodes)
            path_nodes_dict[key] = path_nodes
            support_dict[key] = support
            if not (is_viable_path):
                no_viab_bubble = True
                break
            elif not bubble_end=='t' and path_nodes:
                no_viab_bubble=False
            #print("supported", support)
            print("bubble_end",bubble_end)
            # this denotes the start of the consensus alignment phase: We find the positions for the consensus generation ->we
            # find the nodes which are on either path to be able to tell apart the paths on the bubble
            read_list = get_consensus_positions(value, bubble_start, bubble_end, DG, support)
            pre_consensus_infos[key] = read_list

        if no_viab_bubble:

            supported = Supported(bubble_nodes, viable_bubble)
            bubble_state.append(supported)
            filter_count += 1
            continue
        good_to_pop,cigar,seq_infos,consensus_info_log = align_bubble_nodes(delta_len, all_reads, pre_consensus_infos, work_dir, k_size)
        print("good_to_pop",good_to_pop)
        if good_to_pop:
            # print("Path_nodes type",type(path_nodes))
            Linearization_infos = namedtuple('Linearization_Infos',
                                             'pre_consensus_infos bubble_start bubble_end path_nodes support_dict seq_infos, consensus_log')
            l_infos = Linearization_infos(pre_consensus_infos, bubble_start, bubble_end, path_nodes_dict, support_dict, seq_infos,consensus_info_log)
            Pop_infos = namedtuple('Pop_infos', 'bubble_nodes,lin_infos')
            pop_infos = Pop_infos(bubble_nodes, l_infos)
            popable_bubbles.append(pop_infos)
            viable_bubble = good_to_pop
        else:
            print("bubble_end_non_viable", bubble_end)
            print("bubble not poppable",bubble_nodes)
            print(path_nodes_dict)
            print(support_dict)
            print(pre_consensus_infos)
            no_pop_tuple=(bubble_nodes,cigar)
            no_pop_tup_list.append(no_pop_tuple)
        supported = Supported(bubble_nodes, viable_bubble)
        print("popable",supported)
        bubble_state.append(supported)
    return bubble_state, popable_bubbles,no_pop_tup_list
    # generate_subgraph(DG, bubble_nodes)

"""Helper method that filters the poppable_bubbles to only contain non-overlapping bubbles, as overlapping bubbles yield inconsistencies in our graph
INPUT:      poppable_bubbles:      a list of bubbles that are poppable
OUTPUT:     new_poppable:          a list of bubbles that do not overlap->ready for popping
"""
def find_disjoint_bubbles(poppable_bubbles):
    #print("Filter touched", poppable_bubbles)
    known_nodes = []
    new_poppable = []
    for bubble in poppable_bubbles:
        is_touched = False
        for node in bubble.bubble_nodes:
            if node in known_nodes:
                is_touched = True
        if not is_touched:
            known_nodes.extend(bubble.bubble_nodes)
            new_poppable.append(bubble)
    return new_poppable


"""Actual linearization process of our bubbles
 The main idea of this method is: 1. Get all the nodes which are in both paths and calculate their avg distance to bubble_start
                                  2. Sort the nodes by distance
                                  3. Add the edges connecting the nodes in the indicated order
        INPUT:      DG      Our directed Graph
                consensus_infos:    
                bubble_start        The start node of our bubble
                bubble_end:         The end node of our bubble
                path_nodes:          A list of nodes which make up a path in our bubble
                edges_to_delete:     A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
                support_dict:
"""


def linearize_bubble(DG, pre_consensus_infos, bubble_start, bubble_end, path_nodes, support_dict,seq_infos,consensus_log):
    print("seq_infos",seq_infos)
    print("time for linearization")
    # edges_to_delete: A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
    edges_to_delete,node_dist = remove_edges(DG, pre_consensus_infos, bubble_start, bubble_end, path_nodes,support_dict,consensus_log)
    print("bubbles",bubble_start, ", ",bubble_end)
    #print("Allnodes after sorting", all_nodes)
    #add_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes,node_dist,seq_infos)  # , node_dist)
    prepare_adding_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes, node_dist, seq_infos)
    print("Popped bubble ", path_nodes)

"""Wrapper method that is used for popping a list of bubbles calling linearize_bubble
INPUT:      DG      our graph object
            bubble_pop the list of bubbles that are to be popped

"""
def pop_bubbles(DG, bubble_pop):
    print("Popping_bubbles")
    for bubble in bubble_pop:
        print("bubble_to_pop", bubble)
        linearize_bubble(DG, bubble.lin_infos.pre_consensus_infos, bubble.lin_infos.bubble_start, bubble.lin_infos.bubble_end,
                         bubble.lin_infos.path_nodes, bubble.lin_infos.support_dict,bubble.lin_infos.seq_infos,bubble.lin_infos.consensus_log)#,minimizerdict,all_reads)
        # print("Edges aflin", DG.edges(data=True))
        print("bubble popped", bubble)
        print("Current State of Graph:")
        print(DG.nodes(data=True))
        print(DG.edges(data=True))



"""Helper method used to transform a list of lists into a set of frozenset(makes it easer to compare)
INPUT: input: a list of lists( list of bubbles)
OUTPUT: result: a set of frozensets( set holding a fr
"""


def list_list_to_set_frozenset(input):
    result = set()
    for i in input:
        tup_i = tuple(i)
        set_i = frozenset(tup_i)
        result.add(set_i)
    return result

def eliminate_already_analysed_bubbles(old_bubbles, new_bubbles):
    new_bubbles_set = list_list_to_set_frozenset(new_bubbles)
    old_bubbles_set = list_list_to_set_frozenset(old_bubbles)
    bubbles_to_analyse = new_bubbles_set.difference(old_bubbles_set)
    bubbles_list = list(bubbles_to_analyse)
    bubbles_non_frozen = [list(x) for x in bubbles_list]
    return bubbles_non_frozen



def eliminate_already_analysed_bubbles_inefficient(old_bubbles, new_bubbles):
    bubbles_list = []
    for old in old_bubbles:
        sorted(old)
    for new in new_bubbles:
        sorted(new)
        if new not in old_bubbles:
            bubbles_list.append(new)
    return bubbles_list
"""Helper method used to extend the list of bubbles that has been analysed but not been deemed poppable (old_bubbles)
INPUT:      old_bubbles
            bubble_state
"""
def extend_unpoppable(old_bubbles,poppable_state):
    for newbubble in poppable_state:
        print(newbubble)
        print(newbubble.bubble_supported)
        if not newbubble.bubble_supported:
            print("newbubble",newbubble.bubble_nodes)
            old_bubbles.append(newbubble.bubble_nodes)


"""The backbone of bubble popping: the bubbles are detected and filtered to only yield poppable bubbles. Those are popped then.
INPUT:      DG:         our directed graph
            delta_len:  for all differences longer than delta len we have a structural difference while all differences shorter are deemed to be errors
            all_reads:  list of tuples containing all the reads and (what we need) their sequences
            work_dir:   the current work directory
            k_size:     the length of our k_mers

"""
#TODO:find new name for cycle_basis thingie
def bubble_popping_routine(DG, delta_len, all_reads, work_dir, k_size):
    # find all bubbles present in the graph which are to be popped
    more_to_pop=True
    popped_bubbles = []
    old_bubbles=[]
    no_pop_list = []
    print("Initial state of the graph")
    print(DG.nodes(data=True))
    print(DG.edges(data=True))
    while more_to_pop:
        # find all bubbles in the new state of the graph
        bubbles = find_bubbles(DG)
        bubbles.sort(key=len)
        nr_bubbles = len(bubbles)
        print("Found " + str(nr_bubbles) + " bubbles in our graph")
        print("Bubbles: ", bubbles)
        # find the difference of the new bubbles and the old bubbles (meaning we are only interested in bubbles that have not yet been analysed)
        # TODO: both of the two following lines are working however the more efficient version does yield different order (set)
        # bubbles=eliminate_already_analysed_bubbles(old_bubbles, new_bubbles)
        bubbles = eliminate_already_analysed_bubbles_inefficient(old_bubbles, bubbles)

        bubble_state, poppable_bubbles, no_pop_tup_list = find_poppable_bubbles(DG, delta_len, all_reads, work_dir, k_size,
                                                                         bubbles)
        print("Bub_state",bubble_state)
        # extending old_bubbles: This list of bubbles consists of all bubbles deemed to be not poppable
        extend_unpoppable(old_bubbles,bubble_state)
        print("Bubble_state",bubble_state)
        print("old_bubbles",old_bubbles)
        #this list serves as a debug tool: It stores the cigar string of the bubble alignment so we can see whether all bubbles are correctly identified
        no_pop_list.extend(no_pop_tup_list)
        print("WhileBubbles", bubbles)
        print("poppable_bubbles", poppable_bubbles)
        disjoint_poppable_bubbles = find_disjoint_bubbles(poppable_bubbles)
        print("disjoint_poppable_bubbles", disjoint_poppable_bubbles)
        # pop the poppable bubbles by linearizing them
        pop_bubbles(DG, disjoint_poppable_bubbles)
        popped_bubbles.extend(disjoint_poppable_bubbles)
        print("Actually popped", popped_bubbles)
        print("Current state of the graph")
        print(DG.nodes(data=True))
        print(DG.edges(data=True))
        if not poppable_bubbles:
            more_to_pop=False
    print("NO_POP",no_pop_list)

    # draw_Graph(DG)


"""Overall method used to simplify the graph
During this method: - Bubbles are identified and if possible popped     
                    - Nodes are merged        
INPUT:  DG: Directed Graph before simplifications,
        max_bubblesize: maximum number of elements making up a bubble
        delta_len: Maximum length differences in between two reads
OUTPUT: DG: Graph after simplification took place    """


# TODO: Overall: add relative distances to all edges in the graph/ make sure all edges have relative distances
def simplifyGraph(DG, delta_len, all_reads, work_dir, k_size):
    #TODO: add true minimizers
    print("Simplifying the Graph (Merging nodes, popping bubbles)")
    list_of_cycles = find_repetative_regions(DG)
    print(list_of_cycles)
    #print("Current State of Graph:")
    #print(DG.nodes(data=True))
    #print(DG.edges(data=True))
    #draw_Graph(DG)
    bubble_popping_routine(DG, delta_len, all_reads, work_dir, k_size)
    print("Popping bubbles done")
    merge_nodes(DG)