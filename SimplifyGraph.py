import networkx as nx
from collections import Counter, namedtuple
from consensus import *
import matplotlib.pyplot as plt
from IsoformGeneration import *

"""Helper function used to plot my graph. Taken from GraphGeneration.
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
INPUT:      DG          Directed Graph to plot
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
function to find cycles in the graph (which denote repetitive regions)
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
Helper method used to find all the reads which are in startnode and may be part of a bubble(at least the next node is part of the bubble)
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


def find_bubble_start(DG, listofnodes):
    # print("fbslistofnodes")
    # print(listofnodes)
    still_viable = True
    out = set()
    outedges = []
    for node in listofnodes:
        # print(node)
        out_edges_all = DG.out_edges(node)
        # print("Out edges all",out_edges_all)
        for out_edge in out_edges_all:
            if (out_edge[1] in listofnodes):
                out.add(out_edge[1])
                outedges.append(out_edge)
    set_a = set(listofnodes)
    set_b = out
    subtraction = set_a - set_b
    list_of_strings = [str(s) for s in subtraction]
    # print("LOSstart", list_of_strings)
    if len(list_of_strings) > 1:
        still_viable = False
    bubble_reads = set()
    for out_edge in outedges:
        # print(out_edge)
        edge_infos = DG[out_edge[0]][out_edge[1]]["edge_supp"]
        for entry in edge_infos:
            bubble_reads.add(entry)
    startnode = " ".join(list_of_strings)
    # print("startnode",startnode)
    start_reads = get_start_reads(DG, startnode, listofnodes)
    return still_viable, startnode, start_reads


"""Helper method for get_bubble_start_end: This method finds the maximum node in the bubble
The method iterates through the nodes in a bubble and collects all their in_nodes. For The maximum node there will not be an in_node
INPUT:      listofnodes:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            DG:         the directed graph we want to pop the bubble in
OUTPUT:     endnode: The node deemed to be the end node of the bubble (given as tuple)
            still_viable: Boolean value denoting whether the cycle is still a contender for being a bubble
            end_reads: All reads that are supporting the bubble_end_node as well as one edge from a bubble node
            """


def find_bubble_end(DG, listofnodes):
    real_bubble = True
    in_nodes = set()
    inedges = []
    bubble_end_reads = set()
    for node in listofnodes:
        in_edges_all = DG.in_edges(node)
        for in_edge in in_edges_all:
            if in_edge[0] in listofnodes:
                in_nodes.add(in_edge[0])
                inedges.append(in_edge)
    set_a = set(listofnodes)
    set_b = in_nodes
    subtraction = set_a - set_b
    print(subtraction)
    list_of_strings = [str(s) for s in subtraction]
    # print("LOSend",list_of_strings)
    if len(list_of_strings) > 1:
        real_bubble = False
    print("inedges",inedges)
    for in_edge in inedges:
        print(in_edge)
        print(list_of_strings)
        # if(in_edge[1] in list_of_strings):
        if (in_edge[1] == list_of_strings[0]):
            edge_infos = DG[in_edge[0]][in_edge[1]]["edge_supp"]
            for entry in edge_infos:
                bubble_end_reads.add(entry)
    endnode = " ".join(list_of_strings)
    end_reads = get_end_reads(DG, endnode, listofnodes)
    return real_bubble, endnode, end_reads


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
    print(listofnodes)
    real_bubble_s, startnode, bubble_start_reads = find_bubble_start(DG, listofnodes)
    real_bubble_e, endnode, bubble_end_reads = find_bubble_end(DG, listofnodes)
    real_bubble = real_bubble_e and real_bubble_s
    shared_reads = list(bubble_start_reads.intersection(bubble_end_reads))
    print("shared_reads:", shared_reads)
    return startnode, endnode, real_bubble, shared_reads


def all_same(items):
    return all(x == items[0] for x in items)


# TODO: we get an empty sequence into the consensus generation. WHY?!. It seems we
def generate_consensus_path(work_dir, consensus_attributes, reads, k_size):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    startseqlist = []
    endseqlist = []
    for i, (q_id, pos1, pos2) in enumerate(consensus_attributes, 3):
        print("Printing full seq:", reads[q_id][1])
        if pos2 == 0:
            pos2 = len(reads[q_id][1]) - k_size
        seq = reads[q_id][1][pos1: pos2 + k_size]
        # startseq=reads[q_id][1][pos1:pos1+k_size]
        endseq = reads[q_id][1][pos2:pos2 + k_size]
        # startseqlist.append(startseq)
        endseqlist.append(endseq)
        print(q_id, "from ", pos1, "to", pos2 + k_size, ": ", seq)
        reads_path.write(">{0}\n{1}\n".format(str(q_id) + str(pos1) + str(pos2), seq))
    # print("start",startseqlist)
    print("end", endseqlist)
    # assert all_same(startseqlist)
    # assert all_same(endseqlist)
    reads_path.close()
    # print(reads_path.name)
    # sys.exit()
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
    print("spoa_ref", spoa_ref)
    return spoa_ref


def parse_cigar_differences(cigar_string, delta_len):
    good_to_pop = True
    for elem in cigar_string:
        cig_len = elem[0]
        cig_type = elem[1]
        # all matches are absolutely fine
        if (cig_type != '=') and (cig_type != 'M'):
            # we have a nonmatch, now we have to figure out whether this is due to an exon or indel (everything with len>delta_len is defined as exon)
            if cig_len > delta_len:
                # we found an exon being the difference between the paths->no bubble popping feasible
                good_to_pop = False
    return good_to_pop


def collect_bubble_nodes(path_nodes, consensus_infos, DG, support_dict):
    all_nodes = []
    print("Consensus_infos", consensus_infos)
    print(path_nodes)
    all_shared = []
    for startnode_id, path_node_list in path_nodes.items():
        shared_reads = support_dict[startnode_id]
        print("shared", shared_reads)
        all_shared.extend(shared_reads)
        print(path_node_list)
        for path_node in path_node_list:
            path_node_infos = DG.nodes[path_node]['reads']
            distance = 0
            print("path_node ", path_node)
            for i, s_read in enumerate(shared_reads):
                pos_tuple = path_node_infos[s_read]
                print("pos_tuple")
                print(pos_tuple[0])
                distance += pos_tuple[0]
            final_distance = distance / (i + 1)
            print("final_distance", final_distance)
            all_nodes_tuple = (path_node, final_distance)
            all_nodes.append(all_nodes_tuple)
    print("all_nodes", all_nodes)
    return all_nodes, all_shared


# TODO: fix bug in this method (currently main bug)
"""Helper method which is used to find the distance of a node to the start node. This is done by using the read information for both nodes
    INPUT:          DG      Our directed Graph
                    bubble_start        The start node of our bubble
                    path_node           The node for which we would like to retreive the distance
                    path_reads:         The reads which we are looking at to get the distances
    OUTPUT:         dist_to_start:      mean distance between the two nodes (we only look at the path_reads)
"""


def get_distance_to_start(DG, bubble_start, path_node, path_reads,support_dict):
    print("bubblestart", bubble_start)
    print("path_node", path_node)
    print("Pathnodeinfos", DG.nodes[path_node]['reads'])
    print("support_dict",support_dict)
    #actual_reads = [a_tuple[0] for a_tuple in path_reads]
    #print("actual_reads", actual_reads)
    sum = 0
    # we iterate over actual_reads (the read ids we got from path_reads)
    for i, read in enumerate(support_dict):
        start_infos = DG.nodes[bubble_start]['reads']
        print("startinfos", start_infos)
        start_tuple = start_infos[read]
        print("starttuple", start_tuple)
        start_pos = start_tuple[1]
        print("start_pos", start_pos)
        print("path_node", path_node)
        node_infos = DG.nodes[path_node]['reads']
        # print("node_infos",node_infos)
        end_tuple = node_infos[read]
        end_pos = end_tuple[0]
        print("end_pos", end_pos)
        dist = end_pos - start_pos
        print("dist", dist)
        sum += dist
        print("sum", sum)
    dist_to_start = sum / (i + 1)
    print("disttostart", dist_to_start)
    return dist_to_start


"""
Helper method utilized by linearize bubbles which removes the edges we have to get rid of(all edges which connect 2 nodes which are part of the current bubble
    INPUT:      DG      Our directed Graph
                path_reads:         The reads which we are looking at to get the distances
                bubble_start        The start node of our bubble
                bubble_end:         The end node of our bubble
                path_nodes          A list of nodes which make up a path in our bubble
    OUTPUT:     edges_to_delete         A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
"""


def remove_edges(DG, path_reads, bubble_start, bubble_end, path_nodes,support_dict):
    print("Path nodes", path_nodes)
    print("PAth_reads", path_reads)
    print("Bubble_start", bubble_start)
    print("support_dict",support_dict)
    node_distances = {}
    tup_list = []
    # we store all the infos for the deleted edges into edges_to_delete
    edges_to_delete = {}
    # we have to delete the edges which connect the nodes that are inside the bubble
    for startnode_id, path_node_list in path_nodes.items():
        print("PathNodeList", path_node_list)
        print("startnode_id", startnode_id)
        print("startnode", bubble_start)
        pnl_start = path_node_list[0]
        print("pnl_start", pnl_start)
        # we mark the edge bubble_start->path_node_list[0] and add its infos to edges_to_delete
        # print(DG.edges)
        edges_to_delete[bubble_start, pnl_start] = DG[bubble_start][pnl_start]
        dist = get_distance_to_start(DG, bubble_start, pnl_start, path_reads[startnode_id],support_dict[startnode_id])
        node_distances[pnl_start] = dist
        tup=(pnl_start,dist)
        if not tup in tup_list:
            tup_list.append(tup)
        for index, path_node in enumerate(path_node_list):
            print(index, ", ", path_node)
            dist = get_distance_to_start(DG, bubble_start, path_node, path_reads[startnode_id],support_dict[startnode_id])
            node_distances[path_node] = dist
            tup = (path_node, dist)
            if not tup in tup_list:
                tup_list.append(tup)
            if path_node != path_node_list[-1]:
                edges_to_delete[path_node, path_node_list[index + 1]] = DG[path_node][path_node_list[index + 1]]
            else:
                print("Bubble EndReads", bubble_end)
                print(DG.nodes[bubble_end])
                entry = DG.get_edge_data(path_node, bubble_end)
                edges_to_delete[path_node, bubble_end] = entry

    print("Edges To Delete", edges_to_delete)
    for edge, edge_infos in edges_to_delete.items():
        print(edge)
        DG.remove_edge(edge[0], edge[1])
    print("node_distances ", node_distances)
    return edges_to_delete, node_distances,tup_list


"""Helper method: Adds additional support values for the bubble_nodes(needed to get a consistent graph)
INPUT:      DG:
            nextnode:       the node of which we need the support
OUTPUT:     additional_support:  dictionary of additional support values
"""


def additional_node_support(DG, new_support, this_dist, s_infos, node, nextnode_other_path):
    additional_support = {}
    print("Node_dist", this_dist)
    # print(DG.nodes(data=True))
    # print(node_dist[nextnode])
    # this_dist=node_dist[node]
    for r_id in new_support:
        # TODO: This if/else is not fully functional yet, but helps to improve algo stability
        if r_id not in s_infos.keys():
            other_reads = DG.nodes[nextnode_other_path]['reads']
            # print("other_reads",other_reads)
            # print(r_id," not in in S")
            additional_support[r_id] = other_reads[r_id]
        else:
            # print("r_id",r_id)
            # print("S_Infos", s_infos[r_id])
            start_pos, end_pos = s_infos[r_id]
            additional_support[r_id] = (int(this_dist + end_pos), int(this_dist + end_pos))
    print("Additional node_supp after", additional_support)
    return additional_support


"""Helper method used to merge two dicts with the values being lists 
INPUT: dict1: the first dict
       dict2: the second dict
OUTPUT: merged_dict: dictionary containing all keys in dict1 and dict2
"""


def merge_two_dicts(dict1, dict2):
    merged_dict = {}
    print("dict1", dict1)
    print("dict2", dict2)
    for key, value in dict1.items():
        # if key in dict2:
        # print("mergeval",value)
        # val_list=value
        # print("Key",dict2[key])
        # val_list.extend(dict2[key])
        # merged_dict[key]=val_list
        # else:
        merged_dict[key] = value
    for key2, val2 in dict2.items():
        if key2 not in merged_dict:
            merged_dict[key2] = val2
    return merged_dict


# TODO: For this method we have a bug: Our paths are not sorted accordingly to the allnodes. This results in unfindable nodes
def add_edges(DG, all_nodes, edges_to_delete, bubble_start, bubble_end, all_shared, path_nodes,actual_node_distances):  # ,node_dist):
    counter = 0
    path1 = []
    path2 = []
    print("Bubble_end_node", bubble_end)
    print("EdgestoDelete", edges_to_delete)
    print("Allnodeslen", len(actual_node_distances))
    print("allnodes", actual_node_distances)
    print("pathnodes", path_nodes)
    # we assign both paths to variables to make them easier accessible.
    for id, path in path_nodes.items():
        print(id, path)
        if counter == 0:
            path1 = path
        else:
            path2 = path
        counter = counter + 1
    new_node_supp_dict = {}
    print(path1)
    print(path2)
    s_infos = DG.nodes[bubble_start]['reads']
    print("s INFOS", s_infos, "of node ", bubble_start)
    prevnode1 = bubble_start
    prevnode2 = bubble_start
    nextnode1 = path1[0]
    nextnode2 = path2[0]
    prevnode = bubble_start
    finalprev = None
    # this is the main part of the linearization. We iterate over all_nodes and try to find out which path the nodes belong to.
    # This info is needed as we need the current state ob both paths to add the correct edge_support and node_support to the graph
    for nodetup in actual_node_distances:
        print("Nodetuple", nodetup)
        node = nodetup[0]
        print("Node", node)
        reads_for_node = []
        reads_for_node.extend(all_shared)
        print("P1", path1)
        print("P2", path2)
        new_edge_supp1 = edges_to_delete[prevnode1, nextnode1]['edge_supp']
        print("ES 1 from ", prevnode1, " to ", nextnode1)
        print("NES1", new_edge_supp1)
        new_edge_supp2 = edges_to_delete[prevnode2, nextnode2]['edge_supp']
        full_edge_supp = new_edge_supp1 + new_edge_supp2
        full_edge_supp_final = []
        [full_edge_supp_final.append(x) for x in full_edge_supp if x not in full_edge_supp_final]
        print("full_edge_support_final:", full_edge_supp_final)
        node_supp = []
        real_pos = nodetup[1]
        # if the next node is from path1: pop the node out of path1 and set nextnode1 to the
        if node in path1:
            prevnode1 = path1.pop(0)
            finalprev = prevnode1
            print("P1", path2)
            if len(path1) < 1:
                nextnode1 = bubble_end
                print("Nextnode1", nextnode1)
                if len(path2) < 1:
                    nextno2 = bubble_end
                else:
                    nextno2 = path2[0]
            else:
                if len(path2) < 1:
                    nextno2 = bubble_end
                else:
                    nextno2 = path2[0]
                nextnode1 = path1[0]

                print(nextnode1)
            node_supp = DG.nodes[nextnode1]['reads']
            print("Edge_supp_bubble_pop", new_edge_supp2)

            additional_supp = additional_node_support(DG, new_edge_supp2, real_pos, s_infos, node, nextno2)
            print("additional_support", additional_supp)
            old_node_supp = DG.nodes[node]['reads']
            new_node_supp_dict[node] = merge_two_dicts(additional_supp, old_node_supp)
            # node_supp.update(additional_supp)

        elif node in path2:
            prevnode2 = path2.pop(0)
            finalprev = prevnode2
            print("P2", path2)
            print("prevnode2", prevnode2)
            if len(path2) < 1:
                nextnode2 = bubble_end
                print("Nextnode2", nextnode2)
                if len(path1) < 1:
                    nextno1 = bubble_end
                else:
                    nextno1 = path1[0]
            else:
                if len(path1) < 1:
                    nextno1 = bubble_end
                else:
                    nextno1 = path1[0]
                nextnode2 = path2[0]
                print("Nextnode2", nextnode2)
            # node_supp = DG.nodes[nextnode2]['reads']
            print("new_edge_supp", new_edge_supp1)
            additional_supp = additional_node_support(DG, new_edge_supp1, real_pos, s_infos, node, nextno1)
            old_node_supp = DG.nodes[node]['reads']
            new_node_supp_dict[node] = merge_two_dicts(additional_supp, old_node_supp)
            # node_supp.update(additional_supp)
        else:
            print("Error, ", node, " neither in ", path1, " nor in ", path2)
        print("node_supp", node_supp)
        print("ES 2 from ", prevnode2, " to ", nextnode2)
        print("NES2", new_edge_supp2)
        # print(DG.edges(data=True))
        DG.add_edge(prevnode, node, edge_supp=full_edge_supp_final)
        print("Adding edge from ", prevnode, "to ", node)
        # print(DG.edges(data=True))
        print(node)
        prevnode = node
        print("nodetup", nodetup)
    print("updated nodes")
    print(new_node_supp_dict)
    # print(DG.nodes(data=True))
    nx.set_node_attributes(DG, new_node_supp_dict, "reads")
    print("updated nodes")
    # print(DG.nodes(data=True))
    if len(path2) > 0:
        prevnode2 = path2.pop()
    if len(path1) > 0:
        prevnode1 = path1.pop()
    # end of for loop - in the next lines we add the information for the last edge between all_nodes[-1] and bubble_end (sink_node)
    new_edge_supp1 = edges_to_delete[prevnode1, bubble_end]['edge_supp']
    new_edge_supp2 = edges_to_delete[prevnode2, bubble_end]['edge_supp']
    full_edge_supp = new_edge_supp1 + new_edge_supp2
    full_edge_supp_final = []
    [full_edge_supp_final.append(x) for x in full_edge_supp if x not in full_edge_supp_final]
    DG.add_edge(prevnode, bubble_end, edge_supp=full_edge_supp_final)
    print("Adding node from ", prevnode, "to ", bubble_end)
    # print("FinalEdges",DG.edges(data=True))


""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
           delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """


def align_and_linearize_bubble_nodes(delta_len, all_reads, consensus_infos, work_dir, k_size):
    consensus_list = []
    for path_node, consensus_attributes in consensus_infos.items():
        print("consensus", consensus_attributes)
        if len(consensus_attributes) > 1:
            con = generate_consensus_path(work_dir, consensus_attributes, all_reads, k_size)
        else:
            (q_id, pos1, pos2) = consensus_attributes[0]
            print("consensus_attributes", q_id, ", ", pos1, ", ", pos2)
            con = all_reads[q_id][1][pos1: pos2 + k_size]
        consensus_list.append(con)
    print(consensus_list)
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
    good_to_pop = parse_cigar_differences(cigar_tuples, delta_len)

    return good_to_pop


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
INPUT:          DG:                 The directed graph
                path_start:         the start node of the path
                initial_support:    a list of reads supporting the edge bubble_start,path_start
                cyclce:             List of nodes which make up the "bubble"
                bubble_end:         The node deemed to be the sink of the "bubble"
OUTPUT:         is_viable:          A boolean value telling us whether there is a set of reads that consistently supports the path
                visited_nodes:      a list of nodes in the order they were visited
                intersect supp:    the set of reads that are supporting the full path 

"""


def test_path_viability(DG, path_start, initial_support, cycle, bubble_end):
    curr_support = initial_support
    curr_node = path_start
    visited_nodes = []
    intersect_supp = []
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
                # if we have an overlap, this means that at least one read from our initial_support also supports this edge
                if intersect_supp:
                    # we have found a node to continue
                    further_node = True
                    # set the node to be our current_node
                    curr_node = next_node
                    # update the set of still supporting reads
                    curr_support = intersect_supp
                    print("curr_support")
                    # as we only expect one node to be viable: break the for loop to continue in the next node
                    break
        # we did not find any node we could continue with->our path is not the path of a real bubble
        if not further_node:
            is_viable = False

            return is_viable, visited_nodes, intersect_supp
    # we were able to reach bubble_end->we just confirmed that we have a bubble path
    is_viable = True
    print("Visited_nodes ", visited_nodes)
    print("intersect_supp", intersect_supp)
    return (is_viable, visited_nodes, intersect_supp)


"""Helper method used for the bubble detection. This function finds the nodes which mark the start points of the bubble paths
INPUT:      Cycle:          List of nodes which make up the bubble
            bubble_start:   Source node of the bubble
            bubble_end:     Sink node of the bubble
            DG:             The directed graph
            shared_reads:   list of reads that are supporting the source node as well as the sink node of the bubble
OUTPUT:     path_starts:    

"""


def get_path_starts(cycle, bubble_start, bubble_end, DG, shared_reads):
    # We want to find the nodes, which denote the start points for each path(as we have to find out which reads are in which path)
    min_node_out, max_node_in = get_path_nodes(cycle, bubble_start, bubble_end, DG)
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


def get_path_reads_length(r_ids, bubble_start, bubble_end, DG, shared_reads):
    id_counter = 0
    read_length = 0
    read_list = []
    max_node_infos = DG.nodes[bubble_end]['reads']
    print("max_node_infos")
    print(max_node_infos)
    min_node_infos = DG.nodes[bubble_start]['reads']
    print("min_node_infos")
    print(min_node_infos)
    print("rids", r_ids)
    print("shared_reads", shared_reads)
    for r_id in r_ids:
        if r_id in shared_reads:
            bubble_end_pos = max_node_infos[r_id]
            bubble_start_pos = min_node_infos[r_id]
            if bubble_start_pos == (-1, -1):
                continue
            start_of_bubble = bubble_start_pos[1]
            end_of_bubble = bubble_end_pos[0]
            entry = (r_id, start_of_bubble, end_of_bubble)
            print(entry)
            read_list.append(entry)
            read_length += (end_of_bubble - start_of_bubble)
            id_counter += 1
    path_length = read_length / id_counter
    return path_length, read_list


"""
function which finds the bubbles and if feasible pops them in the graph
INPUT:
    DG:         The directed graph in which we pop the bubbles
    delta_len   Maximum length difference for which the paths are merged ie the bubble is popped
    all_reads   dictionary containing all the reads(string sequence) and their ids


OUTPUT:
    DG:         Directed Graph containing the popped bubbles"""


def filter_bubbles(DG, delta_len, all_reads, work_dir, k_size, bubbles):
    popable_bubbles = []
    bubble_state = []

    # just for debugging and curiosity reasons: We introduce an integer counting the number of pseudobubbles (not true bubbles)
    filter_count = 0
    Popable = namedtuple('Popable', 'bubble_nodes, bubble_poppable')
    # iterate over the different bubbles
    for listint, bubble_nodes in enumerate(bubbles):
        viable_bubble = False
        print("bubble:", bubble_nodes)
        # find the minimum and maximum for each bubble
        bubble_start, bubble_end, real_bubble, shared_reads = get_bubble_start_end(DG, bubble_nodes)
        # we have to figure out whether we are having a real bubble, this is indicated by the boolean value real_bubble
        # If the bubble is not a true bubble: Skip this "bubble"
        if not real_bubble:
            filter_count += 1
            popable = Popable(bubble_nodes, viable_bubble)
            bubble_state.append(popable)
            print("Filtered ", filter_count, " bubbles out")
            continue
        # find the nodes which are on either path to be able to tell apart the paths on the bubble
        path_starts = get_path_starts(bubble_nodes, bubble_start, bubble_end, DG, shared_reads)
        readlen_dict = {}
        consensus_infos = {}
        path_nodes_dict = {}
        no_viab_bubble = False
        support_dict = {}
        print("Pathstarts", path_starts)
        for key, value in path_starts.items():
            (is_viable_path, path_nodes, support) = test_path_viability(DG, key, value, bubble_nodes, bubble_end)
            path_nodes_dict[key] = path_nodes
            support_dict[key] = support
            if not (is_viable_path) or not path_nodes:
                no_viab_bubble = True
                break
            print("supported", support)
            path_length, read_list = get_path_reads_length(value, bubble_start, bubble_end, DG, support)
            print("read_list (consensus:_infos", read_list)
            readlen_dict[key] = path_length
            consensus_infos[key] = read_list
        if no_viab_bubble:
            popable = Popable(bubble_nodes, viable_bubble)
            bubble_state.append(popable)
            filter_count += 1
            print("Filtered ", filter_count, " bubbles out")
            continue
        print("min", bubble_start)
        print("max", bubble_end)
        # print("readlendict",readlen_dict)
        # compare the reads of both paths: If they differ by less than delta_len: Pop the bubble
        # lengthdiff=None
        # for node,length in readlen_dict.items():
        #    if not lengthdiff:
        #        lengthdiff = length
        #        print(lengthdiff)
        #    else:
        #        lengthdiff = abs(lengthdiff-length)
        # if lengthdiff < delta_len:
        print("consensus_infos_before_linearization", consensus_infos)
        print("Linearizing bubble ", bubble_nodes)
        # generate_subgraph(DG, bubble_nodes)
        good_to_pop = align_and_linearize_bubble_nodes(delta_len, all_reads, consensus_infos, work_dir, k_size)
        if good_to_pop:
            # print("Path_nodes type",type(path_nodes))
            Linearization_infos = namedtuple('Linearization_Infos',
                                             'consensus bubble_start bubble_end path_nodes support_dict')
            l_infos = Linearization_infos(consensus_infos, bubble_start, bubble_end, path_nodes_dict, support_dict)
            Pop_infos = namedtuple('Pop_infos', 'bubble_nodes,lin_infos')
            pop_infos = Pop_infos(bubble_nodes, l_infos)
            popable_bubbles.append(pop_infos)
            viable_bubble = good_to_pop

        popable = Popable(bubble_nodes, viable_bubble)
        bubble_state.append(popable)
    return bubble_state, popable_bubbles
    # generate_subgraph(DG, bubble_nodes)


def filter_touched(poppable_bubbles):
    print("Filter touched", poppable_bubbles)
    known_nodes = []
    new_poppable = []
    for bubble in poppable_bubbles:
        print(bubble)
        print(known_nodes)
        is_touched = False
        for node in bubble.bubble_nodes:
            print(node)
            if node in known_nodes:
                is_touched = True
        if not is_touched:
            known_nodes.extend(bubble.bubble_nodes)
            new_poppable.append(bubble)
            print("Added", new_poppable)
    # known_nodes=poppable_bubbles[0].bubble_nodes
    print("New", new_poppable)
    return new_poppable


"""Actual linearization process of our bubbles
        INPUT:      DG      Our directed Graph
                consensus_infos:    
                bubble_start        The start node of our bubble
                bubble_end:         The end node of our bubble
                path_nodes:          A list of nodes which make up a path in our bubble
                edges_to_delete:     A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
                support_dict:
"""


def linearize_bubble(DG, consensus_infos, bubble_start, bubble_end, path_nodes, support_dict):
    # The main idea of this method is: 1. Get all the nodes which are in both paths and calculate their avg distance to bubble_start
    #                                 2. Sort the nodes by distance
    #                                 3. Add the edges connecting the nodes in the indicated order
    print("time for linearization")
    # print("Edges befdel", DG.edges(data=True))
    # edges_to_delete: A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
    edges_to_delete,node_dist,actual_node_distances = remove_edges(DG, consensus_infos, bubble_start, bubble_end, path_nodes,support_dict)

    # print("Edges afdel", DG.edges(data=True))
    all_nodes, all_shared = collect_bubble_nodes(path_nodes, consensus_infos, DG, support_dict)
    print("Allnodes before sorting", all_nodes)
    actual_node_distances.sort(key=lambda tup: tup[1])
    print("Allnodes after sorting", all_nodes)
    add_edges(DG, all_nodes, edges_to_delete, bubble_start, bubble_end, all_shared, path_nodes,actual_node_distances)  # , node_dist)
    print("Popped bubble ", path_nodes)


def pop_bubbles(DG, bubble_pop):
    print("Popping_bubbles")
    for bubble in bubble_pop:
        print("bubble_to_pop", bubble)
        linearize_bubble(DG, bubble.lin_infos.consensus, bubble.lin_infos.bubble_start, bubble.lin_infos.bubble_end,
                         bubble.lin_infos.path_nodes, bubble.lin_infos.support_dict)
        # print("Edges aflin", DG.edges(data=True))
        print("bubble popped", bubble)


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


# TODO: make sure this is correctly implemented
def eliminate_already_analysed_bubbles_inefficient(old_bubbles, new_bubbles):
    bubbles_list = []
    for old in old_bubbles:
        sorted(old)
    for new in new_bubbles:
        sorted(new)
        if new not in old_bubbles:
            bubbles_list.append(new)
    return bubbles_list


def list_list_append(original, new):
    for i in new:
        original.append(i)
    return original


def delete_not_popped_bubbles_from_old_bubbles(bubble_state, old_bubbles):
    print("Poppable_bub", bubble_state)
    print("Old_bub before", old_bubbles)
    for pop_state in bubble_state:
        print("pop_state", pop_state)
        is_poppable = pop_state.bubble_poppable
        if is_poppable:
            old_bubbles.remove(pop_state.bubble_nodes)

    print("Old_bub after", old_bubbles)


def bubble_popping_routine(DG, delta_len, all_reads, work_dir, k_size):
    # find all bubbles present in the graph which are to be popped
    bubbles = find_bubbles(DG)
    bubbles.sort(key=len)
    nr_bubbles = len(bubbles)
    print("Found " + str(nr_bubbles) + " bubbles in our graph")
    print("Bubbles: ", bubbles)
    popped_bubbles = []
    # setting old_bubbles: This list of bubbles consists of all bubbles analysed by the algorithm used to filter out already analysed bubbles
    old_bubbles = bubbles
    # filter out bubbles which are not poppable
    bubble_state, poppable_bubbles = filter_bubbles(DG, delta_len, all_reads, work_dir, k_size, bubbles)
    # as long as we find new bubbles (that we have not yet analysed)
    while poppable_bubbles:
        print(DG.nodes(data=True))
        print(DG.edges(data=True))
        print("WhileBubbles", bubbles)
        #draw_Graph(DG)
        # filter out bubbles which are not poppable
        #bubble_state, poppable_bubbles = filter_bubbles(DG, delta_len, all_reads, work_dir, k_size, bubbles)
        print("Bub_state", bubble_state)
        print("poppable_bubbles", poppable_bubbles)
        # draw_Graph(DG)
        print("old_bubbles before")
        delete_not_popped_bubbles_from_old_bubbles(bubble_state, old_bubbles)
        poppable_bubbles_filtered = filter_touched(poppable_bubbles)
        print("Poppable_bubbles_touched", poppable_bubbles_filtered)
        popped_bubbles.extend(poppable_bubbles_filtered)
        # pop the popable bubbles by linearizing them
        pop_bubbles(DG, poppable_bubbles_filtered)
        # find all bubbles in the new state of the graph
        new_bubbles = find_bubbles(DG)
        new_bubbles.sort(key=len)
        nr_bubbles = len(bubbles)
        print("Found " + str(nr_bubbles) + " bubbles in our graph")
        print("Bubbles (not here): ", bubbles)
        # find the difference of the new bubbles and the old bubbles (meaning we are only interested in bubbles that have not yet been analysed)
        # TODO: both of the two following lines are working however the more efficient version does yield different order (set)
        # bubbles=eliminate_already_analysed_bubbles(old_bubbles, new_bubbles)
        bubbles = eliminate_already_analysed_bubbles_inefficient(old_bubbles, new_bubbles)
        nr_bubbles = len(bubbles)
        print("Found " + str(nr_bubbles) + " bubbles in our graph (Eliminated)")
        print("Bubbles: ", bubbles)
        # print("bnf",bubbles_non_frozen)
        list_list_append(old_bubbles, bubbles)
        print("Popped: ", popped_bubbles)
        # filter out bubbles which are not poppable
        bubble_state, poppable_bubbles = filter_bubbles(DG, delta_len, all_reads, work_dir, k_size, bubbles)
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
    print("Simplifying the Graph (Merging nodes, popping bubbles)")
    list_of_cycles = find_repetative_regions(DG)
    print(list_of_cycles)
    # TODO: rewrite this code: We want to introduce a method which appoints relative edge distances (avg(start8endnode))-avg(end(startnode)) to the whole graph
    # TODO: Use relative distances for all computaations.
    bubble_popping_routine(DG, delta_len, all_reads, work_dir, k_size)
    print("Popping bubbles done")
    merge_nodes(DG)