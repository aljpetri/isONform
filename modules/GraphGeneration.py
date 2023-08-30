import tempfile
import shutil
import networkx as nx
import pickle
import os

from collections import namedtuple


Read_infos = namedtuple('Read_Infos', 'start_mini_end end_mini_start original_support')

"""IsONform script containing the methods used to generate the Directed Graph from the Intervals coming from the Weighted Interval Scheduling Problem
Author: Alexander Petri
The main method in this script was used for debugging therefore is not used during IsONforms actual run.
"""

def isCyclicUtil(DG, nodes_dict, node):
    """
    Used to find out whether adding 'node' led to the graph having a cycle
    Implemented according to https://www.geeksforgeeks.org/detect-cycle-in-a-graph/
    """
    # Mark current node as visited and
    # adds to recursion stack
    CNode = namedtuple('CNode', 'visited recStack')
    cnode = CNode(True, True)
    nodes_dict[node] = cnode
    # Recur for all neighbours
    # if any neighbour is visited and in
    # recStack then graph is cyclic
    for out_edge in DG.out_edges(node):
        neighbour = out_edge[1]
        if neighbour not in nodes_dict:
            cnode = CNode(False, False)
            nodes_dict[neighbour] = cnode
        if not nodes_dict[neighbour].visited:
            if isCyclicUtil(DG, nodes_dict, neighbour):
                return True
        elif nodes_dict[neighbour].recStack:
            return True
    # The node needs to be popped from
    # recursion stack before function ends
    prev_visited = nodes_dict[node].visited
    nodes_dict[node] = CNode(prev_visited, False)
    return False


def dfs(node, path, DG, visited):
    visited.add(node)
    path.add(node)
    for neighbor in DG.neighbors(node):
        if neighbor in visited:
            if neighbor in path:
                return True
        else:
            if dfs(neighbor, path, DG, visited):
                return True
    path.pop()
    return False

def cycle_finder(DG, start_node):
    """
    Method used to detect cycles in our graph structure
    INPUT:
        DG:         the graph structure
        start_node: a node that should be located in a cycle
    OUTPUT:
        bool:       indicator whether a cycle has been found
    """
    nodes_dict = {}
    #CNode = namedtuple('CNode', 'visited recStack', defaults=(False, False))
    #cnode = CNode(False, False)
    #nodes_dict[start_node] = cnode
    if isCyclicUtil(DG, nodes_dict, start_node):
        return True
    return False

def cycle_finder_old(DG, start_node):
    """Method used to detect cycles in our graph structure
    INPUT:      DG:         the graph structure
                start_node: a node that should be located in a cycle
    OUTPUT:     bool:       indicator whether a cycle has been found
    Implemented with the help of chatGPT
    """
    visited = set()
    path = [start_node]
    while len(path) > 0:
        current_node = path.pop(0)
        print("cyc_find_curr_node",current_node)
        visited.add(current_node)
        for neighbor in DG.neighbors(current_node):
            if neighbor in visited: #and neighbor != "t":
                print(neighbor," found again")
                # We have found a cycle!
                return True
            else:
                path.append(neighbor)
    # No cycle was found
    return False


def add_prior_read_infos(inter, r_id, prior_read_infos, name, k):
    """
    Function to add read information of the current interval to prior_read_infos, if they do not exist there yet
    INPUT:      inter:              the current interval holding the information
                r_id:               the read we are looking at
                prior_read_infos:   information about other reads
                name:               name of the node we appoint the information to
                k:                  minimizer length
    OUTPUT:     prior_read_infos:   information about other reads extended by this intervals infos
    """
    read_id = inter[3][slice(0, len(inter[3]),
                             3)]  # recover the read id from the array of instances which was delivered with all_intervals_for_graph
    start_coord = inter[3][slice(1, len(inter[3]),
                                 3)]  # recover the start coordinate of an interval from the array of instances
    end_coord = inter[3][slice(2, len(inter[3]), 3)]
    # iterate through the retreived information and store it in prior_read_infos only for subsequent reads
    for i, r in enumerate(read_id):
        if not r <= r_id:
            start = start_coord[i] + k
            end = end_coord[i]
            tuple_for_data_structure = (r, start, end)
            if tuple_for_data_structure not in prior_read_infos:
                prior_read_infos[tuple_for_data_structure] = name


def convert_array_to_hash(info_array):
    """Helper method to convert the array delivered with all_intervals_for_graph into a hash value to more efficiently look up occurence
    This method additionally deletes the first three entries of the array as they contain the infos about this interval occurence, which changes in between instances
    INPUT:  info_array: The array, which was delivered with the interval to indicate where the interval occurs in other reads
    OUTPUT: hash(tup):  A has value of the tuple the shortened interval was converted into. This hash makes it easy to see whether the interval is already present in the read
    """
    # preprocessing: Delete the first three elements from the array, as they contain the information about this occurrence
    for x in range(0, 3):
        if info_array:
            info_array.pop(0)

    tup = hash(tuple(info_array))
    return tup


def record_cycle(current_read_state, cycle_start):
    """Function to get all nodes which are part of a cycle
        INPUT:  current_read_state: The current state of known_intervals[r_id-1] (known_intervals of the current read)
                cycle_start:        The last occurence of the repeating node before the cycle starts
        OUTPUT: cycle_nodes:        A list of entries indicating which nodes are in the cycle, having the following form: (startpos, node_id, endpos)
    """
    cycle_nodes = []
    indices = [i for i, tupl in enumerate(current_read_state) if tupl[1] == cycle_start]
    index = indices[0]
    for element in range(index, len(current_read_state)):
        cycle_nodes.append(current_read_state[element])

    return cycle_nodes, current_read_state[index]


def find_next_node(thisstartpos, info_array, known_cycles, current_read_state, k, cycle_start, previous_end, DG,
                   delta_len):
    """Helper method for generateGraphfromIntervals
    INPUT:  thisstartpos:          The startposition of the cycle
            info_array:
            known_cycles:
            current_read_state:
            k:
            cycle_start:
            previous_end:
            DG:
    OUTPUT: found:
            found_next_node:

    """
    intermediate_cycles = []
    indices = [i for i, tupl in enumerate(current_read_state) if tupl[1] == cycle_start]
    index = indices[0]
    for key, value in known_cycles.items():
        known_cycle_rid = key[0]
        node_appearance = value[0]
        for i in range(0, len(info_array) - 2):
            if info_array[i] == known_cycle_rid and info_array[i + 1] == node_appearance[0] - k and info_array[i + 2] == \
                    node_appearance[2]:
                intermediate_cycles.append(value)
    found = False
    found_next_node = None
    for cyc in intermediate_cycles:
        possible_cyc = True
        for i, node in enumerate(cyc):
            next_node = node[1]
            if index + i < len(current_read_state):
                if not node[1] == current_read_state[index + i][1]:
                    possible_cyc = False
                    break
                previous_node = node[1]
            else:
                break
        if possible_cyc:
            this_len = thisstartpos - previous_end
            prev_len = DG[previous_node][next_node]["length"]
            len_difference = abs(this_len - prev_len)
            if len_difference < delta_len:
                found_next_node = next_node
                found = True
                break
    return found, found_next_node


def add_edge_support(edge_support, previous_node, name, r_id):
    edge_support[previous_node, name] = []
    edge_support[previous_node, name].append(r_id)


def known_old_node_action(alternatives_filtered, previous_node, this_len, nodes_for_graph, inter, k, seq, node_sequence,
                          r_id, DG, edge_support, alternative_nodes, old_node, alt_info_tuple, name):
    # if we have found a node which this info can be added to
    if alternatives_filtered:
        node_info = alternatives_filtered[0]
        # print("Ninfo",node_info)
        name = node_info[0]
        # update the read information of node name
        prev_nodelist = nodes_for_graph[name]
        r_infos = Read_infos(inter[0], inter[1], True)
        end_mini_seq = seq[inter[1]:inter[1] + k]
        node_sequence[name] = end_mini_seq
        prev_nodelist[r_id] = r_infos
        nodes_for_graph[name] = prev_nodelist
        if not DG.has_edge(previous_node, name):
            no_edge_found_action(DG, previous_node, name, this_len, inter, r_id, seq, node_sequence, nodes_for_graph,
                                 edge_support, k)
        else:
            edge_info = edge_support[previous_node, name]
            if r_id not in edge_info:
                edge_info.append(r_id)
                edge_support[previous_node, name] = edge_info
    # if we have not found a node which this info can be added to
    else:
        no_node_to_add_to_action(alternative_nodes, old_node, alt_info_tuple, inter, seq, k, node_sequence, name, r_id,
                                 nodes_for_graph, DG, this_len, previous_node, edge_support)


def no_edge_found_action(DG, previous_node, name, length, inter, r_id, seq, node_sequence, nodes_for_graph,
                         edge_support, k):
    DG.add_edge(previous_node, name, length=length)
    cycle_added2 = cycle_finder(DG, previous_node)
    if cycle_added2:
        print("cycle added between ",previous_node," and ",name)
        DG.remove_edge(previous_node, name)
        name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
        if not DG.has_node(name):
            DG.add_node(name)
        nodelist = {}
        r_infos = Read_infos(inter[0], inter[1], True)
        end_mini_seq = seq[inter[1]:inter[1] + k]
        node_sequence[name] = end_mini_seq
        nodelist[r_id] = r_infos
        # nodelist[r_id] = (inter[0], inter[1])
        # prev_nodelist[r_id] = r_infos
        nodes_for_graph[name] = nodelist
        DG.add_edge(previous_node, name, length=length)
        if DEBUG:
            print("adding edge no_edge_found_action" , previous_node, ",", name)
        add_edge_support(edge_support, previous_node, name, r_id)
    else:
        if DEBUG:
            print("no_edge_found_action no cycle added", previous_node, ",", name)
        add_edge_support(edge_support, previous_node, name, r_id)


def no_connecting_edge_action(seq, nodes_for_graph, name, inter, k, node_sequence, r_id, this_len, DG, previous_node,
                              edge_support):
    if DEBUG:
        print("no edge")
    # update the read information of node name
    prev_nodelist = nodes_for_graph[name]

    r_infos = Read_infos(inter[0], inter[1], True)
    end_mini_seq = seq[inter[1]:inter[1] + k]
    node_sequence[name] = end_mini_seq
    prev_nodelist[r_id] = r_infos
    nodes_for_graph[name] = prev_nodelist
    length = this_len
    DG.add_edge(previous_node, name, length=length)
    cycle_added2 = cycle_finder(DG, previous_node)
    if cycle_added2:
        DG.remove_edge(previous_node, name)
        if DEBUG:
            print(DG.edges)
        name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
        if not DG.has_node(name):
            DG.add_node(name)
        nodelist = {}
        r_infos = Read_infos(inter[0], inter[1], True)
        end_mini_seq = seq[inter[1]:inter[1] + k]
        node_sequence[name] = end_mini_seq
        nodelist[r_id] = r_infos
        nodes_for_graph[name] = nodelist
        DG.add_edge(previous_node, name, length=length)
        add_edge_support(edge_support, previous_node, name, r_id)
    else:
        add_edge_support(edge_support, previous_node, name, r_id)


def new_interval_action(seq, inter, r_id, node_sequence, nodes_for_graph, known_intervals, k, previous_end,
                        prior_read_infos, DG, previous_node, edge_support, name):
    #print("NIA")
    nodelist = {}
    # add a node into nodes_for_graph

    # add the read information for the node

    r_infos = Read_infos(inter[0], inter[1], True)
    end_mini_seq = seq[inter[1]:inter[1] + k]
    node_sequence[name] = end_mini_seq
    nodelist[r_id] = r_infos
    nodes_for_graph[name] = nodelist
    known_intervals[r_id - 1].append((inter[0], name, inter[1]))
    # get the length between the previous end and this nodes start
    length = inter[0] - previous_end
    add_prior_read_infos(inter, r_id, prior_read_infos, name, k)
    DG.add_node(name)
    # connect the node to the previous one
    DG.add_edge(previous_node, name, length=length)
    if DEBUG:
        print("adding edge clear 1035", previous_node, ",", name)
    add_edge_support(edge_support, previous_node, name, r_id)


def no_node_to_add_to_action(alternative_nodes, old_node, alt_info_tuple, inter, seq, k, node_sequence, name, r_id,
                             nodes_for_graph, DG, this_len, previous_node, edge_support):
    # add a new entry to alternative_nodes[old_node] to enable finding this instance
    alternative_nodes[old_node].append(alt_info_tuple)
    if DEBUG:
        print("oldnodeFurther", old_node)
    # add the read information for the node
    nodelist = {}
    r_infos = Read_infos(inter[0], inter[1], True)
    end_mini_seq = seq[inter[1]:inter[1] + k]
    node_sequence[name] = end_mini_seq
    nodelist[r_id] = r_infos
    nodes_for_graph[name] = nodelist
    DG.add_node(name)
    # get the length between the previous end and this nodes start
    length = this_len
    # connect the node to the previous one
    DG.add_edge(previous_node, name, length=length)
    add_edge_support(edge_support, previous_node, name, r_id)


def cycle_added_action(DG, previous_node, name, inter, r_id, seq, node_sequence, nodes_for_graph, k, length,
                       edge_support):
    DG.remove_edge(previous_node, name)
    name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
    if not DG.has_node(name):
        DG.add_node(name)
    nodelist = {}
    r_infos = Read_infos(inter[0], inter[1], True)
    end_mini_seq = seq[inter[1]:inter[1] + k]
    node_sequence[name] = end_mini_seq
    nodelist[r_id] = r_infos
    nodes_for_graph[name] = nodelist
    DG.add_edge(previous_node, name, length=length)
    add_edge_support(edge_support, previous_node, name, r_id)


def generateGraphfromIntervals(all_intervals_for_graph, k, delta_len, read_len_dict, all_reads):
    """ generates a networkx graph from the intervals given in all_intervals_for_graph.
    # INPUT:    all_intervals_for_graph:    A dictonary holding lists of minimizer intervals.
                k:                          K-mer length
                delta_len:                  integer holding the maximum lenght difference for two reads to still land in the same Isoform
                readlen_dict:               dictionary holding the read_id as key and the length of the read as value
    #OUTPUT:    result                      a tuple of different output values of the algo
                DG                          The Networkx Graph object
                known_intervals             A list of intervals used to check the correctness of the algo
                reads_for_isoforms          A dictionary holding the reads which are to be put in the same isoform
                reads_at_start_dict

    """

    DG = nx.DiGraph()
    # add the read ids to the startend_list
    reads_at_start_dict = {}
    reads_at_end_dict = {}
    reads_for_isoforms = []
    # holds the r_id as key and a list of tuples as value: For identification of reads, also used to ensure correctness of graph
    known_intervals = []
    # the following dictionary is supposed to hold the end minimizer sequence for each node
    node_sequence = {}
    edge_support = {}
    prior_read_infos = {}
    nodes_for_graph = {}
    alternative_nodes = {}
    alt_cyc_nodes={}
    for i in range(1, len(all_intervals_for_graph) + 1):
        reads_at_start_dict[i] = Read_infos(0, 0, True)
        reads_for_isoforms.append(i)
    for i in range(1, len(read_len_dict) + 1):
        reads_at_end_dict[i] = Read_infos(read_len_dict[i], read_len_dict[i], True)
    # a source and a sink node are added to the graph in order to have a well-defined start and end for the paths
    DG.add_node("s", reads=reads_at_start_dict, end_mini_seq='')
    # print("REads at end dict",reads_at_end_dict)
    DG.add_node("t", reads=reads_at_end_dict, end_mini_seq='')
    # adds an empty list for each r_id to known_intervals. To those lists, tuples, representing the intervals are added
    # for _ in itertools.repeat(None, len(all_intervals_for_graph)):
    for i in range(len(all_intervals_for_graph)):
        known_intervals.append([])
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    # intervals_for_read holds all intervals which make up the solution for the WIS of a read
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        if DEBUG:
            print('CURR READ:', r_id)
            print("NR Nodes", len(DG.nodes()))
            print("NR Edges", len(DG.edges()))
        containscycle = False
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"
        previous_end = 0
        # read_hashs is a dictionary storing hash values as keys and the respective node ids as values
        read_hashs = {}
        # the name of each node is defined to be readID, startminimizerpos , endminimizerpos
        # iterate over all intervals, which are in the solution of a certain read
        for pos, inter in enumerate(intervals_for_read):
            if DEBUG:
                print("PRNode", previous_node)
            info_tuple = (r_id, inter[0], inter[1])
            if DEBUG:
                print("INter", inter[0])
                print(info_tuple)
            # generate hash value of the intervals' infos
            curr_hash = convert_array_to_hash(inter[3])
            if curr_hash in read_hashs:
                is_repetative = True
                #cycle_start = read_hashs[curr_hash][-1]
            else:
                is_repetative = False
            # access prior_read_infos, if the same interval was already found in previous reads
            if info_tuple in prior_read_infos:
                if DEBUG:
                    print("in prior read infos")
                # if the interval repeats during this read
                if is_repetative:
                    if DEBUG:
                        print("cycle not close enough")
                    nodelist = {}
                    this_len = inter[0] - previous_end
                    # add a node into nodes_for_graph
                    name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                    seq = all_reads[r_id][1]
                    end_mini_seq = seq[inter[1]:inter[1] + k]
                    node_sequence[name] = end_mini_seq
                    r_infos = Read_infos(inter[0], inter[1], True)
                    nodelist[r_id] = r_infos
                    nodes_for_graph[name] = nodelist
                    DG.add_node(name)
                    # get the length between the previous end and this nodes start
                    length = this_len
                    # connect the node to the previous one
                    if DEBUG:
                        print("Adding edge, newnode no similar cycles connecting prevnode, name")
                    DG.add_edge(previous_node, name, length=length)
                    if DEBUG:
                        print("adding edge clear 553", previous_node, ",", name)
                    add_edge_support(edge_support, previous_node, name, r_id)
                # the interval did not repeat during this read
                else:
                    if DEBUG:
                        print("non repetative interval")
                    # get the name by accessing prior_read_infos
                    name = prior_read_infos[info_tuple]
                    #print("NAME",name)
                    # get length from previous_end to this_start
                    this_len = inter[0] - previous_end
                    # if there is not yet an edge connecting previous_node and name: add edge
                    if not DG.has_edge(previous_node, name):
                        if DEBUG:
                            print("no edge between ", previous_node, " and ", name)
                        seq = all_reads[r_id][1]
                        # no_connecting_edge_action(seq, nodes_for_graph, name, inter, k, node_sequence, r_id, this_len, DG, previous_node, edge_support)

                        # update the read information of node name
                        prev_nodelist = nodes_for_graph[name]
                        r_infos = Read_infos(inter[0], inter[1], True)
                        end_mini_seq = seq[inter[1]:inter[1] + k]
                        node_sequence[name] = end_mini_seq
                        prev_nodelist[r_id] = r_infos
                        nodes_for_graph[name] = prev_nodelist
                        length = this_len
                        #is_cyclic = SimplifyGraph.isCyclic(DG)
                        #print("SimplifyGraph_cyclic ", is_cyclic)
                        DG.add_edge(previous_node, name, length=length)
                        #print("edge ", previous_node, " to ", name, "added")
                        cycle_added2 = cycle_finder(DG, previous_node)
                        #is_cyclic = SimplifyGraph.isCyclic(DG)
                        #try:
                        #    nx_cycle=nx.find_cycle(DG)
                        #    print("Networkx found a cycle ",nx_cycle)
                        #except:
                        #    print("no cycle was found by networkx")
                        #if is_cyclic:
                            #    k_size+=1
                            #    w+=1
                            #print("The graph has a cycle - critical error")
                        #if is_cyclic != cycle_added2:
                        #    print("SimplifyGraph_cyclic ",is_cyclic)
                        #    print("GG_cyclic ", cycle_added2, ", ",previous_node)
                        #    print("ERROR- different outcome")
                        #the edge we added did introduce a cycle in our graph. We remove the edge again and instead generate a new node
                        #if cycle_added2:
                        if cycle_added2:
                            if DEBUG:
                                print("cycle added between ", previous_node, " and ", name)
                                print("CYC")
                                old_name = name
                                # alt_cyc_nodes is just used to keep track of how many nodes we add over the whole graph generation
                                if name not in alt_cyc_nodes.keys():
                                    alt_cyc_list = []
                                    alt_cyc_list.append(name)
                                    alt_cyc_nodes[old_name] = alt_cyc_list
                                else:
                                    alt_cyc_list = alt_cyc_nodes[old_name]
                                    alt_cyc_list.append(name)
                                    alt_cyc_nodes[old_name] = alt_cyc_list
                            DG.remove_edge(previous_node, name)
                            #TODO: add alt_cyc_node instead ofadding a new node to reduce overall nr nodes in our graph

                            name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                            if not DG.has_node(name):
                                DG.add_node(name)
                            nodelist = {}

                            r_infos = Read_infos(inter[0], inter[1], True)
                            end_mini_seq = seq[inter[1]:inter[1] + k]
                            node_sequence[name] = end_mini_seq
                            nodelist[r_id] = r_infos
                            nodes_for_graph[name] = nodelist
                            DG.add_edge(previous_node, name, length=length)
                            add_edge_support(edge_support, previous_node, name, r_id)
                        else:
                            add_edge_support(edge_support, previous_node, name, r_id)
                    # if there is an edge connecting previous_node and name: test if length difference is not higher than delta_len
                    else:
                        if DEBUG:
                            print("Has edge")
                            print(previous_node, name)
                        prev_len = DG[previous_node][name]["length"]
                        len_difference = abs(this_len - prev_len)
                        # if the length difference is <delta_len: Just add the readinfos
                        if len_difference < delta_len:
                            # update the read information of node name
                            prev_nodelist = nodes_for_graph[name]
                            seq = all_reads[r_id][1]
                            r_infos = Read_infos(inter[0], inter[1], True)
                            end_mini_seq = seq[inter[1]:inter[1] + k]
                            node_sequence[name] = end_mini_seq
                            prev_nodelist[r_id] = r_infos
                            nodes_for_graph[name] = prev_nodelist
                            edge_info = edge_support[previous_node, name]
                            if r_id not in edge_info:
                                edge_info.append(r_id)
                                edge_support[previous_node, name] = edge_info
                        # if the length difference is >delta len: generate new node, to be able to tell those nodes apart we use alternative_nodes
                        else:
                            #print("LEN diff > DL")
                            nodelist = {}
                            # inappropriate_node
                            old_node = name
                            name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                            alt_info_tuple = (name, previous_node, prev_len)
                            # add old node as key for alternative nodes: This is to find the alternative nodes for this one easily
                            if not (old_node in alternative_nodes):
                                alternative_nodes[old_node] = []
                                DG.add_node(name)
                                r_infos = Read_infos(inter[0], inter[1], True)
                                end_mini_seq = seq[inter[1]:inter[1] + k]
                                node_sequence[name] = end_mini_seq
                                nodelist[r_id] = r_infos
                                nodes_for_graph[name] = nodelist
                                DG.add_edge(previous_node, name, length=this_len)
                                add_edge_support(edge_support, previous_node, name, r_id)
                            # if we already know old_node (it is a key in alternative_nodes)
                            else:
                                alternative_infos_list = alternative_nodes[old_node]
                                alternatives_filtered = [item for item in alternative_infos_list if
                                                         previous_node == item[1] and abs(
                                                             this_len - item[2]) < delta_len]
                                known_old_node_action(alternatives_filtered, previous_node, this_len,
                                                      nodes_for_graph, inter, k,
                                                      seq, node_sequence, r_id, DG, edge_support, alternative_nodes,
                                                      old_node,
                                                      alt_info_tuple, name)

                # keep known_intervals up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
            # if the information for the interval was not yet found in a previous read (meaning the interval is new)
            else:
                seq = all_reads[r_id][1]
                name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                new_interval_action(seq, inter, r_id, node_sequence, nodes_for_graph, known_intervals, k, previous_end,
                                    prior_read_infos, DG, previous_node, edge_support, name)

            # set the previous node for the next iteration
            previous_node = name
            if DEBUG:
                print("PREVNODE", previous_node)
            previous_end = inter[1]
            # find out whether the current hash has already been added to read_hashs
            if not is_repetative:
                # if the current hash was not yet present in read_hashs: Add it
                read_hashs[curr_hash] = []
                read_hashs[curr_hash].append(name)
        # add an edge from name to "t" as name was the last node in the read
        if not DG.has_edge(name, "t"):
            DG.add_edge(name, "t", length=0)
            if DEBUG:
                print("adding edge clear 1066", name, ",", "t")
            edge_support[name, "t"] = []
            edge_support[name, "t"].append(r_id)
        else:
            edge_info = edge_support[name, "t"]
            if r_id not in edge_info:
                edge_info.append(r_id)
                edge_support[name, "t"] = edge_info
    print("Number of alternative nodes due to cycle",len(alt_cyc_nodes.keys()))
    # set the node attributes to be nodes_for_graph, very convenient way of solving this
    nx.set_node_attributes(DG, nodes_for_graph, name="reads")
    nx.set_node_attributes(DG, node_sequence, name="end_mini_seq")
    nx.set_edge_attributes(DG, edge_support, name='edge_supp')

    # return DG and known_intervals, used to
    result = (DG, reads_for_isoforms)
    return result


"""
Plots the given Directed Graph via Matplotlib
Input: DG   Directed Graph
"""





def get_read_lengths(all_reads):
    """Helper method which extracts the read lengths from all_reads. We will use those during the graph generation to appoint more meaningful information to the node't'
    INPUT: all_reads: dictionary which holds the overall read infos key: r_id, value tuple(readname, sequence, some info i currently don't care about)
    OUTPUT: readlen_dict: dictionary holding the read_id as key and the length of the read as value
    """
    readlen_dict = {}
    for r_id, infos in all_reads.items():
        seq = infos[1]
        seqlen = len(seq)

        readlen_dict[r_id] = seqlen
    return readlen_dict


"""
USED FOR DEBUGGING ONLY-deprecated in IsONform
"""

DEBUG=False
def main():
    # import sys
    # sys.stdout = open('log.txt', 'w')
    print('test')
    w = 10
    reads = 62
    max_seqs_to_spoa = 200
    work_dir = tempfile.mkdtemp()
    # print("Temporary workdirektory:", work_dir)
    k_size = 9
    outfolder = "out_local"
    cwd = os.getcwd()
    # print("CWD",cwd)
    file = open('all_intervals.txt', 'rb')
    all_intervals_for_graph = pickle.load(file)
    file.close()
    # print("All of them", len(all_intervals_for_graph))

    file2 = open('all_reads.txt', 'rb')

    all_reads = pickle.load(file2)
    print("Allreads type")
    # print(type(all_reads))
    # print(all_reads)
    file2.close()
    delta_len = 2 * k_size
    # is_cyclic = True
    # while is_cyclic:
    read_len_dict = get_read_lengths(all_reads)
    print("Generating graph")
    DG, known_intervals, reads_for_isoforms, reads_list = generateGraphfromIntervals(
        all_intervals_for_graph, k_size, delta_len, read_len_dict, all_reads)

    print("edges with attributes:")

    print("Calling the method")

    # simplifyGraph(DG, all_reads,work_dir,k_size,delta_len,known_intervals)
    print("Graph simplified")
    # simplifyGraphOriginal(DG,delta_len, all_reads, work_dir, k_size,known_intervals)
    # print("#Nodes for DG: " + str(DG.number_of_nodes()) + " , #Edges for DG: " + str(DG.number_of_edges()))
    # draw_Graph(DG)
    # print(DG.nodes["s"]["reads"])
    # print(list(DG.nodes(data=True)))
    # print("edges with attributes:")
    # print(DG.edges(data=True))
    # print("ReadNodes")
    # print("all edges for the graph")
    # print([e for e in DG.edges])
    # draw_Graph(DG)
    # print("knownintervals")
    # print(known_intervals)
    # The call for the isoform generation (deprecated during implementation)
    # TODO: invoke isoform_generation again
    # generate_isoforms(DG, all_reads, reads_for_isoforms, work_dir, outfolder, max_seqs_to_spoa)
    # DG.nodes(data=True)

    # print("removing temporary workdir")
    shutil.rmtree(work_dir)


if __name__ == "__main__":
    main()
