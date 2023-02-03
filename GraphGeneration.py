
import networkx as nx
import matplotlib.pyplot as plt
from collections import namedtuple
from SimplifyGraph import *
import tempfile
import shutil



"""IsONform script containing the methods used to generate the Directed Graph from the Intervals coming from the Weighted Interval Scheduling Problem
Author: Alexander Petri
The main method in this script was used for debugging therefore is not used during IsONforms actual run.
"""

"""
Function to convert a list into a string to enable the writing of a graph (Taken from https://www.geeksforgeeks.org/python-program-to-convert-a-list-to-string/)
INPUT: s:   list of elements to be converted into a string
OUTPUT: the string object derived from the list
"""


def listToString(s):
    # initialize an empty string
    str1 = " "

    # return string
    return (str1.join(str(s)))


"""
Function to add read information of the current interval to prior_read_infos, if they do not exist there yet
INPUT:      inter:              the current interval holding the information
            r_id:               the read we are looking at
            prior_read_infos:   information about other reads
            name:               name of the node we appoint the information to
            k:                  minimizer length
OUTPUT:     prior_read_infos:   information about other reads extended by this intervals infos
"""


def add_prior_read_infos(inter, r_id, prior_read_infos, name, k):
    read_id = inter[3][slice(0, len(inter[3]),
                             3)]  # recover the read id from the array of instances which was delivered with all_intervals_for_graph
    start_coord = inter[3][slice(1, len(inter[3]),
                                 3)]  # recover the start coordinate of an interval from the array of instances
    end_coord = inter[3][slice(2, len(inter[3]), 3)]
    # iterate through the retreived information and store it in prior_read_infos only for subsequent reads
    for i, r in enumerate(read_id):
        if not r <= r_id:
            #print(r)
            start = start_coord[i] + k
            end = end_coord[i]
            tuple_for_data_structure = (r, start, end)
            #print(r_id,name,inter)
            #print("TUPLE For DATA STRUCTURE", tuple_for_data_structure)
            if not tuple_for_data_structure in prior_read_infos:
                prior_read_infos[tuple_for_data_structure] = name

    return prior_read_infos


"""Helper method to convert the array delivered with all_intervals_for_graph into a hash value to more efficiently look up occurence
This method additionally deletes the first three entries of the array as they contain the infos about this interval occurence, which changes in between instances
INPUT:  info_array: The array, which was delivered with the interval to indicate where the interval occurs in other reads
OUTPUT: hash(tup):  A has value of the tuple the shortened interval was converted into. This hash makes it easy to see whether the interval is already present in the read
"""


def convert_array_to_hash(info_array):
    # preprocessing: Delete the first three elements from the array, as they contain the information about this occurrence
    for x in range(0, 3):
        if info_array:
            # print("IA",info_array)
            info_array.pop(0)

    tup = tuple(info_array)
    return (hash(tup))


"""Function to get all nodes which are part of a cycle
    INPUT:  current_read_state: The current state of known_intervals[r_id-1] (known_intervals of the current read)
            cycle_start:        The last occurence of the repeating node before the cycle starts
    OUTPUT: cycle_nodes:        A list of entries indicating which nodes are in the cycle, having the following form: (startpos, node_id, endpos)
"""


def record_cycle(current_read_state, cycle_start):
    cycle_nodes = []
    indices = [i for i, tupl in enumerate(current_read_state) if tupl[1] == cycle_start]
    index = indices[0]
    for element in range(index, len(current_read_state)):
        cycle_nodes.append(current_read_state[element])

    return (cycle_nodes, current_read_state[index])


"""Helper method for generateGraphfromIntervals
"""


def find_next_node(thisstartpos, info_array, known_cycles, current_read_state, k, cycle_start, previous_end, DG,
                   delta_len):
    # print(current_read_state)
    feasible_cycles = []
    intermediate_cycles = []
    indices = [i for i, tupl in enumerate(current_read_state) if tupl[1] == cycle_start]
    index = indices[0]
    # print("Known_cycles")
    # print(known_cycles)
    for key, value in known_cycles.items():
        known_cycle_rid = key[0]
        node_appearance = value[0]
        for i in range(0, len(info_array) - 2):
            if info_array[i] == known_cycle_rid and info_array[i + 1] == node_appearance[0] - k and info_array[i + 2] == \
                    node_appearance[2]:
                intermediate_cycles.append(value)
    # print("Intermediate")
    # print(intermediate_cycles)
    found = False
    found_next_node = None
    for cyc in intermediate_cycles:
        possible_cyc = True
        for i, node in enumerate(cyc):
            next_node = node[1]

            if index + i < len(current_read_state):
                # print(current_read_state[index+i][1])
                if not node[1] == current_read_state[index + i][1]:
                    possible_cyc = False
                    break
                previous_node = node[1]
            else:
                break
        if (possible_cyc):
            # TODO figure out whether delta len-(thisstart-previousend)<3. If yes: add this occurence to the cycles' last node,if not:just continue the loop in the hope there are more intermediate cycles
            # print("Possible:")
            # print(cyc)
            this_len = thisstartpos - previous_end
            # print("Previous node"+previous_node)
            # print("Nextnode: "+next_node)
            prev_len = DG[previous_node][next_node]["length"]
            len_difference = abs(this_len - prev_len)
            # print("Len_difference:"+str(len_difference))
            if len_difference < delta_len:
                found_next_node = next_node
                found = True
                break
    # if found:
    # print("Found next node: "+found_next_node)
    return (found, found_next_node)

def find_topo_alternative(topo_graph,previous_node,name,topo_alternatives):
    if name in topo_alternatives:
        #print(topo_alternatives)
        #print(topo_alternatives[name])
        for alternative in topo_alternatives[name]:
            #print(alternative)
            if topo_graph.index(previous_node) < topo_graph.index(alternative):
                return alternative
    return ""
########################################
######## KRISTOFFER START ##############
########################################
def subsequence(seq):
    if not seq:
        return seq
    M = [None] * len(seq)    # offset by 1 (j -> j-1)
    P = [None] * len(seq)

    # Since we have at least one element in our list, we can start by
    # knowing that the there's at least an increasing subsequence of length one:
    # the first element.
    L = 1
    M[0] = 0

    # Looping over the sequence starting from the second element
    for i in range(1, len(seq)):
        # Binary search: we want the largest j <= L
        #  such that seq[M[j]] < seq[i] (default j = 0),
        #  hence we want the lower bound at the end of the search process.
        lower = 0
        upper = L

        # Since the binary search will not look at the upper bound value,
        # we'll have to check that manually
        if seq[M[upper-1]][0] < seq[i][0]:
            j = upper

        else:
            # actual binary search loop
            while upper - lower > 1:
                mid = (upper + lower) // 2
                if seq[M[mid-1]][0] < seq[i][0]:
                    lower = mid
                else:
                    upper = mid

            j = lower    # this will also set the default value to 0

        P[i] = M[j-1]

        if j == L or seq[i][0] < seq[M[j]][0]:
            M[j] = i
            L = max(L, j+1)

    # Building the result: [seq[M[L-1]], seq[P[M[L-1]]], seq[P[P[M[L-1]]]], ...]
    result = []
    pos = M[L-1]
    for _ in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return result[::-1]    # reversing
########################################
######## KRISTOFFER END   ##############
########################################

def find_LIS(DG, r_id, intervals_for_read,k,delta_len):
    top_order = {name: order for order, name in enumerate(nx.topological_sort(DG))}
    anchors = []
    for read_order, (start, stop, weight, supp_int) in enumerate(intervals_for_read):
        _, curr_start, curr_stop = supp_int[:3]
        curr_node_name = "{0}, {1}, {2}".format(curr_start + k, curr_stop, r_id)
        curr_node_distance=curr_stop-(curr_start+k)
        for i in range(len(supp_int) // 3):
            supp_r_id, supp_start, supp_stop = supp_int[3 * i:3 * i + 3]
            if supp_r_id == r_id:
                continue
            supp_node_name = "{0}, {1}, {2}".format(supp_start + k, supp_stop, supp_r_id)
            supp_node_distance=supp_stop-(supp_start+k)
            # print(supp_node_name)
            if supp_node_name in DG:
                if abs(supp_node_distance-curr_node_distance)>delta_len:
                    anchors.append((top_order[supp_node_name], supp_node_name,
                                curr_node_name))  # (order in graph, graph node name, read_node_name)
                # print('Found!')

    #print(anchors, '<-- All anchors')
    ans = subsequence(anchors)
    #print(ans, '<-- Longest increasing subsequence of anchors')

    # or a bit nicer for the eyes
    tmp_1 = [a[0] for a in anchors]  # only the orders
    tmp_2 = [a[0] for a in ans]  # only the orders
    #print(tmp_1)
    #print(tmp_2)
    #print("This above list is the final longer increasing chain of nodes in the read that was found in the graph.")
    #print("The nodes missing from this will have to be added in by forming new nodes or finding them between the LIS (because I use DFG generated topological sort to get the orders")
    # Alex: this is just to check where the LIS was different, i.e., complication _could_ occur. this f statement should be removed/ or sys.exit uncommented
    #if anchors != ans:
        #print('FOUND INTERESTING CASE, stopping')
        #sys.exit()
    return anchors
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
TODO:   add dictionary to store infos about known instances 
        revisit structure of the method and introduce subroutines
"""


#TODO: invoke list_solution to actually merge nodes.
def generateGraphfromIntervals(all_intervals_for_graph, k, delta_len, read_len_dict, all_reads):
    Read_infos = namedtuple('Read_Infos',
                            'start_mini_end end_mini_start original_support')
    DEBUG = False
    DG = nx.DiGraph()
    cycle_nodes = {}
    # add the read ids to the startend_list
    reads_at_start_dict = {}
    reads_at_end_dict = {}
    reads_for_isoforms = []
    for i in range(1, len(all_intervals_for_graph) + 1):
        reads_at_start_dict[i] = Read_infos(0, 0, True)
        reads_for_isoforms.append(i)
    for i in range(1, len(read_len_dict) + 1):
        reads_at_end_dict[i] = Read_infos(read_len_dict[i], read_len_dict[i], True)
        # reads_for_isoforms.append(i)
    # print(reads_at_start_dict)
    # print(reads_at_end_dict)
    # a source and a sink node are added to the graph in order to have a well-defined start and end for the paths
    DG.add_node("s", reads=reads_at_start_dict, end_mini_seq='')
    # print("REads at end dict",reads_at_end_dict)
    DG.add_node("t", reads=reads_at_end_dict, end_mini_seq='')

    # holds the r_id as key and a list of tuples as value: For identification of reads, also used to ensure correctness of graph
    known_intervals = []
    node_overview_read = []
    # adds an empty list for each r_id to known_intervals. To those lists, tuples, representing the intervals are added
    # for _ in itertools.repeat(None, len(all_intervals_for_graph)):
    for i in range(len(all_intervals_for_graph)):
        # known_intervals.append(OrderedSet([]))
        known_intervals.append([])
        node_overview_read.append([])
    # print(known_intervals)
    nodes_for_graph = {}
    alternative_nodes = {}
    topo_alternatives={}
    topo_graph = None
    # the following dictionary is supposed to hold the end minimizer sequence for each node
    node_sequence = {}
    edge_support = {}
    prior_read_infos = {}
    cycle_in_reads = {}
    known_cycles = {}
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    # intervals_for_read holds all intervals which make up the solution for the WIS of a read
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        #if r_id==212 or r_id==251 or r_id==212:# or r_id==198:
        #    DEBUG=True
        #else:
        #    DEBUG=False
        ########################################
        ######## KRISTOFFER START ##############
        ########################################
        #print("READ",r_id)
        #print(r_id, intervals_for_read)
        if DEBUG:
            print('CURR READ:', r_id)
        #print(DG.nodes())  # print(DG.nodes(data=True))
        #print()
        #print(DG.edges(data=True))
        #print(list(nx.topological_sort(DG)))
        #if len(DG) > 2:  # we have added at least one read to the graph already
        #    lis_solution = find_LIS(DG, r_id, intervals_for_read,k,delta_len)
        #else:
        #    lis_solution=[]
        ########################################
        ######## KRISTOFFER END   ##############
        ########################################
        containscycle = False
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"
        previous_end = 0
        # read_hashs is a dictionary storing hash values as keys and the respective node ids as values
        read_hashs = {}

        # the name of each node is defined to be readID, startminimizerpos , endminimizerpos
        # iterate over all intervals, which are in the solution of a certain read
        for pos,inter in enumerate(intervals_for_read):
            #if isCyclic(DG):
            #    print("ERROR-CYCLIC")
            #    print(DG.edges())
            #    sys.exit()
            if DEBUG:
                print("PRNode",previous_node)
            prev_nodelist = {}
            info_tuple = (r_id, inter[0], inter[1])
            if DEBUG:
                print("INter", inter[0])
                print(info_tuple)
            # generate hash value of the intervals' infos
            curr_hash = convert_array_to_hash(inter[3])
            if curr_hash in read_hashs:
                is_repetative = True
                cycle_start = read_hashs[curr_hash][-1]
            else:
                is_repetative = False
            # access prior_read_infos, if the same interval was already found in previous reads
            #lis_sol_occ=[x for x in lis_solution if x[0]==pos]
            #print("LSO",lis_sol_occ)
            if info_tuple in prior_read_infos:
                if DEBUG:
                    #for id,prior_info in prior_read_infos.items():
                        #print(id,prior_info)
                    print("in prior read infos")
                # if the interval repeats during this read
                if is_repetative:
                    if DEBUG:
                        # print(info_tuple)
                        print("Repetative")
                        # print(known_cycles)
                        # print(inter[3])

                    # use find_next_node to access all previously recorded cycles and find the one matching with the current cycle
                    foundnode, merge_address = find_next_node(inter[0], inter[3], known_cycles,
                                                              known_intervals[r_id - 1], k, cycle_start, previous_end,
                                                              DG, delta_len)
                    # if we have found a next node
                    if foundnode:
                        if DEBUG:
                            print("Found node")
                        name = merge_address
                        # only add a new edge if the edge was not present before
                        if not DG.has_edge(previous_node, name):
                            if DEBUG:
                                print("No edge")
                            #if we would harm the topological order in our graph by adding the edge
                            #if True:
                                """if topo_graph.index(previous_node) > topo_graph.index(name):
                                old_name=name
                                #find alternatives for the node in topo_alternatives
                                topo_alt_result=find_topo_alternative(topo_graph,previous_node,name,topo_alternatives)
                                #if we have not found an alternative in topo_alternatives:
                                if not topo_alt_result:
                                    length = this_len
                                    #add new node that either is named 'name' or 'name1..1'
                                    name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                                    if not DG.has_node(name):
                                        DG.add_node(name)
                                    else:
                                        if DEBUG:
                                            print("Node", name, " already present")
                                        while DG.has_node(name):
                                            name = name + 1
                                        if DEBUG:
                                            print("Adding node", name)
                                        DG.add_node(name)
                                    nodelist = {}
                                    r_infos = Read_infos(inter[0], inter[1], True)
                                    end_mini_seq = seq[inter[1]:inter[1] + k]
                                    node_sequence[name] = end_mini_seq
                                    nodelist[r_id] = r_infos
                                    # nodelist[r_id] = (inter[0], inter[1])
                                    #prev_nodelist[r_id] = r_infos
                                    nodes_for_graph[name] = nodelist
                                    if not old_name in topo_alternatives:
                                        # we store the new node as alternative for later possible topological misshaps
                                        alt_list = []
                                        alt_list.append(name)
                                        topo_alternatives[old_name] =alt_list
                                    else:
                                        topo_alternatives[old_name].append(name)
                                    if DEBUG:
                                        print("Adding edge, foundnode, no edge connecting prevnode, name, no topo_alt_result")
                                    if DG.has_edge(previous_node,name):
                                        edge_info = edge_support[previous_node, name]
                                        if not r_id in edge_info:
                                            edge_info.append(r_id)
                                            edge_support[previous_node, name] = edge_info
                                    else:
                                        DG.add_edge(previous_node, name, length=length)
                                        edge_support[previous_node, name] = []
                                        edge_support[previous_node, name].append(r_id)
                                    if r_id > 1:
                                        topo_graph = list(nx.topological_sort(DG))
                                else:
                                    length = this_len
                                    name=topo_alt_result
                                    prev_nodelist = nodes_for_graph[name]
                                    end_mini_seq = seq[inter[1]:inter[1] + k]
                                    r_infos = Read_infos(inter[0], inter[1], True)
                                    node_sequence[name] = end_mini_seq
                                    prev_nodelist[r_id] = r_infos
                                    nodes_for_graph[name] = prev_nodelist
                                    # if DEBUG:
                                    # print("Adding infosd to node", name, " : ", prev_nodelist)
                                    # if name == '601, 652, 20':
                                    #    print("Adding infosf to the node", name, " : ", prev_nodelist)
                                    if not DG.has_edge(previous_node, name):
                                        if DEBUG:
                                            print("Adding edge, foundnode, no edge connecting prevnode, name")
                                        DG.add_edge(previous_node, name, length=length)
                                        edge_support[previous_node, name] = []
                                        edge_support[previous_node, name].append(r_id)
                                    else:
                                        edge_info = edge_support[previous_node, name]
                                        if not r_id in edge_info:
                                            edge_info.append(r_id)
                                            edge_support[previous_node, name] = edge_info
                            else:"""
                                # update the read information of node name
                            prev_nodelist = nodes_for_graph[name]
                            seq = all_reads[r_id][1]
                            r_infos = Read_infos(inter[0], inter[1], True)
                            end_mini_seq = seq[inter[1]:inter[1] + k]
                            if not (end_mini_seq == DG.nodes[name]['end_mini_seq']):
                                print("ERROR: ", end_mini_seq, " not equal to ", DG.nodes[name]['end_mini_seq'])
                            # node_sequence[name]=end_mini_seq
                            prev_nodelist[r_id] = r_infos
                            nodes_for_graph[name] = prev_nodelist
                            # if DEBUG:
                            # print("Adding infos to node",name," : ",prev_nodelist)
                            # if name == '601, 652, 20':
                            #    print("Adding infos to node", name, " : ", prev_nodelist)
                            length = this_len
                            if DEBUG:
                                print("Adding edge, foundnode, topo order not violated  for prevnode, name")
                            DG.add_edge(previous_node, name, length=length)
                            cycle_added = isCyclic(DG)
                            cycle_added2 = cycle_finder(DG,previous_node)
                            cycle_added3 = cycle_finder(DG,name)
                            #if cycle_added != cycle_added2 and cycle_added != cycle_added3:
                            #	print(cycle_added,cycle_added2,cycle_added3)
                            #	sys.exit(1)
                            #cycle_added =isCyclic2(DG,name)
                            #try:
                            #    topo_sort=nx.topological_sort(DG)
                            #    #print(topo_sort[0])
                            #except:
                            if cycle_added2:
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
                                print("adding edge clear509", previous_node, ",", name)
                                edge_support[previous_node, name] = []
                                edge_support[previous_node, name].append(r_id)
                            else:
                                if DEBUG:
                                    print("adding edge hope 501", previous_node, ",", name)
                                edge_support[previous_node, name] = []
                                edge_support[previous_node, name].append(r_id)
                        #we found an edge from previous_node to name. We only have to add the read(-interval) information
                        else:
                            if DEBUG:
                                print("Found edge")
                            # update the read information of node name
                            prev_nodelist = nodes_for_graph[name]
                            seq = all_reads[r_id][1]
                            r_infos = Read_infos(inter[0], inter[1], True)
                            # end_mini_seq = seq[inter[1]:inter[1] + k]

                            # node_sequence[name] = end_mini_seq
                            prev_nodelist[r_id] = r_infos
                            nodes_for_graph[name] = prev_nodelist
                            # if DEBUG:
                            #    print("Adding infos to node",name," : ",prev_nodelist)
                            # if name == '601, 652, 20':
                            #    print("Adding infos to node", name, " : ", prev_nodelist)
                            edge_info = edge_support[previous_node, name]
                            if not r_id in edge_info:
                                edge_info.append(r_id)
                                edge_support[previous_node, name] = edge_info

                    # TODO:Verify that this is sufficient to generate a new node for DG
                    # if the cycles were not close enough to this read: Generate new node and save the cycle( after building up nodes)
                    else:
                        if DEBUG:
                            print("cycle not close enough")
                        nodelist = {}
                        this_len = inter[0] - previous_end
                        # add a node into nodes_for_graph
                        name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                        if name == "166, 197, 9":
                            print("Node ", name)
                        end_mini_seq = seq[inter[1]:inter[1] + k]
                        node_sequence[name] = end_mini_seq
                        r_infos = Read_infos(inter[0], inter[1], True)
                        nodelist[r_id] = r_infos
                        nodes_for_graph[name] = nodelist
                        # if DEBUG:
                        # print("Adding infos to node", name, " : ", nodelist)
                        # if name == '601, 652, 20':
                        #    print("Adding infos to node", name, " : ", nodelist)
                        DG.add_node(name)
                        # get the length between the previous end and this nodes start
                        length = this_len
                        # connect the node to the previous one
                        if DEBUG:
                            print("Adding edge, newnode no similar cycles connecting prevnode, name")
                        DG.add_edge(previous_node, name, length=length)
                        if DEBUG:
                            print("adding edge clear 553", previous_node, ",", name)
                        edge_support[previous_node, name] = []
                        edge_support[previous_node, name].append(r_id)
                        #if r_id > 1:
                            #topo_graph = list(nx.topological_sort(DG))

                # the interval did not repeat during this read
                else:
                    if DEBUG:
                        print("non repetative interval")
                    # get the name by accessing prior_read_infos
                    name = prior_read_infos[info_tuple]
                    # get length from previous_end to this_start
                    this_len = inter[0] - previous_end

                    # len_difference = abs(this_len - prev_len)
                    # if this_len < 0:
                    #    print("Error- length negative.")
                    # if there is not yet an edge connecting previous_node and name: add edge
                    if not DG.has_edge(previous_node, name):
                        if DEBUG:
                            print("no edge")
                        if True:
                            """if topo_graph.index(previous_node) > topo_graph.index(name):
                            old_name = name
                            # find alternatives for the node in topo_alternatives
                            topo_alt_result = find_topo_alternative(topo_graph, previous_node, name, topo_alternatives)
                            # if we have not found an alternative in topo_alternatives:
                            if not topo_alt_result:
                                length = this_len
                                # add new node that either is named 'name' or 'name1..1'
                                name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                                if not DG.has_node(name):
                                    DG.add_node(name)
                                else:
                                    if DEBUG:
                                        print("Node", name, " already present")
                                    while DG.has_node(name):
                                        name = name + 1
                                    if DEBUG:
                                        print("Adding node", name)
                                    DG.add_node(name)
                                nodelist = {}
                                r_infos = Read_infos(inter[0], inter[1], True)
                                end_mini_seq = seq[inter[1]:inter[1] + k]
                                node_sequence[name] = end_mini_seq
                                nodelist[r_id] = r_infos
                                # nodelist[r_id] = (inter[0], inter[1])
                                #prev_nodelist[r_id] = r_infos
                                nodes_for_graph[name] = nodelist
                                if not old_name in topo_alternatives:
                                    # we store the new node as alternative for later possible topological misshaps
                                    alt_list = []
                                    alt_list.append(name)
                                    topo_alternatives[old_name] = alt_list
                                else:
                                    topo_alternatives[old_name].append(name)
                                if DEBUG:
                                    print("Adding edge, no repetition prevnode, name")
                                if DG.has_edge(previous_node, name):
                                    edge_info = edge_support[previous_node, name]
                                    if not r_id in edge_info:
                                        edge_info.append(r_id)
                                        edge_support[previous_node, name] = edge_info
                                else:
                                    DG.add_edge(previous_node, name, length=length)
                                    edge_support[previous_node, name] = []
                                    edge_support[previous_node, name].append(r_id)
                                if r_id > 1:
                                    topo_graph = list(nx.topological_sort(DG))
                            else:
                                length = this_len
                                name = topo_alt_result
                                prev_nodelist = nodes_for_graph[name]
                                end_mini_seq = seq[inter[1]:inter[1] + k]
                                r_infos = Read_infos(inter[0], inter[1], True)
                                node_sequence[name] = end_mini_seq
                                prev_nodelist[r_id] = r_infos
                                nodes_for_graph[name] = prev_nodelist
                                # if DEBUG:
                                # print("Adding infosd to node", name, " : ", prev_nodelist)
                                # if name == '601, 652, 20':
                                #    print("Adding infosf to the node", name, " : ", prev_nodelist)
                                if not DG.has_edge(previous_node, name):
                                    if DEBUG:
                                        print("Adding edge, no repetition , no topo_alt result prevnode, name")
                                    DG.add_edge(previous_node, name, length=length)
                                    edge_support[previous_node, name] = []
                                    edge_support[previous_node, name].append(r_id)
                                else:
                                    edge_info = edge_support[previous_node, name]
                                    if not r_id in edge_info:
                                        edge_info.append(r_id)
                                        edge_support[previous_node, name] = edge_info
                        else:"""
                            # update the read information of node name
                            prev_nodelist = nodes_for_graph[name]

                            seq = all_reads[r_id][1]
                            r_infos = Read_infos(inter[0], inter[1], True)
                            end_mini_seq = seq[inter[1]:inter[1] + k]
                            node_sequence[name] = end_mini_seq
                            prev_nodelist[r_id] = r_infos
                            nodes_for_graph[name] = prev_nodelist
                            # if DEBUG:
                            # print("Adding infos to node", name, " : ", prev_nodelist)
                            # if name == '601, 652, 20':
                            #    print("Adding infos to node", name, " : ", prev_nodelist)
                            # only add a new edge if the edge was not present before
                            length = this_len
                            is_false=False
                            if True:
                                #print("Adding edge, not similar to cycles, topo not violated prevnode, name")
                                DG.add_edge(previous_node, name, length=length)
                                #cycle_added=isCyclic(DG)
                                cycle_added2 = cycle_finder(DG,previous_node)
                                #cycle_added3 = cycle_finder(DG,name)
                                #if cycle_added != cycle_added2 and cycle_added != cycle_added3:
                                #    print(cycle_added,cycle_added2,cycle_added3)
                                #    sys.exit(1)
                                #cycle_added = isCyclic2(DG, name)
                                #if DEBUG:
                                #    print(previous_node,name)
                                    #print(DG.edges)
                                #try:
                                #    print("ADDINGCYCLIC?",)
                                #    topo_sort = nx.topological_sort(DG)
                                #    #print(topo_sort[0])
                                if cycle_added2:
                                    #print("NO TOPO Found")
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
                                    # nodelist[r_id] = (inter[0], inter[1])
                                    # prev_nodelist[r_id] = r_infos
                                    nodes_for_graph[name] = nodelist
                                    #print(nodes_for_graph[name])
                                    #print(previous_node)
                                    DG.add_edge(previous_node, name, length=length)
                                    if DEBUG:
                                        print("adding edge clear 686", previous_node, ",", name)
                                    edge_support[previous_node, name] = []
                                    edge_support[previous_node, name].append(r_id)
                                    #print(edge_support[previous_node, name])
                                    #print(DG.edges)
                                    #print("ISFALSE",is_false)
                                else:
                                    #print("ELSE")
                                    if DEBUG:
                                        print("adding edge hope690", previous_node, ",", name)
                                    edge_support[previous_node, name] = []
                                    edge_support[previous_node, name].append(r_id)
                                #print("ISFALSE",is_false)
                            #DG.add_edge(previous_node, name, length=length)
                            #edge_support[previous_node, name] = []
                            #edge_support[previous_node, name].append(r_id)
                    # TODO: add topological sort and everything after this point
                    # if there is an edge connecting previous_node and name: test if length difference is not higher than delta_len
                    else:
                        if DEBUG:
                            print("Has edge")
                            print(previous_node, name)
                        prev_len = DG[previous_node][name]["length"]
                        # print("Prev_len:"+str(prev_len))
                        len_difference = abs(this_len - prev_len)
                        # print("Len_difference:"+str(len_difference))
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
                            # if DEBUG:
                            # print("Adding infost to mnodes", name, " : ", prev_nodelist)
                            # if name == '601, 652, 20':
                            #    print("Adding infos to the node", name, " : ", prev_nodelist)

                            edge_info = edge_support[previous_node, name]
                            if not r_id in edge_info:
                                edge_info.append(r_id)
                                edge_support[previous_node, name] = edge_info
                        # if the length difference is >delta len: generate new node, to be able to tell those nodes apart we use alternative_nodes
                        else:
                            nodelist = {}
                            # inappropriate_node
                            old_node = name
                            # print("old_node",old_node)
                            # add a node into nodes_for_graph
                            # if DEBUG:
                            #    print("NFG", nodes_for_graph['601, 652, 20'])
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
                                # nodelist[r_id] = (inter[0], inter[1])
                                #prev_nodelist[r_id] = r_infos
                                nodes_for_graph[name] = nodelist
                                if DEBUG:
                                    print("delta len < actual diff, adding new node ")
                                DG.add_edge(previous_node, name, length=this_len)
                                if DEBUG:
                                    print("adding edge clear 750", previous_node, ",", name)
                                edge_support[previous_node, name] = []
                                edge_support[previous_node, name].append(r_id)
                                #if r_id > 1:
                                    #topo_graph = list(nx.topological_sort(DG))
                            # if we already know old_node (it is a key in alternative_nodes)
                            else:
                                alternative_infos_list = alternative_nodes[old_node]
                                alternatives_filtered = [item for item in alternative_infos_list if
                                                         previous_node == item[1] and abs(
                                                             this_len - item[2]) < delta_len]
                                # if we have found a node which this info can be added to
                                if alternatives_filtered:
                                    node_info = alternatives_filtered[0]
                                    # print("Ninfo",node_info)
                                    name = node_info[0]
                                    if True:
                                        """if topo_graph.index(previous_node) > topo_graph.index(name):
                                        old_name = name
                                        # find alternatives for the node in topo_alternatives
                                        topo_alt_result = find_topo_alternative(topo_graph, previous_node, name,
                                                                                topo_alternatives)
                                        # if we have not found an alternative in topo_alternatives:
                                        if not topo_alt_result:
                                            length = this_len
                                            # add new node that either is named 'name' or 'name1..1'
                                            name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                                            if not DG.has_node(name):
                                                DG.add_node(name)
                                            else:
                                                if DEBUG:
                                                    print("Node", name, " already present")
                                                while DG.has_node(name):
                                                    name = name + 1
                                                if DEBUG:
                                                    print("Adding node", name)
                                                DG.add_node(name)
                                            nodelist = {}
                                            r_infos = Read_infos(inter[0], inter[1], True)
                                            end_mini_seq = seq[inter[1]:inter[1] + k]
                                            node_sequence[name] = end_mini_seq
                                            nodelist[r_id] = r_infos
                                            # nodelist[r_id] = (inter[0], inter[1
                                            nodes_for_graph[name] = nodelist
                                            if not old_name in topo_alternatives:
                                                # we store the new node as alternative for later possible topological misshaps
                                                alt_list = []
                                                alt_list.append(name)
                                                topo_alternatives[old_name] = alt_list
                                            else:
                                                topo_alternatives[old_name].append(name)
                                            if DEBUG:
                                                print("oldnode known, found alternative node")
                                            if DG.has_edge(previous_node, name):
                                                edge_info = edge_support[previous_node, name]
                                                if not r_id in edge_info:
                                                    edge_info.append(r_id)
                                                    edge_support[previous_node, name] = edge_info
                                            else:
                                                DG.add_edge(previous_node, name, length=length)
                                                edge_support[previous_node, name] = []
                                                edge_support[previous_node, name].append(r_id)
                                            #if r_id > 1:
                                                #topo_graph = list(nx.topological_sort(DG))
                                        else:
                                            length = this_len
                                            name = topo_alt_result
                                            prev_nodelist = nodes_for_graph[name]
                                            end_mini_seq = seq[inter[1]:inter[1] + k]
                                            r_infos = Read_infos(inter[0], inter[1], True)
                                            node_sequence[name] = end_mini_seq
                                            prev_nodelist[r_id] = r_infos
                                            nodes_for_graph[name] = prev_nodelist
                                            # if DEBUG:
                                            # print("Adding infosd to node", name, " : ", prev_nodelist)
                                            # if name == '601, 652, 20':
                                            #    print("Adding infosf to the node", name, " : ", prev_nodelist)
                                            if not DG.has_edge(previous_node, name):
                                                if DEBUG:
                                                    print("no topo alt results")
                                                DG.add_edge(previous_node, name, length=length)
                                                edge_support[previous_node, name] = []
                                                edge_support[previous_node, name].append(r_id)
                                            else:
                                                edge_info = edge_support[previous_node, name]
                                                if not r_id in edge_info:
                                                    edge_info.append(r_id)
                                                    edge_support[previous_node, name] = edge_info
                                    else:"""
                                        # update the read information of node name
                                        prev_nodelist = nodes_for_graph[name]
                                        r_infos = Read_infos(inter[0], inter[1], True)
                                        end_mini_seq = seq[inter[1]:inter[1] + k]
                                        node_sequence[name] = end_mini_seq
                                        prev_nodelist[r_id] = r_infos
                                        nodes_for_graph[name] = prev_nodelist
                                        # if DEBUG:
                                        # print("Adding infosd to node", name, " : ", prev_nodelist)
                                        # if name == '601, 652, 20':
                                        #    print("Adding infosf to the node", name, " : ", prev_nodelist)
                                        if not DG.has_edge(previous_node, name):
                                            if True:
                                                """if topo_graph.index(previous_node) > topo_graph.index(name):
                                                old_name = name
                                                # find alternatives for the node in topo_alternatives
                                                topo_alt_result = find_topo_alternative(topo_graph, previous_node, name,
                                                                                        topo_alternatives)
                                                # if we have not found an alternative in topo_alternatives:
                                                if not topo_alt_result:
                                                    length = this_len
                                                    # add new node that either is named 'name' or 'name1..1'
                                                    name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                                                    if not DG.has_node(name):
                                                        DG.add_node(name)
                                                    else:
                                                        if DEBUG:
                                                            print("Node", name, " already present")
                                                        while DG.has_node(name):
                                                            name = name + 1
                                                        if DEBUG:
                                                            print("Adding node", name)
                                                        DG.add_node(name)
                                                    nodelist = {}
                                                    r_infos = Read_infos(inter[0], inter[1], True)
                                                    end_mini_seq = seq[inter[1]:inter[1] + k]
                                                    node_sequence[name] = end_mini_seq
                                                    nodelist[r_id] = r_infos
                                                    # nodelist[r_id] = (inter[0], inter[1])
                                                    #prev_nodelist[r_id] = r_infos
                                                    nodes_for_graph[name] = nodelist
                                                    if not old_name in topo_alternatives:
                                                        # we store the new node as alternative for later possible topological misshaps
                                                        alt_list = []
                                                        alt_list.append(name)
                                                        topo_alternatives[old_name] = alt_list
                                                    else:
                                                        topo_alternatives[old_name].append(name)
                                                    if DEBUG:
                                                        print("702")
                                                    if DG.has_edge(previous_node, name):
                                                        edge_info = edge_support[previous_node, name]
                                                        if not r_id in edge_info:
                                                            edge_info.append(r_id)
                                                            edge_support[previous_node, name] = edge_info
                                                    else:
                                                        DG.add_edge(previous_node, name, length=length)
                                                        edge_support[previous_node, name] = []
                                                        edge_support[previous_node, name].append(r_id)
                                                    if r_id > 1:
                                                        topo_graph = list(nx.topological_sort(DG))
                                                else:
                                                    length = this_len
                                                    name = topo_alt_result
                                                    prev_nodelist = nodes_for_graph[name]
                                                    end_mini_seq = seq[inter[1]:inter[1] + k]
                                                    r_infos = Read_infos(inter[0], inter[1], True)
                                                    node_sequence[name] = end_mini_seq
                                                    prev_nodelist[r_id] = r_infos
                                                    nodes_for_graph[name] = prev_nodelist
                                                    # if DEBUG:
                                                    # print("Adding infosd to node", name, " : ", prev_nodelist)
                                                    # if name == '601, 652, 20':
                                                    #    print("Adding infosf to the node", name, " : ", prev_nodelist)
                                                    if not DG.has_edge(previous_node, name):
                                                        if DEBUG:
                                                            print("722")
                                                        DG.add_edge(previous_node, name, length=length)
                                                        edge_support[previous_node, name] = []
                                                        edge_support[previous_node, name].append(r_id)
                                                    else:
                                                        edge_info = edge_support[previous_node, name]
                                                        if not r_id in edge_info:
                                                            edge_info.append(r_id)
                                                            edge_support[previous_node, name] = edge_info
                                            else:"""
                                                if DG.has_edge(previous_node,name):
                                                    edge_info = edge_support[previous_node, name]
                                                    if not r_id in edge_info:
                                                        edge_info.append(r_id)
                                                        edge_support[previous_node, name] = edge_info
                                                else:
                                                    #DG.add_edge(previous_node, name, length=length)
                                                    #edge_support[previous_node, name] = []
                                                    #edge_support[previous_node, name].append(r_id)
                                                    DG.add_edge(previous_node, name, length=length)
                                                    cycle_added = isCyclic(DG)
                                                    cycle_added2 = cycle_finder(DG,previous_node)
                                                    #cycle_added3 = cycle_finder(DG,name)
                                                    #if cycle_added != cycle_added2 and cycle_added != cycle_added3:
                                                    #    print(cycle_added,cycle_added2,cycle_added3)
                                                    #    sys.exit(1)
                                                    #cycle_added = isCyclic2(DG, name)
                                                    #try:
                                                    #    topo_sort=nx.topological_sort(DG)
                                                    #    #print(topo_sort[0])
                                                    #except:
                                                    if cycle_added2:
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
                                                            print("adding edge clear 954", previous_node, ",", name)
                                                        edge_support[previous_node, name] = []
                                                        edge_support[previous_node, name].append(r_id)
                                                    else:
                                                        if DEBUG:
                                                            print("adding edge hope 958", previous_node, ",", name)
                                                        edge_support[previous_node, name] = []
                                                        edge_support[previous_node, name].append(r_id)
                                        else:
                                            edge_info = edge_support[previous_node, name]
                                            if not r_id in edge_info:
                                                edge_info.append(r_id)
                                                edge_support[previous_node, name] = edge_info
                                # if we have not found a node which this info can be added to
                                else:

                                    # add a new entry to alternative_nodes[old_node] to enable finding this instance
                                    alternative_nodes[old_node].append(alt_info_tuple)
                                    if DEBUG:
                                        print("oldnodeFurther", old_node)
                                    # print("name",name)
                                    # add the read information for the node
                                    nodelist = {}
                                    r_infos = Read_infos(inter[0], inter[1], True)
                                    end_mini_seq = seq[inter[1]:inter[1] + k]
                                    node_sequence[name] = end_mini_seq
                                    # if DEBUG:
                                    #    print("NFG", nodes_for_graph['601, 652, 20'])

                                    nodelist[r_id] = r_infos
                                    # if DEBUG:
                                    #    print("NFGbef", nodes_for_graph['601, 652, 20'])
                                    nodes_for_graph[name] = nodelist
                                    # if DEBUG:
                                    #    print("NFGaft", nodes_for_graph['601, 652, 20'])
                                    # print("Adding infosg to node", name, " : ", nodelist)
                                    # if name == '601, 652, 20':
                                    #        print("Adding infosh to the node", name, " : ", nodelist)
                                    DG.add_node(name)
                                    # get the length between the previous end and this nodes start
                                    length = this_len
                                    # connect the node to the previous one
                                    if DEBUG:
                                        print("768")

                                    DG.add_edge(previous_node, name, length=length)
                                    if DEBUG:
                                        print("adding edge clear 999", previous_node, ",", name)
                                    edge_support[previous_node, name] = []
                                    edge_support[previous_node, name].append(r_id)
                                    #if r_id > 1:
                                     #   topo_graph = list(nx.topological_sort(DG))
                # keep known_intervals up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
                node_overview_read[r_id - 1].append(name)
            # if the information for the interval was not yet found in a previous read (meaning the interval is new)
            else:
                nodelist = {}
                # add a node into nodes_for_graph
                name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                # add the read information for the node
                seq = all_reads[r_id][1]
                r_infos = Read_infos(inter[0], inter[1], True)
                end_mini_seq = seq[inter[1]:inter[1] + k]
                node_sequence[name] = end_mini_seq
                nodelist[r_id] = r_infos
                nodes_for_graph[name] = nodelist
                # if DEBUG:
                #    print("Adding infosj to node", name, " : ", nodelist)
                # if name=='601, 652, 20':
                #    print("Adding infosk to node", name, " : ", nodelist)
                # keep known_intervals up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
                node_overview_read[r_id - 1].append(name)
                # get the length between the previous end and this nodes start
                length = inter[0] - previous_end
                prior_read_infos = add_prior_read_infos(inter, r_id, prior_read_infos, name, k)
                DG.add_node(name)
                # connect the node to the previous one
                # print("Adding edge from "+previous_node+" to "+name)
                if DEBUG:
                    print("803")
                DG.add_edge(previous_node, name, length=length)
                if DEBUG:
                    print("adding edge clear 1035", previous_node, ",", name)
                edge_support[previous_node, name] = []
                edge_support[previous_node, name].append(r_id)
                #if r_id > 1:
                #    topo_graph = list(nx.topological_sort(DG))
            # set the previous node for the next iteration
            previous_node = name
            if DEBUG:
                print("PREVNODE",previous_node)
            previous_end = inter[1]
            # find out whether the current hash has already been added to read_hashs
            if not is_repetative:
                # if the current hash was not yet present in read_hashs: Add it
                read_hashs[curr_hash] = []
                read_hashs[curr_hash].append(name)
            # If the current hash was already present in read_hashs: This node is the second occurance of a node, forming a cycle. We have  to record the cycle to make sure other reads can be added if possible
            else:
                current_read_state = known_intervals[r_id - 1]
                newcyc, startnode = record_cycle(current_read_state, cycle_start)
                key_tuple = (r_id, startnode[1])
                known_cycles[key_tuple] = newcyc

                # print("Keytuple")
                # print(key_tuple)
                # print("Current hash already found in the interval.")
                # print(info_tuple)
                read_hashs[curr_hash].append(name)
                containscycle = True
        cycle_in_reads[r_id] = containscycle
        # add an edge from name to "t" as name was the last node in the read
        if not DG.has_edge(name, "t"):
            DG.add_edge(name, "t", length=0)
            if DEBUG:
                print("adding edge clear 1066", name, ",", "t")
            edge_support[name, "t"] = []
            edge_support[name, "t"].append(r_id)
        else:
            edge_info = edge_support[name, "t"]
            if not r_id in edge_info:
                edge_info.append(r_id)
                edge_support[name, "t"] = edge_info
        # print("Finished for read "+str(r_id))
        # print(edge_support[name,"t"])
        #if r_id==1:
        #    topo_graph = list(nx.topological_sort(DG))
    # set the node attributes to be nodes_for_graph, very convenient way of solving this
    # print("NFG",nodes_for_graph)
    #print(nodes_for_graph['59, 102, 50'])
    nx.set_node_attributes(DG, nodes_for_graph, name="reads")
    nx.set_node_attributes(DG, node_sequence, name="end_mini_seq")
    nx.set_edge_attributes(DG, edge_support, name='edge_supp')
    # use the known_intervals data structure to be able to verify the number of nodes appointed to each read
    # check_graph_correctness(known_intervals,all_intervals_for_graph)
    # return DG and known_intervals, used to
    result = (DG, known_intervals, node_overview_read, reads_for_isoforms, reads_at_start_dict)
    return result


def check_graph_correctness(known_intervals, all_intervals_for_graph):
    correct_graph = True
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        if not (len(intervals_for_read) == len(known_intervals[r_id - 1])):
            correct_graph = False
            # print("graph incorrect for read "+str(r_id))
        # print("All_intervals_for_graph at read "+str(r_id)+" has "+str(len(intervals_for_read))+" Intervals")
        for inter in intervals_for_read:
            read_id = inter[3][slice(0, len(inter[3]), 3)]
            number_of_reads_all_intervals = len(read_id)
        # print("Known_intervals at read " + str(r_id) + " has " + str(len(known_intervals[r_id-1])) + " Intervals")

    # if correct_graph:
    #    print("Graph built up correctly")
    # else:
    #    print("ERROR - Incorrect Graph")


"""
Plots the given Directed Graph via Matplotlib
Input: DG   Directed Graph
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
    """Helper method which extracts the read lengths from all_reads. We will use those during the graph generation to appoint more meaningful information to the node't'
    INPUT: all_reads: dictionary which holds the overall read infos key: r_id, value tuple(readname, sequence, some info i currently don't care about)
    OUTPUT: readlen_dict: dictionary holding the read_id as key and the length of the read as value
    """


def get_read_lengths(all_reads):
    readlen_dict = {}
    for r_id, infos in all_reads.items():
        seq = infos[1]
        seqlen = len(seq)

        readlen_dict[r_id] = seqlen
    return readlen_dict


"""
USED FOR DEBUGGING ONLY-deprecated in IsONform
"""


def main():
    #import sys
    #sys.stdout = open('log.txt', 'w')
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
    #is_cyclic = True
    #while is_cyclic:
    read_len_dict = get_read_lengths(all_reads)
    print("Generating graph")
    DG, known_intervals, node_overview_read, reads_for_isoforms, reads_list = generateGraphfromIntervals(
            all_intervals_for_graph, k_size, delta_len, read_len_dict, all_reads)

    print("edges with attributes:")
    # print(DG.edges(data=True))
    # print(DG.nodes(data=True))
    # print(known_intervals)
    # check_graph_correctness(known_intervals,all_intervals_for_graph)
    # print("#Nodes for DG: " + str(DG.number_of_nodes()) + " , #Edges for DG: " + str(DG.number_of_edges()))
    # edgelist = list(DG.edges.data())
    # print(DG.nodes(data=True))
    # print(type(DG.nodes(data=True)))
    # draw_Graph(DG)
    # simplifyGraph(DG, delta_len,all_reads,work_dir,k_size)
    print("Calling the method")

    simplifyGraph(DG, all_reads,work_dir,k_size,delta_len,known_intervals)
    print("Graph simplified")
    # simplifyGraphOriginal(DG,delta_len, all_reads, work_dir, k_size,known_intervals)
    # print("#Nodes for DG: " + str(DG.number_of_nodes()) + " , #Edges for DG: " + str(DG.number_of_edges()))
    # draw_Graph(DG)
    # print(DG.nodes["s"]["reads"])
    # print(list(DG.nodes(data=True)))
    # print("edges with attributes:")
    # print(DG.edges(data=True))
    # print("ReadNodes")
    # print(node_overview_read)
    # print("all edges for the graph")
    # print([e for e in DG.edges])
    # draw_Graph(DG)
    # print("knownintervals")
    # print(known_intervals)
    # The call for the isoform generation (deprecated during implementation)
    # TODO: invoke isoform_generation again
    generate_isoforms(DG, all_reads, reads_for_isoforms, work_dir, outfolder, max_seqs_to_spoa)
    # DG.nodes(data=True)

    # print("removing temporary workdir")
    shutil.rmtree(work_dir)


if __name__ == "__main__":
    main()
