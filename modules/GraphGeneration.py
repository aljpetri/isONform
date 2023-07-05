
import networkx as nx
from collections import namedtuple
from modules import SimplifyGraph


def cycle_finder(DG, start_node):
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
        visited.add(current_node)
        for neighbor in DG.neighbors(current_node):
            if neighbor in visited:
                # We have found a cycle!
                return True
            else:
                path.append(neighbor)
    # No cycle was found
    return False


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


"""Helper method for generateGraphfromIntervals used to identify a cycles that we may have seen before 
    INPUT:  thisstartpos: the start position of the 
            info_array:
            known_cycles:
            current_read_state:
            k:
            cycle_start:
            previous_end:
            DG:
            delta_len:
    OUTPUT: """


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

def generateGraphfromIntervals(all_intervals_for_graph, k, delta_len, read_len_dict, all_reads):
    Read_infos = namedtuple('Read_Infos',
                            'start_mini_end end_mini_start original_support')
    DEBUG = False
    DG = nx.DiGraph()
    # add the read ids to the startend_list
    reads_at_start_dict = {}
    reads_at_end_dict = {}
    reads_for_isoforms = []
    for i in range(1, len(all_intervals_for_graph) + 1):
        reads_at_start_dict[i] = Read_infos(0, 0, True)
        reads_for_isoforms.append(i)
    for i in range(1, len(read_len_dict) + 1):
        reads_at_end_dict[i] = Read_infos(read_len_dict[i], read_len_dict[i], True)

    # a source and a sink node are added to the graph in order to have a well-defined start and end for the paths
    DG.add_node("s", reads=reads_at_start_dict, end_mini_seq='')
    # print("REads at end dict",reads_at_end_dict)
    DG.add_node("t", reads=reads_at_end_dict, end_mini_seq='')

    # holds the r_id as key and a list of tuples as value: For identification of reads, also used to ensure correctness of graph
    known_intervals = []
    # adds an empty list for each r_id to known_intervals. To those lists, tuples, representing the intervals are added
    # for _ in itertools.repeat(None, len(all_intervals_for_graph)):
    for i in range(len(all_intervals_for_graph)):
        # known_intervals.append(OrderedSet([]))
        known_intervals.append([])
    # print(known_intervals)
    nodes_for_graph = {}
    alternative_nodes = {}

    # the following dictionary is supposed to hold the end minimizer sequence for each node
    node_sequence = {}
    edge_support = {}
    prior_read_infos = {}
    cycle_in_reads = {}
    known_cycles = {}
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    # intervals_for_read holds all intervals which make up the solution for the WIS of a read
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        if DEBUG:
            print('CURR READ:', r_id)
        containscycle = False
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"
        previous_end = 0
        # read_hashs is a dictionary storing hash values as keys and the respective node ids as values
        read_hashs = {}

        # the name of each node is defined to be readID, startminimizerpos , endminimizerpos
        # iterate over all intervals, which are in the solution of a certain read
        for pos,inter in enumerate(intervals_for_read):
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
            # access prior_read_infos, if the same interval was already found in previous reads)
            if info_tuple in prior_read_infos:
                if DEBUG:
                    print("in prior read infos")
                # if the interval repeats during this read
                if is_repetative:
                    if DEBUG:
                        print("Repetative")
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
                            prev_nodelist = nodes_for_graph[name]
                            seq = all_reads[r_id][1]
                            r_infos = Read_infos(inter[0], inter[1], True)
                            end_mini_seq = seq[inter[1]:inter[1] + k]
                            if not (end_mini_seq == DG.nodes[name]['end_mini_seq']):
                                print("ERROR: ", end_mini_seq, " not equal to ", DG.nodes[name]['end_mini_seq'])
                            prev_nodelist[r_id] = r_infos
                            nodes_for_graph[name] = prev_nodelist
                            length = this_len
                            if DEBUG:
                                print("Adding edge, foundnode, topo order not violated  for prevnode, name")
                            DG.add_edge(previous_node, name, length=length)
                            cycle_added2 = cycle_finder(DG,previous_node)
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
                            prev_nodelist[r_id] = r_infos
                            nodes_for_graph[name] = prev_nodelist
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
                # the interval did not repeat during this read
                else:
                    if DEBUG:
                        print("non repetative interval")
                    # get the name by accessing prior_read_infos
                    name = prior_read_infos[info_tuple]
                    # get length from previous_end to this_start
                    this_len = inter[0] - previous_end

                    # if there is not yet an edge connecting previous_node and name: add edge
                    if not DG.has_edge(previous_node, name):
                        if DEBUG:
                            print("no edge")
                        # update the read information of node name
                        prev_nodelist = nodes_for_graph[name]
                        seq = all_reads[r_id][1]
                        r_infos = Read_infos(inter[0], inter[1], True)
                        end_mini_seq = seq[inter[1]:inter[1] + k]
                        node_sequence[name] = end_mini_seq
                        prev_nodelist[r_id] = r_infos
                        nodes_for_graph[name] = prev_nodelist
                        # only add a new edge if the edge was not present before
                        length = this_len
                        DG.add_edge(previous_node, name, length=length)
                        cycle_added2 = cycle_finder(DG,previous_node)
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
                            if DEBUG:
                                print("adding edge clear 686", previous_node, ",", name)
                            edge_support[previous_node, name] = []
                            edge_support[previous_node, name].append(r_id)
                        else:
                            if DEBUG:
                                print("adding edge hope690", previous_node, ",", name)
                            edge_support[previous_node, name] = []
                            edge_support[previous_node, name].append(r_id)

                    # TODO: add topological sort and everything after this point
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
                            if not r_id in edge_info:
                                edge_info.append(r_id)
                                edge_support[previous_node, name] = edge_info
                        # if the length difference is >delta len: generate new node, to be able to tell those nodes apart we use alternative_nodes
                        else:
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
                                if DEBUG:
                                    print("delta len < actual diff, adding new node ")
                                DG.add_edge(previous_node, name, length=this_len)
                                if DEBUG:
                                    print("adding edge clear 750", previous_node, ",", name)
                                edge_support[previous_node, name] = []
                                edge_support[previous_node, name].append(r_id)
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

                                    # update the read information of node name
                                    prev_nodelist = nodes_for_graph[name]
                                    r_infos = Read_infos(inter[0], inter[1], True)
                                    end_mini_seq = seq[inter[1]:inter[1] + k]
                                    node_sequence[name] = end_mini_seq
                                    prev_nodelist[r_id] = r_infos
                                    nodes_for_graph[name] = prev_nodelist
                                    if not DG.has_edge(previous_node, name):
                                        if DG.has_edge(previous_node,name):
                                            edge_info = edge_support[previous_node, name]
                                            if not r_id in edge_info:
                                                edge_info.append(r_id)
                                                edge_support[previous_node, name] = edge_info
                                        else:
                                            DG.add_edge(previous_node, name, length=length)
                                            cycle_added2 = cycle_finder(DG,previous_node)
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
                                    if DEBUG:
                                        print("768")

                                    DG.add_edge(previous_node, name, length=length)
                                    if DEBUG:
                                        print("adding edge clear 999", previous_node, ",", name)
                                    edge_support[previous_node, name] = []
                                    edge_support[previous_node, name].append(r_id)
                # keep known_intervals up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
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
                # keep known_intervals up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
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
    nx.set_node_attributes(DG, nodes_for_graph, name="reads")
    nx.set_node_attributes(DG, node_sequence, name="end_mini_seq")
    nx.set_edge_attributes(DG, edge_support, name='edge_supp')
    # return DG and known_intervals, used to
    result = (DG, reads_for_isoforms)
    return result



