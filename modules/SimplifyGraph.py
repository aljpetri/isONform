from collections import namedtuple

import itertools
import networkx as nx
import os



from modules import consensus
from modules import IsoformGeneration

class Readtup:
    def __init__(self, path, supp, non_supp):
        self.path = path
        self.supp = supp
        self.non_supp = non_supp


Read_infos = namedtuple('Read_Infos',
                        'start_mini_end end_mini_start original_support')

CNode = namedtuple('CNode', 'visited recStack')


def isCyclicUtil(DG, nodes_dict, node):
    """
    Used to find out whether adding 'node' led to the graph having a cycle
    Implemented according to https://www.geeksforgeeks.org/detect-cycle-in-a-graph/
    """
    # Mark current node as visited and
    # adds to recursion stack

    cnode = CNode(True, True)
    nodes_dict[node] = cnode
    # Recur for all neighbours
    # if any neighbour is visited and in
    # recStack then graph is cyclic
    for out_edge in DG.out_edges(node):
        neighbour = out_edge[1]
        if not nodes_dict[neighbour].visited:
            if isCyclicUtil(DG, nodes_dict, neighbour):
                #print(neighbour, " revisited")
                return True
        elif nodes_dict[neighbour].recStack:
            return True
    # The node needs to be popped from
    # recursion stack before function ends
    prev_visited = nodes_dict[node].visited
    nodes_dict[node] = CNode(prev_visited, False)
    return False


def isCyclic(DG):
    """Helper method used to detect cycles in our graph:
    INPUT: DG       our directed Graph
    OUTPUT: iscyclic    A boolean value indicating whether a cycle was found
    """
    # Returns true if graph contains cycles else false
    nodes = DG.nodes
    nodes_dict = {}
    #CNode = namedtuple('CNode', 'visited recStack', defaults=(False, False))
    cnode = CNode(False, False)
    for node in nodes:
        nodes_dict[node] = cnode
    for node in nodes:
        if not nodes_dict[node].visited:
            if isCyclicUtil(DG, nodes_dict, node):
                return True
    return False


def find_possible_starts(DG, topo_nodes_dict, possible_starts):
    """Helper method which finds all possible bubble starts in our graph. This is done by collecting all nodes having at
     least 2 out nodes
    INPUT:  DG              Our directed graph object
            TopoNodes       A list of our nodes in topological order
    OUTPUT: possible_starts:    A list holding start_tuples(node_id, start_supp:= all reads in this node)
    """
    # iterate over all nodes in TopoNodes
    for node in topo_nodes_dict.keys():
        # the current node is only a start node if it has an out_degree>1 (due to bubble definition)
        if DG.out_degree(node) > 1:
            # find all reads supporting this node
            start_supp = tuple(DG.nodes[node]['reads'])
            # store node id and support in start_tup
            start_tup = (node, start_supp)
            # append start_tup to the list of possible bubble_starts
            possible_starts.append(start_tup)


def find_possible_ends(DG, topo_nodes_dict, possible_ends):
    """Helper method which finds all possible bubble starts in our graph. This is done by collecting all nodes having at
     least 2 out nodes
    INPUT:  DG              Our directed graph object
            TopoNodes       A list of our nodes in topological order
    OUTPUT: possible_starts:    A list holding start_tuples(node_id, start_supp:= all reads in this node)
    """
    # iterate over all nodes in TopoNodes
    for node in topo_nodes_dict.keys():
        # the current node is only an end node if it has an in_degree>1 (due to bubble definition)
        if DG.in_degree(node) > 1:
            # find all reads supporting this node
            end_supp = tuple(DG.nodes[node]['reads'])
            # store node id and support in end_tup
            end_tup = (node, end_supp)
            # append end_tup to the list of possible bubble_ends
            possible_ends.append(end_tup)


def generate_combinations(possible_starts, possible_ends, topo_nodes_dict, combis):
    """used to generate the combinations of start and end nodes
    """
    for startnode in possible_starts:
        for endnode in possible_ends:
            if topo_nodes_dict[startnode[0]] < topo_nodes_dict[endnode[0]]:
                inter = tuple(sorted(set(startnode[1]).intersection(set(endnode[1]))))
                if len(inter) >= 2:
                    combi = (startnode[0], endnode[0], inter)
                    combis.append(combi)


def filter_combinations(combinations, not_viable, combinations_filtered):
    """we filter out combinations that already have been deemed as not_viable
    """
    # iterate over all combinations
    for combi in combinations:
        # if the combination is viable add it to combinations_filtered
        if combi not in not_viable:
            combinations_filtered.append(combi)

#TODO add marked to stop finding a path as soon as we encounter a marked path(this is due to the path not being usable anymore this iteration)
def find_paths(DG, startnode, endnode, support, all_paths, marked):
    """detect the paths in our bubble
    """
    node_support_left = set(support)
    all_supp = set(support)
    # we iterate as long as still not all support was allocated to a path
    while node_support_left:
        node = startnode
        read = node_support_left.pop()
        current_node_support = node_support_left
        current_node_support.add(read)
        visited_nodes = []
        # As long as we have not visited the bubble end node we continue walking through our graph
        while node != endnode:
            visited_nodes.append(node)
            next_found = False
            ###improvement
            if node in marked:

                break
            ###
            #out_edges = DG.out_edges(node)
            #for edge in out_edges:
            for neighbor in DG.successors(node):
                #edge_supp = DG[edge[0]][edge[1]]['edge_supp']
                edge_supp=DG[node][neighbor]['edge_supp']
                if read in edge_supp:
                    #node = edge[1]
                    node=neighbor
                    current_node_support = current_node_support.intersection(edge_supp)
                    next_found = True
                    break
            if not next_found:
                break
        if current_node_support and next_found:
            node_support_left -= current_node_support
            final_add_support = all_supp - current_node_support
            path_supp_tup = (visited_nodes, tuple(sorted(current_node_support)), final_add_support)
            all_paths.append(path_supp_tup)
        else:
            break


def parse_cigar_diversity(cigar_tuples, delta_perc, delta_len):
    miss_match_length = 0
    alignment_len = 0
    for i, elem in enumerate(cigar_tuples):
        cig_len = elem[0]
        cig_type = elem[1]
        alignment_len += cig_len
        if (cig_type != '=') and (cig_type != 'M'):
            # we want to add up all missmatches to compare to sequence length
            miss_match_length += cig_len
            # we also want to make sure that we do not have too large internal sequence differences
            if cig_len > delta_len:
                return False
    diversity = (miss_match_length / alignment_len)
    max_bp_diff = max(delta_perc * alignment_len, delta_len)
    mod_div_rate = max_bp_diff / alignment_len
    if DEBUG:
        print("diversity", diversity, "mod_div", mod_div_rate)
    if diversity <= mod_div_rate:
        return True
    else:
        return False


def parse_cigar_differences(cigar_string, delta_len):
    """Helper method used to test whether two sequences are close enough to pop their underlying bubble
    INPUT:      cigar_string    a string denoting the cigar output of the alignment
                delta_len       a threshold which we use to distinguish between mutations and errors
    OUTPUT:     good_to_pop     boolean value indicating our finding
    """
    good_to_pop = True
    for i, elem in enumerate(cigar_string):
        cig_len = elem[0]
        cig_type = elem[1]
        # all matches are absolutely fine
        if (cig_type != '=') and (cig_type != 'M'):
            # we have a non match, now we have to figure out whether this is due to an exon or indel (everything with len>delta_len is defined as exon)
            if cig_len > delta_len:
                # we found an exon being the difference between the paths->no bubble popping feasible
                good_to_pop = False
    return good_to_pop


def get_dist_to_prev(DG, prev_node, curr_node):
    prev_supp = DG.nodes[prev_node]['reads']
    curr_supp = DG.nodes[curr_node]['reads']
    intersect_supp = list(set(prev_supp).intersection(curr_supp))
    summation = 0
    avg_dist = 0
    # we did not yet come up with anything better than calculating the average distance for all reads supporting both nodes
    for i, read in enumerate(intersect_supp):
        prev_pos = prev_supp[read].end_mini_start
        curr_pos = curr_supp[read].end_mini_start
        summation += (curr_pos - prev_pos)
        avg_dist = summation / (i + 1)
    return avg_dist


def new_distance_to_start(DG, pnl_start, curr_node, consensus_log):
    seq = consensus_log[pnl_start]
    node_seq = DG.nodes[curr_node]['end_mini_seq']
    possible_pos = seq.find(node_seq)
    return possible_pos


def remove_edges(DG, bubble_start, bubble_end, path_nodes, consensus_log, edges_to_delete, node_distances):
    """
    Helper method utilized by linearize bubbles which removes the edges we have to get rid of(all edges which connect 2 nodes which are part of the current bubble
    Alters edges_to_delete.
        INPUT:      DG                  Our directed Graph
                    bubble_start        The start node of our bubble
                    bubble_end:         The end node of our bubble
                    path_nodes          A list of nodes which make up a path in our bubble
                    consensus_log:
                    edges_to_delete:
                    node_distances:
        OUTPUT:     edges_to_delete         A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
    """
    if DEBUG:
        print("Bubble_start", bubble_start)
        print("Bubble_end", bubble_end)
    # we have to delete the edges which connect the nodes that are inside the bubble
    for one_info in path_nodes:
        path_node_list = one_info[0]
        if path_node_list:
            path_node_list.pop(0)
        if not path_node_list:
            pnl_start = bubble_end
        else:
            pnl_start = path_node_list[0]
        prev_node = bubble_start
        entry = DG.get_edge_data(prev_node, pnl_start)
        edges_to_delete[prev_node, pnl_start] = entry
        if not path_node_list:
            pnl_start = bubble_end
        else:
            pnl_start = path_node_list[0]
        for index, path_node in enumerate(path_node_list):
            inter_dist = new_distance_to_start(DG, pnl_start, path_node, consensus_log)
            if inter_dist == -1:
                dist_to_prev = get_dist_to_prev(DG, prev_node, path_node)
                # we found the distance to the previous node, however we are still missing the distance of the previous node to s
                if prev_node == bubble_start:
                    prev_to_start_dist = 0
                else:
                    prev_to_start_dist = node_distances[prev_node]
                dist = prev_to_start_dist + dist_to_prev
            else:
                dist = inter_dist
            node_distances[path_node] = dist
            if path_node != path_node_list[-1]:
                edges_to_delete[path_node, path_node_list[index + 1]] = DG[path_node][path_node_list[index + 1]]
            else:
                entry = DG.get_edge_data(path_node, bubble_end)
                edges_to_delete[path_node, bubble_end] = entry
    for edge, edge_infos in edges_to_delete.items():
        DG.remove_edge(edge[0], edge[1])
    return edges_to_delete, node_distances


def get_avg_interval_length(DG, node):
    node_support = DG.nodes[node]['reads']
    summation = 0
    i = 0
    for r_id, positions in node_support.items():
        # print(positions)
        i += 1
        summation += (positions.end_mini_start - positions.start_mini_end)
        # print(sum)
    if i == 0:
        finalresult = 0
    else:
        finalresult = int(summation / i)
    return finalresult


def additional_node_support(DG, new_support, node_dist_dict, node, other_prevnode, bubble_start, additional_support):
    """Helper method: Adds additional support values for the bubble_nodes(needed to get a consistent graph)
    INPUT:      DG:
                nextnode:       the node of which we need the support
    OUTPUT:     additional_support:  dictionary of additional support values
    """
    if other_prevnode == bubble_start:
        other_dist = 0
    else:
        other_dist = node_dist_dict[other_prevnode]
    this_dist = node_dist_dict[node]
    # this will contain all reads that are to be added to the node with their respective (virtual) positions
    this_reads = DG.nodes[node]['reads']
    avg_len = get_avg_interval_length(DG, node)
    for r_id in new_support:
        # figure out whether the r_id is part of this bubble_path
        if r_id not in this_reads:
            # what we do if the read meets the bubble and is not in this_reads
            # figure out if r_id is also represented by bubble_start
            # r_id is not in s_infos (meaning the read is not supporting bubble_start), if it would be in s_infos, we do not have to recompute anything
            # get the position of the read in previous node
            # calculate relative distance to start (relative to previous node)
            previous_other_path_reads = DG.nodes[other_prevnode]['reads']
            if r_id in previous_other_path_reads:
                pos_info_tup = previous_other_path_reads[r_id]
                prev_end = pos_info_tup.end_mini_start
                relative_dist = int(this_dist - other_dist)
                newend = prev_end + relative_dist
                newstart = int(newend - avg_len)
                additional_support[r_id] = Read_infos(newstart, newend, False)


def merge_two_dicts(dict1, dict2):
    """Helper method used to merge two dicts with the values being lists
    INPUT: dict1: the first dict
           dict2: the second dict
    OUTPUT: merged_dict: dictionary containing all keys in dict1 and dict2
    """
    merged_dict = {}
    for key, value in dict1.items():
        merged_dict[key] = value
    for key2, val2 in dict2.items():
        if key2 not in merged_dict:
            merged_dict[key2] = val2
    return merged_dict


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


def compare_by_length(nextnode1, nextnode2, bubble_end, nn_contender1, nn_contender2):
    if nextnode1 == bubble_end:
        return nextnode2
    elif nextnode2 == bubble_end:
        return nextnode1
    return nextnode1 if nn_contender1 < nn_contender2 else nextnode2


def find_real_nextnode(nextnode1, nextnode2, bubble_end, DG, topo_nodes_dict):
    if DG.has_edge(nextnode1, nextnode2):
        return nextnode1
    elif DG.has_edge(nextnode2, nextnode1):
        return nextnode2
    else:
        nn_contender1 = topo_nodes_dict[nextnode1]
        nn_contender2 = topo_nodes_dict[nextnode2]
        nextnode = compare_by_length(nextnode1, nextnode2, bubble_end, nn_contender1, nn_contender2)
        return nextnode


def get_next_node(path, bubble_end):
    if len(path) < 2:
        nextnode = bubble_end
    else:
        nextnode = path[1]
    return nextnode


"""
Helper method that finds edges that connect nodes between different paths
INPUT:



OUTPUT:
"""


def find_connecting_edges(path_nodes, DG):
    connecting_edges = set()
    path1 = path_nodes[0][0]
    path2 = path_nodes[1][0]
    for node in path1:
        for onode in path2:
            if DG.has_edge(node, onode):
                connecting_edges.add((node, onode))
            elif DG.has_edge(onode, node):
                connecting_edges.add((onode, node))

    return connecting_edges


def test_conn_end(conn_edges, overall_nextnode):
    conn_list = []
    for conn_edge in conn_edges:
        if overall_nextnode == conn_edge[1]:
            conn_list.append(conn_edge)
    return conn_list


def prepare_adding_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes, node_dist,
                         topo_nodes_dict):  # ,node_dist):
    counter = 0
    path1 = []
    path2 = []
    edge_params = {}
    linearization_order = []
    conn_edges = find_connecting_edges(path_nodes, DG)
    # we assign both paths to variables to make them easier accessible.
    for path_info in path_nodes:
        nodes = path_info[0]
        if counter == 0:
            for p in nodes:
                path1.append(p)
        else:
            for p in nodes:
                path2.append(p)
        counter = counter + 1
    # this dictionary contains the reads with their respected positions that are added to the nodes
    new_node_supp_dict = {}
    linearization_order.append(bubble_start)
    prevnode1 = bubble_start
    prevnode2 = bubble_start
    if path1:
        nextnode1 = path1[0]
    else:
        nextnode1 = bubble_end
    if path2:
        nextnode2 = path2[0]
    else:
        nextnode2 = bubble_end
    prevnode = bubble_start
    while path1 or path2:
        overall_nextnode = find_real_nextnode(nextnode1, nextnode2, bubble_end, DG, topo_nodes_dict)
        conn_list = test_conn_end(conn_edges, overall_nextnode)
        new_edge_supp1 = edges_to_delete[prevnode1, nextnode1]['edge_supp']
        new_edge_supp2 = edges_to_delete[prevnode2, nextnode2]['edge_supp']
        if not conn_list:
            full_edge_supp = new_edge_supp1 + new_edge_supp2
        else:
            this_edge_supp = DG[conn_list[0][0]][overall_nextnode]["edge_supp"]
            full_edge_supp = new_edge_supp1 + new_edge_supp2 + this_edge_supp
        if overall_nextnode in path1:
            nextnode1 = get_next_node(path1, bubble_end)
            additional_supp = {}
            # TOD: we need to add the support of global_prev for both additional_node_support occurrences
            additional_node_support(DG, new_edge_supp2, node_dist, overall_nextnode,
                                    prevnode2, bubble_start, additional_supp)
            path1.remove(overall_nextnode)
            prevnode1 = overall_nextnode
            curr_node = prevnode1
        elif overall_nextnode in path2:
            nextnode2 = get_next_node(path2, bubble_end)
            additional_supp = {}
            additional_node_support(DG, new_edge_supp1, node_dist, overall_nextnode,
                                    prevnode1, bubble_start, additional_supp)
            path2.remove(overall_nextnode)
            prevnode2 = overall_nextnode
            curr_node = prevnode2
        old_node_supp = DG.nodes[overall_nextnode]['reads']
        new_node_supp_dict[overall_nextnode] = merge_two_dicts(additional_supp, old_node_supp)
        linearization_order.append(overall_nextnode)
        nx.set_node_attributes(DG, new_node_supp_dict, "reads")
        edge_params[prevnode, overall_nextnode] = full_edge_supp
        prevnode = curr_node

    new_edge_supp1 = edges_to_delete[prevnode1, bubble_end]['edge_supp']
    new_edge_supp2 = edges_to_delete[prevnode2, bubble_end]['edge_supp']
    full_edge_supp = new_edge_supp1 + new_edge_supp2
    edge_params[prevnode, bubble_end] = full_edge_supp
    linearization_order.append(bubble_end)
    # adding the new edges to the Graph
    for key, value in edge_params.items():
        # if the edge has not been in the graph before
        if not DG.has_edge(key[0], key[1]):
            DG.add_edge(key[0], key[1], edge_supp=value)
        else:
            # The edge was in the graph before, add old support to keep the graph consistent
            old_edge_supp = DG.edges[key[0], key[1]]['edge_supp']
            for old_supp in old_edge_supp:
                if old_supp not in value:
                    value.append(old_supp)
            DG.add_edge(key[0], key[1], edge_supp=value)
    # this is the main part of the linearization. We iterate over all_nodes and try to find out which path the nodes belong to.
    # This info is needed as we need the current state ob both paths to add the correct edge_support and node_support to the graph


def generate_consensus_path(work_dir, consensus_attributes, reads, k_size, spoa_count):
    """Helper method used to generate the consensus sequence for each bubble path
    INPUT:      work_dir:               The work directory we are in
                consensus_attributes:   list of tuples containing the read_id and the positions that make up the limits of our subsequences
                reads:                  dictionary containing the actual reads (the sequences as well as their original ids)
                k_size:                 the length of parameter k needed to extract the correct subsequence length
    OUTPUT:     spoa_ref:               the consensus sequence generated from the reads

    """
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    seq_infos = {}
    endseqlist = []
    reads_path_len = 0
    max_len = 0
    if len(consensus_attributes) > 2:
        for i, (q_id, pos1, pos2) in enumerate(consensus_attributes):
            if pos2 == 0:
                pos2 = len(reads[q_id][1]) - k_size
            if pos1 < (pos2 + k_size):
                if pos1 < 0:
                    pos1 = 0
                seq = reads[q_id][1][pos1: (pos2 + k_size)]
            else:
                seq = ""
            seq_infos[q_id] = (pos1, pos2 + k_size, seq)
            endseq = reads[q_id][1][pos2:pos2 + k_size]
            endseqlist.append(endseq)
            if len(seq) < k_size:
                if len(seq) > max_len:
                    max_len = len(seq)
                elif len(seq) < 1:
                    max_len = 1
            else:
                reads_path_len += 1
                reads_path.write(">{0}\n{1}\n".format(str(q_id) + str(pos1) + str(pos2), seq))
        reads_path.close()
        if reads_path_len > 0:
            spoa_count += 1
            spoa_ref = IsoformGeneration.run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"))
            return spoa_ref, seq_infos, spoa_count
        else:
            string_val = "X" * max_len  # gives you "xxxxxxxxxx"
            return string_val, seq_infos, spoa_count
    else:
        f_id, fstart, fend = consensus_attributes[0]
        e_id, estart, eend = consensus_attributes[1]
        fdist = fend - fstart
        edist = eend - estart
        if fdist > edist:
            consensus = reads[f_id][1][fstart: (fend + k_size)]
            seq_infos[f_id] = (fstart, fend + k_size, consensus)
        else:
            consensus = reads[e_id][1][estart: (eend + k_size)]
            seq_infos[f_id] = (estart, eend + k_size, consensus)
        return consensus, seq_infos, spoa_count


def collect_consensus_reads(consensus_attributes):
    con_reads = set()
    for con_att in consensus_attributes:
        con_reads.add(con_att[0])
    con_reads_fin = frozenset(con_reads)
    return con_reads_fin


def align_bubble_nodes(all_reads, consensus_infos, work_dir, k_size, spoa_count, multi_consensuses, is_megabubble,
                       combination, delta_len):
    consensus_list = []
    consensus_log = {}
    seq_infos = {}
    for path_node, consensus_attributes in consensus_infos.items():
        con_reads = collect_consensus_reads(consensus_attributes)
        combi = (combination[0], combination[1], con_reads)
        if combi in multi_consensuses:
            contains_combi = True
        else:
            contains_combi = False
        if is_megabubble and contains_combi:
            con_infos = multi_consensuses[combi]
            con = con_infos[0]
            spoa_count = con_infos[2]
            consensus_log[path_node] = con
            consensus_list.append(con)
        else:
            if len(consensus_attributes) > 1:
                con, seq_infos_from_fun, spoa_count = generate_consensus_path(work_dir, consensus_attributes, all_reads,
                                                                              k_size, spoa_count)
                if is_megabubble:
                    multi_consensuses[combi] = (con, seq_infos_from_fun, spoa_count)
                if len(con) < 3:
                    consensus_log[path_node] = ""
                    consensus_list.append("")
                else:
                    consensus_log[path_node] = con
                    seq_infos.update(seq_infos_from_fun)
                    consensus_list.append(con)
            else:
                (q_id, pos1, pos2) = consensus_attributes[0]
                if abs(pos2 - pos1) < 3:
                    consensus_log[path_node] = ""
                    consensus_list.append("")
                elif pos2 < pos1:
                    consensus_log[path_node] = ""
                    consensus_list.append("")

                else:
                    con = all_reads[q_id][1][pos1: pos2 + k_size]
                    seq_infos[q_id] = (pos1, pos2 + k_size, con)
                    consensus_log[path_node] = con
                    consensus_list.append(con)

    consensus1 = consensus_list[0]
    consensus2 = consensus_list[1]

    s1_len = len(consensus1)
    s2_len = len(consensus2)
    if s1_len > s2_len:
        longer_len = s1_len
        shorter_len = s2_len
    else:
        longer_len = s2_len
        shorter_len = s1_len
    delta = 0.20
    if shorter_len > delta_len and longer_len > delta_len:
        if (longer_len - shorter_len) / longer_len < delta:
            s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = consensus.parasail_alignment(consensus1, consensus2,
                                                                                               match_score=2,
                                                                                               mismatch_penalty=-8,
                                                                                               # standard: -8
                                                                                               opening_penalty=12,
                                                                                               gap_ext=1)  # opening penalty: standard: 12
            good_to_pop = parse_cigar_diversity(cigar_tuples, delta, delta_len)
            if DEBUG:
                print("GOODTOPOP?", good_to_pop)
            cigar_alignment = (s1_alignment, s2_alignment)
            return good_to_pop, cigar_alignment, seq_infos, consensus_log, spoa_count
        else:
            return False, "", "", consensus_log, spoa_count
    elif shorter_len < delta_len and longer_len < delta_len:
        return True, "", "", consensus_log, spoa_count

    else:
        if (longer_len - shorter_len) < delta_len:
            return True, "", "", consensus_log, spoa_count
        else:
            return False, "", "", consensus_log, spoa_count


def get_consensus_positions(bubble_start, bubble_end, DG, shared_reads):
    read_list = []
    max_node_infos = DG.nodes[bubble_end]['reads']
    min_node_infos = DG.nodes[bubble_start]['reads']
    # print("finding consensus_infos from ",bubble_start,"to ",bubble_end,", shared_reads ", str(shared_reads))
    for r_id in shared_reads:
        bubble_end_pos = max_node_infos[r_id]
        bubble_start_pos = min_node_infos[r_id]
        start_of_bubble = bubble_start_pos.end_mini_start
        end_of_bubble = bubble_end_pos.end_mini_start
        entry = (r_id, start_of_bubble, end_of_bubble)
        read_list.append(entry)
    # print("R_list",read_list)
    return read_list


"""
function which finds the bubbles and if feasible pops them in the graph
INPUT:
    DG:         Our directed graph 
    delta_len   Maximum length difference for which the paths are merged ie the bubble is popped
    all_reads   dictionary containing all the reads(string sequence) and their ids


OUTPUT:
    bubble_state, popable_bubbles,no_pop_tup_list"""


def find_disjoint_bubbles(poppable_bubbles):
    """Helper method that filters the poppable_bubbles to only contain non-overlapping bubbles, as overlapping bubbles yield inconsistencies in our graph
    INPUT:      poppable_bubbles:      a list of bubbles that are poppable
    OUTPUT:     new_poppable:          a list of bubbles that do not overlap->ready for popping
    """
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


def linearize_bubble(DG, bubble_start, bubble_end, path_nodes,
                     consensus_log, topo_nodes_dict):
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
    edges_to_delete = {}
    node_dist = {}
    # edges_to_delete: A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
    edges_to_delete, node_dist = remove_edges(DG, bubble_start, bubble_end, path_nodes, consensus_log, edges_to_delete,
                                              node_dist)
    prepare_adding_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes, node_dist,
                         topo_nodes_dict)


def extend_unpoppable(old_bubbles, poppable_state):
    """Helper method used to extend the list of bubbles that has been analysed but not been deemed poppable (old_bubbles)
    INPUT:      old_bubbles
                bubble_state
    """
    for newbubble in poppable_state:
        if not newbubble.bubble_supported:
            old_bubbles.append(newbubble.bubble_nodes)


def filter_marked_paths(all_paths, marked):
    is_marked = []
    for path_tup in all_paths:
        path = path_tup[0]
        already_marked = False
        for node in path:
            if node in marked:
                already_marked = True
        marked_tuple = (path, already_marked)
        is_marked.append(marked_tuple)
    for tup in is_marked:
        if tup[1]:
            all_paths.remove(tup)


def filter_out_if_marked(all_paths, marked, direct_combis, endnode):
    filter_list = []
    for path in all_paths:
        path_nodes = path[0]
        for node in path_nodes:
            if node in marked:
                if path not in filter_list:
                    filter_list.append(path)
                    break
        for direct_combi in direct_combis:
            # print("DC",direct_combi)
            startnode = direct_combi[0]
            for ident, node in enumerate(path_nodes):
                # make sure the path_nodes list is long enough
                if len(path_nodes) > ident + 1:
                    if node == startnode and path_nodes[ident + 1] == direct_combi[1]:
                        # we add the path to filter_list, as we have already seen the combi as direct combination
                        if path not in filter_list:
                            filter_list.append(path)
                        break
                else:
                    if node == startnode and endnode == direct_combi[1]:
                        if path not in filter_list:
                            filter_list.append(path)
                        break
    # remove all paths that we want to filter out from the initial list of paths
    if filter_list:
        for entry in filter_list:
            all_paths.remove(entry)
    return all_paths


def get_path_length(path_infos, DG, poss_start, poss_end):
    reads = path_infos[1]
    path_len_sum = 0
    start_dict = DG.nodes[poss_start]['reads']
    end_dict = DG.nodes[poss_end]['reads']
    for read in reads:
        start_pos = start_dict[read].end_mini_start
        end_pos = end_dict[read].end_mini_start
        # print("r_id",read)
        path_len = end_pos - start_pos
        path_len_sum += path_len
    avg_path_len = path_len_sum / (len(reads))
    return avg_path_len


def filter_path_if_marked(marked, path):
    # print("filter path if marked for ",marked," and ",path)
    for node in path:
        if node in marked:
            return True
    return False


def find_path(r_id, DG, edge_attr):
    current_node = "s"
    visited_nodes = []
    visited_edges = []

    reached_t = False
    while not reached_t:
        # add current node to the list of visited_nodes
        visited_nodes.append(current_node)
        prev_node = current_node
        # print("CurrnodebefMethod",current_node)
        # print()
        edgelist = list(DG.out_edges(current_node))
        for edge in edgelist:
            edge_reads = edge_attr[edge]
            if r_id in edge_reads:
                current_node = edge[1]

        # print("current node returned by get best supported edge node", current_node)
        edge_tup = (prev_node, current_node)
        # print("edge_tup",edge_tup)
        visited_edges.append(edge_tup)
        if current_node == "t":
            visited_nodes.append("t")
            reached_t = True
        else:
            # print("not there yet")
            if prev_node == current_node:
                print("ERROR")
                break
    return visited_nodes


def find_all_read_paths(DG, all_reads, merged_dict):
    all_read_paths = {}
    all_path_sets = {}
    edge_attr = nx.get_edge_attributes(DG, "edge_supp")
    for r_id in all_reads.keys():
        if r_id not in merged_dict:
            # print(r_id)
            all_read_paths[r_id] = find_path(r_id, DG, edge_attr)
        else:
            all_read_paths[r_id] = all_read_paths[merged_dict[r_id]]
    for r_id, path in all_read_paths.items():
        all_path_sets[r_id] = set(path)
    return all_read_paths, all_path_sets


def update_paths(DG, all_reads, prev_marked, merged_dict, all_paths_s_to_t, all_path_sets):
    new_all_paths_s_to_t = {}
    edge_attr = nx.get_edge_attributes(DG, "edge_supp")
    for r_id in all_reads.keys():
        # if we do not find a node in the readpath that has been marked as changed
        if set(all_path_sets[r_id]).isdisjoint(prev_marked):
            # copy the old path
            new_all_paths_s_to_t[r_id] = all_paths_s_to_t[r_id]
        # we found that the readpath has been changed->Regenerate the readpath
        else:
            if r_id in merged_dict.keys():
                new_all_paths_s_to_t[r_id] = new_all_paths_s_to_t[merged_dict[r_id]]
            else:
                new_all_paths_s_to_t[r_id] = find_path(r_id, DG, edge_attr)
            all_path_sets[r_id].update(new_all_paths_s_to_t[r_id])
    return new_all_paths_s_to_t, all_path_sets
    # TODO write this method to recalculate the paths of reads that were affected by bubble popping


def new_bubble_popping_routine(DG, all_reads, work_dir, k_size, delta_len, slowmode):
    """The backbone of bubble popping: the bubbles are detected and filtered to only yield poppable bubbles. Those are popped afterwards.
    INPUT:      DG:         our directed graph
                all_reads:  list of tuples containing all the reads and (what we need) their sequences
                work_dir:   the current work directory
                k_size:     the length of our k_mers
                delta_len:  for all differences longer than delta len we have a structural difference while all differences shorter are deemed to be errors
    OUTPUT: The simplified graph.
    """
    not_viable_global = set()
    not_viable_multibubble = set()
    has_combinations = True
    old = True
    nr_popped = 0
    overall_pops = 0
    spoa_count = 0
    multi_consensuses = {}
    merged_dict = {}
    iteration_number = 0
    this_it_pops = 0
    initial_edge_nr = len(DG.edges())
    if DEBUG:
        all_paths_s_to_t, all_path_sets = find_all_read_paths(DG, all_reads, merged_dict)
        #print("ALLPATHSINITIAL", all_paths_s_to_t)
    print("SLOW", slowmode)
    if slowmode:
        pop_threshold = 1
    else:
        pop_threshold = max(int(initial_edge_nr / 100), 1)
    # we want to continue the simplification process as long as we find combinations that have not been deemed to be "not viable" to pop
    while has_combinations:
        overall_pops += this_it_pops
        # print("Popthreshold",pop_threshold)

        this_it_pops = 0
        iteration_number += 1
        marked = set()
        direct_combis = []

        print("ITERATION NUMBER " + str(iteration_number))
        print()
        print("GRAPH NR NODES: {0} EDGES: {1} ".format(len(DG.nodes()), len(DG.edges())))
        print()
        # TopoNodes holds the topological order of the nodes in our graph
        TopoNodes = list(nx.topological_sort(DG))
        #we from here on only use topo_nodes_dict which holds the node id as key and the index of the nodeid in TopoNodes (as this improves upon the runtime)
        topo_nodes_dict = {n: i for i, n in enumerate(TopoNodes)}
        poss_starts = []
        # find all possible bubble start nodes in our graph
        find_possible_starts(DG, topo_nodes_dict, poss_starts)
        poss_ends = []
        # find all possible bubble end nodes in our graph
        find_possible_ends(DG, topo_nodes_dict, poss_ends)
        combinations = []
        # generate all combination of bubble start nodes and bubble end nodes in which the poss_starts comes before poss_end in TopoNodes
        generate_combinations(poss_starts, poss_ends, topo_nodes_dict, combinations)
        combinations_filtered = [item for item in combinations if
                                 item not in not_viable_global]  # filters the list and keeps order of elements.
        # if we haven't found any new combinations we successfully finished our bubble popping
        if not combinations_filtered:
            break

        # sort the combinations so that the shortest combinations come first
        sorted_combinations = sorted(combinations_filtered, key=lambda x: topo_nodes_dict[x[1]] - topo_nodes_dict[x[0]])
        if DEBUG:
            print("NOTVIABLE", not_viable_global)
            print("Sorted combis", sorted_combinations)
        # iterate over all combinations
        for combination in sorted_combinations:
            if DEBUG:
                print("this combination", combination)
            is_alignable = True
            all_paths = []
            ###if combination[0] not in marked and combination[1] not in marked:
            find_paths(DG, combination[0], combination[1], combination[2], all_paths, marked)
            ###else:
                ###continue
            initial_all_paths = len(all_paths)
            # if we only did find one viable path from s' to t' we figure the combination was not viable
            if len(all_paths) == 1:
                not_viable_global.add(combination)
            # we cannot touch the same node over and over during one iteration as the bubbles do not display how the graph changed
            all_paths_filtered = filter_out_if_marked(all_paths, marked, direct_combis, combination[1])
            if DEBUG:
                print("all_paths_filtered", all_paths_filtered)
            # consensus_infos stores the positions and read infos for generating the consensus
            consensus_infos = {}

            if len(all_paths_filtered) < 2:
                continue
            # if we found two paths in our bubble
            elif len(all_paths_filtered) == 2:
                if initial_all_paths > 2:
                    is_megabubble = True
                else:
                    is_megabubble = False
                this_combi_reads = tuple(sorted(set(all_paths_filtered[0][1]) | set(all_paths_filtered[1][1])))
                this_combi = (combination[0], combination[1], this_combi_reads)
                if this_combi in not_viable_multibubble:
                    continue
                # we can set the paths by accessing all_paths_filtered as we have a well-defined 2-path bubble
                path1 = all_paths_filtered[0][0]
                path2 = all_paths_filtered[1][0]
                # it can happen that one of the paths only consists of one edge. We then set pathnode to be the endnode of the bubble, the pathnode is the path defining node
                if len(path1) > 1:
                    pathnode1 = path1[1]
                else:
                    pathnode1 = combination[1]
                if len(path2) > 1:
                    pathnode2 = path2[1]
                else:
                    pathnode2 = combination[1]
                # we now would like to know whether the bubble paths have nodes in common
                p_set1 = set(path1[1:])
                p_set2 = set(path2[1:])
                intersect = p_set1.intersection(p_set2)
                # we need to find out whether the paths intersect( have a common node)
                # if they do not have an intersection we can continue
                if not intersect:
                    consensus_infos[pathnode1] = get_consensus_positions(combination[0], combination[1], DG,
                                                                         all_paths_filtered[0][1])
                    consensus_infos[pathnode2] = get_consensus_positions(combination[0], combination[1], DG,
                                                                         all_paths_filtered[1][1])
                # The paths intersect
                else:
                    if initial_all_paths == 2:
                        not_viable_global.add(combination)
                    # if we have more than two paths: we add only this combination as invalid to not_viable_multibubble
                    else:
                        this_combi_reads = tuple(sorted(set(all_paths_filtered[0][1]) | set(all_paths_filtered[1][1])))
                        this_combi = (combination[0], combination[1], this_combi_reads)
                        if DEBUG:
                            print("not_viable_multibubble add", this_combi)
                        # we only know about this combination of paths so we only set not_viable_multibubble
                        not_viable_multibubble.add(this_combi)
                    is_alignable = False
                    if DEBUG:
                        print("Not POPPABLE")
                # print("our consensus infos:",consensus_infos)
                if is_alignable:
                    is_poppable, cigar, seq_infos, consensus_info_log, spoa_count = align_bubble_nodes(all_reads,
                                                                                                       consensus_infos,
                                                                                                       work_dir, k_size,
                                                                                                       spoa_count,
                                                                                                       multi_consensuses,
                                                                                                       is_megabubble,
                                                                                                       this_combi,
                                                                                                       delta_len)

                    if is_poppable:
                        # print("Poppable")
                        linearize_bubble(DG, combination[0], combination[1], all_paths_filtered, consensus_info_log, topo_nodes_dict)
                        this_it_pops += 1
                        nr_popped += 1
                        # add all nodes that have been affected to marked
                        for node in all_paths_filtered[0][0]:
                            marked.add(node)
                        for node in all_paths_filtered[1][0]:
                            marked.add(node)
                        # if we find a directpath from s' to t'
                        if not all_paths_filtered[0][0] or not all_paths_filtered[1][0]:
                            # add the combination to direct_combis
                            combi_list = [combination[0], combination[1]]
                            if combi_list not in direct_combis:
                                direct_combis.append(combi_list)
                    else:
                        if initial_all_paths == 2:
                            not_viable_global.add(combination)
                        else:
                            this_combi_reads = tuple(
                                sorted(set(all_paths_filtered[0][1]) | set(all_paths_filtered[1][1])))
                            this_combi = (combination[0], combination[1], this_combi_reads)
                            if DEBUG:
                                print("not_viable_multibubble add", this_combi)
                            # we only know about this combination of paths so we only set not_viable_multibubble
                            not_viable_multibubble.add(this_combi)
            # we have more than two paths connecting s' and t'. We now want to efficiently compare those paths
            elif len(all_paths_filtered) > 2:
                directpath_marked = False
                if DEBUG:
                    print("more paths in", combination)
                    print("APF", all_paths_filtered)
                # initial_listing is a list holding all possible combinations of paths
                initial_listing = [(p1, p2) for (p1, p2) in itertools.combinations(all_paths_filtered, 2) if (combination[0], combination[1], tuple(sorted(set(p1[1]) | set(p2[1])))) not in not_viable_multibubble]
                if DEBUG:
                    print(initial_listing)
                # if initial_listing is empty: we do not have a viable bubble before us
                if not initial_listing:
                    # print("Not viable now bef :",not_viable_global)
                    not_viable_global.add(combination)
                    # print("Not viable now:",not_viable_global)
                    continue
                # print("APF",all_paths_filtered)
                for path_combi in initial_listing:
                    p1 = path_combi[0]
                    p2 = path_combi[1]
                    p_set1 = set(p1[0][1:])
                    p_set2 = set(p2[0][1:])
                    this_combi_reads = tuple(sorted(set(p1[1]) | set(p2[1])))
                    this_combi = (combination[0], combination[1], this_combi_reads)
                    if p_set1.intersection(p_set2):
                        not_viable_multibubble.add(this_combi)
                        continue
                    consensus_infos = {}
                    if (not p1[0]) or (not p2[0]) and directpath_marked:
                        if DEBUG:
                            print("Marked")
                        continue
                    p1_filtered = filter_path_if_marked(marked, p1[0])
                    p2_filtered = filter_path_if_marked(marked, p2[0])
                    if DEBUG:
                        print()
                    # we con only pop the bubble if both paths have not been affected by previous bubble popping steps
                    if not p1_filtered and not p2_filtered:
                        if not len(p1[0]) > 1:
                            pathnode1 = combination[1]
                        else:
                            pathnode1 = p1[0][1]
                        if not len(p2[0]) > 1:
                            pathnode2 = combination[1]
                        else:
                            pathnode2 = p2[0][1]

                        if DEBUG:
                            print("not filtered")
                        consensus_infos[pathnode1] = get_consensus_positions(combination[0], combination[1], DG, p1[1])
                        consensus_infos[pathnode2] = get_consensus_positions(combination[0], combination[1], DG, p2[1])
                        is_poppable, cigar, seq_infos, consensus_info_log, spoa_count = align_bubble_nodes(all_reads,
                                                                                                           consensus_infos,
                                                                                                           work_dir,
                                                                                                           k_size,
                                                                                                           spoa_count,
                                                                                                           multi_consensuses,
                                                                                                           True,
                                                                                                           this_combi,
                                                                                                           delta_len)
                        if DEBUG:
                            print("Do we pop?", is_poppable)
                        if is_poppable:
                            # print("POPPED_Multi")
                            if DEBUG:
                                print("POPPED_Multi")
                            all_paths_filtered = [p1, p2]
                            if DEBUG:
                                print("init", initial_all_paths)
                            # print("ALL_Paths_filtered",all_paths_filtered)
                            linearize_bubble(DG, combination[0], combination[1], all_paths_filtered, consensus_info_log, topo_nodes_dict)
                            this_it_pops += 1
                            #nr_popped += 1
                            #if (nr_popped % 10) == 0:
                                #print("NR_popped", nr_popped)
                            # add all nodes that were part of the bubble paths to marked
                            for node in all_paths_filtered[0][0]:
                                marked.add(node)
                            for node in all_paths_filtered[1][0]:
                                marked.add(node)
                            if not all_paths_filtered[0][0] or not all_paths_filtered[1][0]:
                                directpath_marked = True
                                combi_list = [combination[0], combination[1]]
                                if combi_list not in direct_combis:
                                    direct_combis.append(combi_list)
                        else:
                            # the combination is not poppable
                            not_viable_multibubble.add(this_combi)

        print("This iterations pops ", this_it_pops)
        if this_it_pops < pop_threshold:
            break
    print("Overall number of bubbles popped", overall_pops)


DEBUG = False


def simplifyGraph(DG, all_reads, work_dir, k_size, delta_len, mode):
    """Overall method used to simplify the graph
    Works as wrapper script for new_bubble_popping_routine
    INPUT:  DG: Directed Graph before simplifications,
            all_reads:Dictionary containing the sequences as well as ids for all reads that we analyze
            work_dir: The current working directory
            k_size: parameter k for our minimizers
    OUTPUT: DG: Graph after simplification took place    """
    print("Simplifying the graph")
    new_bubble_popping_routine(DG, all_reads, work_dir, k_size, delta_len, mode)
