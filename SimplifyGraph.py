import networkx as nx
from collections import Counter, namedtuple
from consensus import *
import matplotlib.pyplot as plt
from IsoformGeneration import *
#from MemoryAnalysis import *
import copy
from recordclass import recordclass
from functools import cmp_to_key

#from pyinstrument import Profiler
import itertools
from pyinstrument import Profiler
"""Helper function used to plot the graph. Taken from GraphGeneration.
    INPUT: DG   Directed Graph to plot
"""
#TODO: Add delta len again!
#Implemented according to https://www.geeksforgeeks.org/detect-cycle-in-a-graph/
def isCyclicUtil(DG, nodes_dict, node):
    # Mark current node as visited and
    # adds to recursion stack
    CNode = namedtuple('CNode', 'visited recStack')
    cnode = CNode(True, True)
    nodes_dict[node]=cnode
    #print(nodes_dict)
    #nodes_dict[node].recStack=True
     #visited[v] = True
    #recStack[v] = True

    # Recur for all neighbours
    # if any neighbour is visited and in
    # recStack then graph is cyclic
    for out_edge in DG.out_edges(node):
        neighbour=out_edge[1]
        if nodes_dict[neighbour].visited == False:
            if isCyclicUtil(DG, nodes_dict,neighbour) == True:
                #print(neighbour)
                return True
        elif nodes_dict[neighbour].recStack == True:
            #print(neighbour)
            return True

    # The node needs to be poped from
    # recursion stack before function ends
    prev_visited=nodes_dict[node].visited
    nodes_dict[node] = CNode(prev_visited,False)
    return False


"""Helper method used to detect cycles in our graph:
INPUT: DG       our directed Graph
OUTPUT: iscyclic    A boolean value indicating whether a cycle was found
"""
def isCyclic(DG):
    # Returns true if graph contains cycles else false
    nodes = DG.nodes
    nodes_dict = {}
    CNode = namedtuple('CNode', 'visited recStack',defaults=(False,False))
    cnode = CNode(False, False)
    for node in nodes:
        nodes_dict[node] = cnode
    for node in nodes:
        if nodes_dict[node].visited == False:
            if isCyclicUtil(DG,nodes_dict,node) == True:
                return True
    return False

"""
Helper method used to visualize the graph in a pdf file
INPUT:      DG          our directed Graph object
            outpath     path where we want to store the visualization
            id          id of the visualization. Wewant to be able to tell different iterations apart from each other
OUTPUT:     .pdf file holding the visualization        
        """
def draw_Graph(DG, outpath = '', id = 0):
    # defines the graph layout to use spectral_layout. Pos gives the position for each node
    pos = nx.spectral_layout(DG)
    # draws the network in the given layout. Only nodes are labeled, the edges are directed but not labeled
    nx.draw_networkx(DG, pos = pos, font_weight='bold', font_size = 2, node_size=50)
    # add labels for the edges and draw them into the graph
    # labels = nx.get_edge_attributes(DG, 'weight')
    # nx.draw_networkx_edge_labels(DG,pos, edge_labels=labels)
    plt.savefig(outpath + "graph_" + str(id) + ".pdf")
    plt.clf()
    # plt.show()





"""Helper method used to generate a subgraph for a list of nodes
INPUT:              DG          Directed Graph to plot
                    bubbles:    list of bubblenodelists to be plotted
OUTPUT:             pdf file holding the visualization of the subgraph
"""
def generate_subgraph(DG, bubble,iter):
    # for bubble in bubbles:
    SG = DG.subgraph(bubble)
    print("SG",SG)
    print(list(SG.edges))
    draw_Graph(SG,id=iter)

"""Helper method which finds all possible bubble starts in our graph. This is done by collecting all nodes having at least 2 out nodes
INPUT:  DG              Our directed graph object
        TopoNodes       A list of our nodes in topological order
OUTPUT: possible_starts:    A list holding start_tuples(node_id, start_supp:= all reads in this node)
"""
def find_possible_starts(DG,TopoNodes,possible_starts):
    #iterate over all nodes in TopoNodes
    for node in TopoNodes:
        #the current node is only a start node if it has an out_degree>1 (due to bubble definition)
        if DG.out_degree(node)>1:
            #find all reads supporting this node
            start_supp=tuple(DG.nodes[node]['reads'])
            #store node id and support in start_tup
            start_tup=(node, start_supp)
            #append start_tup to the list of possible bubble_starts
            possible_starts.append(start_tup)
"""Helper method which finds all possible bubble starts in our graph. This is done by collecting all nodes having at least 2 out nodes
INPUT:  DG              Our directed graph object
        TopoNodes       A list of our nodes in topological order
OUTPUT: possible_starts:    A list holding start_tuples(node_id, start_supp:= all reads in this node)
"""
def find_possible_ends(DG,TopoNodes,possible_ends):
    for node in TopoNodes:
        if DG.in_degree(node)>1:
            end_supp=tuple(DG.nodes[node]['reads'])
            end_tup=(node, end_supp)
            possible_ends.append(end_tup)
"""used to generate the combinations of start and end nodes
"""

def generate_combinations(possible_starts,possible_ends,TopoNodes,combis):
    for startnode in possible_starts:
        for endnode in possible_ends:
            if TopoNodes.index(startnode[0]) < TopoNodes.index(endnode[0]):
                inter=tuple(sorted(set(startnode[1]).intersection(set(endnode[1]))))
                if len(inter)>=2:
                    combi=(startnode[0],endnode[0],inter)
                    combis.append(combi)
"""we filter out combinations that already have been deemed as not_viable
"""
def filter_combinations(combinations,not_viable,combinations_filtered):
    #iterate over all combinations
    for combi in combinations:
        #if the combination is not viable add  it to combinations_filtered
        if combi not in not_viable:
            combinations_filtered.append(combi)
"""detect the paths in our bubble
"""
def find_edges_with_supp(r_id,DG):
    edges=DG.edges(data=True)
    for edge in edges:
        print(edge)
       # if edge["support"]
#TODO: replace combination by single elements-> all tuples as function attributes should be divided into single elements
def find_paths(DG,startnode,endnode,support,all_paths):
    #DEBUG=False
    thebug=False
    if thebug:
        print("FINDPATHS")
        print(DG.edges())
        print("ENDnode",endnode)
        print("Startnode",startnode)
    #print("Support",support)
    #path_and_support will hold the infos concerning the found paths
    node_support_left=set(support)
    all_supp=set(support)
    #we iterate as long as still not all support was allocated to a path
    while node_support_left:
        node = startnode
        #current_node_support = node_support_left
        read=node_support_left.pop()
        if thebug:
            print(read)
        current_node_support = node_support_left
        current_node_support.add(read)
        visited_nodes = []
        if thebug:
            print(node)
        #As long as we have not visited the bubble end node we continue walking through our graph
        while node!=endnode:
            visited_nodes.append(node)
            if thebug:
                print(node)
            out_edges=DG.out_edges(node)
            next_found=False
            for edge in out_edges:
                #if DEBUG:
                #    print("edge",edge)
                edge_supp = DG[edge[0]][edge[1]]['edge_supp']
                if read in edge_supp:
                    node=edge[1]
                    current_node_support=current_node_support.intersection(edge_supp)
                    next_found=True
                    break
            if not next_found:
                    break
        if current_node_support and next_found:
                node_support_left-=current_node_support
                final_add_support=all_supp-current_node_support
                path_supp_tup = (visited_nodes, tuple(sorted(current_node_support)), final_add_support)
                all_paths.append(path_supp_tup)
    if DEBUG:
        print("PATHSFOUND")
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
    # #print("getstartreads")
    for o_edge in out_edges:
        # #print(o_edge)
        nextnode = o_edge[1]
        if nextnode in listofnodes:
            edge_infos = DG[o_edge[0]][nextnode]['edge_supp']
            # #print("edge_infos",edge_infos)
            start_reads.update(edge_infos)
    # #print(start_reads)
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
    # #print("getstartreads")
    for i_edge in in_edges:
        if i_edge[1] in listofnodes:
            end_reads.update(DG[i_edge[0]][i_edge[1]]['edge_supp'])
    # #print(end_reads)
    return end_reads


def parse_cigar_diversity(cigar_tuples,delta_perc,delta_len):
    miss_match_length=0
    alignment_len=0
    #print("Now we are parsing....")
    #print(cigar_tuples)
    too_long_indel = False
    for i, elem in enumerate(cigar_tuples):

        cig_len = elem[0]
        cig_type = elem[1]
        alignment_len += cig_len
        if (cig_type != '=') and (cig_type != 'M'):
            #we want to add up all missmatches to compare to sequence length
            miss_match_length += cig_len
            #we also want to make sure that we do not have too large internal sequence differences
            if cig_len > delta_len:
                #if i>1 and i<(len(cigar_tuples)-1):
                    #print("ELE",elem)
                    #print("No pop due to delta_len")
                    return False
    diversity = (miss_match_length/alignment_len)


    max_bp_diff = max(delta_perc*alignment_len, delta_len)
    mod_div_rate = max_bp_diff/alignment_len
    if DEBUG:
        print("diversity",diversity,"mod_div",mod_div_rate)
    if diversity <= mod_div_rate: #delta_perc:
        return True
    else:
        #print("no pop due to diversity")
        return False
"""Helper method used to test whether two sequences are close enough to pop their underlying bubble
INPUT:      cigar_string    a string denoting the cigar output of the alignment
            delta_len       a threshold which we use to distinguish between mutations and errors
OUTPUT:     good_to_pop     boolean value indicating our finding
"""
def parse_cigar_differences(cigar_string, delta_len):
    good_to_pop = True
    #print("Now we are parsing....")
    #print(cigar_string)
    for i,elem in enumerate(cigar_string):
        cig_len = elem[0]
        cig_type = elem[1]
        # all matches are absolutely fine
        if (cig_type != '=') and (cig_type != 'M'):
            #print(cig_type)
            #print(cig_len)
            # we have a non match, now we have to figure out whether this is due to an exon or indel (everything with len>delta_len is defined as exon)
            if cig_len > delta_len:
                    # we found an exon being the difference between the paths->no bubble popping feasible
                    good_to_pop = False
                    #print(cigar_string)
    #print(good_to_pop)
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
    #print("getting dist to prev")
    #print("prev_reads",prev_node)
    #print("this_reads",curr_node)
    prev_supp=DG.nodes[prev_node]['reads']
    curr_supp=DG.nodes[curr_node]['reads']
    #print(prev_supp)
    #print("curr_supp",curr_supp)
    intersect_supp = list(set(prev_supp).intersection(curr_supp))
    sum=0
    avg_dist=0
    #we did not yet come up with anything better than calculating the average distance for all reads supporting both nodes
    for i,read in enumerate(intersect_supp):
        #print(i+1)
        prev_pos=prev_supp[read].end_mini_start
        curr_pos=curr_supp[read].end_mini_start
        sum+=(curr_pos-prev_pos)
        avg_dist=sum/(i+1)
    return avg_dist


def new_distance_to_start(DG,pnl_start,curr_node,consensus_log):
    #print("consensus_log",consensus_log)
    #print("currnode",curr_node)
    #print(str(pnl_start))
    #print("newdistancetostart")

    seq=consensus_log[pnl_start]
    node_seq=DG.nodes[curr_node]['end_mini_seq']
    #print("Is ",node_seq ," in ",seq," ?")
    possible_pos=seq.find(node_seq)
    #print(possible_pos)
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
def remove_edges(DG, bubble_start, bubble_end, path_nodes,consensus_log,edges_to_delete,node_distances):
    #print("Path nodes", path_nodes)
    ##print("PAth_reads", path_reads)
    if DEBUG:
        print("Bubble_start", bubble_start)
        print("Bubble_end",bubble_end)
    # we have to delete the edges which connect the nodes that are inside the bubble
    for one_info in path_nodes:
        #print("path_nodes", path_nodes)
        path_node_list=one_info[0]
        if path_node_list:
            path_node_list.pop(0)
        #startnode_id=path_node_list[0]
        if not(path_node_list):
            pnl_start=bubble_end
        else:
            pnl_start=path_node_list[0]
        prev_node = bubble_start
        #print("PathNodeList", path_node_list)
        entry = DG.get_edge_data(prev_node, pnl_start)
        edges_to_delete[prev_node,pnl_start]=entry
        if not(path_node_list):
            pnl_start=bubble_end
        else:
            pnl_start=path_node_list[0]
        for index, path_node in enumerate(path_node_list):

            #print("TNI",this_node_infos)
            inter_dist = new_distance_to_start(DG, pnl_start, path_node, consensus_log)
            #print("node_inter:",path_node,",",inter_dist)
            if inter_dist == -1:
                dist_to_prev = get_dist_to_prev(DG, prev_node, path_node)
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
    ##print("Edges To Delete", edges_to_delete)
    for edge, edge_infos in edges_to_delete.items():
        #print("deletingEdge",edge)
        DG.remove_edge(edge[0], edge[1])
    #print("node_distances ", node_distances)
    return edges_to_delete, node_distances



def get_avg_interval_length(DG,node):
    #print("calculating average length for node ", node)
    node_support=DG.nodes[node]['reads']
    sum = 0
    i = 0
    for r_id, positions in node_support.items():
        #print(positions)
        i += 1
        sum += (positions.end_mini_start-positions.start_mini_end)
        #print(sum)
    if i==0:
        finalresult=0
    else:
        finalresult=int(sum/(i))
    #finalresult=100000
    #print("finalresult",finalresult)
    return finalresult



"""Helper method: Adds additional support values for the bubble_nodes(needed to get a consistent graph)
INPUT:      DG:
            nextnode:       the node of which we need the support
OUTPUT:     additional_support:  dictionary of additional support values
"""

def additional_node_support(DG, new_support, node_dist_dict, node,other_prevnode,bubble_start,additional_support):
    Read_infos = namedtuple('Read_Infos',
                            'start_mini_end end_mini_start original_support')
    if other_prevnode==bubble_start:
        other_dist=0
    else:
        other_dist=node_dist_dict[other_prevnode]
    this_dist=node_dist_dict[node]
    #print("other:dist",other_dist)
    #this will contain all reads that are to be added to the node with their respective (virtual) positions
    this_reads = DG.nodes[node]['reads']
    #print("Node_dist_dict", node_dist_dict)
    #print("node",node," thisreads ",this_reads)
    #print("bubble_Start",bubble_start)
    #print("s_infos",s_infos)
    #this_reads=DG.nodes[node]['reads']
    avg_len=get_avg_interval_length(DG,node)
    #a lot of print statements used to debug
    #print("avg_len",avg_len)
    #print("new_support ",new_support)
    #print("prevnode_this_path", prevnode_this_path)
    #print("prevnode_other_path",other_prevnode)
    #start_pos= DG.nodes[bubble_start]['reads']
    for r_id in new_support:
        #print("r_id",r_id)
        #figure out whether the r_id is part of this bubble_path
        if r_id not in this_reads:
            # what we do if the read meets the bubble and is not in this_reads
            # figure out if r_id is also represented by bubble_start
            #if r_id not in s_infos.keys():
                # r_id is not in s_infos (meaning the read is not supporting bubble_start), if it would be in s_infos, we do not have to recompute anything

                #get the position of the read in previous node
                #calculate relative distance to start (relative to previous node)
                #start=position+relative dist
                #end=start+avg_length
                #print("r_id in other path")
                #print("other prevnode",other_prevnode)
                previous_other_path_reads=DG.nodes[other_prevnode]['reads']
                #print("previous_other_path_reads",previous_other_path_reads)
                if r_id in previous_other_path_reads:
                    pos_info_tup=previous_other_path_reads[r_id]
                    prev_end=pos_info_tup.end_mini_start
                    relative_dist= int(this_dist-other_dist)
                    newend=prev_end+relative_dist
                    newstart=int(newend-avg_len)
                    """"#The calculations with relative_dist may be too inaccurate. We want to improve the results a bit
                    #we try to find the read id in the reads of bubble_start
                    if r_id in start_pos:
                        #if r_id is in start_pos, we take its positions there and correct the number for our read
                        this_read_start_infos = start_pos[r_id]
                        this_read_start_bubble = this_read_start_infos.start_mini_end
                        this_read_end_bubble=this_read_start_infos.end_mini_start
                        if newstart < this_read_start_bubble:
                            newstart=this_read_start_bubble
                            if newend< this_read_end_bubble:
                                newend=this_read_end_bubble
                    #we did not find the read_id in our bubble_start-> The only thing we can do is to at least make sure we do not get negative positions
                    else:
                        if newstart<0:
                            newstart=1
                        if newend<0:
                            newend=1"""
                    additional_support[r_id] = Read_infos(newstart, newend, False)

    #print("Additional node_supp after", additional_support)


"""Helper method used to merge two dicts with the values being lists 
INPUT: dict1: the first dict
       dict2: the second dict
OUTPUT: merged_dict: dictionary containing all keys in dict1 and dict2
"""
def merge_two_dicts(dict1, dict2):
    merged_dict = {}
    ##print("dict1", dict1)
    ##print("dict2", dict2)
    for key, value in dict1.items():
        merged_dict[key] = value
    for key2, val2 in dict2.items():
        if key2 not in merged_dict:
            merged_dict[key2] = val2
    return merged_dict

def sort_by_node_label(seq):
    seq_duplicates = sorted(seq, key=lambda x: x[1])
    return seq_duplicates
def find_next_node_of_this_path(actual_node_distances,path1):
    #print("path1",path1)
    #print("actual_node_distances(AND)", actual_node_distances)
    for node in actual_node_distances:
        #print("Node_which_we_found",node)
        if node[0] in path1:
            #print("found ",node, "in both ")
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
def compare_by_length(nextnode1,nextnode2,node_dist,bubble_end, topo_nodes_dict):
    if nextnode1 == bubble_end:
        return nextnode2
    elif nextnode2 == bubble_end:
        return nextnode1
    return nextnode1 if topo_nodes_dict[nextnode1] < topo_nodes_dict[nextnode2] else nextnode2

    dist1 = node_dist[nextnode1]
    dist2 = node_dist[nextnode2]
    #print("nextnode1:",nextnode1,", dist ",dist1,", nextnode2",nextnode2,"dist",dist2)
    if dist1 < dist2-5:
        return nextnode1
    elif dist2 < dist1-5:
        return nextnode2
    else:
        # print(topo_nodes_dict)
        # print(nextnode1, nextnode2)
        # print("We have to come up with something!")
        return nextnode1 if topo_nodes_dict[nextnode1] < topo_nodes_dict[nextnode2] else nextnode2
        #print("is ", nextnode1, " < ", nextnode2, "? As we got ", dist1, "== ", dist2)
    return nextnode1
    # #print("FinalEdges",DG.edges(data=True))


def find_real_nextnode(nextnode1,nextnode2,node_dist,bubble_end,conn_edges,DG, topo_nodes_dict):
    if DG.has_edge(nextnode1,nextnode2):
        return nextnode1
    elif DG.has_edge(nextnode2,nextnode1):
        return nextnode2
    else:
        nextnode = compare_by_length(nextnode1, nextnode2, node_dist, bubble_end, topo_nodes_dict)
        return nextnode


def get_next_node(path, bubble_end):
    #print("get_ndex_node")
    #print(str(type(path)))
    #print(path)
    if len(path)<2:
        nextnode=bubble_end
    else:
        nextnode=path[1]
    return nextnode

"""
Helper method that finds edges that connect nodes between different paths
INPUT:



OUTPUT:
"""
def find_connecting_edges(path_nodes,DG):
    connecting_edges=set()
    path1=path_nodes[0][0]
    path2=path_nodes[1][0]
    for node in path1:
        for onode in path2:
            if DG.has_edge(node,onode):
                connecting_edges.add((node,onode))
            elif DG.has_edge(onode,node):
                connecting_edges.add((onode,node))
    #if connecting_edges:
        #print("connecting_edges",connecting_edges)
    return connecting_edges
def test_conn_end(conn_edges,overall_nextnode):
    conn_list=[]
    #print("connedges",conn_edges)
    for conn_edge in conn_edges:
        #print(conn_edge[1])
        if overall_nextnode == conn_edge[1]:
            conn_list.append(conn_edge)
    return conn_list
#def find_overlap_reads(path1,path2,p1_add_infos,p2_add_infos):
def find_correct_supp(p1_add_infos,p2_add_infos,path1_all_supp,path2_all_supp):
    #as we take the support from the paths including the startnode, we have to clean the supports to only include the reads that are not path_reads for any of the two paths.
    #we first get rid of all path reads for p1
    p2_clean=path2_all_supp-p1_add_infos
    #now we clean out the path reads of p2
    p1_clean=path1_all_supp-p2_add_infos
    #the last and most important step is to get the intersection of reads that are shared between the paths
    path_supp_inter=p1_clean.intersection(p2_clean)
    return path_supp_inter



def find_intersect_nodes(DG,path,intersect_supp):
    p_info_dict = {}
    #print("FIN")
    for node in path:
        # print("node",node)
        node_supp=DG.nodes[node]['reads']
        keys=node_supp.keys()
        #print("keys",keys)
        intersect_reads=intersect_supp.intersection(keys)
        # print("iSupp",intersect_supp)
        # print("IReads",intersect_reads)
        for r in intersect_reads:
            # print("P_info",p_info_dict)
            r_infos=node_supp[r]
            # print("What is ",r," ",r_infos)
            info_tuple = (node, r_infos.end_mini_start)
            if r in p_info_dict:
                already_infos=p_info_dict[r]
                # print("A_infos",already_infos)
                already_infos.append(info_tuple)
                p_info_dict[r]=already_infos
            else:
                empty=[]
                empty.append(info_tuple)
            p_info_dict[r]=empty
    return p_info_dict
def find_intersect_supp_nodes(DG,intersect_supp,path1,path2):
    p1_info_dict=find_intersect_nodes(DG,path1,intersect_supp)
    p2_info_dict=find_intersect_nodes(DG,path2,intersect_supp)
    # print("P1Info",p1_info_dict)
    # print("P2Info", p2_info_dict)
    key_dict={}
    if len(p1_info_dict)>0 and len(p2_info_dict)>0:
        for key in p1_info_dict:
            if key in p2_info_dict:
                p_infos=p1_info_dict[key]+p2_info_dict[key]
                # print("P_infos",p_infos)
                p_infos.sort(key=lambda x:x[1])
                # print("P_infos_s",p_infos)
                key_dict[key]=p_infos
    # print("KD",key_dict)


#TODO: we seem to loose edge infos-> this should happen here
def prepare_adding_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes,node_dist,seq_infos, topo_nodes_dict):  # ,node_dist):
    counter = 0
    path1 = []
    path2 = []
    edge_params={}
    linearization_order=[]
    conn_edges=find_connecting_edges(path_nodes,DG)
    # we assign both paths to variables to make them easier accessible.
    for path_info in path_nodes:
        # print("PI",path_info)
        nodes=path_info[0]
        add_infos=path_info[1]

        ##print("id", nodes[0])
        if counter == 0:
            for p in nodes:
                path1.append(p)

        else:
            for p in nodes:
                path2.append(p)
        counter = counter + 1
    #this dictionary contains the reads with their respected positions that are added to the nodes
    new_node_supp_dict = {}
    #print("path1",path1)
    #print("path2",path2)
    linearization_order.append(bubble_start)
    prevnode1 = bubble_start
    prevnode2 = bubble_start

    if path1:
        nextnode1 = path1[0]
        #print("nextnode1",nextnode1)

    else:
        nextnode1=bubble_end
        ##print("Creepy ands", actual_node_distances)
    if path2:
        nextnode2 = path2[0]
        #print("nextnode2", nextnode2)
    else:
        nextnode2 = bubble_end
        ##print("Creepy ands", actual_node_distances)
    prevnode = bubble_start
    #print("analyzing bubble from", bubble_start, "to", bubble_end)
    while path1 or path2:
        #print("P1",path1)
        #print("P2", path2)

        overall_nextnode=find_real_nextnode(nextnode1,nextnode2,node_dist,bubble_end,conn_edges,DG, topo_nodes_dict)
        #print("overall_nextnode",overall_nextnode)
        conn_list=test_conn_end(conn_edges,overall_nextnode)
        #print("CONNEND?",conn_list)
        new_edge_supp1 = edges_to_delete[prevnode1, nextnode1]['edge_supp']
        #print("NES1 from ",prevnode1," to ",nextnode1,": ",new_edge_supp1)
        new_edge_supp2 = edges_to_delete[prevnode2, nextnode2]['edge_supp']
        #print("NES1 from ", prevnode2, " to ", nextnode2, ": ", new_edge_supp2)

        if not conn_list:
            full_edge_supp = new_edge_supp1 + new_edge_supp2
            #print("FES: ", full_edge_supp)
        else:
            #print("Wanting to add edge support from", prevnode, " to ", overall_nextnode)
            this_edge_supp=DG[conn_list[0][0]][overall_nextnode]["edge_supp"]
            full_edge_supp=new_edge_supp1+new_edge_supp2+this_edge_supp

        if overall_nextnode in path1:
            nextnode1=get_next_node(path1,bubble_end)
            #print("nextnode1",nextnode1)
            additional_supp={}
            #TOD: we need to add the support of global_prev for both additional_node_support occurrences
            additional_node_support(DG, new_edge_supp2, node_dist, overall_nextnode,
                                                       prevnode2, bubble_start,additional_supp)
            #print("additionalNodeSupport1 for",overall_nextnode,":",additional_supp)
            #print("path 1 before remove", path1)
            path1.remove(overall_nextnode)
            #print("path 1 after remove",path1)
            prevnode1=overall_nextnode
            curr_node=prevnode1
        elif overall_nextnode in path2:
            nextnode2=get_next_node(path2,bubble_end)
            additional_supp={}
            additional_node_support(DG, new_edge_supp1, node_dist, overall_nextnode,
                                                       prevnode1, bubble_start,additional_supp)
            #print("additionalNodeSupport2 for", overall_nextnode,":",additional_supp)
            #print("path 2 before remove", path2)
            path2.remove(overall_nextnode)
            #print("path 2 after remove", path2)
            prevnode2=overall_nextnode
            curr_node = prevnode2

            #print("ERRORRRRRR", overall_nextnode," neither in ",nextnode1," nor in ",nextnode2)
        old_node_supp = DG.nodes[overall_nextnode]['reads']
        new_node_supp_dict[overall_nextnode] = merge_two_dicts(additional_supp, old_node_supp)
        linearization_order.append(overall_nextnode)
        nx.set_node_attributes(DG, new_node_supp_dict, "reads")
        edge_params[prevnode,overall_nextnode]=full_edge_supp
        #DG.add_edge(prevnode, overall_nextnode, edge_supp=full_edge_supp)

        ##print("Adding edge from ", prevnode, "to ", overall_nextnode)
        prevnode=curr_node

    new_edge_supp1 = edges_to_delete[prevnode1, bubble_end]['edge_supp']
    new_edge_supp2 = edges_to_delete[prevnode2, bubble_end]['edge_supp']
    full_edge_supp = new_edge_supp1 + new_edge_supp2
    edge_params[prevnode, bubble_end] = full_edge_supp
    linearization_order.append(bubble_end)
    #print("edge_params",edge_params)
    #print("Linearized nodes to be in the following order ",linearization_order)
    #DG.add_edge(prevnode, bubble_end, edge_supp=full_edge_supp)
    #print("Adding edges")
    #adding the new edges to the Graph
    for key,value in edge_params.items():
        #if the edge has not been in the graph before
        if not DG.has_edge(key[0],key[1]):
            DG.add_edge(key[0],key[1],edge_supp=value)
        else:
            #The edge was in the graph before, add old support to keep the graph consistent
            #print("Before:",value)
            #old_edge_supp=DG.edge[key[0]][key[1]]['edge_supp']
            old_edge_supp=DG.edges[key[0], key[1]]['edge_supp']
            #all_edges=DG.edges.data()
            #old_edge_supp=all_edges[key[0]][key[1]]['edge_supp']
            for old_supp in old_edge_supp:
                if not old_supp in value:
                    value.append(old_supp)
            #print("After:",value)
            DG.add_edge(key[0],key[1],edge_supp=value)
            #TODO: make this work


        #print("NEW edges and support for ",key,":",value)
        #print("Adding edge fa from ", key[0], "to ", key[1])
    # this is the main part of the linearization. We iterate over all_nodes and try to find out which path the nodes belong to.
    # This info is needed as we need the current state ob both paths to add the correct edge_support and node_support to the graph
"""Helper method used to generate the consensus sequence for each bubble path
INPUT:      work_dir:               The work directory we are in
            consensus_attributes:   list of tuples containing the read_id and the positions that make up the limits of our subsequences
            reads:                  dictionary containing the actual reads (the sequences as well as their original ids)
            k_size:                 the length of parameter k needed to extract the correct subsequence length
OUTPUT:     spoa_ref:               the consensus sequence generated from the reads

"""
def generate_consensus_path(work_dir, consensus_attributes, reads, k_size,spoa_count):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    seq_infos={}
    endseqlist = []
    reads_path_len=0
    max_len=0
    longest_seq_len=-1
    if len(consensus_attributes)>2:
        #print("More than 2")
        for i, (q_id, pos1, pos2) in enumerate(consensus_attributes):
            #print("consensus_atm:",q_id,", ",pos1,",",pos2)
        #print("read ",q_id)
        #print(pos2)
        #print("Printing full seq:", reads[q_id][1])
            if pos2 == 0:
                #print("TRUE")
                pos2 = len(reads[q_id][1]) - k_size
        #print(pos2)
            if pos1<(pos2+k_size):
                if pos1<0:
                    pos1=0
                seq = reads[q_id][1][pos1: (pos2 + k_size)]
            else:
                seq=""
        #print(seq)
            seq_infos[q_id]=(pos1,pos2+k_size,seq)
            #print("seq_infos",seq_infos)
            endseq = reads[q_id][1][pos2:pos2 + k_size]
            endseqlist.append(endseq)
        #print(q_id, "from ", pos1, "to", pos2 + k_size, ": ", seq)
            if len(seq)<k_size:
                if len(seq)>max_len:
                    max_len=len(seq)
                elif len(seq)<1:
                #print("Seq too short: ",q_id,", ",pos1,",",pos2)
                    max_len=1
            #print("not popping ",q_id)
            else:
                reads_path_len+=1
                reads_path.write(">{0}\n{1}\n".format(str(q_id) + str(pos1) + str(pos2), seq))
        reads_path.close()
    #print("seq_infos",seq_infos)
    #print("RPL",reads_path_len)
        if reads_path_len>0:
            spoa_count+=1
            if (spoa_count%100)==0:
                print("Spoa_count",spoa_count)
            spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir, "spoa_tmp.fa"), "spoa")
        #print("spoa_ref", spoa_ref)
            return spoa_ref,seq_infos,spoa_count
        else:
            string_val = "X" * max_len  # gives you "xxxxxxxxxx"
            return string_val,seq_infos,spoa_count
    else:
        #print("less than two seqs")
        f_id,fstart,fend=consensus_attributes[0]
        e_id,estart,eend=consensus_attributes[1]
        fdist=fend-fstart
        edist=eend-estart
        if fdist>edist:
            consensus=reads[f_id][1][fstart: (fend + k_size)]
            seq_infos[f_id] = (fstart, fend + k_size, consensus)
        else:
            consensus = reads[e_id][1][estart: (eend + k_size)]
            seq_infos[f_id] = (estart, eend + k_size, consensus)
        return consensus, seq_infos, spoa_count
""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
           delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """
def collect_consensus_reads(consensus_attributes):
    con_reads=set()
    for con_att in consensus_attributes:
        con_reads.add(con_att[0])
    con_reads_fin=frozenset(con_reads)
    return con_reads_fin


#TODO: Either we have to change the calculation of the positions or we have to think more about pop threshold
def align_bubble_nodes(all_reads, consensus_infos, work_dir, k_size,spoa_count,multi_consensuses,is_megabubble,combination,delta_len):
    if DEBUG:
        print("aligning")
        print("current consensus_infos",consensus_infos)
    consensus_list = []
    consensus_log= {}
    seq_infos={}
    for path_node, consensus_attributes in consensus_infos.items():
        #print("consensus", consensus_attributes)
        con_reads=collect_consensus_reads(consensus_attributes)
        combi = (combination[0], combination[1], con_reads)
        if combi in multi_consensuses:
            contains_combi=True
        else:
            contains_combi=False
        if is_megabubble and contains_combi:
            con_infos=multi_consensuses[combi]
            con = con_infos[0]
            spoa_count=con_infos[2]
            consensus_log[path_node] = con
            consensus_list.append(con)
        else:
        #if True:
            if len(consensus_attributes) > 1:
                #print("more in this consensus", consensus_attributes)
                con,seq_infos_from_fun,spoa_count = generate_consensus_path(work_dir, consensus_attributes, all_reads, k_size,spoa_count)
                #print(con,seq_infos_from_fun)
                if is_megabubble:
                    multi_consensuses[combi]=(con,seq_infos_from_fun,spoa_count)
                if len(con)<3:
                    consensus_log[path_node] = ""
                    consensus_list.append("")
                else:
                    consensus_log[path_node]=con
                    seq_infos.update(seq_infos_from_fun)
                #print(con)
                    consensus_list.append(con)
            else:
                #print("one elem consensus", consensus_attributes)
                (q_id, pos1, pos2) = consensus_attributes[0]
                if abs(pos2-pos1)<3:
                    consensus_log[path_node] = ""
                    consensus_list.append("")
                elif pos2<pos1:
                    consensus_log[path_node] = ""
                    consensus_list.append("")

                else:
                #print("consensus_attributes", q_id, ", ", pos1, ", ", pos2)
                    con = all_reads[q_id][1][pos1: pos2 + k_size]
                    seq_infos[q_id]=(pos1,pos2+k_size,con)
                #print("single consensus",con)
                    consensus_log[path_node]=con
                    consensus_list.append(con)
    #print("CITEMS",consensus_infos)
    #print(consensus_list)
    #print("consensus_log",consensus_log)
    consensus1 = consensus_list[0]

    consensus2 = consensus_list[1]
    if DEBUG:
        print("consensus1",consensus1)
        print("consensus2",consensus2)
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
        if DEBUG:
            print("long:", longer_len, " , short ", shorter_len, "ratio ", (longer_len - shorter_len) / longer_len)
        if (longer_len - shorter_len) / longer_len < delta:
            #print("popping works")
            s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(consensus1, consensus2,
                                                                                       match_score=2,
                                                                                       mismatch_penalty=-8,#standard: -8
                                                                                       opening_penalty=12, gap_ext=1) #opening penalty: standard: 12
            good_to_pop = parse_cigar_diversity(cigar_tuples, delta, delta_len)
            if DEBUG:
                print("GOODTOPOP?",good_to_pop)
            cigar_alignment = (s1_alignment, s2_alignment)
            if good_to_pop and DEBUG:
                print("deemed to be poppable", cigar_alignment," ",combination,"$$$ ",consensus_infos,"€€",cigar_tuples)
            #else:
                #print("No pop", cigar_alignment, " ", consensus_infos,"€€",cigar_tuples)
            return good_to_pop, cigar_alignment, seq_infos, consensus_log,spoa_count
        else:
            if DEBUG:
                #print("deemed to be poppable", cigar_alignment, "")
                print("No pop-too diverse ", combination)
            return False,"","", consensus_log,spoa_count
    elif shorter_len < delta_len and longer_len < delta_len:
        if DEBUG:
            print("deemed Poppable as short ",consensus_infos)
        return True, "","",consensus_log,spoa_count

    else:
        if (longer_len - shorter_len) < delta_len:
            if DEBUG:
                print("Poppable as short and not too different", consensus_infos)
            return True, "", "", consensus_log,spoa_count
        else:
            if DEBUG:
                print("Not poppable as short but too different", consensus_infos)
            return False, "", "", consensus_log,spoa_count
        #good_to_pop=False

    #print(s1_alignment)
    #print(s2_alignment)
    #print(cigar_string)
    #print(cigar_tuples)
    #print("cigar done")



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
        #print(out_node)
        inter_out_readlist = DG.nodes[out_node]['reads']
        #print(inter_out_readlist)
        #out_readlist = [ for i in inter_out_readlist]
        out_path_reads_list = []
        for read in inter_out_readlist:
            if read in shared_reads:
                out_path_reads_list.append(read)
        path_starts[out_node] = out_path_reads_list
    #print("Pathstarts")
    #print(path_starts)

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
        #print("This is something!")
        #print("visited_nodes", visited_nodes)
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
                ##print("intersect supp",intersect_supp)
                # if we have an overlap, this means that at least one read from our initial_support also supports this edge
                if intersect_supp:
                    # we have found a node to continue
                    further_node = True
                    # set the node to be our current_node
                    curr_node = next_node
                    # update the set of still supporting reads
                    curr_support = intersect_supp
                    ##print("curr_support")
                    # as we only expect one node to be viable: break the for loop to continue in the next node
                    break
        # we did not find any node we could continue with->our path is not the path of a real bubble
        if not further_node:
            is_viable = False
            return is_viable, visited_nodes, intersect_supp
    # we were able to reach bubble_end->we just confirmed that we have a bubble path
    is_viable = True
    #print("Visited_nodes ", visited_nodes)
    ##print("intersect_supp", intersect_supp)
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
    #print("minnodeout",min_node_out)
    # Now we want to get the actual reads for each path
    path_starts = get_path_reads(DG, min_node_out, shared_reads)
    # iterate over all shared reads and get the pathlength for each
    #print("Path_starts ", path_starts)
    return path_starts


"""Helper method for find_bubbles: This method figures out which reads belong to which part from bubble source to bubble sink
INPUT:      cycle:          A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element:    The node deemed to be the starting node of the bubble (given as tuple)
            max_element:    The node deemed to be the end node of the bubble (given as tuple)
            DG:             the directed graph we want to pop the bubble in
            contains_s:     A boolean value denoting whether the cycle contains node "s"
            contains_t:     A boolean value denoting whether the cycle contains node "t"
"""


def get_consensus_positions( bubble_start, bubble_end, DG, shared_reads):
    read_list = []
    max_node_infos = DG.nodes[bubble_end]['reads']
    min_node_infos = DG.nodes[bubble_start]['reads']
    #print("finding consensus_infos from ",bubble_start,"to ",bubble_end,", shared_reads ", str(shared_reads))

    #print(max_node_infos)
    #print(min_node_infos)
    for r_id in shared_reads:
            bubble_end_pos = max_node_infos[r_id]
            bubble_start_pos = min_node_infos[r_id]
            start_of_bubble = bubble_start_pos.end_mini_start
            end_of_bubble = bubble_end_pos.end_mini_start
            entry = (r_id, start_of_bubble, end_of_bubble)
            #if not entry[1]>entry[2]:
            read_list.append(entry)
    #print("R_list",read_list)
    return read_list


"""
function which finds the bubbles and if feasible pops them in the graph
INPUT:
    DG:         Our directed graph 
    delta_len   Maximum length difference for which the paths are merged ie the bubble is popped
    all_reads   dictionary containing all the reads(string sequence) and their ids


OUTPUT:
    bubble_state, popable_bubbles,no_pop_tup_list"""

"""Helper method that filters the poppable_bubbles to only contain non-overlapping bubbles, as overlapping bubbles yield inconsistencies in our graph
INPUT:      poppable_bubbles:      a list of bubbles that are poppable
OUTPUT:     new_poppable:          a list of bubbles that do not overlap->ready for popping
"""
def find_disjoint_bubbles(poppable_bubbles):
    ##print("Filter touched", poppable_bubbles)
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


def linearize_bubble(DG, pre_consensus_infos, bubble_start, bubble_end, path_nodes, support_dict,seq_infos,consensus_log, topo_nodes_dict):
    #print("seq_infos",seq_infos)
    #print("time for linearization")
    edges_to_delete={}
    node_dist={}
    # edges_to_delete: A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)

    edges_to_delete,node_dist = remove_edges(DG, bubble_start, bubble_end, path_nodes,consensus_log,edges_to_delete,node_dist)
    #print("ETD", edges_to_delete)
    #print("bubbles",bubble_start, ", ",bubble_end)
    ##print("Allnodes after sorting", all_nodes)
    #add_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes,node_dist,seq_infos)  # , node_dist)
    #print("PATH nodes", path_nodes)
    prepare_adding_edges(DG, edges_to_delete, bubble_start, bubble_end, path_nodes, node_dist, seq_infos, topo_nodes_dict)
    #print("Popped bubble ", path_nodes)



"""Helper method used to extend the list of bubbles that has been analysed but not been deemed poppable (old_bubbles)
INPUT:      old_bubbles
            bubble_state
"""
def extend_unpoppable(old_bubbles,poppable_state):
    for newbubble in poppable_state:
        #print(newbubble)
        #print(newbubble.bubble_supported)
        if not newbubble.bubble_supported:
            #print("newbubble",newbubble.bubble_nodes)
            old_bubbles.append(newbubble.bubble_nodes)
def filter_marked_paths(all_paths,marked):
    is_marked=[]
    for path_tup in all_paths:
        path=path_tup[0]
        already_marked=False
        for node in path:
            if node in marked:
                already_marked=True
        marked_tuple=(path,already_marked)
        is_marked.append(marked_tuple)
    for tup in is_marked:
        if tup[1]==True:
            all_paths.remove(tup)
def filter_out_if_marked(all_paths,marked,direct_combis,endnode):
    filter_list=[]
    for path in all_paths:
        path_nodes=path[0]
        for node in path_nodes:
            if node in marked:
                if not path in filter_list:
                    filter_list.append(path)
                    break
        for direct_combi in direct_combis:
            #print("DC",direct_combi)
            startnode = direct_combi[0]
            for id,node in enumerate(path_nodes):
                #make sure the path_nodes list is long enough
                if len(path_nodes)>id+1:
                    if node==startnode and path_nodes[id+1]==direct_combi[1]:
                        #we add the path to filter_list, as we have already seen the combi as direct combination
                        if not path in filter_list:
                            filter_list.append(path)
                        break
                else:
                    if node==startnode and endnode==direct_combi[1]:
                        if not path in filter_list:
                            filter_list.append(path)
                        break
    #remove all paths that we want to filter out from the initial list of paths
    if filter_list:
        for entry in filter_list:
            all_paths.remove(entry)
    return all_paths

def get_path_length(path_infos,DG,poss_start,poss_end):
    reads=path_infos[1]
    path_len_sum=0
    start_dict = DG.nodes[poss_start]['reads']
    end_dict = DG.nodes[poss_end]['reads']
    #print(start_dict)
    #print(end_dict)
    #print("number of reads",len(reads))
    for read in reads:
        start_pos=start_dict[read].end_mini_start
        end_pos=end_dict[read].end_mini_start
        #print("r_id",read)
        path_len=end_pos-start_pos
        path_len_sum +=path_len
    avg_path_len=path_len_sum/(len(reads))
    return avg_path_len
def filter_path_if_marked(marked, path):
    #print("filter path if marked for ",marked," and ",path)
    for node in path:
        if node in marked:
            return True
    return False

def subgraph_bubble(DG,itere,all_paths_filtered,bubble_start,bubble_end,error_bubble):
    lista=[]

    for elem in all_paths_filtered:
        lista.append(elem[0])
    flat_list = [item for sublist in lista for item in sublist]
    flat_list.append(bubble_end)
    if not bubble_start in flat_list:
        flat_list.append(bubble_start)
    for error_node in error_bubble:
        if error_node in flat_list:
            #print("APF", all_paths_filtered)
            flat_list.extend(error_bubble)
            generate_subgraph(DG, flat_list, itere)
            break
    #print("FlatList", flat_list)
def find_combi_paths(combination,all_paths):
    #TODO change this function to find the paths for a combination and return the correct value
    Readtup=recordclass('Readtup','path supp non_supp')
    #all_r_ids=set(all_paths.keys())
    startnode=combination[0]
    endnode=combination[1]
    inter=combination[2]
    path_list=[]
    r_ids=[]
    #print("START",startnode, "END",endnode)
    for r_id in inter:
        nodelist=all_paths[r_id]
        #nodelist=list(inter[1] for inter in inter_list)

        start_idx=nodelist.index(startnode)
        end_idx=nodelist.index(endnode)
        path=nodelist[start_idx:end_idx]
        if DEBUG:
            print("PATH",path)
            print("NODELIST", nodelist)
        supp_rid_list=set()
        supp_rid_list.add(r_id)
        diff_set=set(inter).difference(supp_rid_list)
        path_supp_tup = Readtup(path, supp_rid_list, diff_set)
        path_list.append(path_supp_tup)
        r_ids.append(r_id)
    #all_rids=set(r_ids)
    if DEBUG:
        print("PATHLIST",path_list)
    already_merged=[]
    merge_dict={}
    for id,rtup in enumerate(path_list):
        for id2,rtup2 in enumerate(path_list):
            if not rtup2 in already_merged:
                if id<id2:
                    if DEBUG:
                        print(id,", ",id2)
                    if rtup.path==rtup2.path:

                            already_merged.append(rtup2)

                            new_supp=rtup.supp.union(rtup2.supp)
                            """if DEBUG:
                                print("P1", rtup.path)
                                print("P2", rtup2.path)
                                print("Type", str(type(rtup.path)))
                                print("Type", str(type(rtup2.path)))
                                print(new_supp)"""
                            other_supp=set(inter).difference(new_supp)
                            rtup.supp=new_supp
                            rtup.non_supp=other_supp
                            #merge_dict[id]=new_tuple
                            #print(path_list)
    if DEBUG:
        print("Merged_already",already_merged)
        print("MERGEDICT",merge_dict)
    for merged_elem in already_merged:
        path_list.remove(merged_elem)
    new_path_list=[]
    for thispath in path_list:
        new_path=(thispath.path,tuple(thispath.supp),thispath.non_supp)
        new_path_list.append(new_path)
    return path_list
def find_path(r_id,DG,edge_attr):
    current_node="s"
    visited_nodes=[]
    visited_edges=[]

    reached_t=False
    while (not reached_t):
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
    return visited_nodes

def find_all_read_paths(DG,all_reads,merged_dict):
    all_read_paths={}
    all_path_sets={}
    edge_attr = nx.get_edge_attributes(DG, "edge_supp")
    for r_id in all_reads.keys():
        if not r_id in merged_dict:
            print(r_id)
            all_read_paths[r_id]=find_path(r_id,DG,edge_attr)
        else:
            all_read_paths[r_id]=all_read_paths[merged_dict[r_id]]
    for r_id,path in all_read_paths.items():
        all_path_sets[r_id]=set(path)
    return all_read_paths, all_path_sets

def update_merged_reads_dict(merged_dict,all_read_paths):
    print("Hello World")
def generate_equal_reads_dict(all_read_paths):
    equal_reads_dict={}
    for id, path in all_read_paths.items():
        pa1=set(path)
        #print(id)

        for id2, path2 in all_read_paths.items():
            if id<id2:
                if DEBUG:
                    print(id, ":", id2)
                pa2=set(path2)
                if pa1==pa2:
                    equal_reads_dict[id2]=id
    return equal_reads_dict
#TODO: add all_read_sets
def update_paths(DG,all_reads,prev_marked,merged_dict,all_paths_s_to_t,all_path_sets):
    new_all_paths_s_to_t={}
    edge_attr = nx.get_edge_attributes(DG, "edge_supp")
    for r_id in all_reads.keys():
        #if we do not find a node in the readpath that has been marked as changed
        if set(all_path_sets[r_id]).isdisjoint(prev_marked):
            #copy the old path
            new_all_paths_s_to_t[r_id]=all_paths_s_to_t[r_id]
        #we found that the readpath has been changed->Regenerate the readpath
        else:
            if r_id in merged_dict.keys():
                new_all_paths_s_to_t[r_id]=new_all_paths_s_to_t[merged_dict[r_id]]
            else:
                new_all_paths_s_to_t[r_id]=find_path(r_id,DG,edge_attr)
            all_path_sets[r_id].update(new_all_paths_s_to_t[r_id])
    return new_all_paths_s_to_t, all_path_sets
    print("update_read_pahts")
    #TODO write this method to recalculate the paths of reads that were affected by bubble popping
"""The backbone of bubble popping: the bubbles are detected and filtered to only yield poppable bubbles. Those are popped afterwards.
INPUT:      DG:         our directed graph
            delta_len:  for all differences longer than delta len we have a structural difference while all differences shorter are deemed to be errors
            all_reads:  list of tuples containing all the reads and (what we need) their sequences
            work_dir:   the current work directory
            k_size:     the length of our k_mers
OUTPUT: The simplified graph.
"""
def new_bubble_popping_routine(DG, all_reads, work_dir, k_size,delta_len,known_intervals):

    not_viable_global=set()
    not_viable_multibubble = set()
    has_combinations=True
    old=True
    nr_popped=0
    overall_pops=0
    spoa_count=0
    multi_consensuses={}
    merged_dict={}
    iteration_number=0
    this_it_pops=0
    initial_edge_nr=len(DG.edges())
    #all_r_ids=set(all_reads.keys())
    print("finding all paths")
    all_paths_s_to_t, all_path_sets = find_all_read_paths(DG, all_reads, merged_dict)
    print("paths found")
    prev_marked=set()
    #print(all_paths_s_to_t[6])
    #print("ERD",equal_reads_dict)
    #changed_reads=[]
    pop_threshold = int(initial_edge_nr / 200)
    #we want to continue the simplification process as long as we find combinations that have not been deemed to be "not viable" to pop
    while has_combinations:
        overall_pops+=this_it_pops
        print("Popthreshold",pop_threshold)

        this_it_pops = 0
        iteration_number += 1
        all_paths_s_to_t, all_path_sets=update_paths(DG,all_reads,prev_marked,merged_dict,all_paths_s_to_t,all_path_sets)
        marked = set()
        #print(all_paths_s_to_t[6])
        direct_combis = []
        """if iteration_number>5:
            return"""
        #profiler = Profiler()
        #profiler.start()

        print("ITERATION NUMBER "+str(iteration_number))
        print()
        print("GRAPH NR NODES: {0} EDGES: {1} ".format(len(DG.nodes()), len(DG.edges())))
        print()
        #if iteration_number==5:
        #    DEBUG=True
        #else:
        #    DEBUG=False
        numnodes=len(DG.nodes())
        numedges=len(DG.edges())
        #if iteration_number>1:
        #    return
        #if(isCyclic(DG)):
        #    print("Cyclic Graph")
        #    return-1
        #TopoNodes holds the topological order of the nodes in our graph
        TopoNodes = list(nx.topological_sort(DG))
        topo_nodes_dict = {n : i for i,n in enumerate(TopoNodes)}
        poss_starts=[]
        #find all possible bubble start nodes in our graph
        find_possible_starts(DG, TopoNodes,poss_starts)
        poss_ends=[]
        # find all possible bubble end nodes in our graph
        find_possible_ends(DG, TopoNodes,poss_ends)
        #if DEBUG:
            #print("startnodes", poss_starts)
            #print("endnodes", str(poss_ends))
        combinations=[]
        #generate all combination of bubble start nodes and bubble end nodes in which the poss_starts comes before poss_end in TopoNodes
        generate_combinations(poss_starts, poss_ends, TopoNodes,combinations)
        #print("Combinations:", len(combinations))
        combinations_filtered=[]
        #filter out the combinations we already have analysed in a previous iteration
        filter_combinations(combinations, not_viable_global,combinations_filtered)
        #print("Combinations filtered:",len(combinations_filtered))
        #print("not_viable ", not_viable_global)

        #if we haven't found any new combinations we successfully finished our bubble popping
        if not combinations_filtered:
            break
        #if DEBUG:
            #print("combis",combinations_filtered)

        # sort the combinations so that the shortest combinations come first
        sorted_combinations=sorted(combinations_filtered, key=lambda x: TopoNodes.index(x[1])-TopoNodes.index(x[0]))
        #iterate over all combinations
        for combination in sorted_combinations:
            edges=DG.edges()
            """if ('15, 49, 6', '58, 98, 6') in edges:
                print("IT is inside!")
            print(combination[0],combination[1])
            if numnodes==1029 and numedges==1304:
                print("COMBINATION STARTS HERE")
                if 6 in combination[2]:
                    #DEBUG=True
                    print("THIS_COMBI",combination)
                #else:
                    #DEBUG=False
                if combination[0]=='s' and combination[1]=='72, 99, 3':
                    thebug = True

                elif combination[0]=='1673, 1687, 9' and combination[1]=='714, 779, 19':
                    thebug = True
                else: thebug=False
            else:
                thebug=False
            #print(thebug)"""
            thebug=False
            #print(thebug)
            is_alignable=True
            all_paths = []
            if not old:
                all_paths=find_combi_paths(combination,all_paths_s_to_t)
            else:
            #we find all paths from s' to t' via find_paths
                    #if DEBUG:
                    #    cyclic=isCyclic(DG)
                    #    if cyclic:
                    #        print("CYCLIC",cyclic)
                    if thebug:
                        find_edges_with_supp(6, DG)
                        print("FINDPATHS outside")
                    find_paths(DG,combination[0],combination[1],combination[2],all_paths)
                    """if not all_paths_old==all_paths:
                        print("OLD",all_paths_old)
                        print("NEW",all_paths)"""
                    if not all_paths:
                        combiset= {combination[0],combination[1]}
                        #print("NO ALLPATHS found")
                        for r_id in combination[2]:
                            print("Difference",combiset.difference(all_paths_s_to_t[r_id]))
                        #not_viable_global.add(combination)
            #print("FINDPATHS ended")
            #if DEBUG:
                #print("all_paths:", all_paths,"len",len(all_paths))
                #for path in all_paths:
                    #print(path[0], path[1], path[2])
            initial_all_paths=len(all_paths)
            #if DEBUG:
                #print("initial_all_paths", initial_all_paths)
            #print("All paths",all_paths)
            #print("Initial_all_paths",initial_all_paths)
            #if we only did find one viable path from s' to t' we figure the combination was not viable
            if len(all_paths)==1:
                not_viable_global.add(combination)

                #print("No pop- not enough paths")
            all_paths_filtered=[]
            #we cannot touch the same node over and over during one iteration as the bubbles do not display how the graph changed
            all_paths_filtered=filter_out_if_marked(all_paths,marked,direct_combis,combination[1])
            #if DEBUG:
                #print("all_paths_filtered",all_paths_filtered)
            #print("all_paths after:", all_paths, "len", len(all_paths))
            #print("initial_all_paths",initial_all_paths)

            #consensus_infos stores the positions and read infos for generating the consensus
            consensus_infos={}
            #if we found two paths in our bubble
            if len(all_paths_filtered)==2:
                if DEBUG:
                    print("Found exactly 2 paths for this it->filtered")
                if initial_all_paths>2:
                    is_megabubble = True
                    if DEBUG:
                        print("This is a megabubble")
                else:
                    is_megabubble=False

                #print("2-path ", all_paths_filtered)

                this_combi_reads = tuple(sorted(set(all_paths_filtered[0][1]) | set(all_paths_filtered[1][1])))
                this_combi = (combination[0], combination[1], this_combi_reads)


                #print("NVMB",not_viable_multibubble)
                if DEBUG:
                    print("this_combi", this_combi)

                if this_combi in not_viable_multibubble:

                    #print(this_combi ,"in ", not_viable_multibubble)

                    continue
                #we can set the paths by accessing all_paths_filtered as we have a well-defined 2-path bubble
                path1=all_paths_filtered[0][0]
                path2=all_paths_filtered[1][0]

                #for path in all_paths_filtered:
                ##print("current_path",path)
                #it can happen that one of the paths only consists of one edge. We then set pathnode to be the endnode of the bubble, the pathnode is the path defining node
                if len(path1)>1:
                    pathnode1=path1[1]
                else:
                    pathnode1=combination[1]
                if len(path2)>1:
                        pathnode2 = path2[1]
                else:
                        pathnode2 = combination[1]

                #print("pathnode ",pathnode1)
                #print("pathnode2",pathnode2)

                #we now would like to know whether the bubble paths have nodes in common
                p_set1=set(path1[1:])
                p_set2=set(path2[1:])
                intersect=p_set1.intersection(p_set2)
                #print("intersect ",intersect)

                #we need to find out whether the paths intersect( have a common node)
                # if they do not have an intersection we can continue
                if not intersect:
                        consensus_infos[pathnode1]=get_consensus_positions(combination[0], combination[1], DG, all_paths_filtered[0][1])
                        consensus_infos[pathnode2] = get_consensus_positions(combination[0], combination[1], DG, all_paths_filtered[1][1])
                #The paths intersect
                else:
                    if initial_all_paths==2:
                        not_viable_global.add(combination)
                    #if we have more than two paths: we add only this combination as invalid to not_viable_multibubble
                    else:
                        this_combi_reads = tuple(sorted(set(all_paths_filtered[0][1]) | set(all_paths_filtered[1][1])))
                        this_combi = (combination[0], combination[1], this_combi_reads)
                        if DEBUG:
                            print("not_viable_multibubble add", this_combi)
                        # we only know about this combination of paths so we only set not_viable_multibubble
                        not_viable_multibubble.add(this_combi)
                            #print("Not viable now:", not_viable_global)
                    is_alignable=False
                    if DEBUG:
                        print("Not POPPABLE")
                #print("our consensus infos:",consensus_infos)
                if is_alignable:
                    is_poppable,cigar,seq_infos,consensus_info_log,spoa_count=align_bubble_nodes(all_reads,consensus_infos,work_dir,k_size,spoa_count,multi_consensuses,is_megabubble,this_combi,delta_len)

                    if is_poppable:

                        linearize_bubble(DG,consensus_infos,combination[0],combination[1],all_paths_filtered,combination[2],seq_infos,consensus_info_log, topo_nodes_dict)
                        this_it_pops+=1
                        nr_popped += 1
                        #add all nodes that have been affected to marked
                        for node in all_paths_filtered[0][0]:
                            marked.add(node)
                        for node in all_paths_filtered[1][0]:
                            marked.add(node)
                        #if we find a directpath from s' to t'
                        if not all_paths_filtered[0][0] or not all_paths_filtered[1][0]:
                            #add the combination to direct_combis
                            combi_list=[]
                            combi_list.append(combination[0])
                            combi_list.append(combination[1])
                            if not combi_list in direct_combis:
                                direct_combis.append(combi_list)
                    else:
                        #print("not poppable")
                        if initial_all_paths==2:
                            not_viable_global.add(combination)
                            #if DEBUG:
                                #print("Not viable now:", not_viable_global)
                        else:

                            this_combi_reads = tuple(sorted(set(all_paths_filtered[0][1]) | set(all_paths_filtered[1][1])))
                            this_combi = (combination[0], combination[1], this_combi_reads)
                            if DEBUG:
                                print("not_viable_multibubble add",this_combi)
                            # we only know about this combination of paths so we only set not_viable_multibubble
                            not_viable_multibubble.add(this_combi)
            # we have more than two paths connecting s' and t'. We now want to efficiently compare those paths
            elif len(all_paths_filtered)>2:
                directpath_marked=False
                if DEBUG:
                    print("more paths in", combination)
                    print("APF",all_paths_filtered)
                #initial_listing is a list holding all possible combinations of paths
                initial_listing=[(p1,p2) for (p1,p2) in itertools.combinations(all_paths_filtered, 2) if (combination[0],combination[1],tuple(sorted( set(p1[1]) | set(p2[1])))) not in not_viable_multibubble]
                if DEBUG:
                    print(initial_listing)
                #if initial_listing is empty: we do not have a viable bubble before us
                if not initial_listing:
                    #print("Not viable now bef :",not_viable_global)
                    not_viable_global.add(combination)
                    #print("Not viable now:",not_viable_global)
                    continue
                #print("APF",all_paths_filtered)
                for path_combi in initial_listing:
                    #if DEBUG:
                        #print("path_combi", path_combi)
                    p1=path_combi[0]
                    p2=path_combi[1]
                    p_set1 = set(p1[0][1:])
                    p_set2 = set(p2[0][1:])
                    this_combi_reads = tuple(sorted(set(p1[1]) | set(p2[1])))
                    this_combi = (combination[0], combination[1], this_combi_reads)
                    if p_set1.intersection(p_set2):
                        # print("this combi not viable",this_combi)
                        not_viable_multibubble.add(this_combi)
                        continue
                    consensus_infos={}
                    if (not p1[0]) or (not p2[0]) and directpath_marked:
                        if DEBUG:
                            print("Marked")
                        continue
                    p1_filtered = filter_path_if_marked(marked, p1[0])
                    p2_filtered = filter_path_if_marked(marked, p2[0])
                    if DEBUG:
                        print()
                    #we con only pop the bubble if both paths have not been affected by previous bubble popping steps
                    if not p1_filtered and not p2_filtered:
                        if not len(p1[0])>1:
                            pathnode1 = combination[1]
                        else:
                            pathnode1 = p1[0][1]
                        if not len(p2[0])>1:
                            pathnode2 = combination[1]
                        else:
                            pathnode2 = p2[0][1]

                        if DEBUG:
                            print("not filtered")
                        consensus_infos[pathnode1] = get_consensus_positions(combination[0], combination[1], DG, p1[1])
                        consensus_infos[pathnode2] = get_consensus_positions(combination[0], combination[1], DG, p2[1])
                        is_poppable, cigar, seq_infos, consensus_info_log, spoa_count = align_bubble_nodes(all_reads, consensus_infos,
                                                                                           work_dir, k_size,spoa_count, multi_consensuses, True, this_combi, delta_len)
                        if DEBUG:
                            print("Do we pop?",is_poppable)
                        if is_poppable:
                            #print("POPPED_Multi")
                            if DEBUG:
                                print("POPPED_Multi")
                            all_paths_filtered=[]
                            all_paths_filtered.append(p1)
                            all_paths_filtered.append(p2)
                            if DEBUG:
                                print("init",initial_all_paths)
                            #print("ALL_Paths_filtered",all_paths_filtered)
                            linearize_bubble(DG, consensus_infos, combination[0], combination[1], all_paths_filtered,
                                         combination[2], seq_infos, consensus_info_log, topo_nodes_dict)
                            this_it_pops += 1
                            nr_popped+=1
                            if (nr_popped % 10) == 0:
                                print("NR_popped", nr_popped)
                            #add all nodes that were part of the bubble paths to marked
                            for node in all_paths_filtered[0][0]:
                                marked.add(node)
                            for node in all_paths_filtered[1][0]:
                                marked.add(node)
                            if not all_paths_filtered[0][0] or not all_paths_filtered[1][0]:
                                directpath_marked = True
                                combi_list = []
                                combi_list.append(combination[0])
                                combi_list.append(combination[1])
                                if not combi_list in direct_combis:
                                    direct_combis.append(combi_list)
                        else:
                            #the combination is not poppable
                            not_viable_multibubble.add(this_combi)
        print("This iterations pops ", this_it_pops)
        if this_it_pops<pop_threshold:
            break
        prev_marked=marked
        #profiler.stop()

        #profiler.print()
    print("Overall number of bubbles popped", overall_pops)
DEBUG=False

"""Overall method used to simplify the graph
Works as wrapper script for new_bubble_popping_routine       
INPUT:  DG: Directed Graph before simplifications,
        all_reads:Dictionary containing the sequences as well as ids for all reads that we analyze
        work_dir: The current working directory
        k_size: parameter k for our minimizers
OUTPUT: DG: Graph after simplification took place    """
def simplifyGraph(DG, all_reads, work_dir, k_size,delta_len,known_intervals):
    #print("Current State of Graph:")
    ##print(DG.nodes(data=True))
    ##print(DG.edges(data=True))
    #draw_Graph(DG)
    #profiler = Profiler()
    #profiler.start()
    print("Simplifying the graph")
    new_bubble_popping_routine(DG, all_reads, work_dir, k_size,delta_len,known_intervals)
    #profiler.stop()
    #profiler.print()
    #print("Cycles:",list_of_cycles)
    #print("Popping bubbles done")
    #draw_Graph(DG)
