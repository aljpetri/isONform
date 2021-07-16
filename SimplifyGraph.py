import networkx as nx
from collections import Counter
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
def generate_subgraphs(DG,bubbles):
    for bubble in bubbles:
        SG=DG.subgraph(bubble)
        draw_Graph(SG)
"""Helper method used to generate a subgraph for a list of nodes
INPUT:      DG          Directed Graph to plot
                    bubbles:    list of bubblenodelists to be plotted
"""
def generate_subgraph(DG,bubble):
    #for bubble in bubbles:
        SG=DG.subgraph(bubble)
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
    print("fbslistofnodes")
    print(listofnodes)
    real_bubble=True
    out= set()
    outedges=[]
    for node in listofnodes:
        print(node)
        out_edges_all=DG.out_edges(node)
        print(out_edges_all)
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
    print(listofnodes)
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

"""Helper method which is used to find the distance of a node to the start node. This is done by using the read information for both nodes
    INPUT:          DG      Our directed Graph
                    bubble_start        The start node of our bubble
                    path_node           The node for which we would like to retreive the distance
                    path_reads:         The reads which we are looking at to get the distances
    OUTPUT:         dist_to_start:      mean distance between the two nodes (we only look at the path_reads)
"""
def get_distance_to_start(DG,bubble_start, path_node,path_reads):
    actual_reads = [a_tuple[0] for a_tuple in path_reads]
    print("actual_reads",actual_reads)
    sum=0
    #we iterate over actural_reads (the read ids we got from path_reads)
    for i,read in enumerate(actual_reads):
        start_infos=DG.nodes[bubble_start]['reads']
        print("startinfos",start_infos)
        start_tuple=start_infos[read]
        print("starttuple",start_tuple)
        start_pos=start_tuple[1]
        node_infos=DG.nodes[path_node]['reads']
        end_tuple=node_infos[read]
        end_pos=end_tuple[0]
        sum+=end_pos-start_pos
    dist_to_start=sum/(i+1)
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
def remove_edges(DG,path_reads,bubble_start,bubble_end,path_nodes):
    print("Path nodes", path_nodes)
    print("PAth_reads",path_reads)
    tup_list=[]
    #we store all the infos for the deleted edges into edges_to_delete
    edges_to_delete= {}
    #we have to delete the edges which connect the nodes that are inside the bubble
    for startnode_id,path_node_list in path_nodes.items():
        print(path_node_list)
        #we delete the edge bubble_start-path_node_list[0] and add its infos to edges_to_delete
        edges_to_delete[bubble_start,path_node_list[0]]=DG[bubble_start][path_node_list[0]]
        for index,path_node in enumerate(path_node_list):
            print(index,", ",path_node)
            dist=get_distance_to_start(DG,bubble_start, path_node,path_reads[startnode_id])
            tup=dist,path_node
            tup_list.append(tup)
            if path_node!=path_node_list[-1]:
                edges_to_delete[path_node,path_node_list[index+1]]=DG[path_node][path_node_list[index+1]]
            else:
                print("BER",bubble_end)
                print(DG.nodes[bubble_end])
                entry = DG.get_edge_data(path_node,bubble_end)
                edges_to_delete[path_node,bubble_end]=entry

    print("ETD",edges_to_delete)
    for edge,edge_infos in edges_to_delete.items():
        print(edge)
        DG.remove_edge(edge[0],edge[1])
    return edges_to_delete
def generate_consensus_path(work_dir,consensus_attributes,reads,k_size):
    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    for i, (q_id, pos1, pos2) in  enumerate(consensus_attributes, 3):
        if pos2==0:
            pos2=len(reads[q_id][1])-k_size
        seq = reads[q_id][1][pos1: pos2 + k_size]
        reads_path.write(">{0}\n{1}\n".format(str(q_id)+str(pos1)+str(pos2), seq))
    reads_path.close()
    # print(reads_path.name)
    # sys.exit()
    spoa_ref = run_spoa(reads_path.name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
    print("spoa_ref",spoa_ref)
    return spoa_ref
def parse_cigar_differences(cigar_string,delta_len):
    good_to_pop=True
    for elem in cigar_string:
        cig_len=elem[0]
        cig_type=elem[1]
        #all matches are absolutely fine
        if (cig_type!= '=') and (cig_type!='M'):
            #we have a nonmatch, now we have to figure out whether this is due to an exon or indel (everything with len>delta_len is defined as exon)
            if cig_len>delta_len:
                #we found an exon being the difference between the paths->no bubble popping feasible
                good_to_pop=False
    return good_to_pop
def collect_bubble_nodes(path_nodes,consensus_infos,DG,support_dict):
    all_nodes=[]
    print("Consensus_infos",consensus_infos)
    print(path_nodes)
    all_shared=[]
    for startnode_id,path_node_list in path_nodes.items():
        shared_reads=support_dict[startnode_id]
        print("shared",shared_reads)
        all_shared.extend(shared_reads)
        print(path_node_list)
        for path_node in path_node_list:
            path_node_infos=DG.nodes[path_node]['reads']
            distance=0
            for i,s_read in enumerate(shared_reads):
                pos_tuple=path_node_infos[s_read]
                print(pos_tuple)
                distance+=pos_tuple[0]
            final_distance=distance/(i+1)
            print("final_distance",final_distance)
            all_nodes_tuple=(path_node,final_distance)
            all_nodes.append(all_nodes_tuple)
    print("all_nodes",all_nodes)
    return all_nodes,all_shared

def add_edges(DG,all_nodes,edges_to_delete,consensus_infos,bubble_start,bubble_end,all_shared,path_nodes):
    counter=0
    path1=[]
    path2=[]

    print("EdgestoDelete",edges_to_delete)
    prevallnodes=bubble_start
    print("Allnodeslen",len(all_nodes))
    #we assign both paths to variables to make them easier accessible.
    for id,path in path_nodes.items():
        print(id,path)
        if counter==0:
            path1=path
        else:
            path2=path
        counter=counter+1
    print(path1)
    print(path2)
    prevnode1 = bubble_start
    prevnode2 = bubble_start
    nextnode1 = path1[0]
    nextnode2 = path2[0]
    prevnode=bubble_start
    #this is the main part of the linearization. We iterate over all_nodes and try to find out which path the nodes belong to.
    #This info is needed as we need the current state ob both paths to add the correct edge_support to the graph
    for counter,nodetup in enumerate(all_nodes):
        node=nodetup[0]
        print("Node",node)
        reads_for_node=[]
        reads_for_node.extend(all_shared)
        print("P1",path1)
        print("P2", path2)
        if node in path1:
            prevnode1=path1.pop(0)
            if len(path1)<1:
                nextnode1=bubble_end
                print("Nextnode1", nextnode1)
            else:
                nextnode1=path1[0]
        else:
            prevnode2=path2.pop(0)
            if len(path2)<1:
                nextnode2=bubble_end
                print("Nextnode2",nextnode2)
            else:
                nextnode2=path2[0]
        new_edge_supp1=edges_to_delete[prevnode1,nextnode1]['edge_supp']
        print("NES1", new_edge_supp1)
        new_edge_supp2=edges_to_delete[prevnode2,nextnode2]['edge_supp']
        full_edge_supp=new_edge_supp1+new_edge_supp2
        full_edge_supp_final =[]
        [full_edge_supp_final.append(x) for x in full_edge_supp if x not in full_edge_supp_final]
        print("NES2",new_edge_supp2)
        print(DG.edges(data=True))
        DG.add_edge(prevnode,node,edge_supp=full_edge_supp_final)
        print("Adding node from ",prevnode, "to ",node)
        print(DG.edges(data=True))
        prevnode=node
    new_edge_supp1 = edges_to_delete[prevnode1, bubble_end]['edge_supp']
    new_edge_supp2 = edges_to_delete[prevnode2, bubble_end]['edge_supp']
    full_edge_supp = new_edge_supp1 + new_edge_supp2
    full_edge_supp_final = []
    [full_edge_supp_final.append(x) for x in full_edge_supp if x not in full_edge_supp_final]
    DG.add_edge(prevnode, bubble_end, edge_supp=full_edge_supp_final)
    print("Adding node from ", prevnode, "to ", bubble_end)
    print("FinalEdges",DG.edges(data=True))

"""Actual linearization process of our bubbles
        INPUT:      DG      Our directed Graph
                consensus_infos:    
                bubble_start        The start node of our bubble
                bubble_end:         The end node of our bubble
                path_nodes:          A list of nodes which make up a path in our bubble
                edges_to_delete:         A dictionary holding all edges which we want to delete and their respective attributes (needed for adding edges)
                support_dict:
                
"""
def linearize_bubble(DG, consensus_infos, bubble_start, bubble_end, path_nodes,edges_to_delete,support_dict):
    #The main idea of this method is: 1. Get all the nodes which are in both paths and calculate their avg distance to bubble_start
    #                                 2. Sort the nodes by distance
    #                                 3. Add the edges connecting the nodes in the indicated order
    all_nodes,all_shared=collect_bubble_nodes(path_nodes,consensus_infos,DG,support_dict)
    print("Allnodes before sorting",all_nodes)
    all_nodes.sort(key=lambda tup: tup[1])
    print("Allnodes after sorting", all_nodes)
    add_edges(DG,all_nodes,edges_to_delete,consensus_infos,bubble_start,bubble_end,all_shared,path_nodes)
    print("Popped bubble ",path_nodes)
""" Method to simplify the graph. 
    INPUT:  DG  Directed Graph
           delta_len   parameter giving the maximum length difference between two paths to still pop the bubble
    OUTPUT: DG Directed Graph without bubbles
        """
#TODO: We have to generate a consensus of the reads in one bubble path to be able to have a pairwise sequence alignment between the bubble paths
def align_and_linearize_bubble_nodes(DG,bubble_start,bubble_end,delta_len,all_reads,consensus_infos,work_dir,k_size,path_nodes,support_dict):
    #TODO: generate consensus sequences via spoa
    consensus_list=[]
    for path_node,consensus_attributes in consensus_infos.items():
        print("consensus",consensus_attributes)
        if len(consensus_attributes)>1:
            con=generate_consensus_path(work_dir, consensus_attributes, all_reads, k_size)
        else:
            (q_id,pos1,pos2)=consensus_attributes[0]
            con=all_reads[q_id][1][pos1: pos2]
        consensus_list.append(con)
    print(consensus_list)
    consensus1=consensus_list[0]
    consensus2=consensus_list[1]
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score=parasail_alignment(consensus1,consensus2,match_score=2,missmatch_penalty=-2,opening_penalty=3,gap_ext=1)
    print(s1_alignment)
    print(s2_alignment)
    print(cigar_string)
    print(cigar_tuples)
    good_to_pop=parse_cigar_differences(cigar_tuples,delta_len)
    if good_to_pop:
        print("time for linearization")
        print("Edges befdel",DG.edges(data=True))
        edges_to_delete=remove_edges(DG, consensus_infos, bubble_start, bubble_end, path_nodes)
        print("Edges afdel",DG.edges(data=True))
        linearize_bubble(DG, consensus_infos, bubble_start, bubble_end, path_nodes,edges_to_delete,support_dict)
        print("Edges aflin",DG.edges(data=True))
    print(score)
    #TODO implement linearize_bubble to alter DG in order to pop the bubble



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
    visited_nodes=[]
    intersect_supp = []
    #we only iterate until we have reached bubble_end
    while curr_node != bubble_end:
        visited_nodes.append(curr_node)
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
                #if we have an overlap, this means that at least one read from our initial_support also supports this edge
                if intersect_supp:
                    #we have found a node to continue
                    further_node=True
                    #set the node to be our current_node
                    curr_node=next_node
                    #update the set of still supporting reads
                    curr_support=intersect_supp
                    #as we only expect one node to be viable: break the for loop to continue in the next node
                    break
        #we did not find any node we could continue with->our path is not the path of a real bubble
        if not further_node:
            is_viable = False

            return is_viable,visited_nodes,intersect_supp
    #we were able to reach bubble_end->we just confirmed that we have a bubble path
    is_viable=True
    print("Visited_nodes ",visited_nodes)
    return (is_viable, visited_nodes,intersect_supp)
def get_path_starts(cycle,bubble_start,bubble_end,DG,shared_reads):
    #We want to find the nodes, which denote the start points for each path(as we have to find out which reads are in which path)
    min_node_out, max_node_in=get_path_nodes(cycle, bubble_start, bubble_end, DG)
    #Now we want to get the actual reads for each path
    path_starts=get_path_reads(DG,min_node_out,shared_reads)
    #iterate over all shared reads and get the pathlength for each
    print("Path_starts ",path_starts)
    return path_starts
"""Helper method for find_bubbles: This method figures out which reads belong to which part from bubble source to bubble sink
INPUT:      cycle:  A list of nodes which make up the bubble(found to be a bubble by being a cycle in an undirected graph
            min_element: The node deemed to be the starting node of the bubble (given as tuple)
            max_element: The node deemed to be the end node of the bubble (given as tuple)
            DG:         the directed graph we want to pop the bubble in
            contains_s: A boolean value denoting whether the cycle contains node "s"
            contains_t: A boolean value denoting whether the cycle contains node "t"

"""

def get_path_reads_length(r_ids, bubble_start, bubble_end, DG,shared_reads):
    id_counter=0
    read_length=0
    read_list=[]
    for r_id in r_ids:
        if r_id in shared_reads:
                max_node_infos=DG.nodes[bubble_end]['reads']
                #print("max_node_infos")
                print(max_node_infos)
                min_node_infos=DG.nodes[bubble_start]['reads']
                #print("min_node_infos")
                print(min_node_infos)
                bubble_end_pos=max_node_infos[r_id]
                bubble_start_pos=min_node_infos[r_id]
                start_of_bubble=bubble_start_pos[1]
                end_of_bubble=bubble_end_pos[0]
                entry=(r_id,start_of_bubble,end_of_bubble)
                print(entry)
                read_list.append(entry)
                read_length+=(end_of_bubble-start_of_bubble)
                id_counter+=1
    path_length=read_length/id_counter
    return path_length,read_list

"""
function which finds the bubbles and if feasible pops them in the graph

INPUT:
    DG:         The directed graph in which we pop the bubbles
    delta_len   Maximum length difference for which the paths are merged ie the bubble is popped
    all_reads   dictionary containing all the reads(string sequence) and their ids
    
    
OUTPUT:
    DG:         Directed Graph containing the popped bubbles"""

def find_and_pop_bubbles(DG, delta_len,all_reads,work_dir,k_size):
    bubbles=find_bubbles(DG)
    nr_bubbles=len(bubbles)
    print("Found "+str(nr_bubbles)+" bubbles in our graph")
    print("Bubbles: ",bubbles)
    #just for debugging and curiosity reasons: We introduce an integer counting the number of pseudobubbles (not true bubbles)
    filter_count=0
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
        # find the nodes which are on either path to be able to tell apart the paths on the bubble
        #readlen_dict,consensus_infos=get_path_reads_length(bubble_nodes, bubble_start, bubble_end, DG,shared_reads)
        path_starts=get_path_starts(bubble_nodes, bubble_start, bubble_end, DG, shared_reads)
        readlen_dict = {}
        consensus_infos = {}
        path_nodes_dict = {}
        no_viab_bubble=False
        read_len_dict={}
        support_dict={}
        for key, value in path_starts.items():
            (is_viable_path, path_nodes,support) = test_path_viability(DG, key, value, bubble_nodes, bubble_end)
            path_nodes_dict[key]=path_nodes
            support_dict[key]=support
            if not (is_viable_path):
                no_viab_bubble=True
                break
            path_length,read_list=get_path_reads_length(value, bubble_start, bubble_end, DG,shared_reads)
            print("read_list (consensus:_infos",read_list)
            readlen_dict[key]=path_length
            consensus_infos[key]=read_list
        if no_viab_bubble:
            filter_count += 1
            print("Filtered ", filter_count, " bubbles out")
            continue
        #print("Listof normal",listofnormalnodes)
        print("min",bubble_start)
        print("max",bubble_end)
        #print("readlendict",readlen_dict)
        # TODO compare the reads of both paths: If they differ by less than delta_len: Pop the bubble
        lengthdiff=None
        for node,length in readlen_dict.items():
            if not lengthdiff:
                lengthdiff = length
                print(lengthdiff)
            else:
                lengthdiff = abs(lengthdiff-length)
        if lengthdiff < delta_len:
            print("Linearizing bubble ",bubble_nodes)
            generate_subgraph(DG, bubble_nodes)
            align_and_linearize_bubble_nodes(DG,bubble_start,bubble_end,delta_len,all_reads,consensus_infos,work_dir,k_size,path_nodes_dict,support_dict)
            generate_subgraph(DG, bubble_nodes)
"""Overall method used to simplify the graph
During this method: - Bubbles are identified and if possible popped     
                    - Nodes are merged        
INPUT:  DG: Directed Graph before simplifications,
        max_bubblesize: maximum number of elements making up a bubble
        delta_len: Maximum length differences in between two reads
OUTPUT: DG: Graph after simplification took place    """


def simplifyGraph(DG, delta_len,all_reads,work_dir,k_size):
    print("Simplifying the Graph (Merging nodes, popping bubbles)")
    # remove edges which yield self loops, not sure yet whether it makes sense to remove or if needed
    #print("self loops")
    #print(list(nx.selfloop_edges(DG)))
    list_of_cycles = find_repetative_regions(DG)
    print(list_of_cycles)
    #print(DG.edges())
    s_reads = DG.nodes["s"]['reads']

    find_and_pop_bubbles(DG, delta_len,all_reads,work_dir,k_size)
    print("Popping bubbles done")
    merge_nodes(DG)