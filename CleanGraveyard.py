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



def update_known_intervals(known_intervals,x):
    #TODO: write function to update the amount of known intervals (merge nodes) according to the output of get_read_difference
    print("This function updates the graph to show merged nodes instead of real ones")
#TODO: write update function for graph (possible in networkx) as a reference
def find_equal_reads(known_intervals,delta_len,max_bubblesize,max_interval_delta,all_intervals_for_graph):
    #Deprecated file reader for a file containing known_intervals
    # file = open('known_intervals.txt', 'rb')
    # known_intervals = pickle.load(file)
    # file.close()

    # sort the tuples by interval start positions.
    for r_ids, intervals in enumerate(known_intervals):
        # print(type(intervals))
        known_intervals[r_ids] = sorted(intervals, key=lambda x: x[0])
    # for key, value in known_intervals.items():
    #    print(key, value)
    #the final output of the function: A list of lists structure which contains lists of all r_ids of reads which can be merged
    node_ids_all_reads = []
    #a list containing all nodeids (meaning the id's of all minimizer spans) present in a read
    nodes_all_reads = []
    # list to keep track which differences were found between two reads
    differences_overview = []
    #fill the nodes_all_reads with all the intervals found in the reads
    #for r_ids, intervals in enumerate(known_intervals):
    #    nodeids = [x[1] for x in intervals]
    #    nodes_all_reads.append(nodeids)
    #    print(r_ids, nodeids)
    #print("done")
    #list of reads already considered. If a read has been added into a consensus already we do not have to consider it again
    considered_reads = []
    #iterate over all reads' interval lists
    for i in range(0, len(known_intervals)):
        #if a read has not yet been considered before
        if not (i + 1 in considered_reads):
            # list of reads which are equal to the current read. Those reads will be merged with the current read to form a consensus
            reads_to_merge = []
            #add this read to the list of equal reads (+1 as we want to have the read_id, but the iteration starts at 0 (r_ids start with 1)
            reads_to_merge.append(i + 1)
            # add this read to the list of considered reads(+1 as we want to have the read_id, but the iteration starts at 0 (r_ids start with 1). We know you now
            considered_reads.append(i + 1)
            #this is the step in which we compare the read to all other reads to find out whether there are mergeable ones
            for j in range(i + 1, len(known_intervals)):
                # if a read has not yet been considered before
                if not (j + 1 in considered_reads):
                    #differing_nodes=nodes_all_reads[i].symmetric_difference(nodes_all_reads[j])
                    #too_different is a boolean, denoting whether the reads are too different to merge them
                    too_different,differences,bubble_break=get_read_difference(known_intervals,i,j,max_interval_delta,delta_len,max_bubblesize,all_intervals_for_graph)
                    if not too_different:
                        read_to_pop = j + 1
                        reads_to_merge.append(read_to_pop)
                        considered_reads.append(read_to_pop)
                    else:
                        if bubble_break:
                            print("bubble_break")
                        else:
                            print( "no bubble_break")

                        #for elem in different_nodes[i]:
                            #nodes_all_reads.index(elem)
                    #if not different_nodes:
                        #print(different_nodes)
                    #if nodes_all_reads[i] == nodes_all_reads[j]:

            node_ids_all_reads.append(reads_to_merge)
    print("Knownreads:")
    print(str(len(considered_reads)))
    print(considered_reads)
    # for r_ids2,intervals2 in known_intervals.items():
    # do not pop if read is only equal to itself
    #   if not r_ids==r_ids2:
    # do only pop if all intervals of the reads are equal
    # change here to only look at the name(id) and not start/stop anymore
    #      if intervals==intervals2:
    #          popentries=(r_ids,r_ids2)
    #          popitems.append(popentries)
    #          print("deleted read "+str(r_ids2)+"from known_intervals as equal to read "+str(r_ids))
    # for popits in popitems:
    #    if popits[0]<popits[1]:
    #        r_id=popits[1]
    # print(type(r_id))
    #        if r_id in known_intervals.keys():
    #            known_intervals.pop(r_id)
    # for key, value in known_intervals.items():
    #    print(key, value)
    # print("And now for the isoforms")
    for equalreads in node_ids_all_reads:
        print(equalreads)
    print("Found " + str(len(node_ids_all_reads)) + "different isoforms.")
    return node_ids_all_reads
    # for mainid,otherids in isoforms_by_reads.items():
    #    print("Read "+str(mainid) +" is equal to the following reads:")
    #    print(','.join(str(x) for x in otherids))
    # for popits in popitems:
    # if popits[0] < popits[1]:
    #        r_id = popits[1]
    # print(type(r_id))
    #        if r_id in known_intervals2.keys():
    #            known_intervals2.pop(r_id)
    # for key, value in known_intervals2.items():
    #    print(key, value)
#TODO: invoke all_intervals_for_graph to retreive changes in the graph
def get_read_difference(known_intervals,r_id1,r_id2,max_interval_delta,delta_len,max_bubblesize):
    #list_i contains all nodeids for read r_id1
    list_r_id1 = map(lambda item: item[1], known_intervals[r_id1])
    # list_j contains all nodeids for read r_id2
    list_r_id2 = map(lambda item: item[1], known_intervals[r_id2])
    #known_inter_r_id_1 contains the list of tuples saved in known_intervals for read r_id1
    known_inter_r_id_1=known_intervals[r_id1]
    #known_inter_r_id_1 contains the list of tuples saved in known_intervals for read r_id2
    known_inter_r_id_2 = known_intervals[r_id2]
    #a dictionary holding the read id as key and a list of interval ids as a value
    list_difference= {}
    #list of interval ids for read r_id1 (value for list_difference)
    difference1=[]
    # list of interval ids for read r_id2 (value for list_difference)
    difference2 =[]
    #integer value holding the length of the bubble
    bubblesize=0
    #integer value holding the end of the interval previous to the bubble
    prevpos=0
    #list of bubbles for read r_id 1
    bubbles_rid_1_len=[]


    #iterate over the elements of list_r_id1 and find out, if 'item' is also an element of list_r_id_2
    for i,item in enumerate(list_r_id1):
        #if the interval is not in list_r_id2:
        if item not in list_r_id2:
            # list of intervals' ids only found in read r_id_1
            #exclusive_int_r_1 = []
            # list of intervals' ids only found in read r_id_1
            #exclusive_int_r_2 = []
            #if a certain element is in list_r_id1 but not in list_r_id2 it is added to difference1
            difference1.append(item)
            #find out whether difference1 is already bigger than the given threshold
            if len(difference1)>max_interval_delta:
                #if difference1>max_interval_delta, the difference between the reads is too big
                read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
                #The third value in the return statement indicates, whether a too much intervals (or a bubble length) were the reason for telling the reads apart
                return (True, read_ids_for_differences,True)
            #if difference1<max_interval_delta
            else:
                #this is part of a bubble, therefore increase bubblesize
                bubblesize+=1
                #if bubblesize>max_bubblesize do not continue looking at the reads
                if bubblesize>max_bubblesize:
                    read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
                    return (True, read_ids_for_differences,True)
        #if the interval is also in list_r_id2
        else:
            #tupl has the information about the interval(startpos,id,endpos)
            tupl = known_inter_r_id_1[i]
            #if the bubblesize is >0 this means we are at the end of a bubble
            if bubblesize>0:
                #find out how long (Bp) the bubble was
                bubblelen_i=tupl[0]-prevpos
                #append the information about the bubble (interval at end(id),length of bubble) to bubbles_rid_1_len to be able to access it later
                bubbles_rid_1_len.append((item,bubblelen_i))
                #set bubblesize to be 0 as the bubble did end here
                bubblesize=0
            #get the position of the end of this interval(possibly needed for bubblelength calculation)
            prevpos = tupl[2]
    #diff_counter counts the total amount of intervals which are not equal in both reads
    diff_counter = len(difference1)
    bubblesize=0
    prevpos=0
    # iterate over the elements of list_r_id2 and find out which items are not in list_r_id1
    for i,item in enumerate(list_r_id2):
        if item not in list_r_id1:
            difference2.append(item)
            diff_counter+=1
            if diff_counter>max_interval_delta:
                read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
                return (True, read_ids_for_differences,True)
            else:
                bubblesize+=1
                if bubblesize>max_bubblesize:
                    read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
                    return (True,read_ids_for_differences,True)
        else:
            tupl = known_inter_r_id_2[i]
            if bubblesize > 0:
                bubblelen_j = tupl[0] - prevpos
                #bubbles_rid_2_len.append((item, bubblelen_i))
                #find the tuple in r_id_1 which has the same start interval as we need
                suitable_tuple=[item_i for item_i in bubbles_rid_1_len if item_i[0] ==item]
                if suitable_tuple:
                    len_difference=abs(bubblelen_j-suitable_tuple[1])
                else:
                    read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
                    return (True, read_ids_for_differences,False)
                if len_difference>delta_len:
                    read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
                    return (True, read_ids_for_differences,True)
                bubblesize = 0
            prevpos = tupl[2]
            # list_i.remove(item)
            # list_j.remove(item)
    # store the read ids into a tuple to find out which r_id contains which intervals exclusively
    read_ids_for_differences = (r_id1, r_id2, difference1, difference2)
    #list_difference[r_id1]=difference1
    #list_difference[r_id2]=difference2
    return (False,read_ids_for_differences,False)

