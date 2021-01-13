"""function to merge consecutive nodes, if they contain the same reads to simplify the graph
    INPUT: DG Directed Graph
    OUTPUT: DG: Directed Graph with merged nodes
"""
def merge_nodes(DG):
    #iterate over the edges to find all pairs of nodes
    edgesView=DG.edges.data()
    for ed_ge in edgesView:
        startnode = ed_ge[0]
        endnode = ed_ge[1]
        #we only need to know the out degree of the start node and the end degree of the end node
        start_out_degree = DG.out_degree(startnode)
        end_in_degree = DG.in_degree(endnode)
        #if the degrees are both equal to 1 and if none of the nodes is s or t
        if(start_out_degree==end_in_degree==1and startnode!="s"and endnode!="t"):
            #print("Merging nodes "+startnode+" and "+endnode)
            #use the builtin function to merge nodes, prohibiting self_loops decreases the amount of final edges
            DG=nx.contracted_nodes(DG, startnode,endnode,self_loops=False)


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

