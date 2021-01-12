import _pickle as pickle
import collections
from sys import stdout
import networkx as nx
import itertools
import matplotlib.pyplot as plt
from SimplifyGraph import *
from IsoformGeneration import *

def generateSimpleGraphfromIntervals(all_intervals_for_graph):
    G = nx.DiGraph()
    # a source and a sink node are added to the graph in order to have a well-defined start and end for the paths
    G.add_node("s")
    G.add_node("t")
    # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
    previous_node = "s"
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed are not already present
    for r_id, intervals_for_read in all_intervals_for_graph.items():  # intervals_for_read holds all intervals which make up the solution for the WIS of a read
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"
        # iterate over all intervals, which are in the solution of a certain read
        for inter in intervals_for_read:
            # the name of each node is defined to be startminimizerpos , endminimizerpos
            name = str(inter[0]) + ", " + str(inter[1])  # +str(r_id)
            if not G.has_node(name):
                G.add_node(name)

            # add edge between current node and previous node
            G.add_edge(previous_node, name)
            # whatever we did before, we have to set previous_node to the node we looked at in this iteration to keep going
            previous_node = name
        # the last node of a path is connected to the sink node t
        G.add_edge(previous_node, "t")

    return G
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
def add_prior_read_infos(inter,r_id,prior_read_infos,name,k):
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
            if not tuple_for_data_structure in prior_read_infos:
                prior_read_infos[tuple_for_data_structure] = name

    return prior_read_infos

"""Helper method to convert the array delivered with all_intervals_for_graph into a hash value to more efficiently look up occurence
This method additionally deletes the first three entries of the array as they contain the infos about this interval occurence, which changes in between instances
INPUT:  info_array: The array, which was delivered with the interval to indicate where the interval occurs in other reads
OUTPUT: hash(tup):  A has value of the tuple the shortened interval was converted into. This hash makes it easy to see whether the interval is already present in the read
"""

def convert_array_to_hash(info_array):

    #preprocessing: Delete the first three elements from the array, as they contain the information about this occurrence
    for x in range(0, 3):
       info_array.pop(0)

    tup=tuple(info_array)
    return(hash(tup))


"""Function to get all nodes which are part of a cycle
    
    INPUT:  current_read_state: The current state of known_intervals[r_id-1] (known_intervals of the current read)
            cycle_start:        The last occurence of the repeating node before the cycle starts
    OUTPUT: cycle_nodes:        A list of entries indicating which nodes are in the cycle, having the following form: (startpos, node_id, endpos)


"""
def record_cycle(current_read_state,cycle_start):
    cycle_nodes=[]
    indices = [i for i, tupl in enumerate(current_read_state) if tupl[1] == cycle_start]
    index=indices[0]
    for element in range(index,len(current_read_state)):
        cycle_nodes.append(current_read_state[element])

    return(cycle_nodes)


""" generates a networkx graph from the intervals given in all_intervals_for_graph.
# INPUT: all_intervals_for_graph: A dictonary holding lists of minimizer intervals.


TODO:   add dictionary to store infos about known instances 
        revisit structure of the method and introduce subroutines

"""
def generateGraphfromIntervals(all_intervals_for_graph, k,delta_len):
    DG = nx.DiGraph()
    cycle_nodes={}
    #add the read ids to the startend_list
    reads_at_startend_list = []
    reads_for_isoforms=[]
    for i in range(1,len(all_intervals_for_graph)+1):
        reads_at_startend_list.append(i)
        reads_for_isoforms.append(i)
    # a source and a sink node are added to the graph in order to have a well-defined start and end for the paths
    DG.add_node("s",reads=reads_at_startend_list)
    DG.add_node("t",reads=reads_at_startend_list)

    # holds the r_id as key and a list of tuples as value: For identification of reads, also used to ensure correctness of graph
    known_intervals = []
    #adds an empty list for each r_id to known_intervals. To those lists, tuples, representing the intervals are added
    for _ in itertools.repeat(None, len(all_intervals_for_graph)):
        known_intervals.append([])
    #print(known_intervals)
    nodes_for_graph = {}
    prior_read_infos = {}
    #for _ in itertools.repeat(None, len(all_intervals_for_graph)):
    #    prior_read_infos.append([])
    #generate a dictionary key:nodeid,value:[(alternative_node_id, length_difference)]
    alternative_nodes={}
    witnesses = {}
    cycle_in_reads={}
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    for r_id, intervals_for_read in all_intervals_for_graph.items():  # intervals_for_read holds all intervals which make up the solution for the WIS of a read
        print(r_id)
        containscycle=False
        #if not r_id in known_intervals:
        #    known_intervals[r_id] = []
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"
        previous_end=0
        #read_hashs is a dictionary storing hash values as keys and the respective node ids as values
        #TODO:If a node occurs more than twice in a read this, however, could yield issues
        read_hashs={}
        #is_repetative=False
        # the name of each node is defined to be startminimizerpos , endminimizerpos
        # iterate over all intervals, which are in the solution of a certain read
        for inter in intervals_for_read:
            info_tuple = (r_id, inter[0], inter[1])
            #generate hash value of the intervals' infos
            curr_hash=convert_array_to_hash(inter[3])

            #access prior_read_infos, if the same interval was already found in previous reads
            #TODO: use cycle_in_reads to figure out, whether a repeating node can be added to a cycle node. If not create new node.
            #TODO pt2: Think about what happens with self_cycles
            if info_tuple in prior_read_infos:
                name=prior_read_infos[info_tuple]
                #if info_tuple in repetative_infos:
                #    print("This interval is repetative:")
                #    print(info_tuple)

                this_len = inter[0] - previous_end
                if this_len < 0:
                    print("Error- length negative.")
                if not DG.has_edge(previous_node, name):
                    #update the read information of node name
                    prev_nodelist=nodes_for_graph[name]
                    prev_nodelist.append((r_id, inter[0], inter[1]))
                    nodes_for_graph[name]=prev_nodelist
                    #only add a new edge if the edge was not present before
                    length=this_len
                    DG.add_edge(previous_node, name,length=length)
                else:
                    #print(previous_node,name)
                    prev_len = DG[previous_node][name]["length"]
                    #print("Prev_len:"+str(prev_len))
                    len_difference = abs(this_len - prev_len)
                    #print("Len_difference:"+str(len_difference))
                    if len_difference < delta_len:
                        # update the read information of node name
                        prev_nodelist = nodes_for_graph[name]
                        prev_nodelist.append((r_id, inter[0], inter[1]))
                        nodes_for_graph[name] = prev_nodelist
                    else:
                        nodelist = []
                        #inappropriate_node
                        old_node=name

                        # add a node into nodes_for_graph
                        name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                        alt_info_tuple=(name,previous_node,prev_len)
                        if not(old_node in alternative_nodes):
                            alternative_nodes[old_node]=[]
                        else:
                            alternative_infos_list=alternative_nodes[old_node]
                            alternatives_filtered=[item for item in alternative_infos_list if previous_node==item[1]and abs(this_len-item[2])<delta_len]
                            #if we have found a node which this info can be added to
                            if alternatives_filtered:
                                node_info=alternatives_filtered[0]
                                name=node_info[0]
                                # update the read information of node name
                                prev_nodelist = nodes_for_graph[name]
                                prev_nodelist.append((r_id, inter[0], inter[1]))
                                nodes_for_graph[name] = prev_nodelist
                            #if we have not found a node which this info can be added to
                            else:
                                #add a new entry to alternative_nodes[old_node] to enable finding this instance
                                alternative_nodes[old_node].append(alt_info_tuple)
                                # add the read information for the node
                                nodelist.append((r_id, inter[0], inter[1]))
                                nodes_for_graph[name] = nodelist
                                DG.add_node(name)
                                # keep known_intervals up to date
                                #known_intervals[r_id - 1].append((inter[0], name, inter[1]))
                                # get the length between the previous end and this nodes start
                                length = this_len
                                # connect the node to the previous one
                                DG.add_edge(previous_node, name, length=length)
                    # keep known_intervals up to date
                known_intervals[r_id-1].append((inter[0],name,inter[1]))

            # if the information for the interval was not yet found in a previous read (meaning the interval is new)
            else:
                nodelist = []
                # add a node into nodes_for_graph
                name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                #add the read information for the node
                nodelist.append((r_id, inter[0], inter[1]))
                nodes_for_graph[name] = nodelist
                DG.add_node(name)
                #keep known_intervals up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
                # get the length between the previous end and this nodes start
                length = inter[0] - previous_end
                #connect the node to the previous one
                DG.add_edge(previous_node, name,length=length)
                prior_read_infos=add_prior_read_infos(inter,r_id,prior_read_infos,name,k)

            #set the previous node for the next iteration
            previous_node = name
            previous_end=inter[1]
            #find out whether the current hash has already been added to read_hashs
            if curr_hash not in read_hashs:
                #if the current hash was not yet present in read_hashs: Add it
                read_hashs[curr_hash]=[]
                read_hashs[curr_hash].append(name)
            #If the current hash was already present in read_hashs: This node is the second occurance of a node, forming a cycle. We have  to record the cycle to make sure other reads can be added if possible
            else:
                cycle_start=read_hashs[curr_hash][-1]
                current_read_state=known_intervals[r_id-1]
                newcyc=record_cycle(current_read_state,cycle_start)
                #TODO: Think about what to use as key for cycle_in_reads. Is dict even feasible?
                cycle_in_reads
                print("Current hash already found in the interval.")
                print(info_tuple)
                read_hashs[curr_hash].append(name)
                containscycle=True
        cycle_in_reads[r_id]=containscycle
        #add an edge from name to "t" as name was the last node in the read
        if not DG.has_edge(name, "t"):
        # print("I have no idea what I'm doing here!")
            DG.add_edge(name,"t",length=0)
        #delete the entry of prior_read_infos just to save some space (set to empty list to keep coordinates correct
        #prior_read_infos[r_id-1]=[]
    #print("Nodes for graph")
        #print(nodes_for_graph)
    #print(nodes_for_graph)
    #set the node attributes to be nodes_for_graph, very convenient way of solving this
    nx.set_node_attributes(DG,nodes_for_graph,name="reads")
    #print("Repetative Infos")
    #print(repetative_infos)
    #print repetative reads
    #print("Repetative reads")
    #print(repetative_reads)
    #print alternative_nodes
    print("Alternative nodes")
    print(alternative_nodes)
    #return DG and known_intervals as this makes some operations easier
    result = (DG, known_intervals,reads_for_isoforms,reads_at_startend_list)
    return result

def check_graph_correctness(known_intervals,all_intervals_for_graph):
    correct_graph = True
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        if not (len(intervals_for_read)==len(known_intervals[r_id-1])):
            correct_graph=False
            print("graph incorrect for read "+str(r_id))
        print("All_intervals_for_graph at read "+str(r_id)+" has "+str(len(intervals_for_read))+" Intervals")
        for inter in intervals_for_read:
                read_id = inter[3][slice(0, len(inter[3]),3)]
                number_of_reads_all_intervals=len(read_id)
        print("Known_intervals at read " + str(r_id) + " has " + str(len(known_intervals[r_id-1])) + " Intervals")

    if correct_graph:
        print("Graph built up correctly")
    else:
        print("ERROR - Incorrect Graph")
# draws a directed Graph DG
def draw_Graph(DG):
    # defines the graph layout to use spectral_layout. Pos gives the position for each node
    pos = nx.spectral_layout(DG)
    # draws the network in the given layout. Only nodes are labeled, the edges are directed but not labeled
    nx.draw_networkx(DG, pos, font_weight='bold')
    # add labels for the edges and draw them into the graph
    # labels = nx.get_edge_attributes(DG, 'weight')
    # nx.draw_networkx_edge_labels(DG,pos, edge_labels=labels)
    plt.show()

def main():
    k_size=9
    file = open('all_intervals.txt', 'rb')
    all_intervals_for_graph = pickle.load(file)
    file.close()
    delta_len=1
    max_bubblesize=4
    print(all_intervals_for_graph[4])
    DG,known_intervals,reads_for_isoforms,reads_list = generateGraphfromIntervals(all_intervals_for_graph, k_size,delta_len)
    print(known_intervals)
    check_graph_correctness(known_intervals,all_intervals_for_graph)
    print("#Nodes for DG: " + str(DG.number_of_nodes()) + " , #Edges for DG: " + str(DG.number_of_edges()))
    #edgelist = list(DG.edges.data())
    #print(edgelist)
    #simplifyGraph(DG, max_bubblesize, delta_len)
    #draw_Graph(DG)
    #generate_isoforms(DG,reads_for_isoforms)
    #DG.nodes(data=True)
    #print("Number of Nodes for DG:" + str(len(DG)))
    #nodelist = list(DG.nodes)

    #generate_isoforms(DG, reads_list)

    #print("number of edges in DG:" + str(DG.number_of_edges()))
    #for node in nodelist:
    #    print(node)

    #print("Number of Nodes for DG_old:" + str(len(DG_old)))
    #nodelist_old = list(DG_old.nodes)

    #print("number of edges in DG:" + str(DG_old.number_of_edges()))
    #for node in nodelist_old:
    #    print(node)
    #DG2 = generateSimpleGraphfromIntervals(all_intervals_for_graph)
    #add_Nodes(DG,args.delta_len)
    #simplifyGraph(DG,args.max_bubblesize,args.delta_len)
    # att = nx.get_node_attributes(DG, reads)
    # print("749,762 attributes: " + str(att))

    #nodelist = list(DG.nodes(data=True))
    #print(nodelist)
    #print("#Nodes for DG2: " + str(DG2.number_of_nodes()) + " , #Edges for DG: " + str(DG2.number_of_edges()))
    #print()
    #print()
    #draw_Graph(DG)
    # draw_Graph(DG2)
    # writes the graph in GraphML format into a file. Makes it easier to work with the graph later on
    #nx.write_graphml_lxml(DG, "outputgraph.graphml")
    # nx.write_graphml_lxml(DG2, "outputgraph2.graphml")
    #print("finding the reads, which make up the isoforms")

if __name__ == "__main__":
    main()