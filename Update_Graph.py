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
# Function to convert a list into a string to enable the writing of a graph (Taken from https://www.geeksforgeeks.org/python-program-to-convert-a-list-to-string/)
def listToString(s):
    # initialize an empty string
    str1 = " "

    # return string
    return (str1.join(str(s)))
def generateGraph_no_prev_info():
    print("Hello world")
"""function to add read information of the current interval to prior_read_infos, if they do not exist there yet
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
    # print("Node s is added to the graph.")
    DG.add_node("t",reads=reads_at_startend_list)
    # print("Node t is added to the graph.")
    # holds the r_id as key and a list of tuples as value: For identification of reads
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
    repetative_reads = {}
    repetative_infos={}
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    for r_id, intervals_for_read in all_intervals_for_graph.items():  # intervals_for_read holds all intervals which make up the solution for the WIS of a read
        print(r_id)
        #if not r_id in known_intervals:
        #    known_intervals[r_id] = []
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"
        previous_end=0
        # the name of each node is defined to be startminimizerpos , endminimizerpos
        # iterate over all intervals, which are in the solution of a certain read
        for inter in intervals_for_read:
            inter_infos={}
            info_tuple=(r_id,inter[0],inter[1])
            read_id = inter[3][slice(0, len(inter[3]),
                                     3)]  # recover the read id from the array of instances which was delivered with all_intervals_for_graph
            read_id_list= collections.Counter(read_id)
            if not all(value == 1 for value in read_id_list.values()):
                start_coord = inter[3][slice(1, len(inter[3]),
                                             3)]  # recover the start coordinate of an interval from the array of instances
                end_coord = inter[3][slice(2, len(inter[3]), 3)]
                for i,read in enumerate(read_id):
                    if not read in inter_infos:
                        inter_infos[read]=(start_coord[i],end_coord[i])
                    else:
                        rep_tuple=(read,start_coord[i],end_coord[i])
                        repetative_infos[rep_tuple]=[]
                for key,value in read_id_list.items():
                    if value>=2:
                        if not key in repetative_reads:
                            repetative_reads[key]=(r_id)


            #if repetative_reads:
            #    print(repetative_reads)
            #print("Inter")
            #print (inter)
            #read_id = inter[3][slice(0, len(inter[3]),3)]
            #read_count_list = {i:read_id.count(i) for i in read_id}
            #print(read_count_list)

            #access prior_read_infos, if any reads are already known
            if info_tuple in prior_read_infos:
                name=prior_read_infos[info_tuple]
            #if prior_read_infos[r_id-1]:
                #list comprehension magic to find out whether a read with start at inter[0] and end at inter[1] exists in prior_info

                        #witness=repetative_reads.get(r_id)
                        #print("Repetative node")
                        #print(known_intervals[r_id-1])
                        #print(known_tuples)
                #if such an element exists
                #if known_tuples:
                #    if len(known_tuples)>1:
                #        if r_id in repetative_reads:
                #            for i,occurrence in enumerate(known_tuples):
                #                print("Occurrence")
                #                print(occurrence)
                    #known_tuples is a list of tuples
                #    if len(known_tuples)<2:
                #        name = known_tuples[0][0]
                #    else:
                #        for i,tuple in enumerate(known_tuples):
                #            if i>1:
                #                name=tuple[0]
                #                old_name=known_tuples[i-1][0]
                #                start_of_cyc=known_intervals.index(old_name)
                #                old_cycle=known_intervals
                #                known_occurrences=[item for item in known_intervals[r_id-1] if item[1] == name]
                #                if not known_occurrences:
                #                    break

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
    print("Repetative Infos")
    print(repetative_infos)
    #print repetative reads
    print("Repetative reads")
    print(repetative_reads)
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