import _pickle as pickle
from sys import stdout
import networkx as nx
import itertools
import matplotlib.pyplot as plt
from SimplifyGraph import *

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

""" generates a networkx graph from the intervals given in all_intervals_for_graph.
# INPUT: all_intervals_for_graph: A dictonary holding lists of minimizer intervals.
"""
def generateGraphfromIntervals(all_intervals_for_graph, k):
    DG = nx.DiGraph()
    reads_at_startend_list = []
    for i in range(1,len(all_intervals_for_graph)+1):
        reads_at_startend_list.append(i)
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
    prior_read_infos = []
    for _ in itertools.repeat(None, len(all_intervals_for_graph)):
        prior_read_infos.append([])
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    for r_id, intervals_for_read in all_intervals_for_graph.items():  # intervals_for_read holds all intervals which make up the solution for the WIS of a read
        #if not r_id in known_intervals:
        #    known_intervals[r_id] = []
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"

        # the name of each node is defined to be startminimizerpos , endminimizerpos
        #if there are any infos about this read, save t
        if prior_read_infos[r_id - 1]:
            prior_info=prior_read_infos[r_id-1]
        # iterate over all intervals, which are in the solution of a certain read
        for inter in intervals_for_read:
            #access prior_read_infos, if any reads are already known
            if prior_read_infos[r_id-1]:
                #list comprehension magic to find out whether a read with start at inter[0] and end at inter[1] exists in prior_info
                known_tuples=[item for item in prior_info if item[1] == inter[0] and item[2]==inter[1]]
                #if such an element exists
                if known_tuples:
                        #known_tuples is a list of tuples, but should only contain 1 element
                        name = known_tuples[0][0]
                        #update the read information of node name
                        prev_nodelist=nodes_for_graph[name]
                        prev_nodelist.append((r_id, inter[0], inter[1]))
                        nodes_for_graph[name]=prev_nodelist
                        #only add a new edge if the edge was not present before
                        if not DG.has_edge(previous_node,name):
                            DG.add_edge(previous_node, name)
                        #keep known_intervals up to date
                        known_intervals[r_id-1].append((inter[0],name,inter[1]))

                # if the information for the interval was not yet found in a previous read (meaning the interval is new)
                else:
                        #print("else")
                        nodelist = []
                        # add a node into nodes_for_graph
                        name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                        #add the read information for the node
                        nodelist.append((r_id, inter[0], inter[1]))
                        nodes_for_graph[name] = nodelist
                        DG.add_node(name)
                        #keep known_intervas up to date
                        known_intervals[r_id - 1].append((inter[0], name, inter[1]))
                        #connect the node to the previous one
                        DG.add_edge(previous_node, name)
                        read_id = inter[3][slice(0, len(inter[3]),
                                     3)]  # recover the read id from the array of instances which was delivered with all_intervals_for_graph
                        start_coord = inter[3][slice(1, len(inter[3]),
                                         3)]  # recover the start coordinate of an interval from the array of instances
                        end_coord = inter[3][slice(2, len(inter[3]), 3)]
                        #iterate through the retreived information and store it in prior_read_infos only for subsequent reads
                        for i, r in enumerate(read_id):
                            if not r <= r_id:
                                start = start_coord[i] + k
                                end = end_coord[i]
                                tuple_for_data_structure = (name, start, end)
                                prior_read_infos[r-1].append(tuple_for_data_structure)
            #if the information for the interval was not yet found in a previous read(should only happen for read 1)
            else:
                #print("else")
                nodelist = []
                # add a node into nodes_for_graph
                name = str(inter[0]) + ", " + str(inter[1]) + ", " + str(r_id)
                # add the read information for the node
                nodelist.append((r_id, inter[0], inter[1]))
                nodes_for_graph[name] = nodelist
                DG.add_node(name)
                # keep known_intervas up to date
                known_intervals[r_id - 1].append((inter[0], name, inter[1]))
                # connect the node to the previous one
                DG.add_edge(previous_node, name)
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
                        tuple_for_data_structure = (name, start, end)
                        prior_read_infos[r - 1].append(tuple_for_data_structure)
            #set the previous node for the next iteration
            previous_node = name
        #add an edge from name to "t" as name was the last node in the read
        if not DG.has_edge(name, "t"):
        # print("I have no idea what I'm doing here!")
            DG.add_edge(name,"t")
        #delete the entry of prior_read_infos just to save some space (set to empty list to keep coordinates correct
        prior_read_infos[r_id-1]=[]
    #print("Nodes for graph")
        #print(nodes_for_graph)
    #print(nodes_for_graph)
    #set the node attributes to be nodes_for_graph, very convenient way of solving this
    nx.set_node_attributes(DG,nodes_for_graph,name="reads")
    #return DG and known_intervals as this makes some operations easier
    result = (DG, known_intervals)
    return result

def check_graph_correctness(known_intervals,all_intervals_for_graph):
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        correct_graph=True
        if not (len(intervals_for_read)==len(known_intervals[r_id-1])):
            correct_graph=False
            print("graph incorrect for read "+str(r_id))

        #print("All_intervals_for_graph at read "+str(r_id)+" has "+str(len(intervals_for_read))+" Intervals")
        #for inter in intervals_for_read:
        #        read_id = inter[3][slice(0, len(inter[3]),3)]
        #        number_of_reads_all_intervals=len(read_id)
        #print("Known_intervals at read " + str(r_id) + " has " + str(len(known_intervals[r_id-1])) + " Intervals")
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
    delta_len=3
    max_bubblesize=4
    DG,known_intervals = generateGraphfromIntervals(all_intervals_for_graph, k_size)
    check_graph_correctness(known_intervals,all_intervals_for_graph)
    print("#Nodes for DG: " + str(DG.number_of_nodes()) + " , #Edges for DG: " + str(DG.number_of_edges()))
    edgelist = list(DG.edges.data())
    print(edgelist)
    simplifyGraph(DG,max_bubblesize,delta_len)
    #DG.nodes(data=True)
    #print("Number of Nodes for DG:" + str(len(DG)))
    #nodelist = list(DG.nodes)

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