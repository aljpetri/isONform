import _pickle as pickle
from sys import stdout
import networkx as nx
import itertools
import matplotlib.pyplot as plt
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
# generates a networkx graph from the intervals given in all_intervals_for_graph.
# INPUT: all_intervals_for_graph: A dictonary holding lists of minimizer intervals.
def generateGraphfromIntervals(all_intervals_for_graph, k):
    DG = nx.DiGraph()
    reads_at_startend_list = []
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
    # iterate through the different reads stored in all_intervals_for_graph. For each read one path is built up from source to sink if the nodes needed for that are not already present
    for r_id, intervals_for_read in all_intervals_for_graph.items():  # intervals_for_read holds all intervals which make up the solution for the WIS of a read
        #if not r_id in known_intervals:
        #    known_intervals[r_id] = []
        # set previous_node to be s. This is the node all subsequent nodes are going to have edges with
        previous_node = "s"

        # the name of each node is defined to be startminimizerpos , endminimizerpos
        #liste gives a list of all the intervals already present in the read which currently is added to the graph
        liste = known_intervals[r_id-1]
        print("liste:")
        print(liste)
        # iterate over all intervals, which are in the solution of a certain read
        for inter in intervals_for_read:
            #res denotes a list containing all known intervals having startpos
            # do some lambda magic to find the tuples containing inter[0] as first element
            res= []
            for tup in liste:
                if tup[0]==inter[0]:
                    res.append(tup)
            #res = list(filter(lambda x[0]: inter[0] in x, liste))
            #print("inter[0]")
            #print(inter[0])
            #print("res:")
            #print(res)
            if len(res) == 0:
                # only add interval node, if the interval is not already known (meaning it is present in known_intervals)

                # if such an element exists, we do nothing, else we add the interval into the graph and add the information of the interval for each read
                # if not res:
                del res
                # print("adding node "+name)
                name = str(inter[0]) + ", " + str(inter[1])# + ", " + str(r_id)

                read_id = inter[3][slice(0, len(inter[3]),
                                         3)]  # recover the read id from the array of instances which was delivered with all_intervals_for_graph
                start_coord = inter[3][slice(1, len(inter[3]),
                                             3)]  # recover the start coordinate of an interval from the array of instances
                end_coord = inter[3][slice(2, len(inter[3]), 3)]
                reads_at_node_list = []
                #reads_at_node_list.append((r_id,inter[0],inter[1]))
                # adds the instance to each read's overview list
                for i, r in enumerate(read_id):
                    # As not all keys may have been added to the dictionary previously, we add the rest of keys now
                    #if not r in known_intervals:
                    #    known_intervals[r] = []
                    # while the start pos stored in inter[0] has the right position the start positions in the list of instances are at pos-k
                    coord = start_coord[i] + k
                    end = end_coord[i]
                    # generate a tuple having the least amount of information needed to properly build up the graph, denoting one minimizer interval and add it to known_intervals
                    tuple = (coord, name, end)
                    known_intervals[r-1].append(tuple)
                    # print("ReadID "+str(r)+" from "+str(start_coord[i])+" to "+str(end_coord[i]))
                    # DONE: Try to find out what to do with all the edges to be added ->main idea: add edge as soon as node was added.
                    # if node is new: Add edges from previous intervals (Where to get this info?), else:
                    reads_at_node_list.append((r,coord,end))
                if not DG.has_node(name):
                    reads_at_node_string = listToString(reads_at_node_list)
                    DG.add_node(name, reads=reads_at_node_list)
                    # print("Node " + name + " is added to the graph.")
                # add edge between current node and previous node
                if DG.has_node(name) and DG.has_node(previous_node):
                    DG.add_edge(previous_node, name)  # weight=weightval
            # if the node was already added to the graph, we still have to find out, whether more edges need to be added to the graph
            else:
                #name = str(inter[0]) + ", " + str(inter[1])# + ", " + str(r_id)
                #self loops occur, as tup is the same twice(this is due the structure of known intervals i suppose)
                tup = res[0]
                name = tup[1]
                # print("Node "+name+" already present in the graph.")
                # see if previous node and this node are connected by an edge. If not add an edge dedicated to fulfill this need
                if not DG.has_edge(previous_node, name):
                    #Here wrongly self-loops are added to the graph, not sure yet why that happens
                    if DG.has_node(name) and DG.has_node(previous_node):
                        #rule out self loops as they mess up the outcoming graph
                        if not(previous_node == name):
                            #print("wanting to add edge from " +previous_node+" to node "+ name)
                            DG.add_edge(previous_node, name)
            # whatever we did before, we have to set previous_node to the node we looked at in this iteration to keep going
            previous_node = name
            # weightval = r_id
        # the last node of a path is connected to the sink node t
        # if DG.has_node("t") and DG.has_node(previous_node):
        DG.add_edge(previous_node, "t")  # weight=weightval

    #for key, value in known_intervals.items():
    #    print(key, value)
    result = (DG, known_intervals)
    print("Known_intervals")
    print(known_intervals)
    print("Known_intervals done")
    with open('known_intervals.txt', 'wb') as file:
        file.write(pickle.dumps(known_intervals))
    return result
def check_graph_correctness(known_intervals,all_intervals_for_graph):
    for r_id, intervals_for_read in all_intervals_for_graph.items():
        print("All_intervals_for_graph at read "+str(r_id)+" has "+str(len(intervals_for_read))+" Intervals")
        #for inter in intervals_for_read:
        #        read_id = inter[3][slice(0, len(inter[3]),3)]
        #        number_of_reads_all_intervals=len(read_id)
        print("Known_intervals at read " + str(r_id) + " has " + str(len(known_intervals[r_id-1])) + " Intervals")

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
    DG, known_intervals = generateGraphfromIntervals(all_intervals_for_graph, k_size)
    check_graph_correctness(known_intervals,all_intervals_for_graph)
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
    DG2 = generateSimpleGraphfromIntervals(all_intervals_for_graph)
    #add_Nodes(DG,args.delta_len)
    #simplifyGraph(DG,args.max_bubblesize,args.delta_len)
    # att = nx.get_node_attributes(DG, reads)
    # print("749,762 attributes: " + str(att))
    draw_Graph(DG)
    # draw_Graph(DG2)
    # writes the graph in GraphML format into a file. Makes it easier to work with the graph later on
    #nx.write_graphml_lxml(DG, "outputgraph.graphml")
    # nx.write_graphml_lxml(DG2, "outputgraph2.graphml")
    print("finding the reads, which make up the isoforms")

if __name__ == "__main__":
    main()