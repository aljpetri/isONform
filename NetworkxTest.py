
import networkx as nx
import matplotlib.pyplot as plt
def generateGraphfromIntervals(intervals_to_correct):
    DG = nx.DiGraph()
    DG.add_node("s")
    DG.add_node("t")
    #DG.add_node("u")
    #DG.add_edge("s","t",weight=1)
    #DG.add_edge("s", "u", weight=1)
    #DG.add_edges_from([("s", "t",weight=1),("s","u",weight=4)])
    oldname="s"
    for i,inter in enumerate(intervals_to_correct):
        name=str(inter[0])+", "+str(inter[1])
        DG.add_node(name)
        #weightval=str(inter[2])
        DG.add_edge(oldname,name, weight=weightval)
        oldname = name
        #print("i"+str(i))
        #print("Interval from "+str(inter[0])+" to "+str(inter[1]) +", supported by "+ str(inter[2])+" reads.")
    return DG
def draw_Graph(DG):
    pos = nx.planar_layout(DG)
    nx.draw_networkx(DG, pos, font_weight='bold')
    labels = nx.get_edge_attributes(DG, 'weight')
    nx.draw_networkx_edge_labels(DG,pos, edge_labels=labels)
    plt.show()
def main():
    #intervals_to_correct=[(10,47,62),(47,67,20)]
    #DG=generateGraphfromIntervals(intervals_to_correct)
    DG=nx.read_graphml("outputgraph.graphml")
    print(len(DG))
    print(DG.number_of_edges())
    draw_Graph(DG)
    for path in nx.all_simple_paths(DG,"s","t"):
        print(path)

    #G = nx.Graph()
    #i = 1
    #G.add_node(i, pos=(i, i))
    #G.add_node("2", pos=("2", 2))
    #G.add_node(3, pos=(1, 0))
    #G.add_edge(1, 2, weight=0.5)
    #G.add_edge(1, 3, weight=9.8)
    #pos = nx.get_node_attributes(G, 'pos')
    #nx.draw(G, pos)
    #labels = nx.get_edge_attributes(G, 'weight')
    #nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    #plt.show()
if __name__ == '__main__':
    main()