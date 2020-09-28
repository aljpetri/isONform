import networkx as nx
def generateGraphfromIntervals(intervals_to_correct):
    DG = nx.DiGraph()
    DG.add_node("s")
    DG.add_node("t")
    for inter in intervals_to_correct:
        name=str(inter[0])+", "+str(inter[1])
        DG.add_node(name)
        oldname=name
        DG.add_edge(oldname,name,weight=1)
        #print("Interval from "+str(inter[0])+" to "+str(inter[1]) +", supported by "+ str(inter[2])+" reads.")
    return DG
