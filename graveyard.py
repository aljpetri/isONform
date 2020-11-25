def merge_nodes(DG):
    edgesView=DG.edges.data()
    for ed_ge in edgesView:
        startnode = ed_ge[0]
        endnode = ed_ge[1]
        start_out_degree = DG.out_degree(startnode)
        end_out_degree = DG.out_degree(endnode)
        start_in_degree = DG.in_degree(startnode)
        end_in_degree = DG.in_degree(endnode)
        if(start_out_degree==end_in_degree==1and startnode!="s"and endnode!="t"):
            readinfos=nx.get_node_attributes(DG,'reads')
            start_infos=readinfos[startnode]
            end_infos=readinfos[endnode]
            merged_infos=[]

            start_inter = tuple(map(int, startnode.split(', ')))
            start_name=start_inter[0]
            id_name=start_inter[2]
            end_inter=tuple(map(int, endnode.split(', ')))
            end_name = end_inter[1]
            merged_name=str(start_name) + ", " + str(end_name)+", "+str(id_name)
            #print(merged_name)
            for i, start_tuple in enumerate(start_infos):
                end_tuple=end_infos[i]
                if(end_tuple[0]==start_tuple[0]):
                    merged_tuple=(start_tuple[0],start_tuple[1],end_tuple[2])
                    merged_infos.append(merged_tuple)
            in_edges=list(DG.in_edges(startnode))
            out_edges=list(DG.out_edges(endnode))
            #DG.remove_edge(startnode,endnode)
            #DG.remove_node(startnode)
            #DG.remove_node(endnode)
            #DG.add_node(merged_name,reads=merged_infos)
            for i_edge in in_edges:
                print("Incoming")
                print(i_edge)
                #DG.add_edge(i_edge[0],merged_name)
            print("In Edges for "+startnode)
            print(in_edges)
            #readinfos_end=nx.get_node_attributes(DG,endnode)
            print(start_infos)
            print(end_infos)
            print("Merging nodes "+startnode+" and "+endnode)
            DG=nx.contracted_nodes(DG, startnode,endnode)
    return DG