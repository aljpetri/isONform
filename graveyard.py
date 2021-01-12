# if prior_read_infos[r_id-1]:
# list comprehension magic to find out whether a read with start at inter[0] and end at inter[1] exists in prior_info

# witness=repetative_reads.get(r_id)
# print("Repetative node")
# print(known_intervals[r_id-1])
# print(known_tuples)
# if such an element exists
# if known_tuples:
#    if len(known_tuples)>1:
#        if r_id in repetative_reads:
#            for i,occurrence in enumerate(known_tuples):
#                print("Occurrence")
#                print(occurrence)
# known_tuples is a list of tuples
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

"""Preprocessing step during generateGraphfromIntervals:
Inside interval iteration"""
# inter_infos is a dict holding the infos about all interval occurances that already were found
            inter_infos = {}

            read_id = inter[3][slice(0, len(inter[3]),
                                     3)]  # recover the read id from the array of instances which was delivered with all_intervals_for_graph

            #count the occurences of the read ids and store it in read_id_list
            read_id_list= collections.Counter(read_id)

            #PREPROCESSING STEP: the interval information is searched for repeating read id's to identify repetative nodes in other reads
            if not all(value == 1 for value in read_id_list.values()):
                is_repetative=True
                start_coord = inter[3][slice(1, len(inter[3]),
                                                                3)]  # recover the start coordinate of an interval from the array of instances
                end_coord = inter[3][slice(2, len(inter[3]), 3)]# recover the end coordinate of an interval from the array of instances
                #iterate through all the reads mentioned in read_id
                for i,read in enumerate(read_id):
                    #if the read id has not yet been stored in the inter_infos
                    if not read in inter_infos:
                        #add dictionary entry with the new read id and the start and end positions
                        inter_infos[read]=(start_coord[i]+k,end_coord[i])
                    else:
                        first_occ=inter_infos.get(read)
                        rep_tuple=(start_coord[i]+k,end_coord[i])
                        repetative_infos[(read,first_occ,rep_tuple)]=[]
                #print(read_id_list)
                for key,value in read_id_list.items():
                    if value>=2:
                        if not key in witnesses:
                            #appoint this read id as the value of key repeating_r_id in repetative_reads
                            witnesses[key]=(r_id)




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