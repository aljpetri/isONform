import _pickle as pickle

"this script can be used to load known_intervals.txt, which contains the dictionary with the same name generated in main. " \
"The script finds all Isoforms by comparing all reads' minimizer intervals"
def main():
    file = open('known_intervals.txt', 'rb')
    known_intervals = pickle.load(file)
    #print(type(known_intervals))
    #print(known_intervals)
    #the reads which can be deleted, as they are exactly the same as others.

    #known_intervals2=known_intervals
    #sort the tuples by interval start positions.
    for r_ids,intervals in known_intervals.items():
        #print(type(intervals))
        known_intervals[r_ids]=sorted(intervals, key=lambda x: x[0])
    for key, value in known_intervals.items():
        print(key,value)
    node_ids_all_reads=[]
    #find all reads which are completely equal and pop the one with the higher id from the set
    for r_ids,intervals in known_intervals.items():
        nodeids=[x[1] for x in intervals]
        node_ids_all_reads.append(nodeids)
        print(r_ids, nodeids)
    print("done")
    for i in  range(0,len(node_ids_all_reads)-1):
        poppedreads=[]
        poppedreads.append(i)
        for j in range(i+1,len(node_ids_all_reads)):
            if node_ids_all_reads[i]==node_ids_all_reads[j]:
                read_to_pop=j+1
                if read_to_pop in known_intervals.keys():
                    known_intervals.pop(read_to_pop)
                    print("Deleting read "+str(read_to_pop)+" from known_intervals")
                    poppedreads.append(read_to_pop)
        node_ids_all_reads.append(poppedreads)
        #for r_ids2,intervals2 in known_intervals.items():
            #do not pop if read is only equal to itself
         #   if not r_ids==r_ids2:
                #do only pop if all intervals of the reads are equal
                #change here to only look at the name(id) and not start/stop anymore
          #      if intervals==intervals2:
          #          popentries=(r_ids,r_ids2)
          #          popitems.append(popentries)
          #          print("deleted read "+str(r_ids2)+"from known_intervals as equal to read "+str(r_ids))
    #for popits in popitems:
    #    if popits[0]<popits[1]:
    #        r_id=popits[1]
            #print(type(r_id))
    #        if r_id in known_intervals.keys():
    #            known_intervals.pop(r_id)
    for key, value in known_intervals.items():
        print(key,value)
    print("And now for the isoforms")
    for equalreads in node_ids_all_reads:
        print(equalreads)
    #for mainid,otherids in isoforms_by_reads.items():
    #    print("Read "+str(mainid) +" is equal to the following reads:")
    #    print(','.join(str(x) for x in otherids))
    #for popits in popitems:
        #if popits[0] < popits[1]:
    #        r_id = popits[1]
            #print(type(r_id))
    #        if r_id in known_intervals2.keys():
    #            known_intervals2.pop(r_id)
    print("Single occurance reads (no equals found):")
    #for key, value in known_intervals2.items():
    #    print(key, value)
if __name__ == '__main__':
    main()