#! /usr/bin/env python
from __future__ import print_function
import argparse
from array import array
import itertools
import os
import shutil
import sys
import tempfile
import pickle
import operator
from collections import defaultdict, deque


from modules import help_functions, GraphGeneration, batch_merging_parallel, IsoformGeneration, SimplifyGraph

D = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)}


def write_batch(reads,outfolder,batch_pickle):
    this_batch_dict = {}
    for id, (acc, seq, qual) in reads.items():
        this_batch_dict[acc] = seq
    pickle_batch_file = open(os.path.join(outfolder, batch_pickle), 'wb')
    pickle.dump(this_batch_dict, pickle_batch_file)
    pickle_batch_file.close()



def get_read_lengths(all_reads):
    """Helper method which extracts the read lengths from all_reads. We will use those during the graph generation to appoint more meaningful information to the node't'
        INPUT: all_reads: dictionary which holds the overall read infos key: r_id, value tuple(readname, sequence, some info i currently don't care about)
        OUTPUT: readlen_dict: dictionary holding the read_id as key and the length of the read as value
    """
    readlen_dict = {}
    for r_id, infos in all_reads.items():
        seq = infos[1]
        seqlen = len(seq)
        readlen_dict[r_id] = seqlen
    return readlen_dict


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def remove_read_polyA_ends(seq, threshold_len, to_len):
    """
    Funtion that removes the polyA_ends from the reads. We remove all polyAtails longer than threshold_len by transforming them into polyA strings of length to_len.
    INPUT:  seq:            the sequence to be altered
            threshold_len:  the length threshold over which we alter the polyA sequences
            to_len:         the length of the poly_A tails after the alteration
    OUTPUT: seq_mod:        the sequence that has been modified by the function, i.e. the sequence with shortened polyA tails
    """
    #we only want to alter polyA sequences that are located in the end of the read->calculate a window length in which we perform the change
    end_length_window = min(len(seq)//2, 100)
    seq_list = [ seq[:-end_length_window] ]

    for ch, g in itertools.groupby(seq[-end_length_window:]):
        h_len = sum(1 for x in g)
        if h_len > threshold_len and (ch == "A" or ch == "T"):
            seq_list.append(ch*to_len)
        else:
            seq_list.append(ch*h_len)

    seq_mod = "".join([s for s in seq_list])
    return seq_mod




def rindex(lst, value):
    """
        Function that calculates the reverse index. We want to find the last (but still) smallest k-mer in the window not the first
        INPUT:  lst:            the window of kmers as a list
                value:          the smallest kmer
        OUTPUT: minimizer_pos:  the last position of kmer with value in lst (not the first as before)
        """
    return len(lst) - operator.indexOf(reversed(lst), value) - 1

def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    #save the window as a deque and instead of the sequence itself we use its hash value
    window_kmers = deque([hash(seq[i:i+k_size]) for i in range(w +1)])
    #the smallest kmer (or to be precise the smallest hash) we have found in the window
    curr_min = min(window_kmers)
    #we now want the last occurrence of the smallest kmer not the first anymore
    minimizer_pos = rindex(list(window_kmers), curr_min)
    #add the initial minimizer to minimizers
    minimizers = [ (seq[minimizer_pos: minimizer_pos+k_size], minimizer_pos) ] # get the last element if ties in window
    #iterate over the remaining read and find all minimizers therein
    for i in range(w+1, len(seq) - k_size):
        new_kmer = hash(seq[i:i+k_size])
        # updating window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous window's minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer and minimizer_pos < i - w:
            curr_min = min(window_kmers)
            minimizer_pos = rindex(list(window_kmers), curr_min) + i - w
            minimizers.append( (seq[minimizer_pos: minimizer_pos+k_size], minimizer_pos) ) # get the last element if ties in window

        # Previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (seq[i: i+k_size], i) )

    return minimizers

def get_kmer_maximizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = max(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updating window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer:
            curr_min = max(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_kmer > curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers


def get_minimizers_and_positions_compressed(reads, w, k, hash_fcn):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]

        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))

        if hash_fcn == "lex":
            minimizers = get_kmer_minimizers(seq_hpol_comp, k, w)
        elif hash_fcn == "rev_lex":
            minimizers = get_kmer_maximizers(seq_hpol_comp, k, w)
        # indicies we want to take quality values from to get quality string of homopolymer compressed read
        indices = [i for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2]
        indices.append(len(seq) - 1)
        positions_in_non_compressed_sring = [(m, indices[p]) for m, p in minimizers ]
        M[r_id] = positions_in_non_compressed_sring

    return M


def get_minimizers_and_positions(reads, w, k, hash_fcn):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        if hash_fcn == "lex":
            minimizers = get_kmer_minimizers(seq, k, w)
        elif hash_fcn == "rev_lex":
            minimizers = get_kmer_maximizers(seq, k, w)

        M[r_id] = minimizers

    return M


def get_minimizer_combinations_database(M, k, x_low, x_high,reads):
    M2 = defaultdict(lambda: defaultdict(lambda: array("I")))
    tmp_cnt = 0
    forbidden = 'A'*k
    for r_id in M:
        minimizers = M[r_id]
        for (m1,p1), m1_curr_spans in minimizers_comb_iterator(minimizers, k, x_low, x_high):
            for (m2, p2) in m1_curr_spans:
                if m2 == m1 == forbidden:
                    continue

                tmp_cnt +=1
                M2[m1][m2].append(r_id)
                M2[m1][m2].append(p1)
                M2[m1][m2].append(p2)

    print(tmp_cnt, "MINIMIZER COMBINATIONS GENERATED")

    avg_bundance = 0
    singleton_minimzer = 0
    cnt = 1
    for m1 in list(M2.keys()):
        for m2 in list(M2[m1].keys()):
            #if the minimizer_pair is represented at more than one occurrence
            if len(M2[m1][m2]) > 3:
                avg_bundance += len(M2[m1][m2])//3
                cnt += 1
            #the minimizer_combination is only once in the data therefore does not yield viable information for the graph later on
            else:
                del M2[m1][m2]
                singleton_minimzer += 1
            #we also want to filter out minimizer combinations if they are too abundant (more than 3 times per read)
            if len(M2[m1][m2]) // 3 > 3 * len(reads):
                del M2[m1][m2]

    print("Average abundance for non-unique minimizer-combs:", avg_bundance/float(cnt))
    print("Number of singleton minimizer combinations filtered out:", singleton_minimzer)

    return M2


def minimizers_comb_iterator(minimizers, k, x_low, x_high):
    for i, (m1, p1) in enumerate(minimizers[:-1]):
        m1_curr_spans = []
        for j, (m2, p2) in enumerate(minimizers[i+1:]):
            if x_low < p2 - p1 and p2 - p1 <= x_high:
                m1_curr_spans.append((m2, p2))
                # yield (m1,p1), (m2, p2)
            elif p2 - p1 > x_high:
                break
        yield (m1, p1), m1_curr_spans[::-1]


def fill_p2(p, all_intervals_sorted_by_finish):
    stop_to_max_j = {stop: j for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish) if start < stop}
    all_choord_to_max_j = []
    j_max = 0
    for i in range(0, all_intervals_sorted_by_finish[-1][1] + 1):
        if i in stop_to_max_j:
            j_max = stop_to_max_j[i]

        all_choord_to_max_j.append(j_max)
    for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish):
        j_max = all_choord_to_max_j[start]
        p.append(j_max)
    return p


def solve_WIS(all_intervals_sorted_by_finish):
    # Using notation from https://courses.cs.washington.edu/courses/cse521/13wi/slides/06dp-sched.pdf
    p = [None]
    fill_p2(p, all_intervals_sorted_by_finish)
    epsilon = 0.0001
    # w - 1 since the read interval itself is included in the instance
    v = [None] + [(w - 1)*(stop-start + epsilon) for (start, stop, w, _) in all_intervals_sorted_by_finish]
    OPT = [0]
    for j in range(1, len(all_intervals_sorted_by_finish) +1):
        OPT.append( max(v[j] + OPT[ p[j] ], OPT[j-1] ) )
    # Find solution
    opt_indicies = []
    j = len(all_intervals_sorted_by_finish)
    while j >= 0:
        if j == 0:
            break
        if v[j] + OPT[p[j]] > OPT[j-1]:
            opt_indicies.append(j - 1) # we have shifted all indices forward by one so we neew to reduce to j -1 because of indexing in python works
            j = p[j]
        else:
            j -= 1
    return opt_indicies


def get_intervals_to_correct(opt_indicies, all_intervals_sorted_by_finish):
    intervals_to_correct = []
    for j in opt_indicies:
        start, stop, weights, instance = all_intervals_sorted_by_finish[j]
        intervals_to_correct.append((start, stop, weights, instance))

    return intervals_to_correct


def batch(dictionary, size):
    batches = []
    sub_dict = {}
    for i, (acc, seq) in enumerate(dictionary.items()):
        if i > 0 and i % size == 0:
            batches.append(sub_dict)
            sub_dict = {}
            sub_dict[acc] = seq
        else:
            sub_dict[acc] = seq

    if i / size != 0:
        sub_dict[acc] = seq
        batches.append(sub_dict)
    elif len(dictionary) == 1:
        batches.append(sub_dict)

    return batches


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def add_items(seqs, r_id, p1, p2):
    seqs.append(r_id)
    seqs.append(p1)
    seqs.append(p2)


def find_most_supported_span(r_id, m1, p1, m1_curr_spans, minimizer_combinations_database, all_intervals, k_size, delta_len):
    """
    Funtion that detects the most supported spans in the database and adds them to all_intervals.
    INPUT:  r_id:           the id of the read we are workin on
            m1:             the length threshold over which we alter the polyA sequences
            p1:             the length of the poly_A tails after the alteration
            m1_curr_spans:
            minimizer_combinations_database:
            reads:
            all_intervals:  list holding the interval information (empty at this point)
            k_size:
            delta_len:

    OUTPUT: all_intervals:     modified list of all intervals
    """
    #print("Most_supp_span")
    for (m2, p2) in m1_curr_spans:
        relevant_reads = minimizer_combinations_database[m1][m2]
        if len(relevant_reads)//3 >= 3:
            seqs = array("I")
            seqs.append(r_id)
            seqs.append(p1)
            seqs.append(p2)
            for relevant_read_id, pos1, pos2 in grouper(relevant_reads, 3): #relevant_reads:
                if r_id == relevant_read_id:
                    continue
                elif abs((p2-p1)-(pos2-pos1)) < delta_len:
                    seqs.append(relevant_read_id)
                    seqs.append(pos1)
                    seqs.append(pos2)
            all_intervals.append((p1 + k_size, p2,  len(seqs)//3, seqs))


def main(args):
    #Todo: add timestamp
    print("ARGS",args)
    all_batch_reads_dict={}
    # start = time()
    if os.path.exists("mapping.txt"):
        os.remove("mapping.txt")
    outfolder = args.outfolder
    #sys.stdout = open(os.path.join(outfolder,"stdout.txt"), "w")
    # read the file and filter out polyA_ends(via remove_read_polyA_ends)
    all_reads = {i + 1: (acc, remove_read_polyA_ends(seq, 12, 1), qual) for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(args.fastq, 'r')))}
    max_seqs_to_spoa = args.max_seqs_to_spoa
    if len(all_reads) <= args.exact_instance_limit:
        args.exact = True
    if args.set_w_dynamically:
        args.w = args.k + min(7, int(len(all_reads) / 500))
    delta_iso_len_3 = args.delta_iso_len_3
    delta_iso_len_5 = args.delta_iso_len_5
    work_dir = tempfile.mkdtemp()
    print("Temporary workdirektory:", work_dir)

    k_size = args.k
    x_high = args.xmax
    x_low = args.xmin
    if args.parallel:
        filename = args.fastq.split("/")[-1]
        tmp_filename = filename.split("_")
        tmp_lastpart = tmp_filename[-1].split(".")
        p_batch_id = tmp_lastpart[0]
        skipfilename = "skip"+p_batch_id+".fa"

    for batch_id, reads in enumerate(batch(all_reads, args.max_seqs)):


        new_all_reads = {}
        if args.parallel:
            batch_pickle = str(p_batch_id) + "_batch"
        else:
            skipfilename="skip"+str(batch_id)+".fa"
            batch_pickle = str(batch_id) + "batch"
            skipfilename = "skip" + str(batch_id) + ".fa"
        skipfile = open(os.path.join(outfolder, skipfilename), 'w')
        write_batch(reads, outfolder, batch_pickle)
        skipfile=open(os.path.join(outfolder,skipfilename),'w')

        if args.set_w_dynamically:
            # Activates for 'small' clusters with less than 700 reads
            if len(reads) >= 100:
                w = min(args.w, args.k + (
                            len(reads) // 100 + 4))
            elif len(reads) == 1:
                for id, vals in reads.items():
                    (acc, seq, qual) = vals
                    skipfile.write(">{0}\n{1}\n".format(acc, seq))
                print("Not enough reads to work on!")
                continue
            else:
                w = args.k + 1 + len(reads) // 30
        else:
            w = args.w
        #print("Window used for batch:", w)
        iso_abundance = args.iso_abundance
        delta_len = args.delta_len
        graph_id = 1
        print("Working on {0} reads in a batch".format(len(reads)))
        hash_fcn = "lex"
        not_used=0
        #generate all minimizer combinations
        if args.compression:
            minimizer_database = get_minimizers_and_positions_compressed(reads, w, k_size, hash_fcn)
        else:
            minimizer_database = get_minimizers_and_positions(reads, w, k_size, hash_fcn)

        minimizer_combinations_database = get_minimizer_combinations_database(minimizer_database, k_size, x_low,
                                                                              x_high, reads)
        previously_corrected_regions = defaultdict(list)
        all_intervals_for_graph = {}
        for r_id in sorted(reads):
            print(r_id)
            read_min_comb = [((m1, p1), m1_curr_spans) for (m1, p1), m1_curr_spans in
                             minimizers_comb_iterator(minimizer_database[r_id], k_size, x_low, x_high)]

            if args.exact:
                previously_corrected_regions = defaultdict(list)

            read_previously_considered_positions = set(
                [tmp_pos for tmp_p1, tmp_p2, w_tmp, _ in previously_corrected_regions[r_id] for tmp_pos in
                 range(tmp_p1, tmp_p2)])

            if args.verbose:
                if read_previously_considered_positions:
                    eprint("not corrected:", [(p1_, p2_) for p1_, p2_ in
                                              zip(sorted(read_previously_considered_positions)[:-1],
                                                  sorted(read_previously_considered_positions)[1:]) if
                                              p2_ > p1_ + 1])
                else:
                    eprint("not corrected: entire read", )

            if previously_corrected_regions[r_id]:
                read_previously_considered_positions = set(
                    [tmp_pos for tmp_p1, tmp_p2, w_tmp, _ in previously_corrected_regions[r_id] for tmp_pos in
                     range(tmp_p1, tmp_p2)])
                group_id = 0
                pos_group = {}

                if len(read_previously_considered_positions) > 1:
                    sorted_corr_pos = sorted(read_previously_considered_positions)

                    for p1, p2 in zip(sorted_corr_pos[:-1], sorted_corr_pos[1:]):

                        if p2 > p1 + 1:
                            pos_group[p1] = group_id
                            group_id += 1
                            pos_group[p2] = group_id
                        else:
                            pos_group[p1] = group_id

                    if p2 == p1 + 1:
                        pos_group[p2] = group_id

            else:
                read_previously_considered_positions = set()
                pos_group = {}
            all_intervals = []
            prev_visited_intervals = []
            #print("it over read_min_comb")
            for (m1, p1), m1_curr_spans in read_min_comb:
                #print(m1,", ", p1)
                # If any position is not in range of current corrections: then correct, not just start and stop
                not_prev_corrected_spans = [(m2, p2) for (m2, p2) in m1_curr_spans if not (
                        p1 + k_size in read_previously_considered_positions and p2 - 1 in read_previously_considered_positions)]
                set_not_prev = set(not_prev_corrected_spans)
                not_prev_corrected_spans2 = [(m2, p2) for (m2, p2) in m1_curr_spans if
                                             (m2, p2) not in set_not_prev and (
                                                     p1 + k_size in pos_group and p2 - 1 in pos_group and pos_group[
                                                 p1 + k_size] != pos_group[p2 - 1])]
                not_prev_corrected_spans += not_prev_corrected_spans2

                if not_prev_corrected_spans:  # p1 + k_size not in read_previously_considered_positions:
                    find_most_supported_span(r_id, m1, p1, m1_curr_spans, minimizer_combinations_database,
                                              all_intervals, k_size, delta_len)
            # add prev_visited_intervals to intervals to consider
            all_intervals.extend(prev_visited_intervals)

            if previously_corrected_regions[r_id]:  # add previously corrected regions in to the solver

                all_intervals.extend(previously_corrected_regions[r_id])
                del previously_corrected_regions[r_id]

            if not all_intervals:
                not_used+=1
                if DEBUG:
                    eprint("Found no reads to work on")
                vals=reads[r_id]
                (acc, seq, qual) = vals
                skipfile.write(">{0}\n{1}\n".format(acc, seq))
                continue
            else:
                all_intervals.sort(key=lambda x: x[1])
                opt_indicies = solve_WIS(
                    all_intervals)  # solve Weighted Interval Scheduling here to find set of best non overlapping intervals to correct over
                intervals_to_correct = get_intervals_to_correct(opt_indicies[::-1], all_intervals)
                #if we have found intervals in our read: add it to all_intervals_for_graph and give it a graph_id (an internal id)
                if intervals_to_correct:
                    all_intervals_for_graph[graph_id] = intervals_to_correct
                    new_all_reads[graph_id] = reads[r_id]
                    graph_id += 1
                #the read has no intervals in common with other reads-> we add it into skipfile
                else:
                    not_used += 1
                    if DEBUG:
                        eprint("Found no reads to work on")
                    vals = reads[r_id]
                    (acc, seq, qual) = vals
                    skipfile.write(">{0}\n{1}\n".format(acc, seq))

        if not_used>0:
            print("Skipped ",not_used," reads due to not having high enough interval abundance")
        else:
            print("Working on all reads")
        print("Generating the graph")
        all_batch_reads_dict[batch_id] = new_all_reads
        read_len_dict = get_read_lengths(all_reads)
        #for key,value in all_intervals_for_graph.items():
        #    print(key,len(value))
        #print(all_intervals_for_graph)

        #profiler = Profiler()
        #profiler.start()
        #generate the graph from the intervals
        #TODO: add timestamp
        DG,  reads_for_isoforms = GraphGeneration.generateGraphfromIntervals(
            all_intervals_for_graph, k_size, delta_len, read_len_dict,new_all_reads)
        #profiler.stop()

        #profiler.print()
        #test for cyclicity of the graph - a status we cannot continue working on -> if cyclic we get an error
        #is_cyclic = SimplifyGraph.isCyclic(DG)
        #if is_cyclic:
        #    k_size+=1
        #    w+=1
        #    print("The graph has a cycle - critical error")
            #return -1
        #else:
        #    print("No cycle in graph")
        if DEBUG==True:
            print("BATCHID",batch_id)
            for id, value in all_batch_reads_dict.items():
                for other_id,other_val in value.items():
                   print(id,": ",other_id,":",other_val[0],"::",other_val[1])

        mode = args.slow
        #profiler = Profiler()
        #profiler.start()
        #the bubble popping step: We simplify the graph by linearizing all poppable bubbles
        SimplifyGraph.simplifyGraph(DG, new_all_reads, work_dir, k_size, delta_len, mode)
        #profiler.stop()

        #profiler.print()
        #TODO: add delta as user parameter possibly?
        delta = 0.15
        print("Starting to generate Isoforms")

        if args.parallel:
            batch_id = p_batch_id
        #profiler = Profiler()
        #profiler.start()
        #generation of isoforms from the graph structure
        IsoformGeneration.generate_isoforms(DG, new_all_reads, reads_for_isoforms, work_dir, outfolder, batch_id, delta, delta_len, delta_iso_len_3, delta_iso_len_5, max_seqs_to_spoa)
        #profiler.stop()

        print("Isoforms generated-Starting batch merging ")
    if not args.parallel:
            print("Merging the batches with linear strategy")
            #merges the predictions from different batches
            #batch_merging_parallel.join_back_via_batch_merging(args.outfolder, delta, args.delta_len, args.delta_iso_len_3, args.delta_iso_len_5,
            #                                                   args.max_seqs_to_spoa, args.iso_abundance)

    print("removing temporary workdir")
    #sys.stdout.close()
    shutil.rmtree(work_dir)

DEBUG=False

#TODO: add low_output (bool) as well as delta as user parameters and remove slow, merge_sub_isoforms_3, merge_sub_isoforms_5
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="De novo error correction of long-read transcriptome reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 0.3.3')
    parser.add_argument('--fastq', type=str, default=False, help='Path to input fastq file with reads')

    parser.add_argument('--k', type=int, default=20, help='Kmer size')
    parser.add_argument('--w', type=int, default=31, help='Window size')
    parser.add_argument('--xmin', type=int, default=18, help='Upper interval length')
    parser.add_argument('--xmax', type=int, default=80, help='Lower interval length')
    parser.add_argument('--T', type=float, default=0.1, help='Minimum fraction keeping substitution')
    parser.add_argument('--exact', action="store_true", help='Get exact solution for WIS for every read (recalculating weights for each read (much slower but slightly more accuracy,\
                                                                 not to be used for clusters with over ~500 reads)')
    parser.add_argument('--disable_numpy', action="store_true",
                        help='Do not require numpy to be installed, but this version is about 1.5x slower than with numpy.')
    parser.add_argument('--delta_len', type=int, default=5, help='Maximum length difference between two reads intervals for which they would still be merged')
    parser.add_argument('--max_seqs_to_spoa', type=int, default=200, help='Maximum number of seqs to spoa')
    parser.add_argument('--max_seqs', type=int, default=1000,
                        help='Maximum number of seqs to correct at a time (in case of large clusters).')
    parser.add_argument('--exact_instance_limit', type=int, default=0,
                        help='Activates slower exact mode for instance smaller than this limit')
    parser.add_argument('--set_w_dynamically', action="store_true",
                        help='Set w = k + max(2*k, floor(cluster_size/1000)).')
    parser.add_argument('--verbose', action="store_true", help='Print various developer stats.')

    parser.add_argument('--compression', action="store_true", help='Use homopolymer compressed reads. (Deprecated, because we will have fewer \
                                                                        minmimizer combinations to span regions in homopolymenr dense regions. Solution \
                                                                        could be to adjust upper interval length dynamically to guarantee a certain number of spanning intervals.')
    parser.add_argument('--outfolder', type=str, default=None,
                        help='The outfolder of isONform, into which the algorithm will write the results.')
    parser.add_argument('--iso_abundance', type=int,default=5, help='Cutoff parameter: abundance of reads that have to support an isoform to show in results')
    parser.add_argument('--delta_iso_len_3', type=int, default=30,
                        help='Cutoff parameter: maximum length difference at 3prime end, for which subisoforms are still merged into longer isoforms')
    parser.add_argument('--delta_iso_len_5', type=int, default=50,
                        help='Cutoff parameter: maximum length difference at 5prime end, for which subisoforms are still merged into longer isoforms')
    parser.add_argument('--parallel',type=bool,default=False,help='Indicates whether we run the parallelization wrapper script')
    parser.add_argument('--slow',action="store_true", help='EXPERIMENTAL PARAMETER: has high repercussions on run time use the slow mode for the simplification of the graph (bubble popping), slow mode: every bubble gets popped.')
    args = parser.parse_args()

    if args.xmin < 2 * args.k:
        args.xmin = 2 * args.k
        eprint("xmin set to {0}".format(args.xmin))

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    os.environ['PYTHONHASHSEED'] = '0'
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    if 100 < args.w or args.w < args.k:
        eprint('Please specify a window of size larger or equal to k, and smaller than 100.')
        sys.exit(1)

    main(args)
