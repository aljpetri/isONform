#! /usr/bin/env python
"""
This file is a wrapper file to run isONform in parallel (parallelization over batches).
The script was taken from isoncorrect(https://github.com/ksahlin/isoncorrect/blob/master/run_isoncorrect)
by Kristoffer Sahlin and changed by Alexander Petri to be usable with the isONform code base.

"""
# ! /usr/bin/env python

from __future__ import print_function
import os, sys
import argparse
from sys import stdout
import shutil
import subprocess
import glob
import tempfile
import errno
from time import time
from pathlib import Path
import itertools
from Parallelization_side_functions import *
import signal
from multiprocessing import Pool
import multiprocessing as mp
from batch_merging_parallel import *

# import math
import re

def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT
                           ).communicate()[0]
    # print(int(out.split()[0]))
    return int(out.split()[0])


def isONform(data):
    isONform_location, read_fastq_file, outfolder, batch_id, isONform_algorithm_params = data[0], data[1], data[
        2], data[3], data[4]
    mkdir_p(outfolder)
    print("OUT",outfolder)
    isONform_exec = os.path.join(isONform_location, "main.py")
    isONform_error_file = os.path.join(outfolder, "stderr.txt")
    with open(isONform_error_file, "w") as error_file:
        print('Running isONform batch_id:{0}...'.format(batch_id), end=' ')
        stdout.flush()

        #racon_flag = "--use_racon" if isONform_algorithm_params["use_racon"] else ''
        # randstrobes_flag = "--randstrobes" if isONform_algorithm_params["randstrobes"] else ''
        #dyn_flag = "--set_w_dynamically" if isONform_algorithm_params["set_w_dynamically"] else ''
        # null = open("/dev/null", "w")
        isONform_out_file = open(os.path.join(outfolder, "stdout.txt"), "w")
        subprocess.check_call(
                ["python", isONform_exec, "--fastq", read_fastq_file, "--outfolder", outfolder,
                 #"--exact_instance_limit", str(isONform_algorithm_params["exact_instance_limit"]),
                 #"--max_seqs", str(isONform_algorithm_params["max_seqs"]),
                 "--k", str(isONform_algorithm_params["k"]), "--w", str(isONform_algorithm_params["w"]),
                 "--xmin", str(isONform_algorithm_params["xmin"]), "--xmax",
                 str(isONform_algorithm_params["xmax"]),
                 "--exact", "--parallel", "True"
                 #"--T", str(isONform_algorithm_params["T"])
                 ], stderr=error_file, stdout=isONform_out_file)

        print('Done with batch_id:{0}.'.format(batch_id))
        stdout.flush()
    error_file.close()
    isONform_out_file.close()
    return batch_id


def splitfile(indir, tmp_outdir, fname, chunksize,cl_id,ext):
    # from https://stackoverflow.com/a/27641636/2060202
    # fpath, fname = os.path.split(infilepath)
    #cl_id, ext = fname.rsplit('.',1)
    infilepath = os.path.join(indir, fname)
    #infilepath=indir
    # print(fpath, cl_id, ext)
    print("now at splitfile")
    print(indir, tmp_outdir, cl_id, ext)

    i = 0
    written = False
    with open(infilepath) as infile:
        while True:
            outfilepath = os.path.join(tmp_outdir, '{0}_{1}.{2}'.format(cl_id, i, ext) ) #"{}_{}.{}".format(foutpath, fname, i, ext)
            print(outfilepath)
            with open(outfilepath, 'w') as outfile:
                for line in (infile.readline() for _ in range(chunksize)):
                    outfile.write(line)
                written = bool(line)
            # print(os.stat(outfilepath).st_size == 0)
            if os.stat(outfilepath).st_size == 0: # Corner case: Original file is even multiple of max_seqs, hence the last file becomes empty. Remove this
                os.remove(outfilepath)
            if not written:
                break
            i += 1






def symlink_force(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            if not os.path.exists(os.readlink(link_name)):
                print('path %s is a broken symlink' % link_name)
                os.remove(link_name)
                symlink_force(target,link_name)
        else:
            raise e


def split_cluster_in_batches(indir, outdir, tmp_work_dir, max_seqs):
    # create a modified indir
    tmp_work_dir = os.path.join(tmp_work_dir, 'split_in_batches')
    # print(indir)
    mkdir_p(tmp_work_dir)
    smaller_than_max_seqs = False
    # print(sorted(os.listdir(indir), key=lambda x: int(x.split('.')[0])) )
    # sys.exit()
    # add split fiels to this indir
    #for dirName, subdirList, fileList in os.walk(indir):
        #print('Found directory: %s' % dirName)
    #    for fastq_file in fileList:
    #        print('\t%s' % fastq_file)
    pat=Path(indir)
    file_list=list(pat.rglob('*.fastq'))
    print("FLIST",file_list)
    for filepath in file_list:
                print("FPATH",filepath)
                old_fastq_file=str(filepath.resolve())
                path_split=old_fastq_file.split("/")
                folder=path_split[-2]
                print(folder)
                fastq_file=path_split[-1]
                if not folder=="Analysis":
                    print("FQFile",fastq_file)
                    print(indir)
                    #print("PS",path_split)
                    cl_id=path_split[-2]
                    print("CLID",cl_id)
                    #print("FQ",fastq_file)
                    #print(type(fastq_file))
                    #if we have more lines than max_seqs
                    new_indir=os.path.join(indir,folder)
                    print(new_indir)
                    if not smaller_than_max_seqs:
                        num_lines = sum(1 for line in open(os.path.join(new_indir, fastq_file)))
                        print("Number Lines",fastq_file, num_lines)
                        #we reset smaller_than_max_seqs as we now want to see if we really have more than max_seqs reads
                        smaller_than_max_seqs = False if num_lines > 4 * max_seqs else True
                    else:
                        smaller_than_max_seqs = True

                    if not smaller_than_max_seqs:
                        print("Splitting",filepath)
                        ext = fastq_file.rsplit('.', 1)[1]
                        splitfile(new_indir, tmp_work_dir, fastq_file, 4 * max_seqs,cl_id,ext)  # is fastq file
                    else:
                        ext = fastq_file.rsplit('.', 1)[1]
                        #print(fastq_file, "symlinking instead")
                        symlink_force(filepath, os.path.join(tmp_work_dir, '{0}_{1}.{2}'.format(cl_id, 0, ext)))
    return tmp_work_dir


def join_back_corrected_batches_into_cluster(tmp_work_dir, outdir, split_mod, residual):
    print(outdir, tmp_work_dir)
    unique_cl_ids = set()
    for file in os.listdir(tmp_work_dir):
        file = os.fsdecode(file)
        # print(file)
        fname = file.split('_')
        if len(fname) == 2:
            cl_id, batch_id = fname[0], fname[1]  # file.split('_')

            if int(cl_id) % split_mod != residual:
                # print('skipping {0} because args.split_mod:{1} and args.residual:{2} set.'.format(cl_id, args.split_mod, args.residual))
                continue
            unique_cl_ids.add(cl_id)

    for cl_id in unique_cl_ids:
        out_pattern = os.path.join(outdir, cl_id)
        # print(type(tmp_work_dir), type(cl_id))
        batches_pattern = os.path.join(os.fsdecode(outdir), cl_id + '_*')
        # print("joining all", out_pattern, "from", batches_pattern)
        mkdir_p(out_pattern)

        error_file = open(os.path.join(out_pattern, 'cat.stderr'), 'w')
        outfilename = os.path.join(out_pattern, 'isoforms.fastq')
        # print("into outfile", outfilename)

        with open(outfilename, 'wb') as outfile:
            for batch_id in sorted(glob.glob(batches_pattern)):
                # print(batch_id)
                filename = os.path.join(batch_id, 'corrected_reads.fastq')
                if filename == outfilename:
                    # don't want to copy the output into the output
                    continue
                with open(filename, 'rb') as readfile:
                    shutil.copyfileobj(readfile, outfile)
                # print('Removing', batch_id)
                shutil.rmtree(batch_id)
#TODO: Finish implementation of this method and replace call of other merge_function by call to this
"""def join_back_via_batch_merging(tmp_work_dir, outdir, split_mod, residual):
    print("Batch Merging")
    print(outdir, tmp_work_dir)
    unique_cl_ids = set()
    for file in os.listdir(tmp_work_dir):
        file = os.fsdecode(file)
        # print(file)
        fname = file.split('_')
        if len(fname) == 2:
            cl_id, batch_id = fname[0], fname[1]  # file.split('_')

            unique_cl_ids.add(cl_id)
            print(unique_cl_ids)
            for cl_id in unique_cl_ids:
                out_pattern = os.path.join(outdir, cl_id)
                # print(type(tmp_work_dir), type(cl_id))
                batches_pattern = os.path.join(os.fsdecode(outdir), cl_id + '_*')
                # print("joining all", out_pattern, "from", batches_pattern)
                mkdir_p(out_pattern)

                error_file = open(os.path.join(out_pattern, 'cat.stderr'), 'w')
                outfilename = os.path.join(out_pattern, 'isoforms.fastq')
                # print("into outfile", outfilename)"""
def main(args):
    directory = args.fastq_folder  # os.fsencode(args.fastq_folder)
    print(directory)
    isONform_location = os.path.dirname(os.path.realpath(__file__))
    if args.split_wrt_batches:
        tmp_work_dir = tempfile.mkdtemp()
        print("Temporary workdirektory:", tmp_work_dir)
        split_tmp_directory = split_cluster_in_batches(directory, args.outfolder, tmp_work_dir, args.max_seqs)
        split_directory = os.fsencode(split_tmp_directory)
        print("SplitDIR",split_directory)
    else:
        split_directory = os.fsencode(directory)

    instances = []
    for file_ in os.listdir(split_directory):
        print(file_)
        read_fastq_file = os.fsdecode(file_)
        if read_fastq_file.endswith(".fastq"):
            print("True")
            tmp_id= read_fastq_file.split(".")[0]
            snd_tmp_id=tmp_id.split("_")
            cl_id = snd_tmp_id[0]
            batch_id=snd_tmp_id[1]
            outfolder = os.path.join(args.outfolder, cl_id)
            print(batch_id,cl_id)
            print(outfolder)
            #if int(cl_id) % args.split_mod != args.residual:
            #    print('skipping {0} because args.split_mod:{1} and args.residual:{2} set.'.format(batch_id,
            #                                                                                      args.split_mod,
            #                                                                                      args.residual))
            #    continue
            # print(batch_id, outfolder, read_fastq_file, split_directory)
            fastq_file_path = os.path.join(os.fsdecode(split_directory), read_fastq_file)
            print(fastq_file_path)
            compute = True
            if args.keep_old:
                candidate_corrected_file = os.path.join(outfolder, "isoforms.fastq")
                if os.path.isfile(candidate_corrected_file):
                    if wccount(candidate_corrected_file) == wccount(fastq_file_path):
                        print("already computed cluster and complete file", batch_id)
                        compute = False

            if compute:
                print("computing")
                isONform_algorithm_params = {"set_w_dynamically": args.set_w_dynamically,
                                                "exact_instance_limit": args.exact_instance_limit,
                                                "delta_len": args.delta_len,
                                                "k": args.k, "w": args.w, "xmin": args.xmin, "xmax": args.xmax,
                                                "T": args.T, "max_seqs": args.max_seqs, "use_racon": args.use_racon,"parallel": True}
                instances.append(
                    (isONform_location, fastq_file_path, outfolder, batch_id, isONform_algorithm_params))
            # else:
            #     isONform_algorithm_params = {  "set_w_dynamically" : args.set_w_dynamically, "exact_instance_limit" : args.exact_instance_limit, "k": args.k, "w" : args.w, "xmin" : args.xmin, "xmax" :  args.xmax, "T" : args.T }
            #     instances.append((isONform_location, fastq_file_path, outfolder, int(batch_id), isONform_algorithm_params) )
            # print(os.path.join(split_directory, read_fastq_file))
            # continue
        else:
            continue

    # # sys.exit()
    # if args.split_wrt_batches:
    #     instances.sort(key = lambda x: int(x[3].split('_')[0])) # sorting in order of cluster labels
    # else:
    #     instances.sort(key = lambda x: int(x[3])) # sorting in order of cluster size!

    instances.sort(key=lambda x: x[3])  # sorting on batch ids as strings
    print("Printing instances")
    for t in instances:
        print(t)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    mp.set_start_method('spawn')
    print(mp.get_context())
    print("Environment set:", mp.get_context())
    print("Using {0} cores.".format(args.nr_cores))
    start_multi = time()
    pool = Pool(processes=int(args.nr_cores))
    try:
        # res = pool.map_async(isONform, instances, chunksize=1)
        # score_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        start = time()
        for x in pool.imap_unordered(isONform, instances):
            print("{} (Time elapsed: {}s)".format(x, int(time() - start)))
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        sys.exit()
    else:
        pool.close()
    pool.join()

    print("Time elapesd multiprocessing:", time() - start_multi)

    if args.split_wrt_batches:
        file_handling = time()
        #join_back_corrected_batches_into_cluster(split_directory, args.outfolder, args.split_mod, args.residual)
        #join_back_via_batch_merging(split_directory, args.outfolder, args.split_mod, args.residual)
        join_back_via_batch_merging(tmp_work_dir, args.outfolder, args.split_mod, args.residual, args.delta, args.delta_len, args.merge_sub_isoforms_3,
                                    args.merge_sub_isoforms_5, args.delta_iso_len_3, args.delta_iso_len_5, args.max_seqs_to_spoa,args.iso_abundance)
        shutil.rmtree(split_directory)
        print("Joined back batched files in:", time() - file_handling)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="De novo clustering of long-read transcriptome reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.8')
    parser.add_argument('--fastq_folder', type=str, default=False,
                        help='Path to input fastq folder with reads in clusters')
    parser.add_argument('--t', dest="nr_cores", type=int, default=8, help='Number of cores allocated for clustering')
    parser.add_argument('--k', type=int, default=9, help='Kmer size')
    parser.add_argument('--w', type=int, default=20, help='Window size')
    parser.add_argument('--xmin', type=int, default=18, help='Lower interval length')
    parser.add_argument('--xmax', type=int, default=80, help='Upper interval length')
    parser.add_argument('--T', type=float, default=0.1, help='Minimum fraction keeping substitution')
    parser.add_argument('--exact_instance_limit', type=int, default=50,
                        help='Do exact correction for clusters under this threshold')
    # parser.add_argument('--w_equal_k_limit', type=int, default=100,  help='Do not recompute previous results')
    parser.add_argument('--keep_old', action="store_true",
                        help='Do not recompute previous results if corrected_reads.fq is found and has the smae number of reads as input file (i.e., is complete).')
    parser.add_argument('--set_w_dynamically', action="store_true",
                        help='Set w = k + max(2*k, floor(cluster_size/1000)).')
    parser.add_argument('--max_seqs', type=int, default=1000,
                        help='Maximum number of seqs to correct at a time (in case of large clusters).')
    parser.add_argument('--use_racon', action="store_true",
                        help='Use racon to polish consensus after spoa (more time consuming but higher accuracy).')

    parser.add_argument('--split_mod', type=int, default=1, help='Splits cluster ids in n (default=1) partitions by computing residual of cluster_id divided by n.\
                                                                    this parameter needs to be combined with  --residual to take effect.')
    parser.add_argument('--residual', type=int, default=0,
                        help='Run isONform on cluster ids with residual (default 0) of cluster_id divided by --split_mod. ')

    # parser.add_argument('--exact', action="store_true", help='Get exact solution for WIS for evary read (recalculating weights for each read (much slower but slightly more accuracy,\
    #                                                              not to be used for clusters with over ~500 reads)')

    parser.add_argument('--split_wrt_batches', action="store_true",
                        help='Process reads per batch (of max_seqs sequences) instead of per cluster. Significantly decrease runtime when few very large clusters are less than the number of cores used.')

    parser.add_argument('--outfolder', type=str, default=None, help='Outfolder with all corrected reads.')
    parser.add_argument('--randstrobes', action="store_true", help='EXPERIMENTAL PARAMETER: isONform uses paired minimizers (described in isONform paper). This experimental option\
                                                                 uses randstrobes instead of paired minimizers to find shared regions. Randstrobes \
                                                                 reduces memory footprint substantially (and runtime) with only slight increase in post correction quality.')

    parser.add_argument('--layers', type=int, default=argparse.SUPPRESS, help='EXPERIMENTAL PARAMETER: Active when --randstrobes specified.\
                                                                How many "layers" with randstrobes we want per sequence to sample.\
                                                               More layers gives more accureate results but is more memory consuming and slower.\
                                                               It is not reccomended to specify more than 5. ')
    parser.add_argument('--set_layers_manually', action="store_true", help='EXPERIMENTAL PARAMETER: By default isONform sets layers = 1 if nr seqs in batch to be corrected is >= 1000, else layers = 2.\
                                                                            This command will manually pick the number of layers specified with the --layers parameter.')
    parser.add_argument('--delta_len', type=int, default=3,
                        help='Maximum length difference between two reads intervals for which they would still be merged')
    parser.add_argument('--delta',type=float,default=0.1, help='diversity rate used to compare sequences')
    parser.add_argument('--max_seqs_to_spoa', type=int, default=200, help='Maximum number of seqs to spoa')

    parser.add_argument('--verbose', action="store_true", help='Print various developer stats.')

    parser.add_argument('--compression', action="store_true", help='Use homopolymenr compressed reads. (Deprecated, because we will have fewer \
                                                                            minmimizer combinations to span regions in homopolymenr dense regions. Solution \
                                                                            could be to adjust upper interval legnth dynamically to guarantee a certain number of spanning intervals.')
    parser.add_argument('--iso_abundance', type=int, default=1,
                        help='Cutoff parameter: abundance of reads that have to support an isoform to show in results')
    parser.add_argument('--merge_sub_isoforms_3', type=bool, default=True,
                        help='Parameter to determine whether we want to merge sub isoforms (shorter at 3prime end) into bigger isoforms')
    parser.add_argument('--merge_sub_isoforms_5', type=bool, default=True,
                        help='Parameter to determine whether we want to merge sub isoforms (shorter at 5prime end) into bigger isoforms')
    parser.add_argument('--delta_iso_len_3', type=int, default=50,
                        help='Cutoff parameter: maximum length difference at 3prime end, for which subisoforms are still merged into longer isoforms')
    parser.add_argument('--delta_iso_len_5', type=int, default=50,
                        help='Cutoff parameter: maximum length difference at 5prime end, for which subisoforms are still merged into longer isoforms')

    args = parser.parse_args()
    print(len(sys.argv))
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    # if not args.paired_minimizers and 'max_seqs' not in args:
    #     print("max_seqs was not specified and paired_minimizer setting not used. Setting max_seqs to 2000")
    #     args.max_seqs = 2000
    # elif args.paired_minimizers and 'max_seqs' not in args:
    #     print("max_seqs was not specified. Setting max_seqs to 1000")
    #     args.max_seqs = 1000

    if args.set_layers_manually and 'layers' not in args:
        args.layers = 2

    if args.split_mod > 1:
        assert args.residual < args.split_mod

    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
