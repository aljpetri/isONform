#! /usr/bin/env python
"""
This file is a wrapper file to run isONform in parallel (parallelization over batches).
The script was taken from isoncorrect (https://github.com/ksahlin/isoncorrect/blob/master/run_isoncorrect)
by Kristoffer Sahlin and changed by Alexander Petri to be usable with the isONform code base.

"""
# ! /usr/bin/env python

from __future__ import print_function
import argparse
import tempfile
from time import time
from pathlib import Path
import signal
from multiprocessing import Pool
import multiprocessing as mp

import subprocess
import sys
import os
from sys import stdout
import shutil
import errno

from modules import batch_merging_parallel
from modules import help_functions
from modules import Parallelization_side_functions

def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT
                           ).communicate()[0]
    # print(int(out.split()[0]))
    return int(out.split()[0])


def isONform(data):
    isONform_location, read_fastq_file, outfolder, batch_id, isONform_algorithm_params,cl_id = data[0], data[1], data[
        2], data[3], data[4], data[5]
    help_functions.mkdir_p(outfolder)
    #print("OUT",outfolder)
    #print("Algoparams",isONform_algorithm_params)
    isONform_exec = os.path.join(isONform_location, "main.py")
    isONform_error_file = os.path.join(outfolder, "stderr.txt")
    with open(isONform_error_file, "w") as error_file:
        print('Running isONform batch_id:{0}.{1}...'.format(cl_id,batch_id), end=' ')
        stdout.flush()
        isONform_out_file = open(os.path.join(outfolder, "stdout{0}_{1}.txt".format(cl_id,batch_id)), "w")
        subprocess.check_call(
                ["python", isONform_exec, "--fastq", read_fastq_file, "--outfolder", outfolder,
                 "--exact_instance_limit", str(isONform_algorithm_params["exact_instance_limit"]),
                 #"--max_seqs", str(isONform_algorithm_params["max_seqs"]),
                 "--k", str(isONform_algorithm_params["k"]), "--w", str(isONform_algorithm_params["w"]),
                 "--xmin", str(isONform_algorithm_params["xmin"]), "--xmax",
                 str(isONform_algorithm_params["xmax"]),"--delta_len", str(isONform_algorithm_params["delta_len"]),
                 "--exact", "--parallel", "True",  "--delta_iso_len_3", str(isONform_algorithm_params["delta_iso_len_3"]), "--delta_iso_len_5", str(isONform_algorithm_params["delta_iso_len_5"])#, "--slow"
                 #"--T", str(isONform_algorithm_params["T"])
                 ], stderr=error_file, stdout=isONform_out_file)

        print('Done with batch_id:{0}.{1}'.format(cl_id,batch_id))
        stdout.flush()
    error_file.close()
    isONform_out_file.close()
    return batch_id

#splits files containing more than max_seqs reads into smaller files, that can be parallelized upon
def splitfile(indir, tmp_outdir, fname, chunksize,cl_id,ext):
    # from https://stackoverflow.com/a/27641636/2060202
    # fpath, fname = os.path.split(infilepath)
    #cl_id, ext = fname.rsplit('.',1)
    infilepath = os.path.join(indir, fname)
    #infilepath=indir
    # print(fpath, cl_id, ext)
    #print("now at splitfile")
    #print(indir, tmp_outdir, cl_id, ext)

    i = 0
    written = False
    with open(infilepath) as infile:
        while True:
            outfilepath = os.path.join(tmp_outdir, '{0}_{1}.{2}'.format(cl_id, i, ext) ) #"{}_{}.{}".format(foutpath, fname, i, ext)
            #print(outfilepath)
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
    #print("Symlink",link_name)
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

#splits clusters up so that we get smaller batches
def split_cluster_in_batches_corrected(indir, outdir, tmp_work_dir, max_seqs):
    # create a modified indir
    tmp_work_dir = os.path.join(tmp_work_dir, 'split_in_batches')
    # print(indir)
    help_functions.mkdir_p(tmp_work_dir)

    pat = Path(indir)
    #collect all fastq files located in this directory or any subdirectories
    file_list = list(pat.rglob('*.fastq'))
    #print("FLIST",file_list)
    #iterate over the fastq_files
    for filepath in file_list:
        smaller_than_max_seqs = False
        #print("FPATH",filepath)
        old_fastq_file=str(filepath.resolve())
        path_split=old_fastq_file.split("/")
        folder=path_split[-2]
        #print(folder)
        fastq_file=path_split[-1]
        #we do not want to look at the analysis fastq file
        if not folder=="Analysis":
            cl_id=path_split[-2]
            #print("CLID",cl_id)

            #if we have more lines than max_seqs
            new_indir=os.path.join(indir,folder)
            #print(new_indir)
            if not smaller_than_max_seqs:

                num_lines = sum(1 for line in open(os.path.join(new_indir, fastq_file)))
                #print("Number Lines", fastq_file, num_lines)
                #we reset smaller_than_max_seqs as we now want to see if we really have more than max_seqs reads
                smaller_than_max_seqs = False if num_lines > 4 * max_seqs else True
            else:
                smaller_than_max_seqs = True

            if not smaller_than_max_seqs:
                #print("Splitting",filepath)
                ext = fastq_file.rsplit('.', 1)[1]
                splitfile(new_indir, tmp_work_dir, fastq_file, 4 * max_seqs,cl_id,ext)  # is fastq file
            else:
                ext = fastq_file.rsplit('.', 1)[1]
                #print(fastq_file, "symlinking instead")
                symlink_force(filepath, os.path.join(tmp_work_dir, '{0}_{1}.{2}'.format(cl_id, 0, ext)))
    return tmp_work_dir

def split_cluster_in_batches_clust(indir, outdir, tmp_work_dir, max_seqs):
    # create a modified indir
    tmp_work_dir = os.path.join(tmp_work_dir, 'split_in_batches')
    # print(indir)
    #os.mkdir(tmp_work_dir)
    help_functions.mkdir_p(tmp_work_dir)
    #print(tmp_work_dir)
    #print("clust")

    pat = Path(indir)
    file_list = list(pat.rglob('*.fastq'))
    # add split fiels to this indir
    for file_ in file_list:
        smaller_than_max_seqs = False
    #for file_ in sorted(os.listdir(indir), key=lambda x: int(x.split('.')[0])):
        #fastq_path = os.fsdecode(file_)
        old_fastq_file = str(file_.resolve())
        fastq_file = old_fastq_file.split("/")[-1]
        #print("FASTQ",fastq_file)
        if not smaller_than_max_seqs:
            num_lines = sum(1 for line in open(os.path.join(indir, fastq_file)))
            smaller_than_max_seqs = False if num_lines > 4 * max_seqs else True
        else:
            smaller_than_max_seqs = True
        if not smaller_than_max_seqs:
            cl_id, ext = fastq_file.rsplit('.', 1)
            splitfile(indir, tmp_work_dir, fastq_file, 4 * max_seqs, cl_id, ext) # is fastq file
        else:
            cl_id, ext = fastq_file.rsplit('.', 1)
            print("SYMLINK",os.path.join(tmp_work_dir, '{0}_{1}.{2}'.format(cl_id, 0, ext)))
            symlink_force(file_, os.path.join(tmp_work_dir, '{0}_{1}.{2}'.format(cl_id, 0, ext)))
    return tmp_work_dir


#PYTHONHASHSEED = 0
def main(args):
    #print("MERGE?", args.merge_sub_isoforms_3, args.merge_sub_isoforms_5)
    globstart = time()
    directory = args.fastq_folder  # os.fsencode(args.fastq_folder)
    #print(directory)
    #print("ARGS",args)
    isONform_location = os.path.dirname(os.path.realpath(__file__))
    if args.split_wrt_batches:
        if args.tmpdir:
            tmp_work_dir = args.tmpdir
            #curr_work_dir = os.getcwd()
            #os.chdir(tmp_work_dir)
            #try:
                # Create the directory
                #os.makedirs(tmp_work_dir)
            #    print("Directory created successfully.")
            #except FileExistsError:
            #    print("Directory already exists.")
            Parallelization_side_functions.mkdir_p(tmp_work_dir)

            #os.chdir(curr_work_dir)

        else:
            tmp_work_dir = tempfile.mkdtemp()
        #print("SPLITWRTBATCHES")

        print("Temporary workdirektory:", tmp_work_dir)
        if args.clustered:
            split_tmp_directory = split_cluster_in_batches_clust(directory, args.outfolder, tmp_work_dir,
                                                                     args.max_seqs)
        else:
            split_tmp_directory = split_cluster_in_batches_corrected(directory, args.outfolder, tmp_work_dir, args.max_seqs)
        split_directory = os.fsencode(split_tmp_directory)
        #print("SplitDIR",split_directory)
    else:
        split_directory = os.fsencode(directory)

    instances = []
    for file_ in os.listdir(split_directory):
        #print(file_)
        read_fastq_file = os.fsdecode(file_)
        if read_fastq_file.endswith(".fastq"):
            #print("True")
            tmp_id= read_fastq_file.split(".")[0]
            snd_tmp_id=tmp_id.split("_")
            cl_id = snd_tmp_id[0]
            batch_id=snd_tmp_id[1]
            outfolder = os.path.join(args.outfolder, cl_id)
            #print(batch_id,cl_id)
            #print(outfolder)
            fastq_file_path = os.path.join(os.fsdecode(split_directory), read_fastq_file)
            #print(fastq_file_path)
            compute = True
            if args.keep_old:
                candidate_corrected_file = os.path.join(outfolder, "isoforms.fastq")
                if os.path.isfile(candidate_corrected_file):
                    if wccount(candidate_corrected_file) == wccount(fastq_file_path):
                        #print("already computed cluster and complete file", batch_id)
                        compute = False

            if compute:
                #print("computing")
                isONform_algorithm_params = {"set_w_dynamically": args.set_w_dynamically,
                                                "exact_instance_limit": args.exact_instance_limit,
                                                "delta_len": args.delta_len,"--exact": True,
                                                "k": args.k, "w": args.w, "xmin": args.xmin, "xmax": args.xmax,
                                                 "max_seqs": args.max_seqs, "parallel": True, "--slow": True, "delta_iso_len_3": args.delta_iso_len_3,
                                             "delta_iso_len_5": args.delta_iso_len_5}
                instances.append(
                    (isONform_location, fastq_file_path, outfolder, batch_id, isONform_algorithm_params,cl_id))
        else:
            continue

    instances.sort(key=lambda x: x[3])  # sorting on batch ids as strings
    print("Printing instances")
    for t in instances:
        print(t)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass
    #mp.set_start_method('spawn')
    print(mp.get_context())
    print("Environment set:", mp.get_context())
    print("Using {0} cores.".format(args.nr_cores))
    start_multi = time()
    pool = Pool(processes=int(args.nr_cores))
    try:
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
    print("Time elapsed multiprocessing:", time() - start_multi)

    if args.split_wrt_batches:
        print("STILLSPLITWRTBATCHES")
        file_handling = time()
        if args.write_fastq:
            write_fastq = True
        else:
            write_fastq = False
        batch_merging_parallel.join_back_via_batch_merging(args.outfolder, args.delta, args.delta_len, args.delta_iso_len_3, args.delta_iso_len_5, args.max_seqs_to_spoa,args.iso_abundance, write_fastq)
        Parallelization_side_functions.generate_full_output(args.outfolder,write_fastq)
        Parallelization_side_functions.remove_folders(args.outfolder)
        shutil.rmtree(split_directory)
        print("Joined back batched files in:", time() - file_handling)
        print("Finished full algo after :", time() - globstart)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="De novo reconstruction of long-read transcriptome reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 0.3.3')
    parser.add_argument('--fastq_folder', type=str, default=False,
                        help='Path to input fastq folder with reads in clusters')
    parser.add_argument('--t', dest="nr_cores", type=int, default=8, help='Number of cores allocated for clustering')
    parser.add_argument('--k', type=int, default=20, help='Kmer size')
    parser.add_argument('--w', type=int, default=31, help='Window size')
    parser.add_argument('--xmin', type=int, default=18, help='Lower interval length')
    parser.add_argument('--xmax', type=int, default=80, help='Upper interval length')
    parser.add_argument('--exact_instance_limit', type=int, default=50,
                        help='Do exact correction for clusters under this threshold')
    parser.add_argument('--keep_old', action="store_true",
                        help='Do not recompute previous results if corrected_reads.fq is found and has the smae number of reads as input file (i.e., is complete).')
    parser.add_argument('--set_w_dynamically', action="store_true",
                        help='Set w = k + max(2*k, floor(cluster_size/1000)).')
    parser.add_argument('--max_seqs', type=int, default=1000,
                        help='Maximum number of seqs to correct at a time (in case of large clusters).')
    parser.add_argument('--split_wrt_batches', action="store_true",
                        help='Process reads per batch (of max_seqs sequences) instead of per cluster. Significantly decrease runtime when few very large clusters are less than the number of cores used.')
    parser.add_argument('--clustered', action="store_true",
                        help='Indicates whether we use the output of isONclust (i.e. we have uncorrected data)')
    parser.add_argument('--outfolder', type=str, default=None, help='Outfolder with all corrected reads.')
    parser.add_argument('--delta_len', type=int, default=5,
                        help='Maximum length difference between two reads intervals for which they would still be merged')
    parser.add_argument('--delta',type=float,default=0.1, help='diversity rate used to compare sequences')
    parser.add_argument('--max_seqs_to_spoa', type=int, default=200, help='Maximum number of seqs to spoa')
    parser.add_argument('--verbose', action="store_true", help='Print various developer stats.')
    parser.add_argument('--iso_abundance', type=int, default=5,
                        help='Cutoff parameter: abundance of reads that have to support an isoform to show in results')
    parser.add_argument('--delta_iso_len_3', type=int, default=30,
                        help='Cutoff parameter: maximum length difference at 3prime end, for which subisoforms are still merged into longer isoforms')
    parser.add_argument('--delta_iso_len_5', type=int, default=50,
                        help='Cutoff parameter: maximum length difference at 5prime end, for which subisoforms are still merged into longer isoforms')
    parser.add_argument('--tmpdir', type=str,default=None, help='OPTIONAL PARAMETER: Absolute path to custom folder in which to store temporary files. If tmpdir is not specified, isONform will attempt to write the temporary files into the tmp folder on your system. It is advised to only use this parameter if the symlinking does not work on your system.')
    parser.add_argument('--write_fastq', action="store_true", help=' Indicates that we want to ouptut the final output (transcriptome) as fastq file (New standard: fasta)')
    args = parser.parse_args()
    print(len(sys.argv))
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()


    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
