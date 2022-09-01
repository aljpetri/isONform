from __future__ import print_function
import os
import errno
import os.path

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def generate_low_abundance_mapping(outfolder):
    subfolders = [f.path for f in os.scandir(outfolder) if f.is_dir()]
    f = open(os.path.join(outfolder, "transcriptome_mapping_low.txt"), "w")
    for subfolder in subfolders:
        actual_folder = subfolder.split("/")[-1]
        # print(actual_folder)
        if actual_folder.isdigit():
            fname = os.path.join(outfolder, "cluster" + str(actual_folder) + "_mapping_low_abundance.txt")
            # print(fname)

            if os.path.isfile(fname):
                g = open(fname, "r")
                # read content from first file
                for line in g:
                    # append content to second file
                    f.write(line)
def generate_single_mapping(outfolder):
    subfolders = [f.path for f in os.scandir(outfolder) if f.is_dir()]
    f = open(os.path.join(outfolder, "transcriptome_mapping.txt"), "w")
    for subfolder in subfolders:
        actual_folder = subfolder.split("/")[-1]
        # print(actual_folder)
        if actual_folder.isdigit():
            fname = os.path.join(outfolder, "cluster" + str(actual_folder) + "_mapping.txt")
            #print(fname)

            if os.path.isfile(fname):
                g = open(fname, "r")
                # read content from first file
                for line in g:
                    # append content to second file
                    f.write(line)
def generate_single_output(outfolder):
    subfolders = [f.path for f in os.scandir(outfolder) if f.is_dir()]
    f = open(os.path.join(outfolder,"transcriptome.fq"), "w")
    for subfolder in subfolders:
        actual_folder=subfolder.split("/")[-1]
        #print(actual_folder)
        if actual_folder.isdigit():
            fname=os.path.join(outfolder,"cluster"+str(actual_folder)+"_merged.fq")
            #print(fname)

            if os.path.isfile(fname):
                g = open(fname, "r")
                # read content from first file
                for line in g:
                    # append content to second file
                    f.write(line)
            """otherfname= os.path.join(outfolder,"cluster"+str(actual_folder)+"_merged_low_abundance.fq")
            if os.path.isfile(otherfname):
                other_g = open(otherfname, "r")
                # read content from first file
                for other_line in other_g:
                    # append content to second file
                    f.write(other_line)
                #f.write(g.read())"""
def generate_full_output(outfolder):
    generate_single_output(outfolder)
    generate_single_mapping(outfolder)
    generate_low_abundance_mapping(outfolder)