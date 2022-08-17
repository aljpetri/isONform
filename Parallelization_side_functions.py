from __future__ import print_function
import os
import errno

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise



def generate_single_output(outfolder):
    subfolders = [f.path for f in os.scandir(outfolder) if f.is_dir()]
    f = open("transcriptome.fq", "w")
    for subfolder in subfolders:
        actual_folder=subfolder.split("/")[-1]
        print(actual_folder)
        if actual_folder.isdigit():
            openfile=os.path.join(subfolder,"cluster"+str(actual_folder)+"_merged.fa")
            g=open(openfile,"r")
            f.write(g.read())