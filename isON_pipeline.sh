#!/bin/bash
set -e
#the pipeline can be run in different modes:

####ONT data
# ont_with_pychopper: the full pipeline is run in addition to pychopper (pychopper, isONclust,isONcorrect,isONform)
# ont_no_pychopper: only the isONpipeline is run without pychopper (isONclust,isONcorrect, isONform)
# ont_no_pychop_isonclust3
####PACBIO data
# pacbio: for PacBio data runs isONclust and isONform


###### test modes (only for internal use)
# analysis: analysis of ont data:  isONclust,isONcorrect and isONform are run (e.g. analyses on the paper)
# only_isonform: only isONform is run
#!/bin/bash

programname=$0
function usage {
    echo ""
    echo "Runs the full isON pipeline. Please make sure that the input file has been preprocessed with pychopper for ONT data "
    echo ""
    echo "usage: $programname --raw_reads string --outfolder string --num_cores integer --isONform_folder string  --iso_abundance integer --mode string"
    echo ""
    echo "  --raw_reads   	        absolute path to the input file (in fastq format)"
    echo "                          (example: /home/user/Rawdata/raw_reads.fq)"
    echo "  --outfolder             absolute path to the output folder (the folder in which all outputs are stored)"
    echo "                          (example: /home/user/analysis_output)"
    echo "  --num_cores             the number of processors the pipeline may use"
    echo "                          (example: 8)"
    echo "  --isONform_folder       the absolute path to the isONform installation on your machine (leave empty if you have installed isONform via pip)"
    echo "                          (example: /home/user/isONform )"
    echo "  --iso_abundance         threshold which denotes the minimum read support neccessary for an isoform to be called (also minimum number of reads per cluster in isONclust)"
    echo "                          (example: 5)"
    echo "  --mode                  Run mode of the pipeline, possible modes are 'ont_no_pyc' and 'ont_with_pc' for ont data and 'pacbio' for pacbio data"
    echo "                          (example: ont_no_pychopper/ont_with_pychopper/pacbio)"
    echo " For ONT data: use 'ont_no_pychopper' if you want to run the isON pipeline and pychopper, use 'ont_with_pychopper' if you only want to run the isON pipeline. Please run pychopper yourself before running the pipeline."
    echo ""
}

while [ $# -gt 0 ]; do
    # Check if the current argument is "--help"
    if [[ $1 == "--help" ]]; then
        # Call the usage function and exit with status code 0
        usage
        exit 0
    # Check if the current argument is an option starting with "--"
    elif [[ $1 == "--"* ]]; then
        # Extract the option name by removing the leading dashes
        v="${1/--/}"
        # Check if the argument for this option was left empty (then we would have the next argument name as next entry)
        if [[ $2 != "--"* ]]; then
           #The argument was not left empty, therefore we properly set the argument as the value for option
           declare "$v"="$2"
           #We have to shift only in this case
           shift
        fi
    fi
    #This is the shift we have to perform each time
    shift
done

if [[ -z $raw_reads ]]; then
    usage
    die "Missing parameter --raw_reads"
elif [[ -z $outfolder ]]; then
    usage
    die "Missing parameter --outfolder"
elif [[ -z $mode ]]; then
    usage
    die "Missing parameter --mode"
#elif [[ -z $isONform_folder ]]; then
#    isONform_folder=''
    #TODO set isONform folder to '' if not given
fi

echo "Running `basename $0` raw reads: '$raw_reads' outfolder: '$outfolder' num_cores: '$num_cores' isONform_folder:'$isONform_folder' iso_abundance: '$iso_abundance' mode: '$mode'"

mkdir -p $outfolder
##Testing whether the programs are installed properly before attempting to run them
# shellcheck disable=SC1072
if [ $mode != "pacbio" ] && [ $mode != "ont_no_pychop" ]; #we do not want to test for pychopper if we run modes pacbio or ont_no_pychop
then
  pychopper --h
fi
isONclust
#only test for isONcorrect if we want to use it (so not neccessary for Pacbio data)
if [ $mode != "pacbio" ]
then
  run_isoncorrect
fi

if [ -n "$isONform_folder"  ] #the user has given a path to isONform (cloned from github)
  then
  $isONform_folder/isONform_parallel --h
else
  python isONform_parallel --h
fi



if [ $mode == "ont_with_pychopper" ]
then
echo
echo "Will run pychopper (cdna_classifier.py), isONclust, isONcorrect and isONform. Make sure you have these tools installed."
echo "For installation see: https://github.com/ksahlin/isONcorrect#installation and  https://github.com/aljpetri/isONform"
echo

echo
echo "Running pychopper"
echo

pychopper  $raw_reads $outfolder/full_length.fq -t $num_cores

echo
echo "Finished pychopper"
echo

fi

if [ $mode != "only_isonform" ] # this if statement prevents isONclust and isONcorrect from being run
  then
    echo
    echo "Running isONclust"
    echo
  if [ $mode == "ont_with_pychopper" ]
    then
        /usr/bin/time -v isONclust  --t $num_cores  --ont --fastq $outfolder/full_length.fq \
             --outfolder $outfolder/clustering
        /usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $outfolder/full_length.fq --outfolder  $outfolder/clustering/fastq_files
  elif [ $mode == "ont_no_pychopper" ]
   then
       /usr/bin/time -v  isONclust  --t $num_cores  --ont --fastq $raw_reads \
             --outfolder $outfolder/clustering
       /usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $raw_reads --outfolder  $outfolder/clustering/fastq_files
  elif [ $mode == "pacbio" ] #This is the pacbio mode
    then
       /usr/bin/time -v  isONclust  --t $num_cores  --isoseq  --fastq $raw_reads \
             --outfolder $outfolder/clustering
       /usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $raw_reads --outfolder  $outfolder/clustering/fastq_files

  #elif [ $mode == "ont_no_pychop_isonclust3" ]
  #  then
  #    /home/alexanderpetri/Rust/isONclust_rs/target/release/isONclust3  --mode ont --fastq $raw_reads --outfolder $outfolder --n 1
  else #[ $mode != "pacbio" ] && [ $mode != "'ont'" ]
       /usr/bin/time -v  isONclust  --t $num_cores   --fastq $raw_reads \
             --outfolder $outfolder
       /usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $raw_reads --outfolder  $outfolder/clustering/fastq_files


  fi

  echo
  echo "Finished isONclust"
  echo
#conda activate isON311

  if [ $mode != "pacbio" ]
  then
    echo
    echo "Running isONcorrect"
    echo

    /usr/bin/time -v  run_isoncorrect --t $num_cores  --fastq_folder $outfolder/clustering/fastq_files  --outfolder $outfolder/correction/

    echo
    echo "Finished isONcorrect"
    echo
  fi
fi
echo
echo "Merging reads back to single file. Corrected reads per cluster are still stored in: " $outfolder/correction/
echo

echo
echo "Running isONform"
echo
if [ -n "$isONform_folder"  ] #the user has given a path to isONform (cloned from github)
  then
  if [ $mode != "pacbio" ] #i.e. we run in ONT mode
    then
        /usr/bin/time -v  $isONform_folder/isONform_parallel --t $num_cores  --fastq_folder $outfolder/correction --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 10 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance --delta_iso_len_3 30 --delta_iso_len_5 50
  else #we run isONform in pacbio mode (adding the keyword clustered to the command)
        /usr/bin/time -v   $isONform_folder/isONform_parallel --t $num_cores --fastq_folder $outfolder/clustering/fastq_files --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 10 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance  --delta_iso_len_3 30 --delta_iso_len_5 50
  fi


else #the user has not given a path to isONform (pip installation supposed)
  if [ $mode != "pacbio" ] #i.e. we run in ONT mode
    then
        /usr/bin/time -v  isONform_parallel --t $num_cores  --fastq_folder $outfolder/correction/ --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 10 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance  --delta_iso_len_3 30 --delta_iso_len_5 50
  else #we run isONform in pacbio mode (adding the keyword clustered to the command)
        /usr/bin/time -v   isONform_parallel --t $num_cores --fastq_folder $outfolder/clustering/fastq_files --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 10 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance --delta_iso_len_3 30 --delta_iso_len_5 50 
  fi
fi
echo
echo "Finished isONform"
echo

echo
echo "Finished with pipeline and wrote corrected reads into: " $outfolder
echo