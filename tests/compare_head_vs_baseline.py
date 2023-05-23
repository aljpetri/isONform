from modules import help_functions
import argparse
import sys

def main(args):
    baseline_reads = help_functions.readfq(args.baseline)
    head_reads = help_functions.readfq(args.head)
    if baseline_reads == head_reads:
        print("Baseline and head fastq's contain the same reads")
    else:
        print("Baseline and head fastq's contain different reads")
        base_set = set(baseline_reads.items())
        head_set = set(head_reads.items())
        print("reads in baseline but not head:")
        print(base_set - head_set)
        print("reads in head but not baseline:")
        print(head_set - base_set)
        sys.exit(1)
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=str, default=False, help='fastq file generated by the baseline version of the repo')
    parser.add_argument('--head', type=str, default=False, help='fastq file generated by the head version of the repo')
    args = parser.parse_args()
    if __name__ == "__main__":
        main(args)
