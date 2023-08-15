import argparse
import os
import sys
import pandas as pd
from utils import *

def parameter_input():
    parser = argparse.ArgumentParser(description='iVirP (Integrative virome pipeline)')
    parser.add_argument('input1', help='Path to the fastq(.gz) file of read1')
    parser.add_argument('input2', help='Path to the fastq(.gz) file of read2')
    parser.add_argument('--host', help='Name(s) of Bowtie2 index(es) built before, please separate different names with "," when removing more than one hosts (for example: "human,mouse") (default: None, which means this step will be skipped)', default=None)
    parser.add_argument('-a', '--adapter', type=int, help='Choose the adapter file used by trimmomatic to remove adapter sequences (please input the interger code, 0:NexteraPE-PE, 1:TruSeq3-PE-2, 2:TruSeq3-PE, 3:TruSeq2-PE, default:NexteraPE-PE)', default=0)
    parser.add_argument('-o', '--output', help='Path to output directory', default=f'{os.getcwd()}/result')
    parser.add_argument('-t', '--threads', type=int, help='Threads used to run this pipeline (default:1)', default=1)
    parser.add_argument('-r', '--remove_inter_result', action='store_true', help='Remove intermediate results (use carefully, this will only keep the final viral contigs, fastqc and viromeqc reports!)', default=False)
    parser.add_argument('-k', '--keep_log', action='store_true', help='This parameter allows you to continue a killed run, but please make sure you have not changed any files in the output directory.', default=False)
    args = parser.parse_args()
    return args          

if __name__ == '__main__':
    args = parameter_input()
    threads = args.threads
    output = args.output
    sample1 = get_sample_name(args.input1.split('/')[-1])
    sample2 = get_sample_name(args.input2.split('/')[-1])
    sample = sample1[0: -2]

    # create the output directory if it does not exist
    create_output_file(output)

    # check if the user is using "keep_log"
    if args.keep_log is False:
        with open(f"{output}/log.txt", "w") as f:
            f.write("0\n")

    # get log info
    with open(f"{output}/log.txt", "r") as f:
        log = int(f.readline()[0: -1])
    
    # assess quality of sequencing with fastqc
    if log < 1:
        run_fastqc(output, args.input1, args.input2, threads)
        log = 1
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")        
    
    # trim adapters with trimmomatic
    if log < 2:
        if args.adapter == 1:
            adapter = f"{sys.path[0]}/adapters/TruSeq3-PE-2.fa"
        elif args.adapter == 2:
            adapter = f"{sys.path[0]}/adapters/TruSeq3-PE.fa"
        elif args.adapter == 3: 
            adapter = f"{sys.path[0]}/adapters/TruSeq2-PE.fa"      
        else:
            adapter = f"{sys.path[0]}/adapters/NexteraPE-PE.fa"
        run_trim(output, threads, args.input1, args.input2, sample1, sample2, adapter)
        log = 2
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")  

    # trim host sequences with bowtie2
    if log < 3:
        if args.host is not None:
            host_list = args.host.split(',')
        else:
            host_list = []
        index_path = []
        for na in host_list:
            index_path.append(f"{sys.path[0]}/bowtie2_index/{na}/{na}")        
        run_bowtie2(output, threads, sample1, sample2, index_path)
        log = 3
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")  

    # report the data quality after do QCs with viromeQC
    if log < 4:
        run_viromeQC(output, sample1, sample2)
        log = 4
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")  
    
    # assemble contigs with spades
    if log < 5:
        run_spades(output, threads, sample1, sample2)
        log = 5
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 

    # trim short contigs with vsearch
    if log < 6:
        run_vsearch_1(output, sample)
        log = 6
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 

    # find viral contigs with virsorter
    if log < 7:
        run_virsorter(output, threads)
        log = 7
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 
    
    # find viral contigs by comparing contigs to databases with blastn
    if log < 8:
        run_blastn(output, threads)
        log = 8
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 

    # filter blastn results and integrate them with virsorter results
    if log < 9:
        run_combination(output)
        log = 9
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 

    # drop viral contigs with low quality or low completeness with checkv
    if log < 10:
        run_checkv(output, threads)
        log = 10
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 

    if log < 11:
        high_quality_output(output)
        log = 11
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n") 

    # cluster contigs and get final non-redundant contigs
    if log < 12:
        run_vsearch_2(output, threads)
        log = 12
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")  
   
    print("all steps finished")

    if args.remove_inter_result == True:
        remove_inter_result(output)
        
