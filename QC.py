import argparse
import os
import sys
import pandas as pd
from utils import *


def parameter_input():
    parser = argparse.ArgumentParser(description='Integrative viral detection pipeline - QC')
    parser.add_argument('input1', help='Path to the fastq(.gz) file of read1')
    parser.add_argument('input2', help='Path to the fastq(.gz) file of read2')
    parser.add_argument('--host', help='Name(s) of Bowtie2 index(es) built before, please separate different names with "," when removing more than one hosts (for example: "human,mouse") (default: None, which means this step will be skipped)', default=None)
    parser.add_argument('-a', '--adapter', type=int, help='Choose the adapter file used by trimmomatic to remove adapter sequences (please input the interger code, 0:NexteraPE-PE, 1:TruSeq3-PE-2, 2:TruSeq3-PE, 3:TruSeq2-PE, default:NexteraPE-PE)', default=0)
    parser.add_argument('-o', '--output', help='Path to output directory', default=f'{os.getcwd()}/result')
    parser.add_argument('-t', '--threads', type=int, help='Threads used to run this pipeline (default:1)', default=1)
    parser.add_argument('--skip_fastqc', action='store_true', help='Skip the fastqc step', default=False)
    parser.add_argument('--skip_trim', action='store_true', help='Skip the trim adapters step', default=False)
    parser.add_argument('--skip_viromeqc', action='store_true', help='Skip the viromeqc step', default=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parameter_input()
    threads = args.threads
    output = args.output
    sample1 = get_sample_name(args.input1.split('/')[-1])
    sample2 = get_sample_name(args.input2.split('/')[-1])
    sample = sample1[0: -2]

    create_output_file(output)

    # assess quality of sequencing with fastqc
    if args.skip_fastqc is False:
        run_fastqc(output, args.input1, args.input2, threads)
    else:
        print("fastqc skip")
    
    # trim adapters with trimmomatic        
    if os.path.exists(f"{output}/2.trim") is True:
        subprocess.call([f"rm -rf {output}/2.trim"], shell=True)
    subprocess.call([f"mkdir {output}/2.trim"], shell=True)
    if args.skip_trim is False:
        if args.adapter == 1:
            adapter = f"{sys.path[0]}/adapters/TruSeq3-PE-2.fa"
        elif args.adapter == 2:
            adapter = f"{sys.path[0]}/adapters/TruSeq3-PE.fa"
        elif args.adapter == 3: 
            adapter = f"{sys.path[0]}/adapters/TruSeq2-PE.fa"      
        else:
            adapter = f"{sys.path[0]}/adapters/NexteraPE-PE.fa"        
        run_trim(output, threads, args.input1, args.input2, sample1, sample2, adapter)
    else:
        print("trim adapters skip")
        if args.input1[-3:] == '.gz':
            subprocess.call([f"zcat {args.input1} > {output}/2.trim/{sample1}.fastq"], shell=True)
            subprocess.call([f"zcat {args.input2} > {output}/2.trim/{sample2}.fastq"], shell=True)      
        else:
            subprocess.call([f"cp {args.input1} {output}/2.trim/{sample1}.fastq"], shell=True)
            subprocess.call([f"cp {args.input2} {output}/2.trim/{sample1}.fastq"], shell=True)      

    # trim host sequences with bowtie2
    if args.host is not None:
        host_list = args.host.split(',')
    else:
        host_list = []
    index_path = []
    for na in host_list:
        index_path.append(f"{sys.path[0]}/bowtie2_index/{na}/{na}")        
    run_bowtie2(output, threads, sample1, sample2, index_path)          

    # report the data quality after do QCs with viromeQC
    if args.skip_viromeqc is False:
        run_viromeQC(output, sample1, sample2)
    else:
        print("viromeqc skip")

    if args.skip_trim is True:
        subprocess.call([f"rm -rf {output}/2.trim"], shell=True)
    
    print("all steps finished")










