import argparse
import os
import sys
from utils import *

def parameter_input():
    parser = argparse.ArgumentParser(description='Integrative virome pipeline -- Rapid detection of eukaryotic viruses')
    parser.add_argument('input1', help='Path to the fastq(.gz) file of read1')
    parser.add_argument('input2', help='Path to the fastq(.gz) file of read2')
    parser.add_argument('--host', help='Name(s) of Bowtie2 index(es) built before, please separate different names with "," when removing more than one hosts (for example: "human,mouse") (default: None, which means this step will be skipped)', default=None)
    parser.add_argument('-l', '--length', type=int, help='The length of reads', default=100)
    parser.add_argument('-a', '--adapter', type=int, help='Choose the adapter file used by trimmomatic to remove adapter sequences (please input the interger code, 0:NexteraPE-PE, 1:TruSeq3-PE-2, 2:TruSeq3-PE, 3:TruSeq2-PE, default:NexteraPE-PE)', default=0)
    parser.add_argument('-o', '--output', help='Path to output directory', default=f'{os.getcwd()}/result')
    parser.add_argument('-t', '--threads', type=int, help='Threads used to run this pipeline (default:1)', default=1)
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
    
    # Run kraken2
    if log < 5:
        print("Run kraken2")
        if os.path.exists(f"{output}/5.kraken2") is True:
            subprocess.call([f"rm -rf {output}/5.kraken2"], shell=True)
        subprocess.call([f"mkdir {output}/5.kraken2"], shell=True)
        print(f"kraken2 --paired {output}/3.bowtie2/{sample1}.fastq {output}/3.bowtie2/{sample2}.fastq --output {output}/5.kraken2/output.txt --threads {threads} --report {output}/5.kraken2/report.txt --db {sys.path[0]}/kraken2_database")
        ret = subprocess.call([f"kraken2 --paired {output}/3.bowtie2/{sample1}.fastq {output}/3.bowtie2/{sample2}.fastq --output {output}/5.kraken2/output.txt --threads {threads} --report {output}/5.kraken2/report.txt --db {sys.path[0]}/kraken2_database"], shell=True)
        if ret != 0:
            sys.exit("Error: kraken2 error")
        log = 5
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")  

    # Run bracken
    if log < 6:
        print("Run bracken")
        if os.path.exists(f"{output}/6.bracken") is True:
            subprocess.call([f"rm -rf {output}/6.bracken"], shell=True)
        subprocess.call([f"mkdir {output}/6.bracken"], shell=True)
        print(f"bracken -d {sys.path[0]}/kraken2_database -i {output}/5.kraken2/report.txt -o {output}/6.bracken/output.txt -w {output}/6.bracken/report.txt -r {args.length}")
        ret = subprocess.call([f"bracken -d {sys.path[0]}/kraken2_database -i {output}/5.kraken2/report.txt -o {output}/6.bracken/output.txt -w {output}/6.bracken/report.txt -r {args.length}"], shell=True)
        if ret != 0:
            sys.exit("Error: bracken error")
        log = 6
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")  

    # Reform ouput file
    if log < 7:
        print("Run report reformation")
        if os.path.exists(f"{output}/7.final") is True:
            subprocess.call([f"rm -rf {output}/7.final"], shell=True)
        subprocess.call([f"mkdir {output}/7.final"], shell=True)
        print(f"kreport2mpa.py -r {output}/6.bracken/report.txt -o {output}/7.final/report.txt")
        ret = subprocess.call([f"kreport2mpa.py -r {output}/6.bracken/report.txt -o {output}/7.final/report.txt"], shell=True)
        if ret != 0:
            sys.exit("Error: reformation error")
        log = 7
        with open(f"{output}/log.txt", "w") as f:
            f.write(f"{log}\n")       

        
