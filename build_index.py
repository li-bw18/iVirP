import argparse
import os
import sys
import subprocess

def parameter_input():
    parser = argparse.ArgumentParser(description='Build Bowtie2 indexes for host removal')
    parser.add_argument('-f', '--fasta', help='Path to the genome fasta file of host (.fna/.fa/.fasta)')
    parser.add_argument('-n', '--name', help='Name of the Bowtie2 index needed to remove')
    parser.add_argument('-t', '--threads', type=int, help='Threads used to build Bowtie2 index (default:1)', default=1)
    parser.add_argument('-r', '--remove', action='store_true', help='Remove an existing Bowtie2 index', default=False)
    parser.add_argument('-c', '--check', action='store_true', help='Check the built Bowtie2 indexes', default=False)
    args = parser.parse_args()
    return args  

def get_name(input):
    if input[-4:] == '.fna':
        return input[:-4]
    elif input[-3:] == '.fa':
        return input[:-3]
    elif input[-6:] == '.fasta':
        return input[:-6]
    else:
        sys.exit("The input fasta file has a wrong suffix")       

if __name__ == '__main__':
    args = parameter_input()
    if os.path.exists(f"{sys.path[0]}/bowtie2_index") is False:
        subprocess.call([f"mkdir {sys.path[0]}/bowtie2_index"], shell=True)   
    if args.check == True:
        li = os.listdir(f"{sys.path[0]}/bowtie2_index")
        if len(li) == 0:
            print("No existing Bowtie2 index")
        else:
            for na in os.listdir(f"{sys.path[0]}/bowtie2_index"):
                print(na)    
    elif args.remove == True:
        if os.path.exists(f"{sys.path[0]}/bowtie2_index/{args.name}") is True:
            subprocess.call([f"rm -rf {sys.path[0]}/bowtie2_index/{args.name}"], shell=True)  
        else:
            sys.exit("This Bowtie2 index does not exist") 
    else:
        name = get_name(args.fasta.split('/')[-1])
        if os.path.exists(f"{sys.path[0]}/bowtie2_index/{name}") is True:
            subprocess.call([f"rm -rf {sys.path[0]}/bowtie2_index/{name}"], shell=True)  
        subprocess.call([f"mkdir {sys.path[0]}/bowtie2_index/{name}"], shell=True)        
        ret = subprocess.call([f"bowtie2-build --threads {args.threads} {args.fasta} {sys.path[0]}/bowtie2_index/{name}/{name}"], shell=True)
        if ret != 0:
            sys.exit("Index building error")
