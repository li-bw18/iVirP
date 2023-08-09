import argparse
import sys
import subprocess

def parameter_input():
    parser = argparse.ArgumentParser(description='Build the bracken database for rapid detection of eukaryotic viruses')
    parser.add_argument('-l', '--length', type=int, help='The length of reads', default=100)
    parser.add_argument('-t', '--threads', type=int, help='Threads used to run this script (default:1)', default=1)    
    args = parser.parse_args()
    return args          

if __name__ == '__main__':
    args = parameter_input()
    subprocess.call([f'rm {sys.path[0]}/kraken2_database/database*'], shell=True)
    ret = subprocess.call([f'bracken-build -l {args.length} -t {args.threads} -d {sys.path[0]}/kraken2_database'], shell=True)
    if ret != 0:
        sys.exit("install error")
