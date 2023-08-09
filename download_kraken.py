import subprocess
import os
import sys
import argparse
from utils import change_kraken_fasta
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download needed files for rapid detection of eukaryotic viruses')
    parser.add_argument('-t', '--threads', type=int, help='Threads used to run this script', default=1)
    args = parser.parse_args()
    threads = args.threads
    ret = subprocess.call(["wget -c https://cvri.caas.cn/docs/2022-09/e35888f205ca4eafb3807144c0ed5eb5.zip"], shell=True)
    if ret != 0:
        sys.exit("install error")    
    ret = subprocess.call(["unzip *.zip"], shell=True)
    if ret != 0:
        sys.exit("install error")  
    change_kraken_fasta()
    if os.path.exists("kraken2_database") is True:
        subprocess.call(["rm -rf kraken2_database"], shell=True)  
    subprocess.call(["mkdir kraken2_database"], shell=True) 
    ret = subprocess.call([f"kraken2-build --download-taxonomy --threads {threads} --db kraken2_database"], shell=True)
    if ret != 0:
        sys.exit("install error")
    ret = subprocess.call([f"kraken2-build --add-to-library EVRD.fasta --threads {threads} --db kraken2_database"], shell=True)
    if ret != 0:
        sys.exit("install error")
    ret = subprocess.call([f"kraken2-build --build --threads {threads} --db kraken2_database"], shell=True)
    if ret != 0:
        sys.exit("install error")
    subprocess.call(["rm -rf kraken2_info EVRD* *.zip"], shell=True)
    print("Download already finished.")

