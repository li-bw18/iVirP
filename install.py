import subprocess
import os
import argparse
import sys
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Install script for integrative virome pipeline')
    parser.add_argument('-t', '--threads', type=int, help='Threads used to run this install script', default=1)
    args = parser.parse_args()
    his = 0
    with open("install_history.txt", "r") as f:
        his = int(f.readline()[0: -1])
    if his < 1:
        if os.path.exists("viromeqc") is True:
            subprocess.call(["rm -rf viromeqc"], shell=True)    
        ret = subprocess.call(["git clone --recurse-submodules https://github.com/SegataLab/viromeqc.git"], shell=True)
        if ret != 0:
            sys.exit("install error")
        his = 1            
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 2:
        if os.path.exists("VirSorter2") is True:
            subprocess.call(["rm -rf VirSorter2"], shell=True)    
        ret = subprocess.call(["git clone https://github.com/jiarong/VirSorter2.git"], shell=True)
        if ret != 0:
            sys.exit("install error")
        his = 2          
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 3:
        ret = subprocess.call(["pip install -e VirSorter2"], shell=True)
        if ret != 0:
            sys.exit("install error")
        his = 3        
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 4:
        if os.path.exists("viromeqc/index") is True:
            subprocess.call(["rm -rf viromeqc/index"], shell=True)
        ret = subprocess.call(["python viromeqc/viromeQC.py --install --zenodo"], shell=True)
        if ret != 0:
            sys.exit("install error")
        his = 4    
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 5:
        if os.path.exists("VirSorter2/db") is True:
            subprocess.call(["rm -rf VirSorter2/db"], shell=True)
        ret = subprocess.call([f"virsorter setup -d VirSorter2/db -j {args.threads}"], shell=True)
        if ret != 0:
            sys.exit("install error")
        his = 5
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 6:
        ret = subprocess.call(["wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz -P database/"], shell=True)   
        if ret != 0:
            sys.exit("install error")
        if os.path.exists("database/viral.1.1.genomic.fna") is True:
            subprocess.call(["rm database/viral.1.1.genomic.fna"], shell=True)
        ret = subprocess.call(["gunzip -k database/viral.1.1.genomic.fna.gz"], shell=True)
        if ret != 0:
            sys.exit("install error")        
        his = 6
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 7:
        ret = subprocess.call(["makeblastdb -in database/viral.1.1.genomic.fna -dbtype nucl -parse_seqids -title ncbi -out blastn_database/ncbi -logfile blastn_database/ncbi_log"], shell=True)   
        if ret != 0:
            sys.exit("install error")     
        his = 7
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 8:
        if os.path.exists("database/GVD_1.fna") is True:
            subprocess.call(["rm database/GVD_1.fna"], shell=True)
        if os.path.exists("database/GVD_2.fna") is True:
            subprocess.call(["rm database/GVD_2.fna"], shell=True)
        ret = subprocess.call(["gunzip -k database/GVD*.gz"], shell=True)   
        if ret != 0:
            sys.exit("install error")
        his = 8
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 9:
        ret = subprocess.call(["cat database/GVD_1.fna database/GVD_2.fna > database/GVD.fna"], shell=True)
        if ret != 0:
            sys.exit("install error")     
        his = 9
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 10:
        ret = subprocess.call(["makeblastdb -in database/GVD.fna -dbtype nucl -title gvd -out blastn_database/gvd -logfile blastn_database/gvd"], shell=True)   
        if ret != 0:
            sys.exit("install error")     
        his = 10
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 11:
        ret = subprocess.call(["wget -c ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/GPD_sequences.fa.gz -P database/"], shell=True)   
        if ret != 0:
            sys.exit("install error")
        if os.path.exists("database/GPD_sequences.fa") is True:
            subprocess.call(["rm database/GPD_sequences.fa"], shell=True)
        ret = subprocess.call(["gunzip -k database/GPD_sequences.fa.gz"], shell=True)
        if ret != 0:
            sys.exit("install error")        
        his = 11
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 12:
        ret = subprocess.call(["makeblastdb -in database/GPD_sequences.fa -dbtype nucl -parse_seqids -title gpd -out blastn_database/gpd -logfile blastn_database/gpd_log"], shell=True)   
        if ret != 0:
            sys.exit("install error")     
        his = 12
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 13:
        ret = subprocess.call(["wget -c https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08/mgv_contigs.fna -P database/"], shell=True)   
        if ret != 0:
            sys.exit("install error")      
        his = 13
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 14:
        ret = subprocess.call(["makeblastdb -in database/mgv_contigs.fna -dbtype nucl -parse_seqids -title mgv -out blastn_database/mgv -logfile blastn_database/mgv_log"], shell=True)   
        if ret != 0:
            sys.exit("install error")     
        his = 14
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 15:   
        ret = subprocess.call(["checkv download_database ."], shell=True)   
        if ret != 0:
            sys.exit("install error")     
        ret = subprocess.call(["mv checkv* checkv_database"], shell=True)   
        if ret != 0:
            sys.exit("install error")   
        his = 15
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his == 15:
        if os.path.exists("database") is True:
            subprocess.call(["rm -rf database"], shell=True) 
        print("Install already finished.")
