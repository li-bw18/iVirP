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
        ret = subprocess.call(["wget -c https://zenodo.org/record/7700466/files/data.tar.gz"], shell=True)
        if ret != 0:
            sys.exit("install error")
        his = 6
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 7:
        if os.path.exists("data") is True:
            subprocess.call(["rm -rf data"], shell=True)        
        ret = subprocess.call(["tar -xzvf data.tar.gz"], shell=True)
        if ret != 0:
            sys.exit("install error")
        subprocess.call(["rm data.tar.gz"], shell=True)
        his = 7
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his < 8:
        if os.path.exists("blastn_database") is True:
            subprocess.call(["rm -rf blastn_database"], shell=True)
        ret = subprocess.call(["mv data/* ."], shell=True)
        if ret != 0:
            sys.exit("install error")  
        subprocess.call(["rm -rf data"], shell=True)      
        his = 8
        with open("install_history.txt", "w") as f:
            f.write(f"{his}\n")
    if his == 8:
        if os.path.exists("database") is True:
            subprocess.call(["rm -rf database"], shell=True)
        print("Install already finished.")
