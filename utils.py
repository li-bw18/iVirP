import subprocess
import os
import sys
import pandas as pd

def change_kraken_fasta():
    fo = open(f"{sys.path[0]}/EVRD.fasta", "w")
    dat = pd.read_table(f"{sys.path[0]}/kraken2_info/result.txt", header=None)
    i = 0
    with open(f"{sys.path[0]}/EVRD-nt_20220914.fasta", "r") as f:
        line = f.readline()
        while 1:
            if line == '':
                break
            if line[0] == '>':
                sp = line[:-1].split("|")
                if len(sp) >= 2:
                    if dat.iloc[i, 0] == -1:
                        i = i + 1
                        line = f.readline()
                        while line[0] != '>':
                            line = f.readline()
                    else:
                        print(f">kraken:taxid|{int(dat.iloc[i, 0])}|{sp[1]}", file=fo)
                        i = i + 1
                        line = f.readline()
                else:
                    break   
            else:
                print(line[:-1].upper(), file=fo) 
                line = f.readline()

def get_sample_name(file):
    if file[-6: ] == '.fq.gz' or file[-6: ] == '.fastq':
        return file[0: -6]
    elif file[-3: ] == '.fq':
        return file[0: -3]
    elif file[-9: ] == '.fastq.gz':
        return file[0: -9]

def filter_vircontig(output):
    filtered = pd.DataFrame(data=None, columns=['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'database'])
    blastn = []
    virso = []
    if os.path.getsize(f"{output}/8.blastncontigs/crass.out") != 0:
        crass = pd.read_table(f"{output}/8.blastncontigs/crass.out", header=None)
        crass.columns = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'nident', 'qlen', 'slen', 'length', 'mismatch', 'positive', 'ppos', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qcovhsp', 'qcovus', 'qseq', 'sstrand', 'frames']
        crass = crass[crass['pident']>=50]
        crass = crass[crass['evalue']<=1e-10]
        crass = crass[crass['qcovs']>=80]
        for i in range(len(crass)):
            if crass.iloc[i, 0] not in blastn:
                blastn.append(crass.iloc[i, 0])
    crass = crass.iloc[:, :5]
    crass['database'] = ['crass'] * len(crass)
    filtered = pd.concat([filtered, crass], axis=0)
    if os.path.getsize(f"{output}/8.blastncontigs/gpd.out") != 0:
        gpd = pd.read_table(f"{output}/8.blastncontigs/gpd.out", header=None)
        gpd.columns = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'nident', 'qlen', 'slen', 'length', 'mismatch', 'positive', 'ppos', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qcovhsp', 'qcovus', 'qseq', 'sstrand', 'frames']
        gpd = gpd[gpd['pident']>=50]
        gpd = gpd[gpd['evalue']<=1e-10]
        gpd = gpd[gpd['qcovs']>=80]
        for i in range(len(gpd)):
            if gpd.iloc[i, 0] not in blastn:
                blastn.append(gpd.iloc[i, 0])
    gpd = gpd.iloc[:, :5]
    gpd['database'] = ['gpd'] * len(gpd)
    filtered = pd.concat([filtered, gpd], axis=0)
    if os.path.getsize(f"{output}/8.blastncontigs/gvd.out") != 0:
        gvd = pd.read_table(f"{output}/8.blastncontigs/gvd.out", header=None)
        gvd.columns = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'nident', 'qlen', 'slen', 'length', 'mismatch', 'positive', 'ppos', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qcovhsp', 'qcovus', 'qseq', 'sstrand', 'frames']
        gvd = gvd[gvd['pident']>=50]
        gvd = gvd[gvd['evalue']<=1e-10]
        gvd = gvd[gvd['qcovs']>=80]
        for i in range(len(gvd)):
            if gvd.iloc[i, 0] not in blastn:
                blastn.append(gvd.iloc[i, 0])
    gvd = gvd.iloc[:, :5]
    gvd['database'] = ['gvd'] * len(gvd)
    filtered = pd.concat([filtered, gvd], axis=0)
    if os.path.getsize(f"{output}/8.blastncontigs/mgv.out") != 0:
        mgv = pd.read_table(f"{output}/8.blastncontigs/mgv.out", header=None)
        mgv.columns = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'nident', 'qlen', 'slen', 'length', 'mismatch', 'positive', 'ppos', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qcovhsp', 'qcovus', 'qseq', 'sstrand', 'frames']
        mgv = mgv[mgv['pident']>=50]
        mgv = mgv[mgv['evalue']<=1e-10]
        mgv = mgv[mgv['qcovs']>=80]
        for i in range(len(mgv)):
            if mgv.iloc[i, 0] not in blastn:
                blastn.append(mgv.iloc[i, 0])
    mgv = mgv.iloc[:, :5]
    mgv['database'] = ['mgv'] * len(mgv)
    filtered = pd.concat([filtered, mgv], axis=0)
    if os.path.getsize(f"{output}/8.blastncontigs/ncbi.out") != 0:
        ncbi = pd.read_table(f"{output}/8.blastncontigs/ncbi.out", header=None)
        ncbi.columns = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'nident', 'qlen', 'slen', 'length', 'mismatch', 'positive', 'ppos', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qcovhsp', 'qcovus', 'qseq', 'sstrand', 'frames']
        ncbi = ncbi[ncbi['pident']>=50]
        ncbi = ncbi[ncbi['evalue']<=1e-10]
        ncbi = ncbi[ncbi['qcovs']>=80]
        for i in range(len(ncbi)):
            if ncbi.iloc[i, 0] not in blastn:
                blastn.append(ncbi.iloc[i, 0])
    ncbi = ncbi.iloc[:, :5]
    ncbi['database'] = ['ncbi'] * len(ncbi)
    filtered = pd.concat([filtered, ncbi], axis=0)
    if os.path.getsize(f"{output}/7.vircontigs/final-viral-score.tsv") != 0:
        dat = pd.read_table(f"{output}/7.vircontigs/final-viral-score.tsv", header=0)
        for i in range(len(dat)):
            if dat.iloc[i, 0] not in virso:
                virso.append(dat.iloc[i, 0].split('|')[0])
    info = pd.DataFrame(data=None, columns=['contig', 'blastn', 'virsorter'])
    num = 0
    with open(f"{output}/6.filter/contig_1k.fasta") as f:
        line = f.readline()
        if line == '':
            return 1
        f1 = open(f"{output}/9.final-contigs/contigs.fa", 'w')
        while 1:
            contig = line[1: -1]
            out = [contig, 0, 0]
            line = f.readline()
            seq = ''
            while line != '' and line[0] != '>':
                seq = seq + line[0: -1]
                line = f.readline()
            if contig in blastn or contig in virso:
                print(f">{contig}", file=f1)
                print(seq, file=f1)
                if contig in blastn:
                    out[1] = 1
                if contig in virso:
                    out[2] = 1
                info.loc[num] = out
                num += 1
            if line == '':
                break
        f1.close()
    info.to_csv(f"{output}/9.final-contigs/info.txt", header=True, index=False, sep='\t')
    filtered = filtered.sort_values(by=['qcovs', 'pident', 'evalue'],ascending= [False, False, True])
    filtered.drop_duplicates(subset=['qseqid'])
    filtered.to_csv(f"{output}/9.final-contigs/blastn_info.txt", header=True, index=False, sep='\t')
    return 0

def filter_checkv(output):
    dat = pd.read_table(f"{output}/10.checkv/quality_summary.tsv", header=0)
    dat1 = dat[dat["checkv_quality"] == 'Complete']
    dat2 = dat[dat["checkv_quality"] == 'High-quality']
    dat3 = dat[dat["checkv_quality"] == 'Medium-quality']
    checkv = pd.concat([dat1, dat2, dat3])['contig_id'].to_list()
    with open(f"{output}/9.final-contigs/contigs.fa") as f:
        f1 = open(f"{output}/11.high_quality/contigs.fa", 'w')
        while 1:
            line = f.readline()
            if line == '':
                break
            contig = line[1: -1]
            seq = f.readline()[:-1]
            if contig in checkv:
                print(f">{contig}", file=f1)
                print(seq, file=f1)
        f1.close()

def final_info(output):
    checkv = pd.read_table(f"{output}/10.checkv/quality_summary.tsv", header=0, index_col=0)
    blastn = pd.read_table(f"{output}/9.final-contigs/blastn_info.txt", header=0, index_col=0)
    info = pd.read_table(f"{output}/9.final-contigs/info.txt", header=0, index_col=0)
    subprocess.call([f"cat {output}/12.final_non_dup/final.fasta | grep '>' > {output}/12.final_non_dup/temp.txt"], shell=True)
    contig = []
    contig_blastn = []
    with open(f"{output}/12.final_non_dup/temp.txt", "r") as fi:
        while 1:
            line = fi.readline()
            if line == '':
                break
            contig.append(line[1: -1])
    for ct in contig:
        if info.loc[ct, 'blastn'] == 1:
            contig_blastn.append(ct)
    checkv.loc[contig, "checkv_quality"].to_csv(f"{output}/12.final_non_dup/completeness.txt", header=True, index=True, sep='\t')
    info.loc[contig].to_csv(f"{output}/12.final_non_dup/info.txt", header=True, index=True, sep='\t')
    blastn.loc[contig_blastn].to_csv(f"{output}/12.final_non_dup/blastn_info.txt", header=True, index=True, sep='\t')
    subprocess.call([f"rm {output}/12.final_non_dup/temp.txt"], shell=True)

def create_output_file(output):
    if os.path.exists(output) is False:
        subprocess.call([f"mkdir {output}"], shell=True)

def run_fastqc(output, input1, input2, threads):
    print("Run fastqc")
    if os.path.exists(f"{output}/1.fastqc") is True:
        subprocess.call([f"rm -rf {output}/1.fastqc"], shell=True)
    subprocess.call([f"mkdir {output}/1.fastqc"], shell=True)
    print(f"fastqc {input1} {input2} -t {threads} -o {output}/1.fastqc")
    ret = subprocess.call([f"fastqc {input1} {input2} -t {threads} -o {output}/1.fastqc"], shell=True)
    if ret != 0:
        print("Warning: fastqc error")
    
def run_trim(output, threads, input1, input2, sample1, sample2, adapter):
    print("Run trim")
    if os.path.exists(f"{output}/2.trim") is True:
        subprocess.call([f"rm -rf {output}/2.trim"], shell=True)
    subprocess.call([f"mkdir {output}/2.trim"], shell=True)
    print(f"trimmomatic PE -threads {threads} -phred33 {input1} {input2} {output}/2.trim/{sample1}.fastq {output}/2.trim/{sample1}_single.fastq {output}/2.trim/{sample2}.fastq {output}/2.trim/{sample2}_single.fastq ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50")
    ret = subprocess.call([f"trimmomatic PE -threads {threads} -phred33 {input1} {input2} {output}/2.trim/{sample1}.fastq {output}/2.trim/{sample1}_single.fastq {output}/2.trim/{sample2}.fastq {output}/2.trim/{sample2}_single.fastq ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50"], shell=True)
    if ret != 0:
        sys.exit("Error: trimmomatic error")

def run_bowtie2(output, threads, sample1, sample2, index_path):
    if len(index_path) == 0:
        print("No need for running bowtie2")
    else:
        print("Run bowtie2")
    if os.path.exists(f"{output}/3.bowtie2") is True:
        subprocess.call([f"rm -rf {output}/3.bowtie2"], shell=True)
    subprocess.call([f"mkdir {output}/3.bowtie2"], shell=True)
    subprocess.call([f"cp {output}/2.trim/{sample1}.fastq {output}/3.bowtie2/{sample1}.fastq"], shell=True)
    subprocess.call([f"cp {output}/2.trim/{sample2}.fastq {output}/3.bowtie2/{sample2}.fastq"], shell=True)
    for na in index_path:
        print(f"bowtie2 -p {threads} -x {na} -1 {output}/3.bowtie2/{sample1}.fastq -2 {output}/3.bowtie2/{sample2}.fastq --un-conc {output}/3.bowtie2/tmp > {output}/3.bowtie2/tmp.sam")
        ret = subprocess.call([f"bowtie2 -p {threads} -x {na} -1 {output}/3.bowtie2/{sample1}.fastq -2 {output}/3.bowtie2/{sample2}.fastq --un-conc {output}/3.bowtie2/tmp > {output}/3.bowtie2/tmp.sam"], shell=True)
        if ret != 0:
            sys.exit("Error: bowtie2 error")
        subprocess.call([f"mv {output}/3.bowtie2/tmp.1 {output}/3.bowtie2/{sample1}.fastq"], shell=True)
        subprocess.call([f"mv {output}/3.bowtie2/tmp.2 {output}/3.bowtie2/{sample2}.fastq"], shell=True)
    if len(index_path) != 0:
        subprocess.call([f"rm {output}/3.bowtie2/tmp.sam"], shell=True)
    subprocess.call([f"rm -rf {output}/2.trim"], shell=True)

def run_viromeQC(output, sample1, sample2):
    print("Run viromeQC")
    if os.path.exists(f"{output}/4.viromeQC") is True:
        subprocess.call([f"rm -rf {output}/4.viromeQC"], shell=True)
    subprocess.call([f"mkdir {output}/4.viromeQC"], shell=True)
    print(f"python {sys.path[0]}/viromeqc/viromeQC.py -i {output}/3.bowtie2/{sample1}.fastq {output}/3.bowtie2/{sample2}.fastq -o {output}/4.viromeQC/report.txt")
    ret = subprocess.call([f"python {sys.path[0]}/viromeqc/viromeQC.py -i {output}/3.bowtie2/{sample1}.fastq {output}/3.bowtie2/{sample2}.fastq -o {output}/4.viromeQC/report.txt"], shell=True)
    if ret != 0:
        print("Warning: viromeQC error")
    
def run_spades(output, threads, sample1, sample2):
    print("Run spades")
    if os.path.exists(f"{output}/5.spades") is True:
        subprocess.call([f"rm -rf {output}/5.spades"], shell=True)
    subprocess.call([f"mkdir {output}/5.spades"], shell=True)
    print(f"spades.py -1 {output}/3.bowtie2/{sample1}.fastq -2 {output}/3.bowtie2/{sample2}.fastq -t {threads} -o {output}/5.spades")
    ret = subprocess.call([f"spades.py -1 {output}/3.bowtie2/{sample1}.fastq -2 {output}/3.bowtie2/{sample2}.fastq -t {threads} -o {output}/5.spades"], shell=True)
    if ret != 0:
        sys.exit("Error: spades error") 

def run_vsearch_1(output, sample):
    print("Run vsearch (trim short contigs)")
    if os.path.exists(f"{output}/6.filter") is True:
        subprocess.call([f"rm -rf {output}/6.filter"], shell=True)
    subprocess.call([f"mkdir {output}/6.filter"], shell=True)
    print(f"vsearch --sortbylength {output}/5.spades/scaffolds.fasta --minseqlength 1000 --maxseqlength -1 --relabel s{sample}.contig --output {output}/6.filter/contig_1k.fasta")
    ret = subprocess.call([f"vsearch --sortbylength {output}/5.spades/scaffolds.fasta --minseqlength 1000 --maxseqlength -1 --relabel s{sample}.contig --output {output}/6.filter/contig_1k.fasta"], shell=True)
    if ret != 0:
        sys.exit("Error: vsearch error")

def run_virsorter(output, threads):
    print("Run virsorter")
    if os.path.exists(f"{output}/7.vircontigs") is True:
        subprocess.call([f"rm -rf {output}/7.vircontigs"], shell=True)
    subprocess.call([f"mkdir {output}/7.vircontigs"], shell=True)
    print(f"virsorter run -w {output}/7.vircontigs -i {output}/6.filter/contig_1k.fasta -j {threads} all")
    ret = subprocess.call([f"virsorter run -w {output}/7.vircontigs -i {output}/6.filter/contig_1k.fasta -j {threads} all"], shell=True)
    if ret != 0:
        sys.exit("Error: virsorter error")
    
def run_blastn(output, threads):
    print("Run blastn")
    if os.path.exists(f"{output}/8.blastncontigs") is True:
        subprocess.call([f"rm -rf {output}/8.blastncontigs"], shell=True)
    subprocess.call([f"mkdir {output}/8.blastncontigs"], shell=True)
    print(f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/crass -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/crass.out')
    ret = subprocess.call([f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/crass -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/crass.out'], shell=True)
    if ret != 0:
        sys.exit("blastn error")
    print(f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/gpd -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/gpd.out')
    ret = subprocess.call([f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/gpd -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/gpd.out'], shell=True)
    if ret != 0:
        sys.exit("blastn error")
    print(f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/gvd -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/gvd.out')
    ret = subprocess.call([f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/gvd -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/gvd.out'], shell=True)
    if ret != 0:
        sys.exit("blastn error")
    print(f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/mgv -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/mgv.out')
    ret = subprocess.call([f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/mgv -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/mgv.out'], shell=True)
    if ret != 0:
        sys.exit("blastn error")
    print(f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/ncbi -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/ncbi.out')
    ret = subprocess.call([f'blastn -query {output}/6.filter/contig_1k.fasta -db {sys.path[0]}/blastn_database/ncbi -num_threads {threads} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident evalue qcovs nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send bitscore qcovhsp qcovus qseq sstrand frames " -out {output}/8.blastncontigs/ncbi.out'], shell=True)
    if ret != 0:
        sys.exit("Error: blastn error")

def run_combination(output):
    print("Combine virsorter and blastn results")
    if os.path.exists(f"{output}/9.final-contigs") is True:
        subprocess.call([f"rm -rf {output}/9.final-contigs"], shell=True)
    subprocess.call([f"mkdir {output}/9.final-contigs"], shell=True) 
    ret = filter_vircontig(output)
    if ret != 0:
        sys.exit("Error: combine error")

def run_checkv(output, threads):
    print("Run checkv")
    if os.path.exists(f"{output}/10.checkv") is True:
        subprocess.call([f"rm -rf {output}/10.checkv"], shell=True)
    subprocess.call([f"mkdir {output}/10.checkv"], shell=True)
    print(f"checkv end_to_end {output}/9.final-contigs/contigs.fa {output}/10.checkv -d {sys.path[0]}/checkv_database -t {threads}")
    ret = subprocess.call([f"checkv end_to_end {output}/9.final-contigs/contigs.fa {output}/10.checkv -d {sys.path[0]}/checkv_database -t {threads}"], shell=True)
    if ret != 0:
        sys.exit("Error: checkv error")

def high_quality_output(output):
    print("Get final output")
    if os.path.exists(f"{output}/11.high_quality") is True:
        subprocess.call([f"rm -rf {output}/11.high_quality"], shell=True)
    subprocess.call([f"mkdir {output}/11.high_quality"], shell=True) 
    filter_checkv(output)

def run_vsearch_2(output, threads):
    print("Run vsearch (cluster)")
    if os.path.exists(f"{output}/12.final_non_dup") is True:
        subprocess.call([f"rm -rf {output}/12.final_non_dup"], shell=True)
    subprocess.call([f"mkdir {output}/12.final_non_dup"], shell=True)
    print(f"vsearch --cluster_fast {output}/11.high_quality/contigs.fa --id 0.995 --centroids {output}/12.final_non_dup/final.fasta --uc {output}/11.high_quality/clusters.uc --maxseqlength -1 --threads {threads}")
    ret = subprocess.call([f"vsearch --cluster_fast {output}/11.high_quality/contigs.fa --id 0.995 --centroids {output}/12.final_non_dup/final.fasta --uc {output}/11.high_quality/clusters.uc --maxseqlength -1 --threads {threads}"], shell=True)
    if ret != 0:
        sys.exit("Error: vsearch error")
    final_info(output)
  
def remove_inter_result(output):
    subprocess.call([f"rm -rf {output}/3* {output}/5* {output}/6* {output}/7* {output}/8* {output}/9* {output}/10* {output}/11*"], shell=True)

