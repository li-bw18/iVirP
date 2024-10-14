# iVirP (Integrative virome pipeline) #
 
## Description ##

This pipeline can preprocess the raw data (fastq files) of VLP sequencing and find some viral contigs that are reliable with some useful software. In addition, after preprocessing, this pipeline can also estimate the abundance of eukaryotic viruses rapidly with Kraken2.

## Attentions ##

In the current version, for easy usage, this pipeline only supports paired-end fastq files in the phred33 format. Additionaly, the branching function only supports fastq files with 150 reads length. If you need to adapt to other fastq files, feel free to modify the functions in the "utils.py" to avoid errors. We will soon support more other formats in the future.

## Steps to install this pipeline (Internet connection required) ##

### Step 1: prepare the conda environment ###

(1) Download the iVirP software from github

`git clone https://github.com/li-bw18/iVirP.git`

(2) Go to the directory of iVirP

`cd iVirP` (for example)

(3) Run

`conda env create -f iVirP.yaml` (perhaps needs 30-60min)

### Step 2: get into the prepared conda environment ###

`conda activate ivirp`

### Step 3: install some built-in software and download some databases  ###

When you meet any bug in this step, please check the Internet connection and then rerun this step. If the bug is not solved or you meet any bug after finishing this step, please remove the whole directory(iVirP), download a new one from github or somewhere, and then just rerun this script.

We provide two different scripts to run this step:

(1) More sub-steps version:

`python install.py -t 4` (-t determines the threads used to do this step)

(2) Fewer sub-steps version (For those who can download zenodo data fast, we strongly recommend you to use this one):

`python install_easy.py -t 4` (-t determines the threads used to do this step)

<strong>Please use only one version even if you fail to finish this step once. If you really need to change into another version, please remove the whole directory(iVirP), download a new one from github or somewhere, and then just rerun this script.</strong>

### Step 4: build Bowtie2 indexes for removing host sequences  ###

We use Bowtie2 to remove the host sequences, however, the building of host Bowtie2 indexes is both time- and space-consuming. As a result, we did not prepare any available Bowtie2 index during the above installing steps for your personalized analysis.

So this step is designed to build Bowtie2 indexes of host sequences, please provide the genomic fasta file of the host you want to remove, and run the following script:

`python build_index.py -f {directory}/{name}.fasta -t 4` (-t determines the threads used to do this step, file suffixes with '.fa' or '.fna' are also OK)

After this step, the Bowtie2 index of this host will be permanently saved in the "bowtie2_index" directory, and this host is available for host removal. Then, when running this pipeline, you can choose to remove this host with parameter: "--host {name}" ({name} is the name of the fasta file you provided without the file suffix).

For example, after you build a Bowtie2 index of human with the following command:

`python build_index.py -f /home/test/human.fasta -t 4`

you can use '--host human' to remove human host contaminants in this pipeline (find more details in the following instructions).

In addition, this pipeline can save more than one Bowtie2 index, and you can also choose to remove more than one host when running (Make sure you have built their indexes; see -h or the following instructions for more information). As a result, you can add novel or upgrade old Bowtie2 indexes anytime in the future with the above command.

`python build_index.py -r -n {name}` is used to remove an existing Bowtie2 index.

`python build_index.py -c` is used to check the built Bowtie2 indexes.

Do not know how to download the genomic fasta file of a certain host? Here are some recommended web databases for you to use:

(1) NCBI: https://www.ncbi.nlm.nih.gov/genome

(2) GENCODE: https://www.gencodegenes.org (human and mouse)

(3) Ensembl: https://www.ensembl.org

## How to use this pipeline ##

You can use this pipeline by running the "run.py" script in the "iVirP" install directory with several parameters.

Here is a very basic example:

### Step 1: check the host Bowtie2 indexes ###

Run `python build_index.py --check`, and see whether all host Bowtie2 indexes you need are available.

If some of them are not available, see the step 4 of installing above to build them.

### Step 2: run this pipeline ###

`python {iVirP_install_path}/run.py A_1.fq.gz A_2.fq.gz -t 4 --host human,mouse -a 0 -o output`

"A_1.fq.gz" and "A_2.fq.gz" are paired-end sequencing fastq files you want to process.

-o determines the path of the directory to save the intermediate and final result files.

-t determines the threads used to run this pipeline (default:1).

-a determines the adapter file used by trimmomatic to remove adapter sequences (input the interger code, 0:NexteraPE-PE, 1:TruSeq3-PE-2, 2:TruSeq3-PE, 3:TruSeq2-PE, default:NexteraPE-PE).

--host determines name(s) of Bowtie2 index(es) built before, please separate different names with "," when removing more than one hosts (for example: "human,mouse") (default: None, which means this step will be skipped).

The final viral contigs found by this pipeline are saved in "{output_dic}/11.final/contigs.fa".

<strong>If the process of this pipeline is killed for any reason, you can try to continue to run from the killed step instead of the very first step by giving an additional parameter: "-k", please don't change any file or directory in the output directory if you want to use this function.</strong>

Notice: During recent tests, we find that some processes stop because all data has been filtered out and the file has become empty at a certain step, please be careful about this kind of problem when your raw data is quite limited, and we are considering how to deal with these problems in the future version.

For more detailed information about the parameters, please run: 

`python {iVirP_install_path}/run.py -h` or `python {iVirP_install_path}/run.py --help`

## Only run the QC step ##

We also allow only running the QC steps of this pipeline with the "QC.py" script, including the following four steps:

(1) report the quality of sequencing with fastqc

(2) trim adapters with trimmomatic

(3) remove host sequences with bowtie2

(4) evaluate the QC with viromeQC

For example:

`python {iVirP_install_path}/QC.py A_1.fq.gz A_2.fq.gz -t 4 --host human -a 0 -o output`

"A_1.fq.gz" and "A_2.fq.gz" are paired-end sequencing fastq files you want to process.

-o determines the path of the directory to save the intermediate and final result files.

-t determines the threads used to run this pipeline.

-a determines the adapter file used by trimmomatic to remove adapter sequences (input the interger code, 0:NexteraPE-PE, 1:TruSeq3-PE-2, 2:TruSeq3-PE, 3:TruSeq2-PE, default:NexteraPE-PE).

--host determines name(s) of Bowtie2 index(es) built before, please separate different names with "," when removing more than one hosts (for example: "human,mouse") (default: None, which means this step will be skipped).

And you can also skip at least one of the step (1)/(2)/(4) with parameters:

"--skip_fastqc"/"--skip_trim"/"--skip_viromeqc"

Remember: don't forget building bowtie2 indexes before running this script and the cleaned fastq file with be in the directory of "{output_path}/3.bowtie2"

## Discriptions of output files ##

When '-r' parameter is not given, all the intermediate results will be saved in the output directory. Here is a brief introduction of some important files:

### <strong>1.fastqc (step1 output directory)</strong> ###

1.fastqc/{sample}_1(and 2)_fastqc.html(or .zip): results of fastqc

### <strong>3.bowtie2 (step3 output directory)</strong> ###

3.bowtie2/{sample}_1(and 2).fastq: fastq files after removing the sequences of hosts or other contaminations

### <strong>4.viromeQC (step4 output directory)</strong> ###

4.viromeQC/report.txt: results of viromeQC

### <strong>5.spades (step5 output directory)</strong> ###

5.spades/scaffolds.fasta: contigs assembled by spades

### <strong>6.filter (step6 output directory)</strong> ###

6.filter/contig_1k.fasta: contigs with length > 1k nt

### <strong>7.vircontigs (step7 output directory)</strong> ###

7.vircontigs/final-viral-combined.fa: virsorter filtered viral contigs

### <strong>9.final-contigs (step8-9 output directory)</strong> ###

9.final-contigs/contigs.fasta: combined viral contigs from virsorter and blastn

### <strong>11.high_quality (step10 output directory)</strong> ###

11.high_quality/contigs.fa: high-quality viral contigs filtered with checkV

### <strong>12.final_non_dup (step11 output directory)</strong> ###

12.final_non_dup/final.fasta: final non-redundant viral contigs filtered with vsearch (cluster with similarity > 0.995)

12.final_non_dup/completeness.txt: the quality (completeness) of final viral contigs (checkV report)

12.final_non_dup/info.txt: the information about whether a certain contig is found by virsorter or blastn

12.final_non_dup/blastn_info.txt: the detailed information of viral contigs found by blastn (for those contigs with mutiple alignments, only the information of the best alignment will be included)

## Additional function - rapid detection of eukaryotic viruses ##

This function can rapidly give the abundance of some common eukaryotic viruses with Kraken2+Bracken from your raw VLP sequencing fastq files.

<strong>Notice: this function is not that accurate now, please use it carefully!</strong>

To use it, please first do the following step to install this additional function:

(1) Go to the directory of iVirP

`cd iVirP` (for example)

(2) Run

`python download_kraken.py -t 4` (-t determines the threads used to do this step)

Then you can use this function with the following two steps:

### Step 1: build the bracken database ###

Run

`python build_bracken.py -l 150 -t 4` (-l determines the length of your reads; -t determines the threads used to do this step)

Notice: you only need to run this step once before the length of your reads changes; when the length of your reads changes, please rerun this step with new length.

### Step 2: run this function ###

Run

`python run_fast_detection.py -l 150 A_1.fq.gz A_2.fq.gz -t 4 --host human -a 0 -o output` 

Most parameters are the same as those of 'run.py' (main pipeline), because we also run the step1-step4 of the main pipeline for clean reads. -l determines the length of your reads.

Then you can get the final results of the estimated abundance of eukaryotic viruses at different taxonomic levels.

## Just use some software ##

We also welcome those who want to use any software in our pipeline to use our pipeline for installing some software more easily with our simple installing steps, when you have finished the installing steps 1-3, you can use the following software in our "ivirp" conda environment:

checkV: run directly with command `chechv` 

virsorter2: run directly with command `virsorter` 

trimmomatic: run directly with command `trimmomatic` (adapters are at "iVirP/adapters") 

viromeQC: run the script "iVirP/viromeqc/viromeQC.py" with command `python iVirP/viromeqc/viromeQC.py`

## Contact with us ##

(1) Github Issues

(2) Email: lbw22@mails.tsinghua.edu.cn
