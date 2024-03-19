# Introduction to HPC and Bioinformatics CAC 2024

# Class Activity #1 – Installing python packages

Searching for available python packages on Alliance and Frontenac cluster.
```
avail_wheels numpy
avail_wheels numpy --all-version
avail_wheels numpy pandas --all-version
```
DendroPy is a Python library for phylogenetic computing https://dendropy.org/index.html

Find all available versions of `dendropy` and `pyfasta` in the wheelhouse

Install `dendropy` and `pyfasta`

Run `pip freeze` to generate `requirements.txt` file

Check your installation:
```
python -c "import dendropy"
python -c "import dendropy; print(dendropy.__version__)"
```
Create a new virtual environment and install packages using 'requirements.txt'

# Class Activity #2 - Installing software using Apptainer

Search for containers https://hub.docker.com/

If you decide to build an image make sure to set `cache` and `tmp` for `apptainer` in a filesystem where you have read/write permission 

```
mkdir -p /scratch/$USER/apptainer/{cache,tmp}
export APPTAINER_CACHEDIR="/scratch/$USER/apptainer/cache"
export APPTAINER_TMPDIR="/scratch/$USER/apptainer/tmp"
```

Some tools such as `qiime2` are hosted on their own repository https://quay.io/repository/qiime2/core?tab=tags

````
module load apptainer
apptainer build qiime2-2023.5.sif docker://quay.io/qiime2/core:2023.5
````

Recipe file/definition file
```
Bootstrap: docker
From: ubuntu:22.04
Stage: build
%post
    apt-get update && apt-get install -y git
```

Now build a container from the recipe file/definition file

```
apptainer build ubuntu_test_git.sif my_test_def_file.def

apptainer build --sandbox ubuntu_sandbox ubuntu_test_git.sif

```
# Class Activity #3 – Using mdsum

Lets download a file and verify integrity
     
```
wget https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz.md5
```
Verify integrity
```
md5sum -c 18S_fungal_sequences.tar.gz.md5
```
   
Split fasta file using pyfasta (hint: source ~/ENV/bin/activate)
```
pyfasta split -n 6 Triticum_aestivum_subset.IWGSC.cds.all.fa
```
# Class Activity #4 – Text editing

Clone repository `git clone https://github.com/sridhar-ravi/Introduction_to_HPC_Bioinformatics_2024.git`

View fasta header

```
grep "^>" Triticum_aestivum_subset.IWGSC.cds.all.fa
```
Remove spaces from the reference seq
```
sed 's, ,|,g' 
```
Print transcript id and chromosomal position
```
awk -F '[ |,]' '{ if ($0 ~ /^>/) { print $1"|"$3;} else { print $0;}}'
```
# Class Activity #6 – Job submission

Copy *.sam files from `/global/project/Workshop2023/IntroBioInfo/`
