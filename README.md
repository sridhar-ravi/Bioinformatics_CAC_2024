# Introduction to Bioinformatics CAC 2024

# Class Activity #1 – Standard Environment and modules
ssh to your account on Frontenac or one of the Alliance cluster eg: Nibi, Fir, Narval or Rorqual. If you don't have an Alliance account, use the guest or student account. 

```
ssh -X hpcxxxx@login.cac.queensu.ca
ssh -X sa13xxxx@login.cac.queensu.ca
ssh -X username@nibi.alliancecan.ca
ssh -X userxxxx@coss-a.c3.ca
```

Create a folder in your `home` or `scratch` called bioinformatics

Type `$module list` to find default Standard Environment​

`StdEnv/2023` uses GCC 12.3.0, Intel 2023.1, and Open MPI 4.1.5 as defaults

`StdEnv/2020` uses GCC 9.3.0, Intel 2020.1, and Open MPI 4.0.3 as defaults

Type the following command in your terminal to see which versions of `samtools` avaible on out software stack
```
module spider samtools
module load samtools
module list
```
Now switch to `samtools` version `1.18` using the `module load` command.
Try loading module `blast+`
You can also use module spider with wildcard `module -r spider '.*blast.*'`. Now let’s see if we have a module called “Bioconductor”

Let try loading `busco` and use command `module list` to view all loaded modules. Pay attention to `StdEnv`.

Useful resources:

https://docs.alliancecan.ca/wiki/Standard_software_environments

https://docs.alliancecan.ca/wiki/Utiliser_des_modules/en

# Class Activity #2 – screen terminal multiplexer

`screen -S python_install` to start a screen session

`screen -ls` list active screen session

`ctrl` + a + d to minimize and run process in the background

`screen -r python_install` attach screen session

# Class Activity #3 – Installing python packages

Searching for available python packages on Alliance and Frontenac cluster.
```
avail_wheels numpy
avail_wheels numpy --all-version
avail_wheels numpy pandas biopython --all-version
```

Find all available versions of `biopython` and `fastasplit` in the wheelhouse

Install `biopython` and `fastasplit` on python/3.10 with option `-- no-index'

Run `pip freeze` to generate `requirements.txt` file

Check your installation:
```
python -c "import Bio"
python -c "import Bio; print(Bio.__version__)"

```
Create a new virtual environment and install packages using `requirements.txt`

Possible error while resolving dependencies
```
ERROR: Could not find a version that satisfies the requirement torch (from versions: none)
ERROR: No matching distribution found for torch
```

`pip install -r requiremets.txt --no-deps`

https://docs.alliancecan.ca/wiki/Python/en

To make the virtual environment visible on Jupyterhub portal

```
mkdir -p ~/.local/share/jupyter/kernels
pip install --no-index ipykernel
python -m ipykernel install --user --name my_py_env --display-name "pipeline_1 python kernel"
```
JupyterHub portal on Frontenac https://ood-frnt.cac.queensu.ca/pun/sys/dashboard. On the Aliiance clusters use https://docs.alliancecan.ca/wiki/JupyterHub#JupyterHub_on_clusters or https://ondemand.sharcnet.ca/

# Class Activity #4 - Installing R packages

By default, R packages are installed under `$HOME/R/`

```
module load r/4.3.1
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/
R -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")'
```
To install packages after loading R
```
R
> install.packages("ggplot2")
```
Open a R session and view sourced library path

```
> .libPaths()
[1] "/global/teaching-home/sa4898003/.local/R/4.3.1"
[2] "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.3.1/lib64/R/library"
```
Add a custom library path

`> .libPaths(c("/path/to/your/library", .libPaths()))`

RStudio portal on Frontenac `https://ood-frnt.cac.queensu.ca/pun/sys/dashboard`. On the Aliiance clusters use https://docs.alliancecan.ca/wiki/JupyterHub#JupyterHub_on_clusters or https://ondemand.sharcnet.ca/

# Class Activity #5 - Installing software using Apptainer

Search for containers https://hub.docker.com/

If you decide to build an image make sure to set `cache` and `tmp` for `apptainer` in a filesystem where you have read/write permission 

```
mkdir -p /global/scratch/$USER/apptainer/{cache,tmp}
export APPTAINER_CACHEDIR="/global/scratch/$USER/apptainer/cache"
export APPTAINER_TMPDIR="/global/scratch/$USER/apptainer/tmp"
```

Some tools such as `qiime2` are hosted on their own repository https://quay.io/repository/qiime2/core?tab=tags

````
module load apptainer
apptainer build qiime2-2023.5.sif docker://quay.io/qiime2/core:2023.5
````
Exceptions include licensed software that does not allow distribution eg https://docs.nvidia.com/clara/parabricks/latest/index.html

````
module load apptainer
apptainer build parabrick-4.6.0.sif docker://nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
````

Recipe file/definition file
```
Bootstrap: docker
From: ubuntu:22.04
Stage: build
%post

  apt-get update && apt-get -y upgrade
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
  export DEBIAN_FRONTEND=noninteractive
  apt-get -y install fastani

```

Now build a container from the recipe file/definition file

```
apptainer build fastani.sif fastani.def

apptainer shell fastani.sif
Apptainer> fastANI --version

```
For more information https://docs.alliancecan.ca/wiki/Apptainer/en
# Class Activity #6 – Using mdsum

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
git clone https://github.com/sridhar-ravi/Bioinformatics_CAC_2025.git
cd Bioinformatics_CAC_2025/
pyfasta split -n 6 Triticum_aestivum_subset.IWGSC.cds.all.fa
md5sum *.fastq > md5.txt​
md5sum –c md5.txt​
sha256sum *.fastq > sha256.txt​
sha256sum –c sha256.txt​
```

# Class Activity #7 – Text editing

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
# Class Activity #8 – Job submission

Copy *.sam files from `/global/teaching-project/sg4898000/IntroBioinformtaics`
