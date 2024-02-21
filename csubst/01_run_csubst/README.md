## 1. setup conda env for csubst 
```bash
module load rev/23.12 anaconda3/2023.09-0
conda create --name env_csubst39_14 python=3.9 
conda activate env_csubst39_14 

# Installation dependencies and csubst
conda install iqtree # 2.2.6
conda install -c conda-forge pymol-open-source 
conda install bioconda::seqkit
conda install bioconda::mafft # v7.520; 

pip install git+https://github.com/kfuku52/csubst.git@v1.4.0 
# pip install numpy cython ete3 pyvolve matplotlib pypdb biopython # dependencies, but should be all installed by installing csubst. 
pip install https://github.com/etetoolkit/ete/archive/ete4.zip
```
## 2. run
```bash
cd /home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst/nextflow
screen -S csubst
module load rev/23.12 anaconda3/2023.09-0
conda activate env_csubst39_14 
module load nextflow/23.10.0
nextflow csubst_mk.nf -profile gwdg -c nextflow.config 
```