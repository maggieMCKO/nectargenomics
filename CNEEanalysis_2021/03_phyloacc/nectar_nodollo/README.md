## set up phyloacc env
```bash
module load rev/23.12 anaconda3/2023.09-0 
conda env list
conda create -n phyloacc220
conda activate phyloacc220
conda install phyloacc
phyloacc.py --version  # PhyloAcc version 2.2.0 released on April 20, 2023
phyloacc.py --depcheck # All dependencies PASSED.
```

## setup phyloacc and snakemake batches
```bash
module load rev/23.12 anaconda3/2023.09-0 
conda activate phyloacc220

export sub_proj="nectar_nodollo"                            ## <----- sub-proj
cd "$HOME/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/${sub_proj}/"

export scratch_DIR="/scratch/users/$USER/phyloacc2024/"
mkdir -p ${scratch_DIR}

export neutral_model_named="../inputs/mod/nonconserved_4d_named_sp44.mod"
export input_cnee_fa="../inputs/mod/galloseq_gapFixed_sp44.fa"
export input_part_bed="../inputs/mod/galloseq.part.bed"
export input_sp_table="../inputs/mod/commonN.txt"

export MYSCRATCH_input="${scratch_DIR}${sub_proj}/input/"     ## <----- sub-proj
export MYSCRATCH_output="${scratch_DIR}${sub_proj}/output/"   ## <----- sub-proj
mkdir -p ${MYSCRATCH_input} ${MYSCRATCH_output}
export scratch_neutral_model_named="${MYSCRATCH_input}nonconserved_4d_named_sp44.mod"
export scratch_input_cnee_fa="${MYSCRATCH_input}galloseq_gapFixed_sp44.fa"
export scratch_input_part_bed="${MYSCRATCH_input}galloseq.part.bed"

# put named.mod, fasta, cnee bed here
Input_DIR="input/"
mkdir -p ${Input_DIR}
cp ${neutral_model_named} ${Input_DIR}
cp ${input_cnee_fa} ${Input_DIR}
cp ${input_part_bed} ${Input_DIR}
cp ${input_sp_table} ${Input_DIR}

# copy to scratch
cp -r ${Input_DIR} "${scratch_DIR}${sub_proj}"               ## <----- sub-proj

phyloacc.py -a ${scratch_input_cnee_fa} -b ${scratch_input_part_bed} -m ${scratch_neutral_model_named} -t "HLlepAsp1;HLdicExi1;HLlicPen1;HLlicMelCas1;HLphyNov1;HLgraPic1;HLlorGal1;HLtriMol2;HLaraSol1;HLamaAes1;HLphaSup1;HLfloFus1;HLcalAnn5" -r st -burnin 1000 -mcmc 7000 -n 1 -p 12 -j 24 -batch 100 -mem 72 -part "medium" -time 3
```

## run snakemake
```bash
screen -S nect_without

snakemake -k -p -s /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/phyloacc-out-02-07-2024.11-55-54/phyloacc-job-files/snakemake/run_phyloacc.smk --configfile /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/phyloacc-out-02-07-2024.11-55-54/phyloacc-job-files/snakemake/phyloacc-config.yaml --profile /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/phyloacc-out-02-07-2024.11-55-54/phyloacc-job-files/snakemake/profiles/slurm_profile --dryrun
```

## gather outputs after all snakemake jobs are done
```bash
phyloacc_post.py -i phyloacc-out-02-07-2024.11-55-54
```

