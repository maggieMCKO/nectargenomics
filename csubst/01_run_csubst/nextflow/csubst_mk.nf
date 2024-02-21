#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

params.conda_env_path = "/usr/users/mko/.conda/envs/env_csubst39_14"

params.foreground = "$HOME/Nectar/analysis/csubst/00_inputs/foreground.tsv"
params.tree_file="$HOME/Nectar/analysis/csubst/00_inputs/new_sp_tree.nw"
params.sp_ls="$HOME/Nectar/analysis/csubst/00_inputs/new_sp_list.txt"
params.aln_input="$HOME/Nectar/analysis/csubst/00_inputs/nt_align_renamed"
params.transcript = "$HOME/Nectar/analysis/csubst/00_inputs/checked.transcripts.lst"
// params.transcript = "$HOME/Nectar/analysis/csubst/01_run_csubst_2024_proto/sunbrid_affected.txt"

params.prune_tree_script = "$projectDir/prune_ete4.py"
params.outdir = "results"

params.size = 1000

workflow {

    query_ch = Channel.fromPath( params.transcript ).splitText() { it.strip() }.buffer(size: params.size, remainder: true)
    
    // run 1
    ch1 = rm_sunbird( query_ch ) 
    ch1.flatten().set{ ch_aln }

    // run 2
    ch2 = prune_tree( ch1 ) 
    ch2.flatten().set{ ch_tre }

    // prep for 3
    ch_aln
        .map { [it.toString().split("/").last().split("rm_HLcinPul1_")[1].split(".fa")[0], it] }
        .set { ch_aln }
    ch_tre
        .map { [it.toString().split("/").last().split("pruned_")[1].split(".nw")[0], it] }
        .set { ch_tre }
    // ch_aln.view()
    // ch_tre.view()
    
    // run 3
    ch_aln.combine(ch_tre, by:0) | runcsubst
}

process rm_sunbird {
    tag "rm_sunbird $transcript"
    publishDir params.outdir, mode : 'copy' 

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        val id
    
    output:
        // tuple val ("${INDEX}"), path("aln/rm_HLcinPul1_${INDEX}.fa"), emit: aln // doesnt work
        path("aln/rm_HLcinPul1_*.fa"), emit: aln

    """
    for INDEX in ${id.join(' ')}; do
        echo "task \${INDEX}"
        seqkit grep -f ${params.sp_ls} "${params.aln_input}/\${INDEX}/\${INDEX}_renamed.fa" -o "aln/rm_HLcinPul1_\${INDEX}.fa"
    done
    """
}

process prune_tree {
    tag "prune_tree $transcript"
    publishDir params.outdir, mode : 'copy' 

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        val alns
    
    output:
        // path "tmp.txt", emit: pruned_tree // for test
        path ("tree/pruned_*.nw"), emit: pruned_tree
        // tuple val (stem), path ("tree/pruned_*.nw"), emit:pruned_tree // doesnt work

    """    
    for aln in ${alns.join(' ')}; do
        stem=\$(basename \${aln} | sed "s/rm_HLcinPul1_//")
        stem=\${stem%%.fa}
        mkdir -p tree
        python ${params.prune_tree_script} ${params.tree_file} \$aln "tree/pruned_\${stem}.nw"
    done
    """
}

process runcsubst {
    tag "csubst $transcript"
    publishDir params.outdir, mode : 'copy',
    saveAs: {
        filename -> "csubst_out/$id/$filename"
    }

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        tuple val(id), path(aln), path(tree)
    
    output:
        path "csubst_*"
        
    """
    csubst analyze \
    --alignment_file $aln \
    --rooted_tree_file $tree \
    --foreground ${params.foreground} \
    --fg_exclude_wg no \
    --max_arity 10 \
    --exhaustive_until 2 &> csubst_${id}.log
    """
}
