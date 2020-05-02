
#### Phylofit
http://compgen.cshl.edu/phast/phyloFit-tutorial.php
1. prepare input
    - maf
    - tree: NectarTree01.nw
    - chicken gff
        - download
        ```bash
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
        gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz
        ```
        - prepare acc
          - https://www.ncbi.nlm.nih.gov/assembly/GCF_000002315.6
          - create galgal6a_2315.6_acc.tsv
            - col 1: from
            - col 2: to
        - run ```fix_galgalanno.sh```, this script does:
            - run Tim's perl script 'replace_chrs.pl'
                - replace NC_xxx to chr
            - fix six rows with 'Curated  Genomic' which messed up the gff, replace it to 'Curated_Genomic'
2. run ```run_phylofit_de.sh```
