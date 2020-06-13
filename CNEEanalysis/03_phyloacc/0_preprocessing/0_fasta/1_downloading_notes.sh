cd /home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/preprocessing/0_fasta
mkdir -p assembly_report

# 1. HLacaPus1
cp /home/mpg08/mko/Nectar/Results2/acaPus/Genome/ver1/BRTHfinal.assembly.fasta HLacaPus1.fasta &

# 2. HLcalAnn5
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_assembly_report.txt
mv assembly_report/GCF_003957555.1_bCalAnn1_v1.p_assembly_report.txt assembly_report/HLcalAnn5_assembly_report.txt

# 3. HLcalPug1
cp /home/mpg08/mko/Nectar/Results2/calPug/Genome/GCF_001431845.1/GCF_001431845.1_ASM143184v1_genomic.fna HLcalPug1.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/431/845/GCF_001431845.1_ASM143184v1/GCF_001431845.1_ASM143184v1_assembly_report.txt
mv assembly_report/GCF_001431845.1_ASM143184v1_assembly_report.txt assembly_report/HLcalPug1_assembly_report.txt

# 4. HLchaPel1
cp /home/mpg08/mko/Nectar/Results2/chaPel/Genome/GCF_000747805.1/GCF_000747805.1_ChaPel_1.0_genomic.fna HLchaPel1.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/747/805/GCF_000747805.1_ChaPel_1.0/GCF_000747805.1_ChaPel_1.0_assembly_report.txt
mv assembly_report/GCF_000747805.1_ChaPel_1.0_assembly_report.txt assembly_report/HLchaPel1_assembly_report.txt

# 5. HLcolLiv2
cp /home/mpg08/mko/Nectar/Results2/colLiv/Genome/GCA_000337935.2/GCA_000337935.2_Cliv_2.1_genomic.fna HLcolLiv2.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/337/935/GCA_000337935.2_Cliv_2.1/GCA_000337935.2_Cliv_2.1_assembly_report.txt
mv assembly_report/GCA_000337935.2_Cliv_2.1_assembly_report.txt assembly_report/HLcolLiv2_assembly_report.txt

# 6. HLcorCor3
cp /home/mpg08/mko/Nectar/Results2/corCor/Genome/GCF_000738735.2/GCF_000738735.2_ASM73873v2_genomic.fna HLcorCor3.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.2_ASM73873v2/GCF_000738735.2_ASM73873v2_assembly_report.txt
mv assembly_report/GCF_000738735.2_ASM73873v2_assembly_report.txt assembly_report/HLcorCor3_assembly_report.txt

# 7. HLempTra1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/031/625/GCF_003031625.1_ASM303162v1/GCF_003031625.1_ASM303162v1_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/031/625/GCF_003031625.1_ASM303162v1/GCF_003031625.1_ASM303162v1_assembly_report.txt
mv assembly_report/GCF_003031625.1_ASM303162v1_assembly_report.txt assembly_report/HLempTra1_assembly_report.txt

# 8. HLfalTin1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/332/995/GCA_010332995.1_FalTin1.0/GCA_010332995.1_FalTin1.0_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/332/995/GCA_010332995.1_FalTin1.0/GCA_010332995.1_FalTin1.0_assembly_report.txt
mv assembly_report/GCA_010332995.1_FalTin1.0_assembly_report.txt assembly_report/HLfalTin1_assembly_report.txt

# 9. HLfloFus1
cp /home/mpg08/mko/Nectar/Results2/floFus/Genome/floFus3/floFus3.fa HLfloFus1.fasta &

# 10. HLfurRuf1
cp /home/mpg08/mko/Nectar/Results2/furRuf/Genome/10x/furRuf.fasta HLfurRuf1.fasta &

# 11. HLlicCas1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/360/975/GCA_008360975.1_HeHo_1.0/GCA_008360975.1_HeHo_1.0_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/360/975/GCA_008360975.1_HeHo_1.0/GCA_008360975.1_HeHo_1.0_assembly_report.txt
mv assembly_report/GCA_008360975.1_HeHo_1.0_assembly_report.txt assembly_report/HLlicCas1_assembly_report.txt

# 12. HLmalCya1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/741/485/GCA_009741485.1_mCya_1.0/GCA_009741485.1_mCya_1.0_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/741/485/GCA_009741485.1_mCya_1.0/GCA_009741485.1_mCya_1.0_assembly_report.txt
mv assembly_report/GCA_009741485.1_mCya_1.0_assembly_report.txt assembly_report/HLmalCya1_assembly_report.txt

# 13. HLparMaj1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_assembly_report.txt
mv assembly_report/GCF_001522545.3_Parus_major1.1_assembly_report.txt assembly_report/HLparMaj1_assembly_report.txt

# 14. HLphyNov1
cp /home/mpg08/mko/Nectar/Results2/phyNov/Genome/ver1/NHHOfinal.assembly.fasta HLphyNov1.fasta &

# 15. HLserCan1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/115/625/GCF_007115625.1_cibio_Scana_2019/GCF_007115625.1_cibio_Scana_2019_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/115/625/GCF_007115625.1_cibio_Scana_2019/GCF_007115625.1_cibio_Scana_2019_assembly_report.txt
mv assembly_report/GCF_007115625.1_cibio_Scana_2019_assembly_report.txt assembly_report/HLserCan1_assembly_report.txt

# 16. HLstrHab1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/027/225/GCF_004027225.2_bStrHab1.2.pri/GCF_004027225.2_bStrHab1.2.pri_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/027/225/GCF_004027225.2_bStrHab1.2.pri/GCF_004027225.2_bStrHab1.2.pri_assembly_report.txt
mv assembly_report/GCF_004027225.2_bStrHab1.2.pri_assembly_report.txt assembly_report/HLstrHab1_assembly_report.txt

# 17. HLtaeGut4: not GCF_008822105.2, Katya uses GCA_003957565.2
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna.gz &
# wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_assembly_report.txt
# mv assembly_report/GCF_008822105.2_bTaeGut2.pat.W.v2_assembly_report.txt assembly_report/HLtaeGut4_assembly_report.txt
# Katya uses GCA_003957565.2
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/565/GCA_003957565.2_bTaeGut1_v1.p/GCA_003957565.2_bTaeGut1_v1.p_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/565/GCA_003957565.2_bTaeGut1_v1.p/GCA_003957565.2_bTaeGut1_v1.p_assembly_report.txt
mv assembly_report/GCA_003957565.2_bTaeGut1_v1.p_assembly_report.txt assembly_report/HLtaeGut4_assembly_report.txt

# 18. HLtriMol2
cp /home/mpg08/mko/Nectar/Results2/triMol/Genome/HiC/triMol_3DDNA.fasta HLtriMol2.fasta &

# 19.HLtytAlb2 ? GCF_000687205.1 in shared sheet, but representative genome is GCF_902150015.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/150/015/GCF_902150015.1_genome_assembly_l500/GCF_902150015.1_genome_assembly_l500_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/150/015/GCF_902150015.1_genome_assembly_l500/GCF_902150015.1_genome_assembly_l500_assembly_report.txt
mv assembly_report/GCF_902150015.1_genome_assembly_l500_assembly_report.txt assembly_report/HLtytAlb2_assembly_report.txt

# 20. HLzosLat1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/735/GCA_001281735.1_ASM128173v1/GCA_001281735.1_ASM128173v1_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/735/GCA_001281735.1_ASM128173v1/GCA_001281735.1_ASM128173v1_assembly_report.txt
mv assembly_report/GCA_001281735.1_ASM128173v1_assembly_report.txt assembly_report/HLzosLat1_assembly_report.txt

# 21. aptFor1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/145/GCF_000699145.1_ASM69914v1/GCF_000699145.1_ASM69914v1_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/145/GCF_000699145.1_ASM69914v1/GCF_000699145.1_ASM69914v1_assembly_report.txt
mv assembly_report/GCF_000699145.1_ASM69914v1_assembly_report.txt assembly_report/aptFor1_assembly_report.txt

# 22. cucCan1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/709/325/GCF_000709325.1_ASM70932v1/GCF_000709325.1_ASM70932v1_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/709/325/GCF_000709325.1_ASM70932v1/GCF_000709325.1_ASM70932v1_assembly_report.txt
mv assembly_report/GCF_000709325.1_ASM70932v1_assembly_report.txt assembly_report/cucCan1_assembly_report.txt

# 23. falChe1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/975/GCF_000337975.1_F_cherrug_v1.0/GCF_000337975.1_F_cherrug_v1.0_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/975/GCF_000337975.1_F_cherrug_v1.0/GCF_000337975.1_F_cherrug_v1.0_assembly_report.txt
mv assembly_report/GCF_000337975.1_F_cherrug_v1.0_assembly_report.txt assembly_report/falChe1_assembly_report.txt

# 24. falPer1
cp /home/mpg08/mko/Nectar/Results2/falPer/Genome/GCF_000337955.1/GCF_000337955.1_F_peregrinus_v1.0_genomic.fna falPer1.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/955/GCF_000337955.1_F_peregrinus_v1.0/GCF_000337955.1_F_peregrinus_v1.0_assembly_report.txt
mv assembly_report/GCF_000337955.1_F_peregrinus_v1.0_assembly_report.txt assembly_report/falPer1_assembly_report.txt

# 25. ficAlb2
cp /home/mpg08/mko/Nectar/Results2/ficAlb/Genome/GCF_000247815.1/GCF_000247815.1_FicAlb1.5_genomic.fna ficAlb2.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/815/GCF_000247815.1_FicAlb1.5/GCF_000247815.1_FicAlb1.5_assembly_report.txt
mv assembly_report/GCF_000247815.1_FicAlb1.5_assembly_report.txt assembly_report/ficAlb2_assembly_report.txt

# 26. galGal6
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_assembly_report.txt
mv assembly_report/GCF_000002315.6_GRCg6a_assembly_report.txt assembly_report/galGal6_assembly_report.txt

# 27. halLeu1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/465/GCF_000737465.1_Haliaeetus_leucocephalus-4.0/GCF_000737465.1_Haliaeetus_leucocephalus-4.0_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/465/GCF_000737465.1_Haliaeetus_leucocephalus-4.0/GCF_000737465.1_Haliaeetus_leucocephalus-4.0_assembly_report.txt
mv assembly_report/GCF_000737465.1_Haliaeetus_leucocephalus-4.0_assembly_report.txt assembly_report/halLeu1_assembly_report.txt

# 28. melUnd1
cp /home/mpg08/mko/Nectar/Results2/melUnd/Genome/GCF_000238935.1/GCF_000238935.1_Melopsittacus_undulatus_6.3_genomic.fna melUnd1.fasta &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/238/935/GCA_000238935.1_Melopsittacus_undulatus_6.3/GCA_000238935.1_Melopsittacus_undulatus_6.3_assembly_report.txt
mv assembly_report/GCA_000238935.1_Melopsittacus_undulatus_6.3_assembly_report.txt assembly_report/melUnd1_assembly_report.txt

# 29. opiHoa1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/692/075/GCF_000692075.1_ASM69207v1/GCF_000692075.1_ASM69207v1_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/692/075/GCF_000692075.1_ASM69207v1/GCF_000692075.1_ASM69207v1_assembly_report.txt
mv assembly_report/GCF_000692075.1_ASM69207v1_assembly_report.txt assembly_report/opiHoa1_assembly_report.txt

# 30. pseHum1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/331/425/GCF_000331425.1_PseHum1.0/GCF_000331425.1_PseHum1.0_genomic.fna.gz &
wget -P assembly_report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/331/425/GCF_000331425.1_PseHum1.0/GCF_000331425.1_PseHum1.0_assembly_report.txt
mv assembly_report/GCF_000331425.1_PseHum1.0_assembly_report.txt assembly_report/pseHum1_assembly_report.txt


gunzip *.gz

# 2. HLcalAnn5
mv GCF_003957555.1_bCalAnn1_v1.p_genomic.fna HLcalAnn5.fasta &

# 7. HLempTra1
mv GCF_003031625.1_ASM303162v1_genomic.fna HLempTra1.fasta &

# 8. HLfalTin1[][][][][][]
mv GCA_010332995.1_FalTin1.0_genomic.fna HLfalTin1.fasta &

# 11. HLlicCas1
mv GCA_008360975.1_HeHo_1.0_genomic.fna HLlicCas1.fasta &

# 12. HLmalCya1
mv GCA_009741485.1_mCya_1.0_genomic.fna HLmalCya1.fasta &

# 13. HLparMaj1
mv GCF_001522545.3_Parus_major1.1_genomic.fna HLparMaj1.fasta &

# 15. HLserCan1
mv GCF_007115625.1_cibio_Scana_2019_genomic.fna HLserCan1.fasta &

# 16. HLstrHab1
mv GCF_004027225.2_bStrHab1.2.pri_genomic.fna HLstrHab1.fasta &

# 17. HLtaeGut4
# mv GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna HLtaeGut4.fasta &
mv GCA_003957565.2_bTaeGut1_v1.p_genomic.fna.gz HLtaeGut4.fasta &

# 19.HLtytAlb2 ? GCF_000687205.1 in shared sheet, but representative genome is GCF_902150015.1
mv GCF_902150015.1_genome_assembly_l500_genomic.fna HLtytAlb2.fasta &

# 20. HLzosLat1
mv GCA_001281735.1_ASM128173v1_genomic.fna  HLzosLat1.fasta &

# 21. aptFor1
mv GCF_000699145.1_ASM69914v1_genomic.fna aptFor1.fasta &

# 22. cucCan1
mv GCF_000709325.1_ASM70932v1_genomic.fna cucCan1.fasta &

# 23. falChe1
mv GCF_000337975.1_F_cherrug_v1.0_genomic.fna falChe1.fasta &

# 26. galGal6
mv GCF_000002315.6_GRCg6a_genomic.fna galGal6.fasta &

# 27. halLeu1
mv GCF_000737465.1_Haliaeetus_leucocephalus-4.0_genomic.fna halLeu1.fasta &

# 29. opiHoa1
mv GCF_000692075.1_ASM69207v1_genomic.fna opiHoa1.fasta &

# 30. pseHum1
mv GCF_000331425.1_PseHum1.0_genomic.fna pseHum1.fasta &

mkdir ori
mv *.fasta ori
