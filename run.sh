#!/bin/bash
#rm run-m_ref_g_mmseq_a_diamond_c_mcl-test-main.txt
# #/usr/bin/time -v  panta main -t 12 -m g_mmseq_a_diamond_c_mcl_rand  --dont-split  -o out/test/g_mmseq_a_diamond_c_mcl_rand_all -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/all/*.gff >> run-g_mmseq_a_diamond_c_mcl_rand-test-all.txt 2>&1
#rm run-g_mmseq_a_diamond_c_mcl-test-main-apoa-cons.txt
pip install .
#rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-main3-123.txt 2>&1
#rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add3-123_4.txt
#rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-main3-1234.txt
#rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-main3-4.txt
#rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add3-4_123.txt 
#rm -r out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_123_4
#rm -r out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_1234
#rm -r  out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_4_123
#/usr/bin/time -v  panta main -m m_ref_g_mmseq_a_diamond_c_mcl -a protein --poa --dont-split -o  out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_123_4 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_gff_123/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-main3-123.txt 2>&1
#/usr/bin/time -v  panta add -a protein --poa --dont-split -c  out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_123_4 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_4/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add3-123_4.txt 2>&1
#/usr/bin/time -v  panta main -m m_ref_g_mmseq_a_diamond_c_mcl -a protein --poa --dont-split -o out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_1234 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_gff_1234/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-main3-1234.txt 2>&1
#/usr/bin/time -v  panta main -m m_ref_g_mmseq_a_diamond_c_mcl -a protein --poa --dont-split -o out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_4_123 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_4/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-main3-4.txt 2>&1
#/usr/bin/time -v  panta add -a protein --poa --dont-split -c out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_nearest_group_4_123 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_gff_123/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add3-4_123.txt 2>&1
#/usr/bin/time -v ./diamond blastp -q out/Sp600/g_mmseq_a_diamond_c_mcl_apoa_cons_123_4/temp/combined.faa -d out/Sp600/g_mmseq_a_diamond_c_mcl_apoa_cons_123_4/consenus_diamond_db.dmnd -p 16 --evalue 1e-06 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> out/Sp600/g_mmseq_a_diamond_c_mcl_apoa_cons_123_4/diamond.tsv >> run-diamond-sp600-4.txt 2>&1
#for i in {70..99}
#do
    #let d = $i/100.0
#    /usr/bin/time -v panta build -i uniprotkb_reviewed_true_AND_taxonomy_id_2024_07_24.fasta -t 8 -d 0.$i -c 0.$i >> run-build-ref.txt 2>&1
#done
#rm -r hmmer_gene_families
#rm run-build-ref-hmmer.txt
#/usr/bin/time -v panta build -i uniprotkb_reviewed_true_AND_taxonomy_id_2024_07_24.fasta -t 8 -d 0.98 -c 0.98 >> run-build-ref-hmmer.txt 2>&1
# pip install .
# rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-123.txt 2>&1
# rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add-123_4.txt
# rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-1234.txt
# rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-4.txt
# rm run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add-4_123.txt 
# rm -r out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_123_4
# rm -r out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_1234
# rm -r  out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_4_123
# /usr/bin/time -v  panta main -m m_ref_g_mmseq_a_diamond_c_mcl -a protein --poa --dont-split -o  out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_123_4 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_gff_123/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-123.txt 2>&1
# /usr/bin/time -v  panta add -a protein --poa --dont-split -c  out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_123_4 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_4/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add-123_4.txt 2>&1
# /usr/bin/time -v  panta main -m m_ref_g_mmseq_a_diamond_c_mcl -a protein --poa --dont-split -o out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_1234 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_gff_1234/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-1234.txt 2>&1
# /usr/bin/time -v  panta main -m m_ref_g_mmseq_a_diamond_c_mcl -a protein --poa --dont-split -o out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_4_123 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_4/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-4.txt 2>&1
# /usr/bin/time -v  panta add -a protein --poa --dont-split -c out/Sp600/m_ref_g_mmseq_a_diamond_c_mcl_4_123 -g /media/ktht/Store/Quang/bio/CollectionPanta/TestInQualifiedGFF/Sp600_batches/Sp600/batch_gff_123/*.gff >> run-m_ref_g_mmseq_a_diamond_c_mcl-Sp600-add-4_123.txt 2>&1
#/usr/bin/time -v  panta main  --dont-split -o  out/Kp1500/g_cdhit_a_diamond_c_mcl_Kp1500_12345 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_1/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_2/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_3/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_4/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_5/*.gff>> run-g_cdhit_a_diamond_c_mcl-Kp1500-12345.txt 2>&1
#pip install .
# rm run-g_mmseq_a_diamond_c_mcl_unique3.1-Kp1500-1.txt 2>&1
# rm run-g_mmseq_a_diamond_c_mcl_unique3.1-Kp1500-add-12.txt
# rm run-g_mmseq_a_diamond_c_mcl_unique3.1-Kp1500-add-123.txt
# # # rm run-g_mmseq_a_diamond_c_mcl_unique-Kp1500-add-1234.txt
# # # rm run-g_mmseq_a_diamond_c_mcl_unique2-Kp1500-add-12345.txt
#rm run-g_mmseq_a_diamond_c_mcl_unique_cc-Kp1500-12345.txt

#rm -r out/Kp1500/g_mmseq_a_diamond_c_mcl_unique_cc_123_45
#rm -r out/Kp1500/g_mmseq_a_diamond_c_mcl_unique_cc_1234_5
# rm -r out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5

#/usr/bin/time -v  panta main -m g_faiss_a_faiss_c_mcl   --dont-split -o  out/Kp1500/g_faiss_a_faiss_c_mcl -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_1/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_2/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_3/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique_cc-Kp1500-123.txt 2>&1
# /usr/bin/time -v  panta main -m g_mmseq_a_diamond_c_mcl_cc   --dont-split -o  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_1/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique3.1-Kp1500-1.txt 2>&1

#/usr/bin/time -v  panta add --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_2/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique3.1-Kp1500-add-12.txt 2>&1
# # # #/usr/bin/time -v  panta add --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_2/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique3-Kp1500-add-12.txt 2>&1
# /usr/bin/time -v  panta add --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5 -g  /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_3/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique3.1-Kp1500-add-123.txt 2>&1
#/usr/bin/time -v  panta add  --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique_cc_123_45 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_4/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_5/*.gff>> run-g_mmseq_a_diamond_c_mcl_unique_cc-Kp1500-add-123_45.txt 2>&1
#/usr/bin/time -v  panta add  --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_4/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique3-Kp1500-add-1234.txt 2>&1

#/usr/bin/time -v  panta add  --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique3_1_2_3_4_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_5/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique3-Kp1500-add-12345.txt 2>&1
#/usr/bin/time -v  panta main -m g_faiss_a_faiss_c_mcl  --dont-split -o out/Kp1500/g_faiss_a_faiss_c_mcl_12345 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_1/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_2/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_3/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_4/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_5/*.gff >> run-g_faiss_a_faiss_c_mcl-Kp1500-12345.txt 2>&1
#/usr/bin/time -v  panta main -m g_mmseq_a_diamond_c_mcl_cc --dont-split -o out/Kp1500/g_mmseq_a_diamond_c_mcl_unique_cc_1234_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_1/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_2/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_3/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_4/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique_cc-Kp1500-1234.txt 2>&1
#/usr/bin/time -v  panta add  --dont-split -c  out/Kp1500/g_mmseq_a_diamond_c_mcl_unique_cc_1234_5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Kp1500/batch_5/*.gff>> run-g_mmseq_a_diamond_c_mcl_unique_cc-Kp1500-add-1234_5.txt 2>&1
#panaroo --merge_paralogs -t 20 --clean-mode strict
#rm run-g_mmseq_a_ems_c_mcl_unique-test-main.txt
rm run-g_mmseq_a_ems_c_mcl_unique-test-all.txt
#rm run-g_mmseq_a_diamond_c_mcl_unique5-test-add-1.txt
#rm run-g_mmseq_a_diamond_c_mcl_unique5-test-add-2.txt
#rm run-g_mmseq_a_diamond_c_mcl_unique5-test-add-3.txt
rm -r out/test/g_mmseq_a_ems_c_mcl_unique
#rm -r out/test/g_mmseq_a_diamond_c_mcl_unique5_all
#rm -r out/test/g_mmseq_a_diamond_c_mcl_unique5_12
#/usr/bin/time -v  panta main -t 12 -m g_mmseq_a_ems_c_mcl  --dont-split  -o out/test/g_mmseq_a_ems_c_mcl_unique -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/main/*.gff >> run-g_mmseq_a_ems_c_mcl_unique-test-main.txt 2>&1
/usr/bin/time -v  panta main -t 12 -m g_mmseq_a_ems_c_mcl  --dont-split  -o out/test/g_mmseq_a_ems_c_mcl_unique -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/main/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add2/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add3/*.gff >> run-g_mmseq_a_ems_c_mcl_unique-test-all.txt 2>&1

#/usr/bin/time -v  panta add --dont-split -c  out/test/g_mmseq_a_diamond_c_mcl_unique5 -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique5-test-add-1.txt 2>&1
#/usr/bin/time -v  panta add --dont-split -c  out/test/g_mmseq_a_diamond_c_mcl_unique5 -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add2/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique5-test-add-2.txt 2>&1
#/usr/bin/time -v  panta add --dont-split -c  out/test/g_mmseq_a_diamond_c_mcl_unique5 -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add3/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique5-test-add-3.txt 2>&1
#/usr/bin/time -v  panta main -t 12 -m g_mmseq_a_diamond_c_mcl_cc  --dont-split  -o out/test/g_mmseq_a_diamond_c_mcl_unique5_all -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/main/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add2/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add3/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique5-test-all.txt 2>&1
#/usr/bin/time -v  panta main -t 12 -m g_mmseq_a_diamond_c_mcl_cc  --dont-split  -o out/test/g_mmseq_a_diamond_c_mcl_unique5_12 -g /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/main/*.gff /media/ktht/ExtraUbuntu/amromics/panta2/panta/examples/test/add/*.gff >> run-g_mmseq_a_diamond_c_mcl_unique5-test-12.txt 2>&1

# rm -r out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5
# /usr/bin/time -v panta main -t 12 -m g_mmseq_a_diamond_c_mcl_cc  --dont-split -o out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Ecoli/2020-Q1/*gff.gz >> run-unique5_ecoli_2020_q1.txt 2>&1
# for i in 2020 2021 2022
# do
#     for j in 1 2 3 4
#     do
#         if [ $i == 2020 ]; then
#             if [ $j == 1 ]; then
#                 continue
#             fi
#         fi
#         rm run-unique5_ecoli_${i}_q${j}.txt
#         /usr/bin/time -v panta add  -t 12 --dont-split -c out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Ecoli/${i}-Q${j}/*gff.gz >> run-unique5_ecoli_${i}_q${j}.txt 2>&1
#     done
# done
# /usr/bin/time -v  mmseqs easy-cluster out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/temp/combined.faa out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/temp/mmseq out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/temp/tmp --min-seq-id 1 -c 1 --threads 12 >> run-mmseq_easyscluster_combinedfaa.txt 2>&1
# /usr/bin/time -v  mmseqs easy-cluster out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/temp/combined_unique_seq.faa out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/temp/ident/mmseq out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/temp/ident/tmp --min-seq-id 1 -c 1 --threads 12 >> run-mmseq_easyscluster_combinedfaa.txt 2>&1
#/usr/bin/time -v ./diamond blastp -q out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/similar_seqs.fasta -d out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/diamond_db -p 12 --evalue 0.000001 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs 2000 > out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique3/diamond.tsv >> run-diamondblast_combinedfaa.txt 2>&1
#/usr/bin/time -v panta add -t 12 --dont-split -c out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5 -g /media/ktht/ExtraUbuntu/amromics/panta2/Ecoli/S13/*gff.gz >> run-unique5.1_ecoli_s13.txt 2>&1
#/usr/bin/time -v ./diamond blastp -q out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5/temp/reduced_new_unique_seqs.fasta -d out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5/representative_clusters.fasta -p 12 --evalue 1e-06 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 > out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5/temp/matching_new_sequences_diamond.tsv  >> run-diamondblast_reduced_new_unique_seqs.txt 2>&1
#/usr/bin/time -v ./diamond blastp -q out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5/representative_clusters.fasta -d out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5/representative_clusters.fasta -p 4 --evalue 1e-6 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> out/Ecoli30k/g_mmseq_a_diamond_c_mcl_unique5/pairwise_rep_clusters_diamond_test.tsv >> run-diamondblast_representative.txt 2>&1