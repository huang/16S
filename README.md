# 16S



## 1. Run FastQC to allow manual inspection of the quality of sequences
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/

## 2. Rename the files
cd raw_data
for file in *.fastq.gz; do mv $file $(echo $file | cut -d'_' -f1 | cut -d'-' -f1-1)_$(echo $file | cut -d'_' -f4).fastq.gz; done
cd ..


paste -d" " stool_f1.txt stool_f5.txt > stool_f1_f5.sh


./*.sh






## 3.1. trim paired-end reads
mkdir trim_data trimmed_unpaired
cd raw_data




for file in 3-9141vag_R1.fastq.gz 16-9148stuhl_R1.fastq.gz 4-9140vag_R1.fastq.gz 15-9136stuhl_R1.fastq.gz 2-9133vag_R1.fastq.gz 11-9133stuhl_R1.fastq.gz 7-9148vag_R1.fastq.gz 10-9135stuhl_R1.fastq.gz 12-9141stuhl_R1.fastq.gz 5-9161vag_R1.fastq.gz ; do java -jar /home/jhuang/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 16 $file ${file/_R1/_R2} ../trim_data/$file ../trimmed_unpaired/$file ../trim_data/${file/_R1/_R2} ../trimmed_unpaired/${file/_R1/_R2} ILLUMINACLIP:/home/jhuang/Tools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 AVGQUAL:20; done 2> trimmomatic_pe.log
for file in 19-neg_R1.fastq.gz 1-9135vag_R1.fastq.gz 6-9136vag_R1.fastq.gz 9-9143vag_R1.fastq.gz 8-9147vag_R1.fastq.gz 18-9143stuhl_R1.fastq.gz 14-9161stuhl_R1.fastq.gz 17-9147stuhl_R1.fastq.gz 13-9140stuhl_R1.fastq.gz; do java -jar /home/jhuang/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 16 $file ${file/_R1/_R2} ../trim_data/$file ../trimmed_unpaired/$file ../trim_data/${file/_R1/_R2} ../trimmed_unpaired/${file/_R1/_R2} ILLUMINACLIP:/home/jhuang/Tools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 AVGQUAL:20; done 2> trimmomatic_pe2.log

#for file in p8_stool_d00_R1.fastq.gz p8_stool_d00r_R1.fastq.gz  p8_stool_d20_R1.fastq.gz p8_stool_d20r_R1.fastq.gz; do java -jar /home/jhuang/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 16 $file ${file/_R1/_R2} ../trim_data/$file ../trimmed_unpaired/$file ../trim_data/${file/_R1/_R2} ../trimmed_unpaired/${file/_R1/_R2} ILLUMINACLIP:/home/jhuang/Tools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 AVGQUAL:20; done 2> trimmomatic_pe.log


## NOTE that step 4 (pandaseq) is better than 3.2, since it removes the primers instead of matching the primers. --> spring into step 4


#3.2 
mkdir trim_contigs;
for file in trim_data/*_1P.fastq.gz; do mkdir trim_contigs/$(echo $file | cut -d '/' -f2 | cut -d'_' -f1); done
for file in trim_data/*.fastq.gz; do cp $file trim_contigs/$(echo $file | cut -d '/' -f2 | cut -d'_' -f1); done
cd trim_contigs;
for dir in *; do \
  mothur "#make.file(inputdir=./$dir, type=gz, prefix=stability); make.contigs(file=./$dir/stability.files, processors=15);"; \
done
#cat all log-files
#rm mothur.*.logfile
for dir in *; do \
  mothur "#summary.seqs(fasta=./$dir/stability.trim.contigs.fasta, processors=15);"; \
done
#cat all log-files2
#rm mothur.*.logfile
for dir in *; do \
  mothur "#screen.seqs(fasta=./$dir/stability.trim.contigs.fasta, group=./$dir/stability.contigs.groups, maxambig=0, minlength=380, maxlength=460, maxhomop=6);"; \
done

#WARNING: tolerant version
#for dir in 56-488ST 57-514ST; do \
#mothur "#screen.seqs(fasta=./$dir/stability.trim.contigs.fasta, group=./$dir/stability.contigs.groups, maxambig=6, minlength=200, maxlength=502, maxhomop=10);"; \
#done


#1246
mkdir filtered_reads_with_primer1
for dir in trim_contigs/*/; do \
  echo "bbduk.sh -Xmx1g in=$dir/stability.trim.contigs.good.fasta outm=filtered_reads_with_primer1/$(echo $dir | cut -d'/' -f2).assembled_filtered.fasta k=17 literal=CCTACGGGNGGCWGCAG mm=f rcomp=t copyundefined 2>>filter1_log.txt;" \
done

mkdir filtered_reads_with_primer2
for dir in filtered_reads_with_primer1/*.12.assembled_filtered.fasta; do \
  bbduk.sh -Xmx1g in=$dir outm=filtered_reads_with_primer2/$(echo $dir | cut -d '/' -f2 | cut -d'.' -f1-2).assembled_filtered.fasta k=21 literal=GACTACHVGGGTATCTAATCC mm=f rcomp=t copyundefined 2>>filter2_log.txt; \
done



1-2
literal=GYTCAGRDYKAACGCTGGCGG

4-6

#Here are the primers used for the 16S V3- V4 sequencing. Underlined are the ILLUMINA adapter sequences.
3-4
-p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CCTACGGGNGGCWGCAG
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GACTACHVGGGTATCTAATCC

CCTACGGGNGGCWGCAG
CCTACGGGAGGCAGCAGTGGGGAATATTGCGCAATGGGGGGAACCCTGACGGAGCCATGCCGCGTGAATGAAGAAGGCCTTCGGGTTGTAAAGTTCTTTGGGTGATGAGGAAGGGGTTAAGATTAATACCCTTTTTACTTACGGTTAGCTAAACCGGAACCACCGGCTACCTCGTCCCAGCCGGCCCGGGAATCATGGGGGTGGTATGCGTTATGCGGAANACTNGGGGNAANNGGGCCCNAGGNNGTTTTCTGAAGAGGGTGTAAGAGCAGCGCTCGACGCGCGGCGTGGATCATCGATCTGAGAGCTTCAGGCCATATGGAGGGAGTGAAATTCAAGGGTGACCTGTAAAGGCGAAATGACATTCACGAGTCCTCGAGGCGTATGCCACCTGTTGGGCATTACCTGTCACTTATGAGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCGTGTAGTC
                                  GACTACHVGGGTATCTAATCC

CCTACGGGNGGCWGCAG
CCTACGGGGGGCTGCCAAGGCGGCCTGGTTCTGTCAGATTGTCACCGCTTGCCGTGAGACTTGTGACCCTTACCCCGTCACTGGTAGTCAGGGTTAGCCGCGCGGAGCCCAAGGGATGCTCCTTGGGGCAAACTAATTCGGCCCAGCTCTTGTAGGCCTTCACCTGGGGTTCCCACTTTGTTCTCTGGGCCGCCTGGAACCACCCCGTGGGCCCCGGGTCTCTTTTTCCCCCCTCTGCGGATTACCCGGTTCCAAGCACGCCATACGATGGCTCCCAGTCCCTCGCCGCGCCCGGGCTACGATGTCGAGGACCTTCGATGAGCTGCTAACCGGAGTGCAAGTCACGCTCGAGGAGTGGAAATTTTATACTGTTTAGATCGCAGAAAAGTCCATCTACCTCACGAGCCAGCACTCGTCATAAGCTGATAATAGAGGGATCAGCCCTTGAACAACCAGAAGAGCCTTAGAGATGGAAAGCTAGATACCCCGGTAGTC
                                                                    GACTACHVGGGTATCTAATCC
                                                                    
     180088      --> 90044
#-p CCTAC GGGNG GCWGC AG -q GACTA CHVGG GTATC TAATC C 
## 4. stitch the reads together; and filter the reads by quality, length, and primer (see len-dis-gg_v3v4.pdf, cutoff length could be 350, 360 or 380)
mkdir pandaseq.out
#-p primer       Forward primer sequence or number of bases to be removed.
#-q primer       Reverse primer sequence or number of bases to be removed.

#@M03701:126:000000000-B6DW9:1:1101:18034:1879 1:N:0:9
#CCTACGGGGGGCTGCAG #TAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGGTTTCGGCTCGTAAAGCTCTGTTGGTAGTGAAGAAAGATAGAGGTAGTAACTGGCCTTTATTTGACGGTAATTACCTAGAAAGTCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGGT
#@M03701:126:000000000-B6DW9:1:1101:18034:1879 2:N:0:9
#GACTACCCGGGTATCTAATCC #TGTTCGCTACCCATGCTTTCGAGCCTCAGCGTCAGTTGCAGACCAGAGAGCCGCCTTCGCCACTGGTGTTCTTCCATATATCTACGCATTCCACCGCTACACATGGAGTTCCACTCTCCTCTTCTGCACTCAAGTTCAACAGTTTCTGATGCAATTCTCCGGTTGAGCCGAAGGCTTTCACATCAGACTTATTGAACCGCATGCACTCGCTTTACGCCCAATAAATCCGG
#-->
#>M03701:126:000000000-B6DW9:1:1101:18034:1879:9
#TAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGGTTTCGGCTCGTAAAGCTCTGTTGGTAGTGAAGAAAGATAGAGGTAGTAACTGGCCTTTATTTGACGGTAATTACCTAGAAAGTCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGTTCAATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTTGAACTTGAGTGCAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGCAACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACA
mkdir pandaseq.out
#for file in trim_data/*_R1.fastq.gz; do pandaseq -f ${file} -r ${file/_R1.fastq.gz/_R2.fastq.gz} -l 300 -p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC  -w pandaseq.out/$(echo $file | cut -d'/' -f2 | cut -d'_' -f1)_merged.fasta >> LOG_pandaseq; done

for file in trim_data/p*_R1.fastq.gz; do pandaseq -f ${file} -r ${file/_R1.fastq.gz/_R2.fastq.gz} -l 300 -p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC  -w pandaseq.out/$(echo $file | cut -d'/' -f2 | cut -d'_' -f1-3)_merged.fasta >> LOG_pandaseq; done

## 5. Create two QIIME mapping files
validate_mapping_file.py -m map.txt


#for sample in SRS048083 SRS020119 SRS054488 SRS022480 SRS013543 SRS056505 SRS023488 SRS021541 SRS022414 SRS057122 SRS012191 SRS042290 SRS022987 SRS013762 SRS020176 SRS023047 SRS023422 SRS021910 SRS024511 SRS022924 SRS056656 SRS023851; do
#seqtk seq -A ${sample}.fastq | seqkit seq -w 60 - > ${sample}.fasta
#done
#cp processed_controls/*.fasta pandaseq.out/
## 6. Combine files into a labeled file
add_qiime_labels.py -i pandaseq.out -m map_corrected.txt -c FileInput -o combined_fasta
add_qiime_labels.py -i pandaseq.out -m map_corrected_swab.txt -c FileInput -o combined_fasta_swab
add_qiime_labels.py -i pandaseq.out -m map_corrected_stool.txt -c FileInput -o combined_fasta_stool


## 7. Remove chimeric sequences using usearch
cd combined_fasta
pyfasta split -n 100 combined_seqs.fna
for i in {000..099}; do echo "identify_chimeric_seqs.py -i combined_fasta/combined_seqs.fna.${i} -m usearch61 -o usearch_checked_combined.${i}/ -r ~/REFs/gg_97_otus_4feb2011_fw_rc.fasta --threads=14;" >> uchime_commands.sh; done
mv uchime_commands.sh ..
./uchime_commands.sh
cat usearch_checked_combined.000/chimeras.txt usearch_checked_combined.001/chimeras.txt usearch_checked_combined.002/chimeras.txt usearch_checked_combined.003/chimeras.txt usearch_checked_combined.004/chimeras.txt usearch_checked_combined.005/chimeras.txt usearch_checked_combined.006/chimeras.txt usearch_checked_combined.007/chimeras.txt usearch_checked_combined.008/chimeras.txt usearch_checked_combined.009/chimeras.txt usearch_checked_combined.010/chimeras.txt usearch_checked_combined.011/chimeras.txt usearch_checked_combined.012/chimeras.txt usearch_checked_combined.013/chimeras.txt usearch_checked_combined.014/chimeras.txt usearch_checked_combined.015/chimeras.txt usearch_checked_combined.016/chimeras.txt usearch_checked_combined.017/chimeras.txt usearch_checked_combined.018/chimeras.txt usearch_checked_combined.019/chimeras.txt usearch_checked_combined.020/chimeras.txt usearch_checked_combined.021/chimeras.txt usearch_checked_combined.022/chimeras.txt usearch_checked_combined.023/chimeras.txt usearch_checked_combined.024/chimeras.txt usearch_checked_combined.025/chimeras.txt usearch_checked_combined.026/chimeras.txt usearch_checked_combined.027/chimeras.txt usearch_checked_combined.028/chimeras.txt usearch_checked_combined.029/chimeras.txt usearch_checked_combined.030/chimeras.txt usearch_checked_combined.031/chimeras.txt usearch_checked_combined.032/chimeras.txt usearch_checked_combined.033/chimeras.txt usearch_checked_combined.034/chimeras.txt usearch_checked_combined.035/chimeras.txt usearch_checked_combined.036/chimeras.txt usearch_checked_combined.037/chimeras.txt usearch_checked_combined.038/chimeras.txt usearch_checked_combined.039/chimeras.txt usearch_checked_combined.040/chimeras.txt usearch_checked_combined.041/chimeras.txt usearch_checked_combined.042/chimeras.txt usearch_checked_combined.043/chimeras.txt usearch_checked_combined.044/chimeras.txt usearch_checked_combined.045/chimeras.txt usearch_checked_combined.046/chimeras.txt usearch_checked_combined.047/chimeras.txt usearch_checked_combined.048/chimeras.txt usearch_checked_combined.049/chimeras.txt usearch_checked_combined.050/chimeras.txt usearch_checked_combined.051/chimeras.txt usearch_checked_combined.052/chimeras.txt usearch_checked_combined.053/chimeras.txt usearch_checked_combined.054/chimeras.txt usearch_checked_combined.055/chimeras.txt usearch_checked_combined.056/chimeras.txt usearch_checked_combined.057/chimeras.txt usearch_checked_combined.058/chimeras.txt usearch_checked_combined.059/chimeras.txt usearch_checked_combined.060/chimeras.txt usearch_checked_combined.061/chimeras.txt usearch_checked_combined.062/chimeras.txt usearch_checked_combined.063/chimeras.txt usearch_checked_combined.064/chimeras.txt usearch_checked_combined.065/chimeras.txt usearch_checked_combined.066/chimeras.txt usearch_checked_combined.067/chimeras.txt usearch_checked_combined.068/chimeras.txt usearch_checked_combined.069/chimeras.txt usearch_checked_combined.070/chimeras.txt usearch_checked_combined.071/chimeras.txt usearch_checked_combined.072/chimeras.txt usearch_checked_combined.073/chimeras.txt usearch_checked_combined.074/chimeras.txt usearch_checked_combined.075/chimeras.txt usearch_checked_combined.076/chimeras.txt usearch_checked_combined.077/chimeras.txt usearch_checked_combined.078/chimeras.txt usearch_checked_combined.079/chimeras.txt usearch_checked_combined.080/chimeras.txt usearch_checked_combined.081/chimeras.txt usearch_checked_combined.082/chimeras.txt usearch_checked_combined.083/chimeras.txt usearch_checked_combined.084/chimeras.txt usearch_checked_combined.085/chimeras.txt usearch_checked_combined.086/chimeras.txt usearch_checked_combined.087/chimeras.txt usearch_checked_combined.088/chimeras.txt usearch_checked_combined.089/chimeras.txt usearch_checked_combined.090/chimeras.txt usearch_checked_combined.091/chimeras.txt usearch_checked_combined.092/chimeras.txt usearch_checked_combined.093/chimeras.txt usearch_checked_combined.094/chimeras.txt usearch_checked_combined.095/chimeras.txt usearch_checked_combined.096/chimeras.txt usearch_checked_combined.097/chimeras.txt usearch_checked_combined.098/chimeras.txt usearch_checked_combined.099/chimeras.txt > chimeras.txt
filter_fasta.py -f combined_fasta/combined_seqs.fna -o combined_fasta/combined_nonchimera_seqs.fna -s chimeras.txt -n;
rm -rf usearch_checked_combined.0*
#mkdir usearch_checked_combined_DEL
#mv usearch_checked_combined.0* usearch_checked_combined_DEL




cd combined_fasta_swab
pyfasta split -n 100 combined_seqs.fna
for i in {000..099}; do echo "identify_chimeric_seqs.py -i combined_fasta_swab/combined_seqs.fna.${i} -m usearch61 -o usearch_checked_combined_swab.${i}/ -r ~/REFs/gg_97_otus_4feb2011_fw_rc.fasta --threads=14;" >> uchime_commands_swab.sh; done
mv uchime_commands_swab.sh ..
./uchime_commands_swab.sh
cat usearch_checked_combined_swab.000/chimeras.txt usearch_checked_combined_swab.001/chimeras.txt usearch_checked_combined_swab.002/chimeras.txt usearch_checked_combined_swab.003/chimeras.txt usearch_checked_combined_swab.004/chimeras.txt usearch_checked_combined_swab.005/chimeras.txt usearch_checked_combined_swab.006/chimeras.txt usearch_checked_combined_swab.007/chimeras.txt usearch_checked_combined_swab.008/chimeras.txt usearch_checked_combined_swab.009/chimeras.txt usearch_checked_combined_swab.010/chimeras.txt usearch_checked_combined_swab.011/chimeras.txt usearch_checked_combined_swab.012/chimeras.txt usearch_checked_combined_swab.013/chimeras.txt usearch_checked_combined_swab.014/chimeras.txt usearch_checked_combined_swab.015/chimeras.txt usearch_checked_combined_swab.016/chimeras.txt usearch_checked_combined_swab.017/chimeras.txt usearch_checked_combined_swab.018/chimeras.txt usearch_checked_combined_swab.019/chimeras.txt usearch_checked_combined_swab.020/chimeras.txt usearch_checked_combined_swab.021/chimeras.txt usearch_checked_combined_swab.022/chimeras.txt usearch_checked_combined_swab.023/chimeras.txt usearch_checked_combined_swab.024/chimeras.txt usearch_checked_combined_swab.025/chimeras.txt usearch_checked_combined_swab.026/chimeras.txt usearch_checked_combined_swab.027/chimeras.txt usearch_checked_combined_swab.028/chimeras.txt usearch_checked_combined_swab.029/chimeras.txt usearch_checked_combined_swab.030/chimeras.txt usearch_checked_combined_swab.031/chimeras.txt usearch_checked_combined_swab.032/chimeras.txt usearch_checked_combined_swab.033/chimeras.txt usearch_checked_combined_swab.034/chimeras.txt usearch_checked_combined_swab.035/chimeras.txt usearch_checked_combined_swab.036/chimeras.txt usearch_checked_combined_swab.037/chimeras.txt usearch_checked_combined_swab.038/chimeras.txt usearch_checked_combined_swab.039/chimeras.txt usearch_checked_combined_swab.040/chimeras.txt usearch_checked_combined_swab.041/chimeras.txt usearch_checked_combined_swab.042/chimeras.txt usearch_checked_combined_swab.043/chimeras.txt usearch_checked_combined_swab.044/chimeras.txt usearch_checked_combined_swab.045/chimeras.txt usearch_checked_combined_swab.046/chimeras.txt usearch_checked_combined_swab.047/chimeras.txt usearch_checked_combined_swab.048/chimeras.txt usearch_checked_combined_swab.049/chimeras.txt usearch_checked_combined_swab.050/chimeras.txt usearch_checked_combined_swab.051/chimeras.txt usearch_checked_combined_swab.052/chimeras.txt usearch_checked_combined_swab.053/chimeras.txt usearch_checked_combined_swab.054/chimeras.txt usearch_checked_combined_swab.055/chimeras.txt usearch_checked_combined_swab.056/chimeras.txt usearch_checked_combined_swab.057/chimeras.txt usearch_checked_combined_swab.058/chimeras.txt usearch_checked_combined_swab.059/chimeras.txt usearch_checked_combined_swab.060/chimeras.txt usearch_checked_combined_swab.061/chimeras.txt usearch_checked_combined_swab.062/chimeras.txt usearch_checked_combined_swab.063/chimeras.txt usearch_checked_combined_swab.064/chimeras.txt usearch_checked_combined_swab.065/chimeras.txt usearch_checked_combined_swab.066/chimeras.txt usearch_checked_combined_swab.067/chimeras.txt usearch_checked_combined_swab.068/chimeras.txt usearch_checked_combined_swab.069/chimeras.txt usearch_checked_combined_swab.070/chimeras.txt usearch_checked_combined_swab.071/chimeras.txt usearch_checked_combined_swab.072/chimeras.txt usearch_checked_combined_swab.073/chimeras.txt usearch_checked_combined_swab.074/chimeras.txt usearch_checked_combined_swab.075/chimeras.txt usearch_checked_combined_swab.076/chimeras.txt usearch_checked_combined_swab.077/chimeras.txt usearch_checked_combined_swab.078/chimeras.txt usearch_checked_combined_swab.079/chimeras.txt usearch_checked_combined_swab.080/chimeras.txt usearch_checked_combined_swab.081/chimeras.txt usearch_checked_combined_swab.082/chimeras.txt usearch_checked_combined_swab.083/chimeras.txt usearch_checked_combined_swab.084/chimeras.txt usearch_checked_combined_swab.085/chimeras.txt usearch_checked_combined_swab.086/chimeras.txt usearch_checked_combined_swab.087/chimeras.txt usearch_checked_combined_swab.088/chimeras.txt usearch_checked_combined_swab.089/chimeras.txt usearch_checked_combined_swab.090/chimeras.txt usearch_checked_combined_swab.091/chimeras.txt usearch_checked_combined_swab.092/chimeras.txt usearch_checked_combined_swab.093/chimeras.txt usearch_checked_combined_swab.094/chimeras.txt usearch_checked_combined_swab.095/chimeras.txt usearch_checked_combined_swab.096/chimeras.txt usearch_checked_combined_swab.097/chimeras.txt usearch_checked_combined_swab.098/chimeras.txt usearch_checked_combined_swab.099/chimeras.txt > chimeras_swab.txt
filter_fasta.py -f combined_fasta_swab/combined_seqs.fna -o combined_fasta_swab/combined_nonchimera_seqs.fna -s chimeras_swab.txt -n;
rm -rf usearch_checked_combined_*.0*








cd combined_fasta_stool
pyfasta split -n 100 combined_seqs.fna
for i in {000..099}; do echo "identify_chimeric_seqs.py -i combined_fasta_stool/combined_seqs.fna.${i} -m usearch61 -o usearch_checked_combined_stool.${i}/ -r ~/REFs/gg_97_otus_4feb2011_fw_rc.fasta --threads=14;" >> uchime_commands_stool.sh; done
mv uchime_commands_stool.sh ..
./uchime_commands_stool.sh
cat usearch_checked_combined_stool.000/chimeras.txt usearch_checked_combined_stool.001/chimeras.txt usearch_checked_combined_stool.002/chimeras.txt usearch_checked_combined_stool.003/chimeras.txt usearch_checked_combined_stool.004/chimeras.txt usearch_checked_combined_stool.005/chimeras.txt usearch_checked_combined_stool.006/chimeras.txt usearch_checked_combined_stool.007/chimeras.txt usearch_checked_combined_stool.008/chimeras.txt usearch_checked_combined_stool.009/chimeras.txt usearch_checked_combined_stool.010/chimeras.txt usearch_checked_combined_stool.011/chimeras.txt usearch_checked_combined_stool.012/chimeras.txt usearch_checked_combined_stool.013/chimeras.txt usearch_checked_combined_stool.014/chimeras.txt usearch_checked_combined_stool.015/chimeras.txt usearch_checked_combined_stool.016/chimeras.txt usearch_checked_combined_stool.017/chimeras.txt usearch_checked_combined_stool.018/chimeras.txt usearch_checked_combined_stool.019/chimeras.txt usearch_checked_combined_stool.020/chimeras.txt usearch_checked_combined_stool.021/chimeras.txt usearch_checked_combined_stool.022/chimeras.txt usearch_checked_combined_stool.023/chimeras.txt usearch_checked_combined_stool.024/chimeras.txt usearch_checked_combined_stool.025/chimeras.txt usearch_checked_combined_stool.026/chimeras.txt usearch_checked_combined_stool.027/chimeras.txt usearch_checked_combined_stool.028/chimeras.txt usearch_checked_combined_stool.029/chimeras.txt usearch_checked_combined_stool.030/chimeras.txt usearch_checked_combined_stool.031/chimeras.txt usearch_checked_combined_stool.032/chimeras.txt usearch_checked_combined_stool.033/chimeras.txt usearch_checked_combined_stool.034/chimeras.txt usearch_checked_combined_stool.035/chimeras.txt usearch_checked_combined_stool.036/chimeras.txt usearch_checked_combined_stool.037/chimeras.txt usearch_checked_combined_stool.038/chimeras.txt usearch_checked_combined_stool.039/chimeras.txt usearch_checked_combined_stool.040/chimeras.txt usearch_checked_combined_stool.041/chimeras.txt usearch_checked_combined_stool.042/chimeras.txt usearch_checked_combined_stool.043/chimeras.txt usearch_checked_combined_stool.044/chimeras.txt usearch_checked_combined_stool.045/chimeras.txt usearch_checked_combined_stool.046/chimeras.txt usearch_checked_combined_stool.047/chimeras.txt usearch_checked_combined_stool.048/chimeras.txt usearch_checked_combined_stool.049/chimeras.txt usearch_checked_combined_stool.050/chimeras.txt usearch_checked_combined_stool.051/chimeras.txt usearch_checked_combined_stool.052/chimeras.txt usearch_checked_combined_stool.053/chimeras.txt usearch_checked_combined_stool.054/chimeras.txt usearch_checked_combined_stool.055/chimeras.txt usearch_checked_combined_stool.056/chimeras.txt usearch_checked_combined_stool.057/chimeras.txt usearch_checked_combined_stool.058/chimeras.txt usearch_checked_combined_stool.059/chimeras.txt usearch_checked_combined_stool.060/chimeras.txt usearch_checked_combined_stool.061/chimeras.txt usearch_checked_combined_stool.062/chimeras.txt usearch_checked_combined_stool.063/chimeras.txt usearch_checked_combined_stool.064/chimeras.txt usearch_checked_combined_stool.065/chimeras.txt usearch_checked_combined_stool.066/chimeras.txt usearch_checked_combined_stool.067/chimeras.txt usearch_checked_combined_stool.068/chimeras.txt usearch_checked_combined_stool.069/chimeras.txt usearch_checked_combined_stool.070/chimeras.txt usearch_checked_combined_stool.071/chimeras.txt usearch_checked_combined_stool.072/chimeras.txt usearch_checked_combined_stool.073/chimeras.txt usearch_checked_combined_stool.074/chimeras.txt usearch_checked_combined_stool.075/chimeras.txt usearch_checked_combined_stool.076/chimeras.txt usearch_checked_combined_stool.077/chimeras.txt usearch_checked_combined_stool.078/chimeras.txt usearch_checked_combined_stool.079/chimeras.txt usearch_checked_combined_stool.080/chimeras.txt usearch_checked_combined_stool.081/chimeras.txt usearch_checked_combined_stool.082/chimeras.txt usearch_checked_combined_stool.083/chimeras.txt usearch_checked_combined_stool.084/chimeras.txt usearch_checked_combined_stool.085/chimeras.txt usearch_checked_combined_stool.086/chimeras.txt usearch_checked_combined_stool.087/chimeras.txt usearch_checked_combined_stool.088/chimeras.txt usearch_checked_combined_stool.089/chimeras.txt usearch_checked_combined_stool.090/chimeras.txt usearch_checked_combined_stool.091/chimeras.txt usearch_checked_combined_stool.092/chimeras.txt usearch_checked_combined_stool.093/chimeras.txt usearch_checked_combined_stool.094/chimeras.txt usearch_checked_combined_stool.095/chimeras.txt usearch_checked_combined_stool.096/chimeras.txt usearch_checked_combined_stool.097/chimeras.txt usearch_checked_combined_stool.098/chimeras.txt usearch_checked_combined_stool.099/chimeras.txt > chimeras_stool.txt
filter_fasta.py -f combined_fasta_stool/combined_seqs.fna -o combined_fasta_stool/combined_nonchimera_seqs.fna -s chimeras_stool.txt -n;
rm -rf usearch_checked_combined_*.0*




## 8. Create OTU picking parameter file, and run the entire QIIME open reference picking pipeline with usearch61 for reference picking and usearch61_ref for de novo OTU picking
echo "pick_otus:similarity 0.97" > clustering_params.txt
echo "assign_taxonomy:similarity 0.97" >> clustering_params.txt
echo "parallel_align_seqs_pynast:template_fp /home/jhuang/REFs/SILVA_132_QIIME_release/core_alignment/80_core_alignment.fna" >> clustering_params.txt
echo "assign_taxonomy:reference_seqs_fp /home/jhuang/REFs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna" >> clustering_params.txt
echo "assign_taxonomy:id_to_taxonomy_fp /home/jhuang/REFs/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt" >> clustering_params.txt
echo "alpha_diversity:metrics chao1,observed_otus,shannon,PD_whole_tree" >> clustering_params.txt
pick_open_reference_otus.py -r/home/jhuang/REFs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -i combined_fasta/combined_nonchimera_seqs.fna -o clustering/ -p clustering_params.txt --parallel
pick_open_reference_otus.py -r/home/jhuang/REFs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -i combined_fasta_swab/combined_nonchimera_seqs.fna -o clustering_swab/ -p clustering_params.txt --parallel
pick_open_reference_otus.py -r/home/jhuang/REFs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -i combined_fasta_stool/combined_nonchimera_seqs.fna -o clustering_stool/ -p clustering_params.txt --parallel






#9.1 for control data
summarize_taxa_through_plots.py -i clustering34/otu_table_mc2_w_tax_no_pynast_failures.biom -o plots/taxa_summary34 -s
mv usearch_checked_combined usearch_checked_combined_ctrl
mv combined34_fasta combined34_fasta_ctrl
mv clustering34 clustering34_ctrl
mv plots plots_ctrl


## 9.2 for other data: core diversity analyses
core_diversity_analyses.py -o./core_diversity_e4753 -i./clustering/otu_table_mc2_w_tax_no_pynast_failures.biom -m./map_corrected.txt -t./clustering/rep_set.tre -e4753 -p./clustering_params.txt
core_diversity_analyses.py -o./core_diversity_e4753_swab -i./clustering_swab/otu_table_mc2_w_tax_no_pynast_failures.biom -m./map_corrected_swab.txt -t./clustering_swab/rep_set.tre -e4753 -p./clustering_params.txt
core_diversity_analyses.py -o./core_diversity_e4753_stool -i./clustering_stool/otu_table_mc2_w_tax_no_pynast_failures.biom -m./map_corrected_stool.txt -t./clustering_stool/rep_set.tre -e4753 -p./clustering_params.txt
## --> STEP 12 if using Phyloseq.Rmd


## 10. supplements of core diversity analyses
-----------------------
---- by SampleType ----
gunzip ./core_diversity_e6899_g6/table_mc6899.biom.gz
mkdir ./core_diversity_e6899_g6/taxa_plots_SampleType
collapse_samples.py -m ./map_corrected.txt -b./core_diversity_e6899_g6/table_mc6899.biom --output_biom_fp ./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table.biom --output_mapping_fp ./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_map_corrected.txt --collapse_fields "SampleType"
gzip ./core_diversity_e6899_g6/table_mc6899.biom

sort_otu_table.py -i./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table.biom -o./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted.biom
summarize_taxa.py -i./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted.biom -o./core_diversity_e6899_g6/taxa_plots_SampleType/

plot_taxa_summary.py -i./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted_L2.txt,./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted_L3.txt,./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted_L4.txt,./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted_L5.txt,./core_diversity_e6899_g6/taxa_plots_SampleType/SampleType_otu_table_sorted_L6.txt -o./core_diversity_e6899_g6/taxa_plots_SampleType/taxa_summary_plots/

##alpha diversity##
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/PD_whole_tree.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_PD_whole_tree -n 9999
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/chao1.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_chao1 -n 9999
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/observed_otus.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_observed_otus -n 9999
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/shannon.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_shannon -n 9999
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/PD_whole_tree.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_PD_whole_tree_tt -t parametric
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/chao1.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_chao1_tt -t parametric
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/observed_otus.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_observed_otus_tt -t parametric
compare_alpha_diversity.py -i./core_diversity_e6899_g6/arare_max6899/alpha_div_collated/shannon.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/arare_max6899_SampleType/compare_shannon_tt -t parametric

##beta diversity statistics##
make_distance_boxplots.py -d./core_diversity_e6899_g6/bdiv_even6899/weighted_unifrac_dm.txt -f"SampleType" -o./core_diversity_e6899_g6/bdiv_even6899_SampleType/weighted_unifrac_boxplots/ -m ./map_corrected.txt --save_raw_data -n 999
make_distance_boxplots.py -d./core_diversity_e6899_g6/bdiv_even6899/unweighted_unifrac_dm.txt -f"SampleType" -o./core_diversity_e6899_g6/bdiv_even6899_SampleType/unweighted_unifrac_boxplots/ -m ./map_corrected.txt --save_raw_data -n 999
#make_distance_boxplots.py -d./core_diversity_e6899_g6/bdiv_even6899/unweighted_unifrac_dm.txt -f"SampleType" -o./core_diversity_e6899_g6/bdiv_even6899_SampleType/unweighted_unifrac_boxplots/ -m ./map_corrected.txt -g png
compare_categories.py --method adonis -i./core_diversity_e6899_g6/bdiv_even6899/unweighted_unifrac_dm.txt -m./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/bdiv_even6899_SampleType/adonis_out -n 999
compare_categories.py --method anosim -i./core_diversity_e6899_g6/bdiv_even6899/unweighted_unifrac_dm.txt -m./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/bdiv_even6899_SampleType/unweighted_anosim_out -n 999
compare_categories.py --method anosim -i./core_diversity_e6899_g6/bdiv_even6899/weighted_unifrac_dm.txt -m./map_corrected.txt -c "SampleType" -o./core_diversity_e6899_g6/bdiv_even6899_SampleType/weighted_anosim_out -n 999

##using even.biom file to generate group significance##
gunzip ./core_diversity_e6899_g6/table_even6899.biom.gz
group_significance.py -i./core_diversity_e6899_g6/table_even6899.biom -m./map_corrected.txt -c "SampleType" -s kruskal_wallis -o./core_diversity_e6899_g6/group_significance_SampleType_kw_ocs.txt --biom_samples_are_superset --print_non_overlap
group_significance.py -i./core_diversity_e6899_g6/table_even6899.biom -m./map_corrected.txt -c "SampleType" -s g_test -o./core_diversity_e6899_g6/group_significance_SampleType_gtest_ocs.txt
gzip ./core_diversity_e6899_g6/table_even6899.biom



## 11. mv index.html .index.html
cp /media/jhuang/Elements/Data_16S_arckNov_re/core_diversity_e10967/index.html ./


## 12. run Phyloseq.Rmd to get Phyloseq.html (under qiime1-env)
#https://github.com/vaulot/R_tutorials/tree/master/phyloseq
cp ~/DATA/Data_Laura_16S/core_diversity_e5347/Phyloseq.Rmd ./

R -e "rmarkdown::render('Phyloseq.Rmd',output_file='Phyloseq.html')"
#http://girke.bioinformatics.ucr.edu/CSHL_RNAseq/mydoc/mydoc_Rbasics_13/
install.packages("rmarkdown")
install.packages("rmdformats")
install.packages("kableExtra")
library(readxl)
library(dplyr)
library(kableExtra)
library(knitr)
library(rmarkdown)
setwd("~/DATA/Data_Laura_16S_2/core_diversity_e4769")
#"/home/jhuang/DATA/Data_Fran_16S_Exp27/core_diversity_e44117_Exp27"
rmarkdown::render('Phyloseq.Rmd',output_file='Phyloseq.html')


#-Outline
#-Data import
#-Normalization and Filtering
#-Relative Abundance on Phylum and Class level

#-Alpha Diversity
#-Beta Diversity
#-Differential abundance analysis by diet
#-Differential abundance analysis by flora

#https://mikheyevlab.github.io/drosophila-microbiome-selection/
#https://github.com/vaulot/R_tutorials/blob/master/phyloseq/Phyloseq_tutorial.Rmd
