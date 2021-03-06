#' Scripts to process WES data
#' These snipetts were provided by John Zhang lab at MDACC.


Somatic SNPs are called using   path/mutect/1.1.4/muTect-1.1.4.jar
( --reference_sequence Homo_sapiens_assembly19.fasta
--dbsnp dbsnp_132_b37.leftAligned.vcf
 --cosmic /path/hg19/annotations/cosmic/hg19_cosmic_v54_120711.vcf)



Pindel calls are using path/pindel/current/pindel  (version 0.2.5b9)
(path/pindel/0.2.5b9/pindel  -f Homo_sapiens_assembly19.fasta)


Platypus
$python Platypus.py  callVariants  --refFile=Homo_sapiens_assembly19.fasta --bamFiles=in.bam --output=o.vcf  --minBaseQual=1

Steps to preprcess -

java -Xmx8g -jar MarkDuplicates.jar I=in.bwa.rg.bam O=o.markbam M=o.marktsv  VALIDATION_STRINGENCY=SILENT  REMOVE_DUPLICATES=false

 java -Xmx8g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R Homo_sapiens_assembly19.fasta -I o.markbam -o o.realignlist -nt 16 \
 -known Mills_and_1000G_gold_standard.indels.b37.vcf -known 1000g/1000G_phase1.indels.b37.vcf --allow_potentially_misencoded_quality_scores

 java -Xmx8g -jar GenomeAnalysisTK.jar -T IndelRealigner  -R Homo_sapiens_assembly19.fasta  -I $omarkbam -targetIntervals $orealignlist \
 -known Mills_and_1000G_gold_standard.indels.b37.vcf -known 1000g/1000G_phase1.indels.b37.vcf   -o o.realignbam --allow_potentially_misencoded_quality_scores

 java -Xmx8g -jar GenomeAnalysisTK.jar -T BaseRecalibrator  -R Homo_sapiens_assembly19.fasta -I o.realignbam \
  -knownSites dbsnp_132.b37.vcf -knownSites Mills_and_1000G_gold_standard.indels.b37.vcf  -knownSites 1000G_phase1.indels.b37.vcf \
  -allowPotentiallyMisencodedQuals -o o.recalibcsv

 java -Xmx8g -jar GenomeAnalysisTK.jar -T PrintReads   -R Homo_sapiens_assembly19.fasta -I  o.realignbam -BQSR o.recalibcsv -o o.bwa_recalibed.bam";



##
picard/1.119/MarkDuplicates.jar
gatk/3.3-0/GenomeAnalysisTK.jar
