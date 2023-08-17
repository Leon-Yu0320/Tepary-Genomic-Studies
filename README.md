
# **Genomic analysis of Tepary bean diversity panel**

**Author: Li'ang Yu, Andrew Nelson**\
**Date: May 17th, 2023**

Tepary bean diversity panel includes **342 samples** for population genomic studies. The sequencing was performed with Illumina short reads platform (PE 150 bp) along with the depth ranged from 7.78 to 16.69 (calculated based on reference genome size: 512.6 Mb). Please refer [**summary**](https://docs.google.com/spreadsheets/d/1NlpSLPAWTAu47ZDhKmWScodSrdvKM34cjbY-7dVPOvA/edit#gid=0) for details, including the reads number, reads trimming ratio, mapping rate, unique mapping rate, and etc. All raw reads of these accessions are stored on in-house server workspace (see the following workflow for details). 


- [**Genomic analysis of Tepary bean diversity panel**](#genomic-analysis-of-tepary-bean-diversity-panel)
  - [Varaints calling data process](#varaints-calling-data-process)
    - [Download read from common beans as outgroup samples (speceis) for phylogeny](#download-read-from-common-beans-as-outgroup-samples-speceis-for-phylogeny)
    - [Download reads of the 13 unique accessions from publications for phylogeny](#download-reads-of-the-13-unique-accessions-from-publications-for-phylogeny)
    - [Download and pre-process all samples](#download-and-pre-process-all-samples)
    - [Merge fastqs of accessions with subreads](#merge-fastqs-of-accessions-with-subreads)
    - [Mapping reads to reference genome](#mapping-reads-to-reference-genome)
    - [Remove duplicates within the bam files (parameter 1)](#remove-duplicates-within-the-bam-files-parameter-1)
    - [Remove duplicates within the bam files (parameter 2)](#remove-duplicates-within-the-bam-files-parameter-2)
    - [SetNmMdAndUqTags for bam files](#setnmmdanduqtags-for-bam-files)
    - [Mark Duplicates and index bam files](#mark-duplicates-and-index-bam-files)
    - [Using haplotypecaller of gatk4 to generate variants](#using-haplotypecaller-of-gatk4-to-generate-variants)
    - [combine gvcfs based on list (The two list were comprised of 364 and 363 samples)](#combine-gvcfs-based-on-list-the-two-list-were-comprised-of-364-and-363-samples)
    - [Genotype GCVCFs using combined panel file](#genotype-gcvcfs-using-combined-panel-file)
  - [Filter Variants for PCA and phylogeny (outgourp speceis included)](#filter-variants-for-pca-and-phylogeny-outgourp-speceis-included)
      - [Generate contig list for split purpose](#generate-contig-list-for-split-purpose)
    - [Split VCFs by chromosomes for filtering](#split-vcfs-by-chromosomes-for-filtering)
    - [GATK VariantFiltration using fixed parameter](#gatk-variantfiltration-using-fixed-parameter)
    - [Select the bi-allilic SNPs and INDELs](#select-the-bi-allilic-snps-and-indels)
    - [Reformat the VCF with index information added for the SNP filter](#reformat-the-vcf-with-index-information-added-for-the-snp-filter)
    - [Filter the 4DTV loci using reference genome annotation](#filter-the-4dtv-loci-using-reference-genome-annotation)
    - [Perform the LD-prune to reduce the SNPS with strong linkage](#perform-the-ld-prune-to-reduce-the-snps-with-strong-linkage)
  - [Filter Variants for GWAS and other population genomics studies (outgourp speceis excluded)](#filter-variants-for-gwas-and-other-population-genomics-studies-outgourp-speceis-excluded)
    - [Split VCFs to chromosomes](#split-vcfs-to-chromosomes)
    - [GATK hard filteration](#gatk-hard-filteration)
    - [Select bi-allilic variants and filter based on MAF score](#select-bi-allilic-variants-and-filter-based-on-maf-score)
    - [Rename variants ID](#rename-variants-id)
    - [Merge chromosomal-level SNPs and INDELs data for annotation](#merge-chromosomal-level-snps-and-indels-data-for-annotation)
  - [Annotate variants using VEP (Variant Effect Predictor)](#annotate-variants-using-vep-variant-effect-predictor)
  - [Build the phylogeny using IQtree](#build-the-phylogeny-using-iqtree)
  - [Build the PCA of panel](#build-the-pca-of-panel)

## Varaints calling data process
### Download read from common beans as outgroup samples (speceis) for phylogeny

[**SRR**](https://acsess.onlinelibrary.wiley.com/doi/10.3835/plantgenome2017.08.0068) ID and [**publication**](https://acsess.onlinelibrary.wiley.com/doi/10.3835/plantgenome2017.08.0068) for common bean samples were highlighted. For sample SRR5807699, using the ENA web browser to download

```
/vol1/fastq/SRR580/009/SRR5807699/SRR5807699_1.fastq.gz
/vol1/fastq/SRR580/009/SRR5807699/SRR5807699_2.fastq.gz
```
**Download** from [**ENA page**](https://www.ebi.ac.uk/ena/browser/view/SRR5807699)
```bash

 ### forward reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR580/009/SRR5807699/SRR5807699_1.fastq.gz \
    	/data/home/nelsonlab/tepary_DNA/01_CommonBean/

 ### reverse reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR580/009/SRR5807699/SRR5807699_2.fastq.gz \
    	/data/home/nelsonlab/tepary_DNA/01_CommonBean/

```

### Download reads of the 13 unique accessions from publications for phylogeny 
See punblication [**Bio project info**](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=607288)

The 13 accession list
```
/vol1/fastq/SRR114/078/SRR11455478/SRR11455478_1.fastq.gz /vol1/fastq/SRR114/078/SRR11455478/SRR11455478_2.fastq.gz
/vol1/fastq/SRR114/087/SRR11455487/SRR11455487_1.fastq.gz /vol1/fastq/SRR114/087/SRR11455487/SRR11455487_2.fastq.gz
/vol1/fastq/SRR114/040/SRR11455540/SRR11455540_1.fastq.gz /vol1/fastq/SRR114/040/SRR11455540/SRR11455540_2.fastq.gz
/vol1/fastq/SRR114/051/SRR11455551/SRR11455551_1.fastq.gz /vol1/fastq/SRR114/051/SRR11455551/SRR11455551_2.fastq.gz
/vol1/fastq/SRR114/055/SRR11455555/SRR11455555_1.fastq.gz /vol1/fastq/SRR114/055/SRR11455555/SRR11455555_2.fastq.gz
/vol1/fastq/SRR114/086/SRR11455486/SRR11455486_1.fastq.gz /vol1/fastq/SRR114/086/SRR11455486/SRR11455486_2.fastq.gz
/vol1/fastq/SRR114/060/SRR11455560/SRR11455560_1.fastq.gz /vol1/fastq/SRR114/060/SRR11455560/SRR11455560_2.fastq.gz
/vol1/fastq/SRR114/083/SRR11455483/SRR11455483_1.fastq.gz /vol1/fastq/SRR114/083/SRR11455483/SRR11455483_2.fastq.gz
/vol1/fastq/SRR114/092/SRR11455492/SRR11455492_1.fastq.gz /vol1/fastq/SRR114/092/SRR11455492/SRR11455492_2.fastq.gz
/vol1/fastq/SRR114/095/SRR11455495/SRR11455495_1.fastq.gz /vol1/fastq/SRR114/095/SRR11455495/SRR11455495_2.fastq.gz
/vol1/fastq/SRR114/033/SRR11455533/SRR11455533_1.fastq.gz /vol1/fastq/SRR114/033/SRR11455533/SRR11455533_2.fastq.gz
/vol1/fastq/SRR114/041/SRR11455541/SRR11455541_1.fastq.gz /vol1/fastq/SRR114/041/SRR11455541/SRR11455541_2.fastq.gz
/vol1/fastq/SRR114/059/SRR11455559/SRR11455559_1.fastq.gz /vol1/fastq/SRR114/059/SRR11455559/SRR11455559_2.fastq.gz
```

```bash

work_dir="/data/home/nelsonlab/tepary_DNA/02_NC-panel"

IFS=$'\n';
for LINE in $(cat $work_dir/sample.list);
do
	pair1=$(echo ${LINE} | awk '{ print $1}')
	pair2=$(echo ${LINE} | awk '{ print $2}')

  ### forward reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:$pair1 \
    	/data/home/nelsonlab/tepary_DNA/02_NC-panel/

  ### reverse reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:$pair2 \
    	/data/home/nelsonlab/tepary_DNA/02_NC-panel/
done
```


### Download and pre-process all samples
**for the first batch**
```bash
### Trim reads

raw_dir="/data/home/nelsonlab/tepary_DNA/first_batch/01.RawData/01_fastq"
trimmed_dir="/data/home/nelsonlab/tepary_DNA/first_batch/02_CleanReads"
work_dir="/data/home/nelsonlab/tepary_DNA/first_batch"

### Create acession list 
ls $raw_dir | sed 's/1.fq/2.fq/g' list
cat $raw_dir/list | sort | uniq | sed s'/_2.fq.gz/_/g' | grep -v 'list' > $work_dir/Batch1.list

for i in $(cat $work_dir/Batch1.list);
do
  fastp \
    --in1 ${raw_dir}/${i}1.fq.gz --out1 ${trimmed_dir}/${i}_trim_R1.fastq.gz \
    --in2 ${raw_dir}/${i}2.fq.gz --out2 ${trimmed_dir}/${i}_trim_R2.fastq.gz \
    --thread 100 \
    --trim_front1 5 --trim_tail1 5 \
    --cut_front --cut_front_window_size 3 \
    --cut_front_mean_quality 15 --cut_tail \
    --cut_tail_window_size 3 --cut_tail_mean_quality 15 \
    --n_base_limit 5 --average_qual 10 --length_required 75 \
    --dedup --dup_calc_accuracy 6 --reads_to_process 450000000
done
```
**for the second batch**
```bash
raw_dir="/data/home/nelsonlab/tepary_DNA/second_batch/01.RawData/01_fastq"
trimmed_dir="/data/home/nelsonlab/tepary_DNA/second_batch/02_CleanReads"
work_dir="/data/home/nelsonlab/tepary_DNA/second_batch"

### Create trimming list 
ls $raw_dir | sed 's/1.fq/2.fq/g' list
cat $raw_dir/list | sort | uniq | sed s'/_2.fq.gz/_/g' | grep -v 'list' > $work_dir/Batch2.list

for i in $(cat $work_dir/Batch2.list);
do
  fastp \
    --in1 ${raw_dir}/${i}1.fq.gz --out1 ${trimmed_dir}/${i}trim_R1.fastq.gz \
    --in2 ${raw_dir}/${i}2.fq.gz --out2 ${trimmed_dir}/${i}trim_R2.fastq.gz \
    --thread 16 \
    --trim_front1 5 --trim_tail1 5 \
    --cut_front --cut_front_window_size 3 \
    --cut_front_mean_quality 15 --cut_tail \
    --cut_tail_window_size 3 --cut_tail_mean_quality 15 \
    --n_base_limit 5 --average_qual 10 --length_required 75 \
    --dedup --dup_calc_accuracy 6 --reads_to_process 450000000
done
```
### Merge fastqs of accessions with subreads
There are six accessions with multiple reads
```
4	T_137_ZKDN230004067	
4	T_139_ZKDN230004101	
4	T_223_ZKDN230003988	
4	T_434_ZKDN230003985	
4	T_442_ZKDN230004055	
4	T_92_ZKDN230004173	
```

```bash
cat T_137_ZKDN230004067*_trim_R1.fastq.gz > T_137_ZKDN230004067-merge_trim_R1.fastq.gz
cat T_137_ZKDN230004067*_trim_R2.fastq.gz > T_137_ZKDN230004067-merge_trim_R2.fastq.gz

cat T_139_ZKDN230004101*_trim_R1.fastq.gz > T_139_ZKDN230004101-merge_trim_R1.fastq.gz
cat T_139_ZKDN230004101*_trim_R2.fastq.gz > T_139_ZKDN230004101-merge_trim_R2.fastq.gz

cat T_223_ZKDN230003988*_trim_R1.fastq.gz > T_223_ZKDN230003988-merge_trim_R1.fastq.gz
cat T_223_ZKDN230003988*_trim_R2.fastq.gz > T_223_ZKDN230003988-merge_trim_R2.fastq.gz

cat T_434_ZKDN230003985*_trim_R1.fastq.gz > T_434_ZKDN230003985-merge_trim_R1.fastq.gz
cat T_434_ZKDN230003985*_trim_R2.fastq.gz > T_434_ZKDN230003985-merge_trim_R2.fastq.gz

cat T_442_ZKDN230004055*_trim_R1.fastq.gz > T_442_ZKDN230004055-merge_trim_R1.fastq.gz
cat T_442_ZKDN230004055*_trim_R2.fastq.gz > T_442_ZKDN230004055-merge_trim_R2.fastq.gz

cat T_92_ZKDN230004173*_trim_R1.fastq.gz > T_92_ZKDN230004173-merge_trim_R1.fastq.gz
cat T_92_ZKDN230004173*_trim_R2.fastq.gz > T_92_ZKDN230004173-merge_trim_R2.fastq.gz
```

### Mapping reads to reference genome 
```bash
trimmed_dir1="/data/home/nelsonlab/tepary_DNA/first_batch/02_CleanReads"
trimmed_dir2="/data/home/nelsonlab/tepary_DNA/second_batch/02_CleanReads"
Ref_dir="/data/home/nelsonlab/tepary_DNA"
bam_dir="/data/home/nelsonlab/tepary_DNA/03_Alignment"

### Generate the Batch1-mapping file for loop
ls $trimmed_dir1 | cut -d "_" -f1,2 | sort | uniq  | wc -l > $Ref_dir/Batch1-mapping.list

### index the reference genome
bwa index $Ref_dir/Pacutifolius_580_v1.0.fa

### align reads by accessions
for i in $(cat $Ref_dir/Batch1-mapping.list);
do
  bwa mem -M -R "@RG\tID:$i\tSM:$i\tPL:ILLUMINA" \
    -t 4  \
    $Ref_dir/Pacutifolius_580_v1.0.fa \
    ${trimmed_dir1}/${i}*trim_R1.fastq.gz \
    ${trimmed_dir1}/${i}*trim_R2.fastq.gz | samtools view -@10 -bS - -o ${bam_dir}/${i}.bam &

    ### check numbers of reads within each bam file
    samtools flagstat ${bam_dir}/${i}.bam > ${bam_dir}/${i}_flagstat.txt
    
    ### summarize mappeds reads into logStats file
    grep 'mapped' ${bam_dir}/${i}_flagstat.txt >> bwaMapReads.stats
done

### align reads by accessions
for i in $(cat $Ref_dir/Batch2-mapping.list);
do
  bwa mem -M -R "@RG\tID:$i\tSM:$i\tPL:ILLUMINA" \
    -t 4  \
    $Ref_dir/Pacutifolius_580_v1.0.fa \
    ${trimmed_dir1}/${i}_*trim_R1.fastq.gz \
    ${trimmed_dir1}/${i}_*trim_R2.fastq.gz | samtools view -@10 -bS - -o ${bam_dir}/${i}.bam &

    ### check numbers of reads within each bam file
    samtools flagstat ${bam_dir}/${i}.bam > ${bam_dir}/${i}_flagstat.txt
    
    ### summarize mappeds reads into logStats file
    grep 'mapped' ${bam_dir}/${i}_flagstat.txt >> bwaMapReads.stats
done
```

### Remove duplicates within the bam files (parameter 1)
```bash
bam_dir="/data/home/nelsonlab/tepary_DNA/03_Alignment"
uniq_dir="/data/home/nelsonlab/tepary_DNA/04_UniqMap"
Ref_dir="/data/home/nelsonlab/tepary_DNA"

for i in $(cat $Ref_dir/Tepary-sample);
do
  ### remove duplicates

  ### Filterin based on mapping quality
  sambamba view \
    -t 2 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment or mate_is_unmapped or chimeric) and not ([XA] != null or [SA] != null) and proper_pair" \
    ${bam_dir}/${i}.bam \
    -o ${uniq_dir}/${i}.uniq.bam &

  ### Filterin based on mapping quality
  sambamba view \
    -t 2 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
    ${bam_dir}/${i}.bam \
    -o ${uniq_dir}/${i}.uniq.bam &

  ### check numbers of reads within each bam file
  samtools flagstat ${uniq_dir}/${i}.uniq.bam > ${bam_dir}/${i}.uniq.flagstat.txt &

  ### summarize mappeds reads into logStats file 
  grep 'mapped' ${uniq_dir}/${i}.uniq.flagstat.txt >> UniqMapReads.stats

  ### sort bam files
  samtools sort -m 4G ${uniq_dir}/${i}.uniq.bam > ${uniq_dir}/${i}.uniq.sort.bam &

  ### check reads number after sort
  samtools flagstat ${uniq_dir}/${i}.uniq.sort.bam > ${uniq_dir}/${i}.uniqSort.flagstat.txt

done
```

### Remove duplicates within the bam files (parameter 2)
```bash
bam_dir="/data/home/nelsonlab/tepary_DNA/03_Alignment"
uniq_dir="/data/home/nelsonlab/tepary_DNA/05_UniqMap_Para2"
Ref_dir="/data/home/nelsonlab/tepary_DNA"

for i in $(cat $Ref_dir/Tepary-sample);
do

  ### Filterin based on mapping quality and uniq map
  sambamba view \
    -t 4 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
    ${bam_dir}/${i}.bam \
    -o ${uniq_dir}/${i}.uniq.bam &

  ### check numbers of reads within each bam file
  samtools flagstat ${uniq_dir}/${i}.uniq.bam > ${uniq_dir}/${i}.uniq.flagstat.txt &

  ### summarize mappeds reads into logStats file 
  grep 'mapped' ${uniq_dir}/${i}.uniq.flagstat.txt >> UniqMapReads.stats

  ### sort bam files
  samtools sort -m 4G ${uniq_dir}/${i}.uniq.bam > ${uniq_dir}/${i}.uniq.sort.bam &

  ### check reads number after sort
  samtools flagstat ${uniq_dir}/${i}.uniq.sort.bam > ${uniq_dir}/${i}.uniqSort.flagstat.txt
    
done
```
### SetNmMdAndUqTags for bam files
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
uniq_dir="/data/home/nelsonlab/tepary_DNA/04_UniqMap"
fix_dir="/data/home/nelsonlab/tepary_DNA/06_Fix"
picard_dir="/data/home/lyu/mambaforge/share/picard-2.18.23-0"

for i in $(cat $Ref_dir/Tepary-sample);
do
java -Xms1g -Xmx5g -jar \
    ${picard_dir}/picard.jar SetNmMdAndUqTags \
      I=${uniq_dir}/${i}.uniq.sort.bam \
      O=${fix_dir}/${i}.fixed_nm.bam \
      R=$Ref_dir/Pacutifolius_580_v1.0.fa &
done
```

### Mark Duplicates and index bam files
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
uniq_dir="/data/home/nelsonlab/tepary_DNA/04_UniqMap"
MarkDup_dir="/data/home/nelsonlab/tepary_DNA/07_MarkDup"
picard_dir="/data/home/lyu/mambaforge/share/picard-2.18.23-0"

for i in $(cat $Ref_dir/Tepary-sample);
do
  ### markduplicates
  java -Xms100g -Xmx500g -jar \
    ${picard_dir}/picard.jar MarkDuplicates \
      I=${fix_dir}/${i}.fixed_nm.bam \
      O=${MarkDup_dir}/${i}.mark_dup.bam \
      METRICS_FILE=${MarkDup_dir}/${i}.markdup.metric.txt

  samtools index -@ 80 ${MarkDup_dir}/${i}.mark_dup.bam

  samtools flagstat ${MarkDup_dir}/${i}.mark_dup.bam > ${MarkDup_dir}/${i}.mark_dup.flagstat.txt
done
```

### Using haplotypecaller of gatk4 to generate variants 
**gatk4 Version information**

```
Using GATK jar /data/home/lyu/mambaforge/envs/gatk4/share/gatk4-4.2.5.0-0/gatk-package-4.2.5.0-local.jar
Running:
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /data/home/lyu/mambaforge/envs/gatk4/share/gatk4-4.2.5.0-0/gatk-package-4.2.5.0-local.jar --version
The Genome Analysis Toolkit (GATK) v4.2.5.0
HTSJDK Version: 2.24.1
Picard Version: 2.25.4
```

```bash
vcf_dir="/data/home/nelsonlab/tepary_DNA/08_gVCFs"
MarkDup_dir="/data/home/nelsonlab/tepary_DNA/07_MarkDup"
Ref_dir="/data/home/nelsonlab/tepary_DNA"
work_dir="/data/home/nelsonlab/tepary_DNA"

### Create the dict file for reference genome
picard CreateSequenceDictionary \
  R=$Ref_dir/Pacutifolius_580_v1.0.fa \
  O=$Ref_dir/Pacutifolius_580_v1.0.dict

### activate gatk4 environment for mamba 
mamba activate gatk4

### gatk variants calling
for i in $(cat $Ref_dir/Tepary-sample);
do
  gatk HaplotypeCaller \
      --java-options  "-Xmx4g -XX:ParallelGCThreads=1" \
      -R $Ref_dir/Pacutifolius_580_v1.0.fa \
      -I ${MarkDup_dir}/${i}.mark_dup.bam \
      -O ${vcf_dir}/${i}.g.vcf.gz \
      -ERC GVCF \
      -ploidy 2 \
      --max-alternate-alleles 6 \
      --max-genotype-count 1024
done
```
**NOTE:**
Here we merged the VCFs in two types" \
  All samples (commonBean included for phylogeney) \
  Only tepary beans included for other population genomics studies

### combine gvcfs based on list (The two list were comprised of 364 and 363 samples)
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
work_dir="/data/home/nelsonlab/tepary_DNA"
sample_dir="/data/home/nelsonlab/tepary_DNA/08_gVCFs"

gatk CombineGVCFs \
 --java-options  "-Xmx512g -XX:ParallelGCThreads=50" \
 -R $Ref_dir/Pacutifolius_580_v1.0.fa \
 --variant ${sample_dir}/TeparySample.list \
 -O ${work_dir}/TeparyDNA.gvcf
```

### Genotype GCVCFs using combined panel file
```bash
work_dir="/data/home/nelsonlab/tepary_DNA"
Ref_dir="/data/home/nelsonlab/tepary_DNA"

gatk GenotypeGVCFs \
 --java-options  "-Xmx512g -XX:ParallelGCThreads=100" \
 -R $Ref_dir/Pacutifolius_580_v1.0.fa \
 -V $work_dir/TeparyDNA.gvcf \
 -o $work_dir/TeparyDNA.vcf
```

## Filter Variants for PCA and phylogeny (outgourp speceis included)
Check link for defination of [**parameter**](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QUAL_QD_GQ_Formulation_fDG.htm) 

**DP is determined by average total depth (total depth/population size = individual sample)** \
**minGQ=20**: 1% chance that the call is incorrect** \
**minQUAL=20**: 1 % chance that there is no variant at the site**

#### Generate contig list for split purpose
See format below
```
scaffold_15	125718
scaffold_16	91946
scaffold_17	91617
scaffold_21	82945
scaffold_23	79196
scaffold_25	78192
scaffold_26	76070
scaffold_28	73412
scaffold_30	70361
scaffold_32	69992
scaffold_37	67828
```

### Split VCFs by chromosomes for filtering

**The low-missing rate  variants will be retained during the chromosomes Split process**
```bash
vcf_dir="/data/home/nelsonlab/tepary_DNA"
work_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"

### For chromosomes
for i in $(seq -w 01 11); do
    bcftools view ${vcf_dir}/TeparyDNA-tree.gvcf.gz \
    --regions Chr$i \
    -i 'F_MISSING<0.6' | bgzip -@ 100 > ${work_dir}/Chr$i.vcf.gz

    tabix -p vcf ${work_dir}/Chr$i.vcf.gz
done

### For scaffolds
bcftools view ${vcf_dir}/TeparyDNA-tree.gvcf.gz \
  --regions-file ${work_dir}/Scaffold.list \
  -i 'F_MISSING<0.6' | bgzip -@ 100 > ${work_dir}/Chr00.vcf.gz

tabix -p vcf ${work_dir}/Chr00.vcf.gz
```

### GATK VariantFiltration using fixed parameter
**Variants hard filtering**

```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/SNPs"
INDELs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/INDELs"

### Filter SNPs and INDELs using GATK
for i in $(seq -w 00 11); do

    ### Hard filter for SNPs
    gatk VariantFiltration \
        --java-options  "-Xmx800g -XX:ParallelGCThreads=100" \
        -R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $VCF_dir/Chr$i.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "my_snp_filter" \
        -O ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz

    tabix -p vcf ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz

    ### Hard filter for INDELs
      gatk VariantFiltration \
      --java-options  "-Xmx800g -XX:ParallelGCThreads=100" \
      -R $Ref_dir/Pacutifolius_580_v1.0.fa \
      -V $VCF_dir/Chr$i.vcf.gz \
      --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
      --filter-name "my_snp_filter" \
      -O ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz

    tabix -p vcf ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz

done
```
### Select the bi-allilic SNPs and INDELs
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/SNPs"
INDELs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/INDELs"

for i in $(seq -w 00 11); do

    ### Extract the hard filtered SNPs and take only the Bi-allilic SNPs
    bcftools view ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz \
        --apply-filters PASS \
        --types snps -m 2 -M 2  | bgzip -@ 100 > ${SNPs_dir}/Chr$i.filter.SNPs.vcf.gz

    tabix -p vcf ${SNPs_dir}/Chr$i.filter.SNPs.vcf.gz

  ### Extract the hard filtered INDELs and take only the Bi-allilic INDELs
   bcftools view ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz \
        --apply-filters PASS \
        --types indels -m 2 -M 2  | bgzip -@ 100 > ${INDELs_dir}/Chr$i.filter.INDELs.vcf.gz

    tabix -p vcf ${INDELs_dir}/Chr$i.filter.INDELs.vcf.gz

    # Remove loci with missing <0.1 and maf < 0.05
    #vcftools --gzvcf $SNPs_dir/Chr$i.filtered.vcf.gz --maf 0.05 --max-missing 0.7 --recode --stdout | bgzip -@ 100 > $SNPs_dir/Chr$i.filtered2.vcf.gz
    #tabix -p vcf $SNPs_dir/Chr$i.filtered2.vcf.gz

done
  ```

### Reformat the VCF with index information added for the SNP filter

use the Bcftools to annotate the coordiantes information (works for the compressed .gz file)
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/SNPs"
INDELs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/INDELs"

for i in $(seq -w 00 11); do

    bcftools annotate \
      --set-id +'%CHROM\_%POS' \
      --threads 100 \
      ${SNPs_dir}/Chr$i.filter.SNPs.vcf.gz \
      --output-type z  > ${SNPs_dir}/Chr$i.clean.SNPs.vcf.gz
done
```
See below as a example before and after:
**before the filtering**
```
Chr01   27      .        C       G       144.78  .       AC=4;AF=0.018;AN=220;DP=176;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.33>
Chr01   40      .        T       C       139.21  .       AC=4;AF=0.017;AN=240;DP=197;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.33>
Chr01   76      .        A       C       208.64  .       AC=6;AF=0.021;AN=288;DP=280;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.33>
Chr01   80      .        G       A       74.02   .       AC=2;AF=5.988e-03;AN=334;DP=282;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=>      
```
**after the filtering**
```
Chr01   27      Chr01_27        C       G       144.78  .       AC=4;AF=0.018;AN=220;DP=176;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.33>
Chr01   40      Chr01_40        T       C       139.21  .       AC=4;AF=0.017;AN=240;DP=197;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.33>
Chr01   76      Chr01_76        A       C       208.64  .       AC=6;AF=0.021;AN=288;DP=280;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.33>
Chr01   80      Chr01_80        G       A       74.02   .       AC=2;AF=5.988e-03;AN=334;DP=282;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=>   
```

### Filter the 4DTV loci using reference genome annotation

```bash
anno_dir="/data/home/nelsonlab/tepary_DNA"
results_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/04_4DTV"
code_dir="/data/home/lyu/02_script"

# Call 4DTV loci
perl $code_dir/4DTV_scan.pl -g $anno_dir/Pacutifolius_580_v1.0.fa \
    -f $anno_dir/Pacutifolius_580_v1.0.gene_exons.gff3 \
    -o $results_dir/4DTV_list
```
Numbers of 4DTV loci based on annotation
```
10,068,611 
```
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/SNPs"
4DTV_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/04_4DTV"

for i in $(seq -w 01 11); do

    bcftools view \
      -R ${4DTV_dir}/Pacutifolius_4dtv.tsv \
      --threads 100 \
      ${SNPs_dir}/Chr$i.clean.SNPs.vcf.gz \
      --output-type z  > ${4DTV_dir}/Chr$i.4DTV.vcf.gz
done
```

### Perform the LD-prune to reduce the SNPS with strong linkage
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/SNPs"
LD_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/03_LD_prune"


for i in $(seq -w 01 11); do
    ### build the library
    plink --vcf $SNPs_dir/Chr$i.SNP.cleanVCF.gz \
      --make-bed \
      --allow-extra-chr \
      --out $LD_dir/Chr$i.SNP_bfile

    ### Use LD and distance to prune
    plink --bfile $LD_dir/Chr$i.SNP_bfile \
      --indep-pairwise 1000 500 0.2 \
      --allow-extra-chr \
      --out $LD_dir/Chr$i.SNP_pruned

    ### Extract the pruned loci 
    bcftools view \
      --include ID==@Chr$i.SNP_pruned.prune.in \
      $SNPs_dir/Chr$i.SNP.cleanVCF.gz > $LD_dir/Chr$i.prune.vcf
done

```
Summary of the LD-prone results
```
Pruned 32127 variants from chromosome 1, leaving 465.
Pruned 39672 variants from chromosome 2, leaving 522.
Pruned 35088 variants from chromosome 3, leaving 566.
Pruned 21404 variants from chromosome 4, leaving 412.
Pruned 21598 variants from chromosome 5, leaving 384.
Pruned 26563 variants from chromosome 6, leaving 452.
Pruned 33865 variants from chromosome 7, leaving 511.
Pruned 35968 variants from chromosome 8, leaving 521.
Pruned 31452 variants from chromosome 9, leaving 470.
Pruned 18141 variants from chromosome 10, leaving 454.
Pruned 25067 variants from chromosome 11, leaving 513.
```

## Filter Variants for GWAS and other population genomics studies (outgourp speceis excluded)
### Split VCFs to chromosomes
```bash
vcf_dir="/data/home/nelsonlab/tepary_DNA"
work_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/01_filter"

### For chromosomes
for i in $(seq -w 01 11); do
    bcftools view ${vcf_dir}/TeparyDNA-GWAS-genotype.gvcf.gz \
    --regions Chr$i \
    -i 'F_MISSING<0.7' | bgzip -@ 100 > ${work_dir}/Chr$i.vcf.gz

    tabix -p vcf ${work_dir}/Chr$i.vcf.gz
done

### For scaffolds
bcftools view ${work_dir}/TeparyDNA-GWAS-genotype.gvcf.gz \
  --regions-file ${work_dir}/Scaffold.list \
  -i 'F_MISSING<0.7' | bgzip -@ 100 > ${work_dir}/Chr00.vcf.gz

tabix -p vcf ${work_dir}/Chr00.vcf.gz

```
### GATK hard filteration 
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/01_filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/02_Hardfilter/SNPs"
INDELs_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/02_Hardfilter/INDELs"

### Filter SNPs and INDELs using GATK
for i in $(seq -w 00 11); do

    ### Hard filter for SNPs
    gatk VariantFiltration \
        --java-options  "-Xmx800g -XX:ParallelGCThreads=100" \
        -R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $VCF_dir/Chr$i.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "my_snp_filter" \
        -O ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz

    tabix -p vcf ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz

    ### Hard filter for INDELs
      gatk VariantFiltration \
      	--java-options  "-Xmx800g -XX:ParallelGCThreads=100" \
      	-R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $VCF_dir/Chr$i.vcf.gz \
      	--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
      	--filter-name "my_indel_filter" \
      	-O ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz

    tabix -p vcf ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz

done
```

### Select bi-allilic variants and filter based on MAF score
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/02_Hardfilter/SNPs"
INDELs_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/02_Hardfilter/INDELs"

for i in $(seq -w 01 11); do

    ### Extract the hard filtered SNPs and take only the Bi-allilic SNPs
    tabix -p vcf ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz

    bcftools view ${SNPs_dir}/Chr$i.hard.SNPs.vcf.gz \
        --apply-filters PASS \
        --types snps -m 2 -M 2  | bgzip -@ 100 > ${SNPs_dir}/Chr$i.filter.SNPs.vcf.gz

    tabix -p vcf ${SNPs_dir}/Chr$i.filter.SNPs.vcf.gz

  ### Extract the hard filtered INDELs and take only the Bi-allilic INDELs

   tabix -p vcf ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz

   bcftools view ${INDELs_dir}/Chr$i.hard.INDELs.vcf.gz \
        --apply-filters PASS \
        --types indels -m 2 -M 2  | bgzip -@ 100 > ${INDELs_dir}/Chr$i.filter.INDELs.vcf.gz

    tabix -p vcf ${INDELs_dir}/Chr$i.filter.INDELs.vcf.gz

    # Remove loci with missing <0.1 and maf < 0.05
    vcftools --gzvcf $SNPs_dir/Chr$i.filter.SNPs.vcf.gz --maf 0.05 --max-missing 0.7 --recode --stdout | bgzip -@ 100 > $SNPs_dir/Chr$i.MAF-filter.vcf.gz
    tabix -p vcf $SNPs_dir/Chr$i.MAF-filter.vcf.gz

    vcftools --gzvcf $INDELs_dir/Chr$i.filter.INDELs.vcf.gz --maf 0.05 --max-missing 0.7 --recode --stdout | bgzip -@ 100 > $INDELs_dir/Chr$i.MAF-INDELs.vcf.gz
    tabix -p vcf $INDELs_dir/Chr$i.MAF-INDELs.vcf.gz
done

```
### Rename variants ID 
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/10_GWAS/03_MAF"

for i in $(seq -w 01 11); do

    /data/home/nelsonlab/tepary_DNA/10_GWAS/bcftools/bcftools annotate \
      --set-id +'%CHROM\_%POS' \
      --threads 100 \
      ${SNPs_dir}/Chr$i.MAF-filter.vcf.gz \
      --output-type z  > ${SNPs_dir}/Chr$i.MAF-clean.vcf.gz
done

```

### Merge chromosomal-level SNPs and INDELs data for annotation
```bash
cd /data/home/nelsonlab/tepary_DNA/10_GWAS/03_MAF

bcftools concat Chr01.MAF-clean.vcf.gz \
	Chr02.MAF-clean.vcf.gz \
	Chr03.MAF-clean.vcf.gz \
	Chr04.MAF-clean.vcf.gz \
	Chr05.MAF-clean.vcf.gz \
	Chr06.MAF-clean.vcf.gz \
	Chr07.MAF-clean.vcf.gz \
	Chr08.MAF-clean.vcf.gz \
	Chr09.MAF-clean.vcf.gz \
	Chr10.MAF-clean.vcf.gz \
	Chr11.MAF-clean.vcf.gz --output-type z > TeparyBean-GWAS-MAF.vcf.gz

bcftools concat Chr01.MAF-INDELs.vcf.gz \
	Chr02.MAF-INDELs.vcf.gz \
	Chr03.MAF-INDELs.vcf.gz \
	Chr04.MAF-INDELs.vcf.gz \
	Chr05.MAF-INDELs.vcf.gz \
	Chr06.MAF-INDELs.vcf.gz \
	Chr07.MAF-INDELs.vcf.gz \
	Chr08.MAF-INDELs.vcf.gz \
	Chr09.MAF-INDELs.vcf.gz \
	Chr10.MAF-INDELs.vcf.gz \
	Chr11.MAF-INDELs.vcf.gz --output-type z > TeparyBean-MAF-INDELs.vcf.gz

```

## Annotate variants using VEP (Variant Effect Predictor) 
**NOTE**: (Have to use the gff3 with the exon information)
This version is only applicable on vash server due to configuration issue

```bash
genome_dir="/home/liangyu/vep_data"

### prepare tehe VCF and ref genome file
mkdir $HOME/vep_data
ln -s $genome_dir/Pacutifolius_580_v1.0.fa $HOME/vep_data
scp -P 22 lyu@ghibli.bti.cornell.edu:/data/home/nelsonlab/tepary_DNA/10_GWAS/03_MAF/TeparyBean-GWAS-MAF.vcf.gz $HOME/vep_data/
scp -P 22 lyu@ghibli.bti.cornell.edu:/data/home/nelsonlab/tepary_DNA/10_GWAS/03_MAF/TeparyBean-MAF-INDELs.vcf.gz $HOME/vep_data/

### Reformat the gff3
grep -v "#" Pacutifolius_580_v1.0.gene_exons.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Pacutifolius.gff.gz
tabix -p gff Pacutifolius.gff.gz
```
**Using Singularity**
```bash
### Perform the annotation for SNPs
singularity exec $HOME/6_Docker/vep.sif \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/TeparyBean-GWAS-MAF.vcf.gz \
    --gff $HOME/vep_data/Pacutifolius.gff.gz \
    --fasta $HOME/vep_data/Pacutifolius_580_v1.0.fa

### Perform the annotation for INDELs
singularity exec $HOME/6_Docker/vep.sif \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/TeparyBean-MAF-INDELs.vcf.gz \
    --gff $HOME/vep_data/Pacutifolius.gff.gz \
    --fasta $HOME/vep_data/Pacutifolius_580_v1.0.fa
```
**Using docker**
```bash
docker pull ensemblorg/ensembl-vep:release_106.1

sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ensemblorg/ensembl-vep \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/TeparyBean-MAF-INDELs.vcf.gz \
    --gff $HOME/vep_data/Pacutifolius.gff.gz \
    --fasta $HOME/vep_data/Pacutifolius_580_v1.0.fa
```

## Build the phylogeny using IQtree
**NOTE** THE first row of phylip file reprsesents the outgroup of TREE!
```bash
#!/bin/bash
script_dir="/data/home/lyu/03_software/vcf2phylip"
vcf_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/05_tree"

### Convert vcf to phylip
python $script_dir/vcf2phylip.py \
  -i $vcf_dir/Tepary-tree.vcf

### REPLACE "*" into "N" also replace the sequnece ID into certain sample name
sed -i 's/*/N/g' Tepary-tree.min4.phy 

### IQ-tree construction 
iqtree -s Tepary-tree.min4.phy -nt 100 -m MFP
```

## Build the PCA of panel
```bash
PCA_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/06_PCA"

### build the chromosome map
bcftools view -H $PCA_dir/Tepary-tree.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > $PCA_dir/Tepary-tree.chrom-map.txt
vcftools --vcf $PCA_dir/Tepary-tree.vcf --plink --chrom-map $PCA_dir/Tepary-tree.chrom-map.txt --out Tepary-tree

### build the plink index
plink --noweb --file $PCA_dir/Tepary-tree --make-bed --out $PCA_dir/Tepary-tree.vcf.plink

### Build the PCA matrix 
gcta --bfile $PCA_dir/Tepary-tree.vcf.plink --make-grm --out $PCA_dir/Tepary-tree_grm

### perform PCA by gcta program 
gcta --grm $PCA_dir/Tepary-tree_grm --pca 3 --out $PCA_dir/Tepary-tree_grm.PCA
```
