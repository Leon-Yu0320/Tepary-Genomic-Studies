# Tepary-Genomic-Studies
This repo includes the pipeline for tepary resequencing data and population genomic studies

**Author: Li'ang Yu, Magda Julkowska, Andrew Nelson**\
**Date: May 17th, 2023**

Tepary bean diversity panel includes **342 samples** for population genomic studies. The sequencing was performed with Illumina short reads platform (PE 150 bp) along with the depth ranged from 7.78 to 16.69 (calculated based on reference genome size: 512.6 Mb). Please refer [**summary**](https://docs.google.com/spreadsheets/d/1NlpSLPAWTAu47ZDhKmWScodSrdvKM34cjbY-7dVPOvA/edit#gid=0) for details, including the reads number, reads trimming ratio, mapping rate, unique mapping rate, and etc. All raw reads of these accessions are stored on in-house server workspace (see the following workflow for details). 

- [**Genomic analysis of Tepary bean diversity panel**](#genomic-analysis-of-tepary-bean-diversity-panel)
  - [**Sequence alignment and variants calling for the first batch**](#sequence-alignment-and-variants-calling-for-the-first-batch)
    - [Download and pre-process all samples](#download-and-pre-process-all-samples)
    - [Merge fastqs of accessions with subreads](#merge-fastqs-of-accessions-with-subreads)
    - [Mapping reads to reference genome](#mapping-reads-to-reference-genome)
    - [Remove duplicates within the bam files (parameter 1)](#remove-duplicates-within-the-bam-files-parameter-1)
    - [Remove duplicates within the bam files (parameter 2)](#remove-duplicates-within-the-bam-files-parameter-2)
    - [SetNmMdAndUqTags for bam files](#setnmmdanduqtags-for-bam-files)
    - [Mark Duplicates and index bam files](#mark-duplicates-and-index-bam-files)
    - [Using haplotypecaller of gatk4 to generate variants](#using-haplotypecaller-of-gatk4-to-generate-variants)
    - [Genotype VCFs](#genotype-vcfs)
    - [Combine VCFs](#combine-vcfs)

## **Sequence alignment and variants calling for the first batch**
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
      I=${uniq_dir}/${i}.uniq.bam \
      O=${MarkDup_dir}/${i}.mark_dup.bam \
      METRICS_FILE=${MarkDup_dir}/${i}.markdup.metric.txt

  samtools index -@ 80 ${MarkDup_dir}/${i}.mark_dup.bam
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

### Genotype VCFs
```bash
vcf_dir="/data/home/lyu/01_project/03_DIV/07_VCF"
MarkDup_dir="/data/home/nelsonlab/tepary_DNA/07_MarkDup"
Ref_dir="/data/home/lyu/01_project/03_DIV/02_Genome"

for i in $(cat $Ref_dir/Tepary-sample);
do
  gatk HaplotypeCaller \
      --java-options  "-Xmx4g -XX:ParallelGCThreads=1" \
      -R $Ref_dir/GhirsutumCoker_698_v1.0.fa \
      -I ${MarkDup_dir}/${i}.mark_dup.bam \
      -O ${vcf_dir}/${i}.g.vcf.gz \
      -ERC GVCF \
      -ploidy 2 \
      --max-alternate-alleles 6 \
      --max-genotype-count 1024 & 
done

```

### Combine VCFs
```bash
vcf_dir="/data/home/lyu/01_project/03_DIV/07_VCF"
MarkDup_dir="/home/lyu/01_project/03_DIV/06_markDup"
Ref_dir="/data/home/lyu/01_project/03_DIV/02_Genome"
work_dir="/data/home/lyu/01_project/03_DIV"

gatk CombineGVCFs \
  --java-options  "-Xmx512g -XX:ParallelGCThreads=80" \
  -R $Ref_dir/GhirsutumCoker_698_v1.0.fa \
  -V ${work_dir}/TAMU.gvcf.list \
  -O ${vcf_dir}/TAMU-combined.gvcf \
```
