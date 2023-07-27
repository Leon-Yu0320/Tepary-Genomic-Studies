
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
    - [combine gvcfs based on list](#combine-gvcfs-based-on-list)
    - [Genotype GCVCFs using combined panel file](#genotype-gcvcfs-using-combined-panel-file)
  - [Filter variants](#filter-variants)
    - [Reformat the VCF with index information added for the SNP filter](#reformat-the-vcf-with-index-information-added-for-the-snp-filter)
    - [Filter the variants of multiple samples](#filter-the-variants-of-multiple-samples)
      - [Generate contig list for split purpose](#generate-contig-list-for-split-purpose)
      - [Split VCFs by chromosomes](#split-vcfs-by-chromosomes)
      - [GATK VariantFiltration using fixed parameter](#gatk-variantfiltration-using-fixed-parameter)
    - [Filter the 4DTV loci using reference genome annotation](#filter-the-4dtv-loci-using-reference-genome-annotation)
    - [Perform the LD-prune to reduce the SNPS with strong linkage](#perform-the-ld-prune-to-reduce-the-snps-with-strong-linkage)
    - [Annotate variants using SnpEff](#annotate-variants-using-snpeff)
    - [Annotate variants using VEP (Variant Effect Predictor)](#annotate-variants-using-vep-variant-effect-predictor)
    - [Extract the variants annotation of interested genes](#extract-the-variants-annotation-of-interested-genes)
    - [Build the phylogeny using IQtree](#build-the-phylogeny-using-iqtree)
    - [Extract the SNPs/INDELs data for IGV visulization purpose](#extract-the-snpsindels-data-for-igv-visulization-purpose)

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

### combine gvcfs based on list
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
 --java-options  "-Xmx512g -XX:ParallelGCThreads=50" \
 -R $Ref_dir/Pacutifolius_580_v1.0.fa \
 -V $work_dir/TeparyDNA.gvcf \
 -o $work_dir/TeparyDNA.vcf  \
```
## Filter variants
### Reformat the VCF with index information added for the SNP filter

```bash
work_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
vcf_dir="/data/home/nelsonlab/tepary_DNA"

### Extract the header information 
grep "#" $vcf_dir/TeparyDNA.vcf  > $workdir/TeparyVCF-header

### reformat the index and combine with the header
grep -v "#" $vcf_dir/TeparyDNA.vcf  | awk '{print $1,$2,$1"_"$2,$0}' | sed 's/ /\t/g' | cut -f1-3,7- | cat $workdir/TeparyVCF-header - > Tepary-clean.vcf

bgzip -@100 Tepary-clean.vcf
tabix -p vcf Tepary-clean.vcf.gz
```
Or use the Bcftools to edit:
```bash
bcftools annotate --set-id +'%CHROM\_%POS' $vcf_dir/TeparyDNA-tree.gvcf.gz > $work_dir/TeparyDNA-tree-clean.gvcf.gz
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

### Filter the variants of multiple samples
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

#### Split VCFs by chromosomes 
```bash
work_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
hard_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter"

cd $work_dir

for i in $(seq -w 01 11); do
    bcftools view Tepary-clean.vcf --regions Chr$i | bgzip -@ 100 > Chr$i.vcf.gz;
    tabix -p vcf Chr$i.vcf.gz
done

# For scaffolds
bcftools view 500.vcf.gz --regions-file Scaffold.list | bgzip -@ 20 > Chr00.vcf.gz
tabix -p vcf Chr00.vcf.gz
```

#### GATK VariantFiltration using fixed parameter
**SNPs filtering**

```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
SNPs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/SNPs"

for i in $(seq -w 01 11); do
    # Select SNPs only 
    gatk SelectVariants \
        --java-options  "-Xmx512g -XX:ParallelGCThreads=100" \
        -R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $VCF_dir/Chr$i.vcf.gz \
        -select-type SNP \
        -O $SNPs_dir/Chr$i.SNP.vcf.gz

    # Hard filter for SNPs
    gatk VariantFiltration \
        --java-options  "-Xmx512g -XX:ParallelGCThreads=50" \
        -R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $SNPs_dir/Chr$i.SNP.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "my_snp_filter" \
        -O $SNPs_dir/Chr$i.mark.vcf.gz

    # Extract the hard filtered SNPs
    bcftools view $SNPs_dir/Chr$i.mark.vcf.gz \
        --apply-filters PASS | \
        bgzip -@ 100 \
        > $SNPs_dir/Chr$i.filtered.vcf.gz

        tabix -p vcf $SNPs_dir/Chr$i.filtered.vcf.gz

    # Remove loci with missing <0.1 and maf < 0.05
    vcftools --gzvcf $SNPs_dir/Chr$i.filtered.vcf.gz --maf 0.05 --max-missing 0.7 --recode --stdout | bgzip -@ 100 > $SNPs_dir/Chr$i.filtered2.vcf.gz
    tabix -p vcf $SNPs_dir/Chr$i.filtered2.vcf.gz

    # Extract the bi-allilic SNPs
    bcftools view $SNPs_dir/Chr$i.filtered2.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 100 > $SNPs_dir/Chr$i.biSNP.vcf.gz
    tabix -p vcf $SNPs_dir/Chr$i.biSNP.vcf.gz

    # Final SNP set with missing cut of 0.3
    bcftools view -i 'F_MISSING<0.3' $SNPs_dir/Chr$i.biSNP.vcf.gz | bgzip -@ 100 > $SNPs_dir/Chr$i.SNP.cleanVCF.gz
done

```

**INDELs filtering**
```bash
Ref_dir="/data/home/nelsonlab/tepary_DNA"
VCF_dir="/data/home/nelsonlab/tepary_DNA/09_Filter"
INDELs_dir="/data/home/nelsonlab/tepary_DNA/09_Filter/02_HardFilter/INDELs"

for i in $(seq -w 00 11); do
    # Select Indel
    gatk SelectVariants \
        -R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $VCF_dir/Chr$i.vcf.gz \
        -select-type INDEL \
        -O $INDELs_dir/Chr$i.INDEL.all.vcf.gz
    # GATK hard filter

    gatk VariantFiltration \
        -R $Ref_dir/Pacutifolius_580_v1.0.fa \
        -V $INDELs_dir/Chr$i.INDEL.all.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "my_snp_filter" \
        -O $INDELs_dir/Chr$i.INDEL.mark.vcf.gz

    # Extract the Indels passing the hard filter and further filter
    bcftools view $INDELs_dir/Chr$i.INDEL.mark.vcf.gz --apply-filters PASS \
        --types indels -m 2 -M 2 -i 'F_MISSING<0.3' --min-af 0.05 | \
        bgzip -@ 80 \
        > $INDELs_dir/Chr$i.INDEL.filtere.vcf.gz
done
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
Numbers of 4DTV loci
```
10,068,611 
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



### Annotate variants using SnpEff
**Configure and install Snpeff** \
**NOTE:** modify the configuration file: snpeff.config
```
# SnpEff configuration file
# Databases are stored here
# E.g.: Information for 'hg19' is stored in data.dir/hg19/
# You can use tilde ('~') as first character to refer to your home directory. 
# Also, a non-absolute path will be relative to config's file dir

data.dir = /home/lyu/03_software/snpEff/data
```
**Run Snpeff**
```bash
#!/bin/bash

snp_dir="/home/lyu/03_software/snpEff"
genome_dir="/data/home/lyu/01_project/02_Cotton/01_Genome"
cd $snp_dir

### linke reference genome
mkdir $snp_dir/data/Cotton
ln -s $genome_dir/Ghirsutum_527_v2.1.gene.gff3 $snp_dir/data/Cotton/genes.gff

### build libraries
snpEff -Xmx512g build -gff3 Cotton -c snpEff.config

### perform annotation
snpEff  -Xmx512g \
    Cotton ${vcf_dir}/cotton_clean.vcf.recode.vcf > ${snp_dir}/Cotton.ann.vcf
```

### Annotate variants using VEP (Variant Effect Predictor) 
**NOTE**: (Have to use the gff3 with the exon information)
This version is only applicable on vash server due to configuration issue

```bash
genome_dir="/home/liangyu/1_Project/5_Cotton/0_data"

### prepare tehe VCF and ref genome file
mkdir $HOME/vep_data
ln -s $genome_dir/Ghirsutum_527_v2.1.gene.gff3 $HOME/vep_data
scp lyu@ghibli.bti.cornell.edu:/home/lyu/01_project/02_Cotton/02_Variants/01_Filter/cotton_clean.vcf.recode.vcf $HOME/vep_data/

### Reformat the gff3
grep -v "#" Ghirsutum_527_v2.1.gene.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Gh.gff.gz
tabix -p gff Gh.gff.gz
```
**Using Singularity**
```bash
### Perform the annotation 
singularity exec $HOME/6_Docker/vep.sif \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/cotton_clean.vcf.recode.vcf \
    --gff $HOME/vep_data/Gh.gff.gz \
    --fasta $HOME/vep_data/Ghirsutum_527_v2.0.fa
```
**Using docker**
```bash
docker pull ensemblorg/ensembl-vep:release_106.1

sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ensemblorg/ensembl-vep \
    vep --dir $HOME/vep_data \
    -i $HOME/vep_data/cotton_clean.vcf.recode.vcf \
    --gff $HOME/vep_data/Gh.gff.gz \
    --fasta $HOME/vep_data/Ghirsutum_527_v2.0.fa
```

### Extract the variants annotation of interested genes
```bash
 for i in $(cat /home/liangyu/vep_data/01_Results/GeneList);
 do 
  ### grep genes
  grep $i Cotton_effect_output.txt >> Select_Varaints;

done
```

### Build the phylogeny using IQtree
**NOTE** THE first row of phylip file reprsesents the outgroup of TREE!
```bash
#!/bin/bash

script_dir="/data/home/lyu/03_software/vcf2phylip"
vcf_dir="/data/home/lyu/01_project/02_Cotton/04_Phylogeny"

### Convert vcf to phylip
python $script_dir/vcf2phylip.py -i $vcf_dir/cotton.prune.vcf

### REPLACE "*" into "N" also replace the sequnece ID into certain sample name
IFS=$'\n';
for LINE in $(cat $work_dir/ID_table | grep -v 'ID');do

    #READ DIRECTORIES
    ID=$(echo ${LINE} | awk '{ print $1}')
    Accession=$(echo ${LINE} | awk '{ print $2 }')

    #Repalace ID as accessions names
    sed -i "s/$ID/$Accession/g" $work_dir/cotton.prune.min4.phy.treefile 
done

sed -i 's/*/N/g' cotton.prune.min4.phy 

### IQ-tree construction 
iqtree -s cotton.prune.min4.phy -nt 125 -m MFP
```

### Extract the SNPs/INDELs data for IGV visulization purpose
```bash
### modify the header (sample name)
sed -i 's/SRR6311566/AC_134_CB_4029/g' cotton_setID.vcf
sed -i 's/SRR6311571/CB_4012_KING_KARAJAZSY/g' cotton_setID.vcf
sed -i 's/SRR6311717/FELISTANA_UA_7_18/g' cotton_setID.vcf
sed -i 's/SRR6311732/FUNTUA_FT_5/g' cotton_setID.vcf
sed -i 's/SRR6311780/LOCKETT_BXL/g' cotton_setID.vcf
sed -i 's/SRR6311500/TAM_86III_26/g' cotton_setID.vcf
sed -i 's/SRR6311857/TAM_91C_34/g' cotton_setID.vcf
sed -i 's/SRR6311856/TAM_94WE_37S/g' cotton_setID.vcf
sed -i 's/SRR6311570/COKER_310/g' cotton_setID.vcf
sed -i 's/SRR6311744/PD3/g' cotton_setID.vcf
sed -i 's/SRR6311842/TIPO_CHACO_UA_4_4/g' cotton_setID.vcf
sed -i 's/SRR6311538/WESTERN_STOOMPROOF/g' cotton_setID.vcf
sed -i 's/SRR6311572/MEXICO_910/g' cotton_setID.vcf
sed -i 's/SRR6311867/VIR_7094_COKER_310/g' cotton_setID.vcf
sed -i 's/SRR6311863/VIR_7223/g' cotton_setID.vcf
sed -i 's/SRR6311864/VIR_7153_D_10/g' cotton_setID.vcf
sed -i 's/SRR6311820/PLAINS/g' cotton_setID.vcf
sed -i 's/SRR6311853/VIR_6615_MCU_5/g' cotton_setID.vcf

### compress files
bgzip -c cotton_setID.vcf > cotton_setID.vcf.gz
 
### Index the vcf.gz
bcftools index cotton_setID.vcf.gz

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13;
do
 
  ### Extract the regions of interests
  bcftools view cotton_setID.vcf.gz --regions D$i > D$i.vcf 
  bcftools view cotton_setID.vcf.gz --regions D$i > D$i.vcf 

  ### Index the VCF file
  bcftools index D$i.vcf.gz
  bcftools index A$i.vcf.gz 

done
```


