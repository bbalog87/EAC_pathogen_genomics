
# Hands-on Lab: Pathogen Genomic Data Manipulation with Linux Tools

This practical session guides trainees through handling and manipulating bioinformatics data formats using standard Unix/Linux tools and essential bioinformatics utilities (SeqKit, Samtools, BCFtools). It is based on real data folders provided under:

```
./bams/         # Contains BAM files (aligned reads)
./fastas/       # Contains FASTA files (reference genomes/genes)
./rawReads/     # Contains gzipped FASTQ files (raw sequencing reads)
./vcf/          # Contains VCF files (variant calls)
```

## 🎯 Learning Objectives

- Familiarize with standard file types used in pathogen genomics
- Practice core Linux commands (`ls`, `grep`, `cut`, `awk`, `sort`, `wc`, `head`, `tail`)
- Use **SeqKit** to process FASTA/FASTQ
- Use **Samtools** to view and summarize BAM
- Use **BCFtools** to query and filter VCF
- Chain Unix tools to solve real-world genomics problems

---

## 1. 📁 FASTA File Operations

```bash
# Count number of sequences
grep -c '^>' fastas/Bacillus.fasta

# View first 10 lines
head -n 10 fastas/Bacillus.fasta

# List all headers (sequence IDs)
grep '^>' fastas/Bacillus.fasta

# Summary statistics
seqkit stats fastas/*.fasta

# Filter entries by pattern
seqkit grep -p "16S" fastas/16s_rRNA.fasta -o fastas/filtered_16s.fasta

# Extract first 100bp of each sequence
seqkit subseq -r 1:100 fastas/Bacillus.fasta -o fastas/Bacillus.100bp.fasta

# Convert FASTA to table
seqkit fx2tab -n -l fastas/Bacillus.fasta | cut -f1,2
```

---

## 2. 📦 FASTQ File Operations

```bash
# Show first 2 reads (8 lines)
zcat rawReads/Monkeypox_PT0747_2023_R1_001.fastq.gz | head -n 8

# Count number of reads
zgrep -c '^@' rawReads/Monkeypox_PT0747_2023_R1_001.fastq.gz

# Convert to FASTA
seqkit fq2fa rawReads/2140_S1_L001_R1_001.fastq.gz -o 2140_S1_L001_R1_001.fasta

# Subsample 5000 reads
seqkit sample -n 5000 rawReads/2140_S1_L001_R1_001.fastq.gz -o subsample_R1.fastq.gz

# Filter reads containing ACGT
seqkit grep -p "ACGT" rawReads/2140_S1_L001_R1_001.fastq.gz -o filtered_ACGT.fastq.gz
```

---

## 3. 🧬 BAM File Analysis with Samtools

```bash
# View alignments
samtools view bams/S8_mpxv.bam | head

# Mapping statistics
samtools flagstat bams/S8_mpxv.bam

# Sort and index
samtools sort bams/S8_mpxv.bam -o bams/S8_mpxv.sorted.bam
samtools index bams/S8_mpxv.sorted.bam

# Read counts per contig
samtools idxstats bams/S8_mpxv.sorted.bam

# Extract unmapped reads
samtools view -f 4 bams/S8_mpxv.sorted.bam > bams/unmapped.sam

# High-quality reads (MAPQ ≥ 30)
samtools view -b -q 30 bams/S8_mpxv.sorted.bam > bams/highMQ.bam

# Convert BAM to FASTQ
samtools fastq bams/S8_mpxv.sorted.bam -1 reads1.fastq -2 reads2.fastq
```

---

## 4. 🧬 VCF File Analysis with grep/awk/bcftools

```bash
# Count variants (exclude header lines)
grep -vc '^#' vcf/S8_mpxv.vcf

# Extract basic variant details
grep -v '^#' vcf/S8_mpxv.vcf | cut -f1,2,4,5 | head

# Filter high-quality variants (QUAL > 30)
bcftools filter -i 'QUAL>30' -O v -o vcf/S8_mpxv.filtered.vcf vcf/S8_mpxv.vcf

# Query selected fields
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' vcf/S8_mpxv.vcf > vcf/variants_summary.txt

# Variant statistics
bcftools stats vcf/S8_mpxv.vcf > vcf/variant_stats.txt
```

---

## 5. ⚙️ Advanced Subsetting & Pipelining

```bash
# Extract first 200 lines of FASTA
head -n 200 fastas/Bacillus.fasta > fastas/subset.fasta

# Extract contigs with "ATG" in header
grep -A1 'ATG' fastas/Bacillus.fasta > ATG_contigs.fa

# Unique chromosomes in VCF
cut -f1 vcf/S8_mpxv.vcf | grep -v '^#' | sort | uniq -c

# Count bases in all FASTA files
wc -m fastas/*.fasta

# Chained filtering: variants on chr1 with DP > 20
grep -v '^#' vcf/S8_mpxv.vcf | awk '$1=="chr1" && $8~/DP=([2-9][0-9]|[1-9][0-9]{2,})/' | wc -l
```

---

## ✅ Expected Outcomes

By the end of this practicum, you should be able to:

- Manipulate common bioinformatics file formats using shell tools
- Use seqkit, samtools, and bcftools for file inspection, filtering, and conversion
- Extract meaningful summaries from sequencing data
- Chain Unix tools into reproducible, scalable workflows
