
# Hands-on Lab: Pathogen Genomic Data Manipulation with Linux Tools

This lab guides you through handling common genomics data formats using Linux tools: `grep`, `cut`, `sort`, `awk`, `wc`, `seqkit`, `samtools`, and `bcftools`. Use the dataset structured in:

```
./bams/         # BAM alignment files
./fastas/       # Reference sequences (FASTA)
./rawReads/     # Raw sequencing reads (FASTQ)
./vcf/          # Called variants (VCF)
```

---

## 1. 📁 FASTA File Operations

**Q1: Count how many sequences are in `Bacillus.fasta`.**

<details>
<summary>Show Solution</summary>

```bash
grep -c '^>' fastas/Bacillus.fasta
```

</details>

**Q2: Preview the first 10 lines of the file.**

<details>
<summary>Show Solution</summary>

```bash
head -n 10 fastas/Bacillus.fasta
```

</details>

**Q3: Print only sequence headers.**

<details>
<summary>Show Solution</summary>

```bash
grep '^>' fastas/Bacillus.fasta
```

</details>

**Q4: Summarize sequence stats across all FASTA files.**

<details>
<summary>Show Solution</summary>

```bash
seqkit stats fastas/*.fasta
```

</details>

**Q5: Extract sequences matching “16S”.**

<details>
<summary>Show Solution</summary>

```bash
seqkit grep -p "16S" fastas/16s_rRNA.fasta -o fastas/filtered_16s.fasta
```

</details>

**Q6: Extract the first 100bp from each Bacillus contig.**

<details>
<summary>Show Solution</summary>

```bash
seqkit subseq -r 1:100 fastas/Bacillus.fasta -o fastas/Bacillus.100bp.fasta
```

</details>

---

## 2. 📦 FASTQ File Operations

**Q7: Show the first two reads from the Monkeypox R1 FASTQ file.**

<details>
<summary>Show Solution</summary>

```bash
zcat rawReads/Monkeypox_PT0747_2023_R1_001.fastq.gz | head -n 8
```

</details>

**Q8: Count total reads.**

<details>
<summary>Show Solution</summary>

```bash
zgrep -c '^@' rawReads/Monkeypox_PT0747_2023_R1_001.fastq.gz
```

</details>

**Q9: Convert R1 reads to FASTA.**

<details>
<summary>Show Solution</summary>

```bash
seqkit fq2fa rawReads/2140_S1_L001_R1_001.fastq.gz -o fastas/2140_S1_L001_R1_001.fasta
```

</details>

**Q10: Subsample 5000 reads.**

<details>
<summary>Show Solution</summary>

```bash
seqkit sample -n 5000 rawReads/2140_S1_L001_R1_001.fastq.gz -o subsample_R1.fastq.gz
```

</details>

**Q11: Filter reads containing “ACGT”.**

<details>
<summary>Show Solution</summary>

```bash
seqkit grep -p "ACGT" rawReads/2140_S1_L001_R1_001.fastq.gz -o filtered_ACGT.fastq.gz
```

</details>

---

## 3. 🧬 BAM File Analysis with Samtools

**Q12: View the first few alignments.**

<details>
<summary>Show Solution</summary>

```bash
samtools view bams/S8_mpxv.bam | head
```

</details>

**Q13: Check alignment statistics.**

<details>
<summary>Show Solution</summary>

```bash
samtools flagstat bams/S8_mpxv.bam
```

</details>

**Q14: Sort and index the BAM file.**

<details>
<summary>Show Solution</summary>

```bash
samtools sort bams/S8_mpxv.bam -o bams/S8_mpxv.sorted.bam
samtools index bams/S8_mpxv.sorted.bam
```

</details>

**Q15: How many reads map to each contig?**

<details>
<summary>Show Solution</summary>

```bash
samtools idxstats bams/S8_mpxv.sorted.bam
```

</details>

**Q16: Extract only unmapped reads.**

<details>
<summary>Show Solution</summary>

```bash
samtools view -f 4 bams/S8_mpxv.sorted.bam > bams/unmapped.sam
```

</details>

---

## 4. 🧬 VCF File Analysis

**Q17: Count the number of called variants.**

<details>
<summary>Show Solution</summary>

```bash
grep -vc '^#' vcf/S8_mpxv.vcf
```

</details>

**Q18: Extract CHROM, POS, REF, ALT from top 10 records.**

<details>
<summary>Show Solution</summary>

```bash
grep -v '^#' vcf/S8_mpxv.vcf | cut -f1,2,4,5 | head
```

</details>

**Q19: Filter variants with QUAL > 30.**

<details>
<summary>Show Solution</summary>

```bash
bcftools filter -i 'QUAL>30' -O v -o vcf/S8_mpxv.filtered.vcf vcf/S8_mpxv.vcf
```

</details>

**Q20: Export CHROM, POS, QUAL to table.**

<details>
<summary>Show Solution</summary>

```bash
bcftools query -f '%CHROM\t%POS\t%QUAL\n' vcf/S8_mpxv.vcf > variants_summary.txt
```

</details>

---

## ✅ Tips

- Use `man <command>` or `<command> --help` for tool documentation.
- Use `<tab>` key for autocomplete in terminal.
- Pipe commands together to build powerful one-liners.
