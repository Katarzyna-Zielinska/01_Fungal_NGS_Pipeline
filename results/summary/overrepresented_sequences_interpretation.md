# Interpretation of Overrepresented Sequences

## 1. Objective

Investigate overrepresented sequences reported by FastQC in the trimmed Candida auris sequencing reads and assess whether their origin can be established using the available reference genome and annotation.

## 2. Methods

- Inspected the Overrepresented sequences module in existing FastQC ZIP reports.
- Examined 20 selected sequences using BLASTN-short against the Candida auris B8441 V2 reference genome.
- Compared the best-hit coordinates with annotated 5S rRNA features in the GFF annotation.
- Compared the selected motifs with the Illumina adapter sequences used in the Cutadapt log for technical sample SRR26185932.
- Inspected existing FastQC reports for trimmed R1 and R2 reads from susceptible samples SRR10059687, SRR10059688, and SRR10059689.

## 3. Results

### 3.1 Persistence after trimming

Overrepresented sequences remained detectable in trimmed reads from all three inspected susceptible samples. Similar motifs were observed in both R1 and R2.

The most frequent selected motif, CCTATTCTGGGGTCCCGGGGGCCAAATACCAAAAGGGTGTTGGCGGTCAA, represented approximately:
- 0.357% of reads in R1 and 0.351% in R2 for SRR10059687;
- 0.305% in R1 and 0.312% in R2 for SRR10059688;
- 0.505% in R1 and 0.502% in R2 for SRR10059689.

These percentages refer to the exact motif reported by FastQC, not to the sum of all related motifs.

### 3.2 Reference-genome matches

Most selected sequences had full-length 50-nt matches at 100% identity in the best-hit summary. Other reported alignments had 98.039% identity and an alignment length of 51 nt.

The hits clustered in two reference regions:
- PEKT02000002.1, approximately 1.236–1.279 Mb;
- PEKT02000005.1, approximately 438–1,339 bp.

The selected sequences overlap or resemble one another; therefore, the 20 entries should not be interpreted as 20 independent repeat families.

### 3.3 Relationship to annotated 5S rRNA

Annotated 5S rRNA features occur near some hits on PEKT02000002.1. The coordinate comparison did not identify direct overlap between the selected best-hit intervals and annotated 5S rRNA features.

Hits on PEKT02000005.1 were several kilobases from the nearest annotated 5S rRNA feature. Proximity alone does not establish a functional relationship.

### 3.4 Adapter comparison

No exact shared fragment of at least 12 nucleotides was detected between the 20 selected sequences and the two adapter sequences examined. This screening test does not exclude shorter or imperfect matches, nor does it establish the origin of the motifs.

The Cutadapt trimming counts documented in this project refer to technical sample SRR26185932. They should not be treated as adapter-trimming statistics for all six cohort samples.

## 4. Interpretation and limitations

The available results establish that selected overrepresented motifs persist after trimming and that similar motifs occur in both read directions in the three inspected susceptible samples. Several selected motifs match regions of the reference genome.

The origin and biological function of these sequences remain unresolved. The available evidence does not justify classifying them as 5S rRNA, adapter contamination, or contamination from an external organism.

This inspection covered trimmed R1 and R2 reports from three susceptible samples. It does not by itself establish the same FastQC pattern for all six cohort samples. Further investigation would be needed to identify the motifs conclusively.

## 5. Conclusion

The overrepresented-sequence warning should be documented as an unresolved QC observation. It should not be treated as proof of sample contamination or as a reason to remove sequences without additional evidence.

## 6. Supporting files

- `top_overrepresented_sequences.fasta`
- `top_overrepresented_sequences_blast.tsv`
- `top_overrepresented_sequences_blast_summary.tsv`
- `top_overrepresented_sequences_exact_counts.tsv`
