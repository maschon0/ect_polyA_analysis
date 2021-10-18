## Scripts for Poly(A) site analysis  
### Arribas-Hernández et al. (2021) eLife  

For labeling polyadenylation site clusters from Smart-seq2 data,  
scripts included in this folder were used in conjunction with the nanoPARE pipeline
(https://github.com/Gregor-Mendel-Institute/nanoPARE).

Scripts were executed in the following order:  
1. **polyA_label_fastq_files.sh**; Trims adapter sequences and records which reads had an oligo-d(T) adapter removed by appending a "_TAG=E" suffix to each FASTQ read name.  
2. Align FASTQ files to the TAIR10 genome using **endMap** from the nanoPARE pipeline.
3. **extract_polyA.sh**; After alignment, BAM reads containing an "end tag" are extracted and coverted to a BEDGRAPH file. Read 3' ends are compared to the TAIR10 genome, and those upstream of a purine-rich region are rerouted as "false" ends to a separate BEDGRAPH file.
4. Identify poly(A) site clusters using **endGraph** from the nanoPARE pipeline.
5. **annotate_polyA_sites.sh**;  Combines poly(A) site clusters from the three samples of each genotype to define a subset of poly(A) clusters that exist in at least two replicates of either genotype.
6.  **polyA_features.R**; Refines poly(A) clusters by removing clusters overlapping more "false" 3' end reads than reads that passed filtering.
7. **analyze_polyA_data.R**; Comparative analysis of poly(A) cluster distribution across wild-type and mutant genotypes. Used to generate plots for Figure 5 F-G.

