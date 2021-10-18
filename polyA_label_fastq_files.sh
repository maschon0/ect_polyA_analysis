#!/bin/bash
# source activate bookend_env
folder=/groups/nodine/lab/general/sequencing_files/2020/200217/M9684
outfolder=$folder/end_labeled_fastq
INDEX=/groups/nodine/lab/members/Schon/lab_files/Annotations/arabidopsis_ensembl/Ensembl.40.ERCC.GFP.kallisto
cd $folder
mkdir -p $outfolder
mkdir -p kallisto_quant

samplenames=( pb1 pb2 pb3 pb4 pb5 pb6 )
for sample in ${samplenames[@]}
do
    echo "Processing $sample..."
    echo "    Trimming adapters..."
    ADAPTER1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    ADAPTER2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
    cutadapt -j 2 -a $ADAPTER1 -A $ADAPTER2 -o $sample.cutadapt.1.fastq -p $sample.cutadapt.2.fastq $folder/fastq/$sample.1.fastq.gz $folder/fastq/$sample.2.fastq.gz
    echo "    End labeling..."
    labelS="AAGCAGTGGTATCAACGCAGAGTACGGG"
    labelE="GAGTACTTTTTTTTTTTTTTTTTTTT+"
    OUT1=$outfolder/$sample.labeled.1.fastq
    OUT2=$outfolder/$sample.labeled.2.fastq
    OUTSINGLE=$outfolder/$sample.labeled.single.fastq
    python fastq_end_label.py -S $labelS -E $labelE --minlen 20 --qualmask 16 --out1 $OUT1 --out2 $OUT2 --single_out $OUTSINGLE $sample.cutadapt.1.fastq $sample.cutadapt.2.fastq
    echo "    Zipping end-labeled files..."
    gzip $OUT1
    gzip $OUT2
    gzip $OUTSINGLE
    echo "    Deleting temp files..."
    rm $sample.cutadapt.1.fastq $sample.cutadapt.2.fastq
    echo "    Quantifying with kallisto"
    kallisto quant -i $INDEX $OUT1.gz $OUT2.gz -o $sample
    mv $sample/abundance.tsv kallisto_quant/$sample.kallisto.tsv
    rm -R $sample
done

