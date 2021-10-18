#!/bin/bash

genome=TAIR10_ERCC.fas
sources=( pb1 pb2 pb3 pb4 pb5 pb6 )

for i in {1..6}
do
    s=${sources[$i]}
    echo "$i ($s)..."
    echo "    ...extracting polyA reads"
    bam=$i/star/Aligned.out.bam
    python sam_to_bed.py --source $s --bed -F $genome $bam > $s.unsorted.bed
    grep -v -P '^#' $s.unsorted.bed | grep -P '([Ee]$|[Ee][^\t]*$)' > $s.E.unsorted.bed    
    echo "    ...sorting bed file"
    sh bed12_sort.sh $s.E.unsorted.bed | python $py_dir/bed_collapse.py --bed12 > $s.E.bed
    
    echo "    ...making bedgraphs"
    python bed_to_bedgraph.py --stranded -D polyA_bedgraphs -N $s --label_includes E --type 3 --bed $s.E.bed
    python bed_to_bedgraph.py --stranded -D polyA_bedgraphs -N $s.false --label_includes e --type 3 --bed $s.E.bed
    echo "$s complete."
done

