#!/bin/bash

# cd /net/gmi.oeaw.ac.at/nodine/lab/members/Schon/nanoPARE_collab/peter_brodersen/protoplast/scripts
#cd /groups/nodine/lab/members/Schon/nanoPARE_collab/peter_brodersen/protoplast/scripts


for sample in ECT_wt ect234
do
  for type in median consensus
  do
    echo $sample $type
    python get_polyA_consensus.py \
      -G /groups/nodine/lab/members/Schon/lab_files/Annotations/TAIR10/TAIR10_ERCC.fas \
      --type $type \
      --features Araport11_terminal_exons.bed \
      -I ../polyA_bedgraph/$sample.E.plus.bedgraph ../polyA_bedgraph/$sample.E.minus.bedgraph \
      > ../polyA_annotation/$sample.$type.bed
  done
done


