#!/bin/zsh

PREFIX=$1
BAM=$2

echo '=============================================================='
echo Processing $BAM
echo '=============================================================='

echo Calculating coverage: bedtools genomecov -ibam $BAM -d -du -strand + -split \> ${PREFIX}_minus.bed

time bedtools genomecov -ibam $BAM -d -du -strand + -split > ${PREFIX}_minus.bed
cat ${PREFIX}_minus.bed | perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (-1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_minus.bedgraph

#echo Sorting bedgraph: sort -k1,1 -k2,2n ${PREFIX}_minus.bedgraph \| bgzip \> ${PREFIX}_minus.bedgraph.gz

#time sort -k1,1 -k2,2n ${PREFIX}_minus.bedgraph | bgzip > ${PREFIX}_minus.bedgraph.gz

#echo Indexing bedgraph: tabix -s 1 -b 2 -e 3 -0 -f ${PREFIX}_minus.bedgraph.gz

#time tabix -s 1 -b 2 -e 3 -0 -f ${PREFIX}_minus.bedgraph.gz

echo Calculating coverage: bedtools genomecov -ibam $BAM -d -du -strand - -split \> ${PREFIX}_plus.bed

time bedtools genomecov -ibam $BAM -d -du -strand - -split > ${PREFIX}_plus.bed
cat ${PREFIX}_plus.bed | perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_plus.bedgraph

#echo Sorting bedgraph: sort -k1,1 -k2,2n ${PREFIX}_plus.bedgraph \| bgzip \> ${PREFIX}_plus.bedgraph.gz

#time sort -k1,1 -k2,2n ${PREFIX}_plus.bedgraph | bgzip > ${PREFIX}_plus.bedgraph.gz

#echo Indexing bedgraph: tabix -s 1 -b 2 -e 3 -0 -f ${PREFIX}_plus.bedgraph.gz

#time tabix -s 1 -b 2 -e 3 -0 -f ${PREFIX}_plus.bedgraph.gz

echo
echo '============================== DONE =============================='
echo
