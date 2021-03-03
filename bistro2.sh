#!/bin/bash

SS=samplesheet.tsv
FF=fastq
REF=ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
AMPL=amplicons_hg38.bed
SAF=$AMPL.saf

NCPU=$(nproc)
NCPU=$((NCPU-1))

awk 'OFS="\t" {print $4, $1, $2, $3, "."}' $AMPL > $SAF

for BASE in $(cut -f1 $SS | sed 1d) ; do

  FQZ1=$FF/$(grep -w $BASE $SS | cut -f2)
  FQZ2=$FF/$(grep -w $BASE $SS | cut -f3)
  BAM=$FF/$BASE.bam
  VCF=$BAM.vcf
  FQ1=$(echo $FQZ1 | sed 's/.gz/-trimmed-pair1.fastq/')
  FQ2=$(echo $FQZ1 | sed 's/.gz/-trimmed-pair2.fastq/')

  if [ ! -r $BAM ] ; then

    skewer -q 20 $FQZ1 $FQZ2

    biscuit align -@ $NCPU $REF $FQ1 $FQ2 \
    | samblaster | samtools sort -o $BAM -O BAM -

    rm $FQ1 $FQ2

    samtools index $BAM

    samtools flagstat $BAM > $BAM.flagstat.txt

    biscuit pileup -u -o $VCF $REF $BAM

    bgzip -f $VCF

    tabix -p vcf $VCF.gz

  fi

  biscuit vcf2bed -e -t cg $VCF.gz > $VCF.cg.bed

  awk '{OFS="\t"}{print $1"."$2,$1,$2,"+",$NF,$(NF-1), 1-$(NF-1)}' $VCF.cg.bed > $FF/$BASE.methylkit.tsv

  biscuit vcf2bed -e -t ch $VCF.gz > $VCF.ch.bed

  featureCounts -F SAF -a $AMPL.saf -o $BAM.tsv $BAM

  sed 1,2d $BAM.tsv | cut -f1,7 | sed "s/^/$BASE\t/" > tmp && mv tmp $BAM.tsv

  sed "s/^/$BASE\t/" $BAM.tsv.summary | sed 1d > tmp && mv tmp $BAM.tsv.summary

done

cat $FF/*bam.tsv > $FF/3col_ontarget.tsv

cat $FF/*bam.tsv.summary > $FF/3col_offtarget.tsv
