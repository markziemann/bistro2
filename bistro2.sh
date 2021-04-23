#!/bin/bash

SS=samplesheet.tsv
FF=fastq
REF=ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
GENEBED=ref/Homo_sapiens.GRCh38.103.gtf.bed
AMPL=amplicons_hg38.bed
SAF=$AMPL.saf
AMPL_FA=$AMPL.fa
AMPL_FA_STARTS=$AMPL.starts.fa
AMPL_FA_ENDS=$AMPL.ends.fa

NCPU=$(nproc)
NCPU=$((NCPU-1))

awk 'OFS="\t" {print $4, $1, $2, $3, "."}' $AMPL > $SAF
bedtools getfasta -fi $REF -bed $AMPL > $AMPL_FA

cut -c-20 $AMPL_FA \
| tr '[:lower:]' '[:upper:]' \
| sed 's/CG/YG/g' | sed 's/C/T/g' > $AMPL_FA_STARTS

revseq -sequence $AMPL_FA -outseq /dev/stdout \
| perl unwrap_fasta.pl - - | cut -d ' ' -f1 | cut -c-20 \
| tr '[:lower:]' '[:upper:]' \
| sed 's/CG/YG/g' | sed 's/C/T/g' > $AMPL_FA_ENDS

for BASE in $(cut -f1 $SS | sed 1d) ; do

  FQZ1=$FF/$(grep -w $BASE $SS | cut -f2)
  FQZ2=$FF/$(grep -w $BASE $SS | cut -f3)
  BAM=$FF/$BASE.bam
  VCF=$BAM.vcf
  FQ1=$(echo $FQZ1 | sed 's/.gz/-trimmed-pair1.fastq/')
  FQ2=$(echo $FQZ1 | sed 's/.gz/-trimmed-pair2.fastq/')

  if [ ! -r $BAM ] ; then

    skewer -q 20 $FQZ1 $FQZ2

    # have removed samblaster because it was excluding reads with out proper mates
    # for WGBS run samplblaster after biscuit and before samtools sort
    biscuit align -@ $NCPU $REF $FQ1 $FQ2 \
    | samtools sort -o $BAM -O BAM -


    rm $FQ1 $FQ2

    samtools index $BAM

    samtools flagstat $BAM > $BAM.flagstat.txt

    # -p means no filtering out unproper pairs
    biscuit pileup -p -u -o $VCF $REF $BAM

    bgzip -f $VCF

    tabix -p vcf $VCF.gz

  fi

  biscuit vcf2bed -e -t cg $VCF.gz > $VCF.cg.bed

  awk '{OFS="\t"}{print $1"."$2,$1,$2,"+", $NF, $(NF-1)*100 , 100-($(NF-1)*100) }' $VCF.cg.bed > $FF/$BASE.methylkit.tsv

  echo "Y.1 Y 1 + 2 1 1" | tr ' ' '\t' >> $FF/$BASE.methylkit.tsv

  biscuit vcf2bed -e -t ch $VCF.gz > $VCF.ch.bed

  featureCounts -F SAF -a $AMPL.saf -o $BAM.tsv $BAM

  sed 1,2d $BAM.tsv | cut -f1,7 | sed "s/^/$BASE\t/" > tmp && mv tmp $BAM.tsv

  sed "s/^/$BASE\t/" $BAM.tsv.summary | sed 1d > tmp && mv tmp $BAM.tsv.summary

  bedtools genomecov -bg -ibam $BAM > $BAM.bedgraph

  echo "chr start end n_reads chr start end gene overlap" | tr ' ' '\t' > $BAM.bedgraph.intersect

  awk '$4>=100' $BAM.bedgraph | bedtools merge -c 4 -o max \
  | bedtools intersect -wo -a - -b $AMPL >> $BAM.bedgraph.intersect

  echo "chr start end n_reads length_bp" | tr ' ' '\t' > $BAM.bedgraph.offtarget

  awk '$4>=100' $BAM.bedgraph | bedtools merge -c 4 -o max \
  | bedtools intersect -v -a - -b $AMPL \
  | awk '{OFS="\t"} {print $0,$3-$2}' >> $BAM.bedgraph.offtarget

  awk '$4>=100' $BAM.bedgraph | bedtools merge -c 4 -o max \
  | bedtools intersect -v -a - -b $AMPL \
  | awk '{OFS="\t"}{print $0,$3-$2}' \
  | sed 's/\t/_/4' \
  | bedtools getfasta -name+ -fi $REF -bed - \
  | paste - - | tr ':' '\t' | tr -s '\t' \
  | awk '{FS="\t"} {print ">",$2,$3,$1,"bp",$4}' \
  | sed 's/ //' |  sed 's/ /:/' | sed 's/>//2' \
  | sed 's/ //2' | sed 's/ /\n/2' > $BAM.bedgraph.offtarget.fa

  cat $BAM.bedgraph.offtarget.fa | paste - - \
  | awk 'length($3)>60'  | awk '{ print $1,$2,substr($3,1,20) }' \
  | sed 's/ /\n/2' \
  | tr '[:lower:]' '[:upper:]' \
  | sed 's/CG/YG/g' | sed 's/C/T/g' > $BAM.bedgraph.offtarget.starts.fa

  revseq -sequence $BAM.bedgraph.offtarget.fa -outseq /dev/stdout \
  | perl unwrap_fasta.pl - -  \
  | sed 's/ Reversed: / /' | paste - - | awk 'length($3)>60'  \
  | awk '{ print $1,$2,substr($3,1,20) }' \
  | sed 's/ /\n/2' \
  | tr '[:lower:]' '[:upper:]' \
  | sed 's/CG/YG/g' | sed 's/C/T/g' > $BAM.bedgraph.offtarget.ends.fa

  echo "chr start end n_reads chr start end gene dist_bp" | tr ' ' '\t' > $BAM.bedgraph.offtarget.genes

  awk '$4>=100' $BAM.bedgraph | bedtools merge -c 4 -o max \
  | bedtools intersect -v -a - -b $AMPL \
  | bedtools closest -d -a - -b $GENEBED >> $BAM.bedgraph.offtarget.genes

done

cat $FF/*bam.tsv > $FF/3col_ontarget.tsv

cat $FF/*bam.tsv.summary > $FF/3col_offtarget.tsv
