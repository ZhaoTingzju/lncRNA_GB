#### Project lncRNA_GB

updated reference: 

**G.arboreum_CRI-A2_assembly_v1.0.fasta (CAS)

**TM-1_V2.1.fa (ZJU)

**Graimondii_221.fa


software used in this study:

Hisat2 version 2.1.0

samtools version 1.6

stringtie 2.0

gffcompare v0.10.4

CPC,py 0.1

HMMER 3.1b2

pfam_scan

TransDecoder V5.3.0

bedtools v2.26.0

the pipline

#### step1 : building genome index

gffread -T -o TM-1_V2.1.gene.gtf TM-1_V2.1.gene.gff3

extract_splice_sites.py TM-1_V2.1.gene.gtf > TM-1_V2.1.ss 

extract_exons.py TM-1_V2.1.gene.gtf > TM-1_V2.1.exon

hisat2-build --ss TM-1_V2.1.ss --exon TM-1_V2.1.exon TM-1_V2.1.fa TM-1_V2.1

#### step 2: mapping

ls *gz|perl -pi -e 's/.R[12].clean.fastq.gz//' > sample.txt

cat sample.txt|while read id; do hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 ${id}.R1.clean.fastq.gz -2 
${id}.R2.clean.fastq.gz -S ${id}.sam

or

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V001-V002.R1.clean.fastq.gz -2 V001-V002.R2.clean.fastq.gz -S V001-V002.sam

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V001-V002-1.R1.clean.fastq.gz -2 V001-V002-1.R2.clean.fastq.gz -S V001-V002-1.sam

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V003-V004.R1.clean.fastq.gz -2 V003-V004.R2.clean.fastq.gz -S V003-V004.sam

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V003-V004-1.R1.clean.fastq.gz -2 V003-V004-1.R2.clean.fastq.gz -S V003-V004-1.sam

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V007-V008.R1.clean.fastq.gz -2 V007-V008.R2.clean.fastq.gz -S V007-V008.sam

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V007-V008-1.R1.clean.fastq.gz -2 V007-V008-1.R2.clean.fastq.gz -S V007-V008-1.sam

hisat2 -p 8 --dta -x ../../ref/TM-1_V2.1 -1 V009-V010.R1.clean.fastq.gz -2 V009-V010.R2.clean.fastq.gz -S V009-V010.sam

If its single ended data and forward stranded you need to set: -rna-strandness F If its paired end data and forward stranded you need to set: -rna-strandness FR

Similar, if its reverse stranded for SE data: -rna-strandness R or (Paired End, reverse stranded) -rna-strandness RF

most strand-specific data are fr-firststrand , use -rna-strandness RF

 ls *sam|while read id; do /public/home/zhaoting/biosoftware/samtools-1.6/samtools sort -@ 1 -o ${id}.sorted.bam $id 
 
or

/public/home/zhaoting/biosoftware/samtools-1.6/samtools sort -@ 1 -o V001-V002-1.sam.sorted.bam V001-V002-1.sam && rm V001-V002-1.sam

/public/home/zhaoting/biosoftware/samtools-1.6/samtools sort -@ 1 -o V001-V002.sam.sorted.bam V001-V002.sam && rm V001-V002.sam

/public/home/zhaoting/biosoftware/samtools-1.6/samtools sort -@ 1 -o V003-V004-1.sam.sorted.bam V003-V004-1.sam && rm V003-V004-1.sam

/public/home/zhaoting/biosoftware/samtools-1.6/samtools sort -@ 1 -o V003-V004.sam.sorted.bam V003-V004.sam && rm V003-V004.sam

#### step 3 :assemble

ls *bam|while read id; do stringtie -G ../../ref/TM-1_V2.1.gene.gtf -p 2 -o ${id}.gtf $id;done

or

stringtie -G ../../ref/TM-1_V2.1.gene.gtf -p 2 -o V001-V002-1.sam.sorted.bam.gtf V001-V002-1.sam.sorted.bam

stringtie -G ../../ref/TM-1_V2.1.gene.gtf -p 2 -o V001-V002.sam.sorted.bam.gtf V001-V002.sam.sorted.bam

stringtie -G ../../ref/TM-1_V2.1.gene.gtf -p 2 -o V003-V004-1.sam.sorted.bam.gtf V003-V004-1.sam.sorted.bam

stringtie -G ../../ref/TM-1_V2.1.gene.gtf -p 2 -o V003-V004.sam.sorted.bam.gtf V003-V004.sam.sorted.bam

add --rf if this data is fr-firstrand


#### step 4: lncRNA prediction

ls *gtf > mergelist.txt

stringtie --merge -o merged.gtf -c 3 ./mergelist.txt

gffcompare -r ../ref/TM-1_V2.1.gene.gtf -p 4 merged.gtf -o merged_lncRNA

awk '$3 == "x"|| $3 == "u"|| $3 == "i" {print $0}' merged_lncRNA.merged.gtf.tmap > novel.gtf.tmap

awk '$11 >200 {print}' novel.gtf.tmap > novel.longRNA.gtf.tmap

awk '{print $5}' novel.longRNA.gtf.tmap | perl ~/zt_script/extract_gtf_by_name.pl merged.gtf - > novel.longRNA.gtf

gffread -g ~/eGWAS/ref/TM-1_V2.1.fa -w exon.fa ./novel.longRNA.gtf

TransDecoder.LongOrfs -t exon.fa # This step generated a file named longest_orfs.ped

pfam_scan.pl -cpu 8 -fasta ./exon.fa.transdecoder_dir/longest_orfs.pep -dir ~/biosoftware/PfamScan/base/ > pfam_scan_eQTL.txt

perl -ne 'print if /Domain/' pfam_scan_eQTL.txt |perl -ne 'print $1."\n" if /::(.*?)::/' > Novel.transcript_with_domain.txt

~/biosoftware/CPC2-beta/bin/CPC2.py -i exon.fa -o cpc_output.txt

perl -ne 'print if /noncoding/' cpc_output.txt |cut -f 1 > Novel_transcript_cpc_nocoding.txt

cat Novel_transcript_cpc_nocoding.txt Novel.transcript_with_domain.txt |sort|uniq -d > intersection.txt

sort Novel_transcript_cpc_nocoding.txt intersection.txt intersection.txt |uniq -u > lncRNA_list.txt

cat lncRNA_list.txt| perl ~/zt_script/extract_gtf_by_name.pl merged.gtf - > LncRNA.gtf

cat  LncRNA.gtf TM-1_V2.1.gene.gtf > lncRNA_mRNA.gtf

rm intersection.txt

sortBed -i lncRNA_mRNA.gtf > 1 && mv 1 lncRNA_mRNA.gtf

rm Novel_transcript_cpc_nocoding.txt

rm Novel.transcript_with_domain.txt

rm pfam_scan_eQTL.txt

#### step5: Estimate transcript abundances

ls *bam|while read id; do stringtie -A ${id}_fpkm.txt --rf -p 8 -G ./lncRNA_mRNA.gtf -o temp.gtf $id;done

or

stringtie -A V001-V002-1.sam.sorted.bam_fpkm.txt -p 8 --rf -G ../../ref/TM-1_V2.1.gene_lncRNA.gtf -o temp.gtf V001-V002-1.sam.sorted.bam

stringtie -A V001-V002.sam.sorted.bam_fpkm.txt -p 8 --rf -G ../../ref/TM-1_V2.1.gene_lncRNA.gtf -o temp.gtf V001-V002.sam.sorted.bam

stringtie -A V003-V004-1.sam.sorted.bam_fpkm.txt -p 8 --rf -G ../../ref/TM-1_V2.1.gene_lncRNA.gtf -o temp.gtf V003-V004-1.sam.sorted.bam

stringtie -A V003-V004.sam.sorted.bam_fpkm.txt -p 8 --rf -G ../../ref/TM-1_V2.1.gene_lncRNA.gtf -o temp.gtf V003-V004.sam.sorted.bam

stringtie -A V007-V008-1.sam.sorted.bam_fpkm.txt -p 8 --rf -G ../../ref/TM-1_V2.1.gene_lncRNA.gtf -o temp.gtf V007-V008-1.sam.sorted.bam

ls *fpkm.txt|while read id; do perl -ne 'print unless /^STRG/' $id > 1 && mv 1 $id;done

#### step5: gtf2bed

perl -ne 'print if /\ttranscript\t/' lncRNA_mRNA.gtf|perl -pi -e 's/gene_id.*?transcript_id "//'|perl -pi -e 's/";//'|awk '{print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' > Ga_lncRNA_mRNA.bed



