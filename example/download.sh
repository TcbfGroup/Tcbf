hg19=http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gtf1=http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.ncbiRefSeq.gtf.gz
rheMac8=https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.gz
gtf2=https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/genes/rheMac8.ncbiRefSeq.gtf.gz
mm10=http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz
gtf3=http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
for line in $hg19 $gtf1 $rheMac8 $gtf2 $mm10 $gtf3;do
   wget $line;
  done

pigz -d *gz

for line in *gtf;do
  gffread -E $line -o-> $line.gff3;
done
