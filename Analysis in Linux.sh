
# STAR for sequencing data mapping
STAR --runThreadN 18 --genomeDir $index --readFilesCommand zcat --readFilesIn $fq1 $fq2 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --quantMode TranscriptomeSAM GeneCounts --runMode alignReads --twopassMode Basic --outSAMunmapped None
# "$index" : Index built with STAR
# "$fq1 $fq2" : R1 and R2 of double ended sequencing data

# Transcriptional quantification using kallisto
kallisto quant -i ${index} -o ${sample} --bias $fq1 $fq2
# "${index}" : Index built with kallisto
# "${sample}" : Output folder name
# "$fq1 $fq2" : R1 and R2 of double ended sequencing data