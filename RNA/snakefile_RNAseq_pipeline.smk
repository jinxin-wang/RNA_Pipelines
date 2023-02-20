import os

## Author : Ismael Padioleau
## May 2018, revision august 2019 for flamingo cluster
## This pipeline will be used to produce RNA-seq analysis.
## The pipeline starts with reads files (FASTQ(.gz)), alignement with STAR (BAM)
## and quantification with QTLtools.

## To run  activate conda env : conda activate pipeline_RNA_seq
## /mnt/beegfs/software/conda/bin/snakemake -s Snakefile_RNA-seq_flamingo.py -c 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}mb -p {params.queue}' --configfile configuration_b37_RNASeq_flamingo.json -j 25 --latency-wait 120 -n
## Get DAG /mnt/beegfs/software/conda/bin/snakemake -s Snakefile_RNA-seq_flamingo.py -c 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}mb -p {params.queue}' --configfile configuration_b37_RNASeq_flamingo.json -j 25 --latency-wait 120 "MYGOAL" --dag   | dot -Tsvg > dag.svg
## Where MYGOAL is the file you want to produce

## Fisrt we define the list of file to produce
## Two possibilities :
##      - Run in a directory nammed fastq, with all fastq, in the working directory
##      or
##      - Specify files to produce, in a file called 'files_to_produce.tsv', with one file per line
##          The pipeline detect this file and run all steps requiered to produce asked outputs

## Enable DEBUG MODE
# sys.stdin = open('/dev/stdin')
# import pdb; pdb.set_trace()

## Single-end mapping must be named {SAMPLE_ID}_0.fastq.gz
## Pair-end mapping must be named {SAMPLE_ID}_1.fastq.gz and {SAMPLE_ID}_2.fastq.gz
SAMPLES = []
if os.path.isfile("files_to_produce.tsv") :
    SAMPLE_INPUT_LIST = open("files_to_produce.tsv",'r')
    for line in SAMPLE_INPUT_LIST :
        SAMPLES.append(line.strip())
else :
    SAMPLES,PAIRED_OR_SINGLE_END = glob_wildcards("fastq/{samples,.+}_{paired_or_single_end,[01]}.fastq.gz")

## Get all fastq to generate the QC
FASTQ_SAMPLES, = glob_wildcards("fastq/{name}.fastq.gz")

if os.path.isfile("files_to_produce.tsv"):
    rule all:
        input:
            expand(SAMPLES)
else:
    rule all:
            input : 
                expand("bam/{sample}.bam", sample = SAMPLES),
                expand('fastq_QC/{fastq_sample}_fastqc.html', fastq_sample=FASTQ_SAMPLES),
                expand('fastq_QC_clean/{fastq_sample}_fastqc.html', fastq_sample=FASTQ_SAMPLES),
                expand("mapping_QC/{sample}_flagstat.txt", sample=SAMPLES),
                # expand('QTLtools/{sample}.exon.count.bed', sample=SAMPLES),
                expand('HTSeq_geneNAME_count/{sample}_geneNAME_count.table', sample=SAMPLES),
                expand('HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.table', sample=SAMPLES),
                expand('HTSeq_geneID_count/{sample}_geneID_count.table', sample=SAMPLES),
                expand('HTSeq_transcriptID_count/{sample}_transcriptID_count.table', sample=SAMPLES),
                expand("rseqc_geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf", sample=SAMPLES),
                expand('rseqc_tin/{sample}.tin.xls', sample=SAMPLES),
                expand('rseqc_read_duplication/{sample}.DupRate_plot.pdf', sample=SAMPLES)
            
## A rule to generate fastq quality control
## [J. WANG] add -t {threads} in shell
##           raise threads 1 to 16
##           raise mem_mb 2000 to 51200
##           change mediumq to shortq
rule fastqc:
    input:
        fastq='fastq/{fastq_sample}.fastq.gz'
    output:
        'fastq_QC/{fastq_sample}_fastqc.html',
        'fastq_QC/{fastq_sample}_fastqc.zip'
    log:
        "logs/fastq_QC/{fastq_sample}_fastqc.html.log"
    threads : 16
    resources:
        mem_mb = 51200
    params:
        queue = "shortq",
        fastqc = config["APP_FASTQC"],
        adapters = config["FASTQC_ADAPTERS"]
    shell:
        '{params.fastqc} -t {threads} -a {params.adapters} -o fastq_QC/ {input.fastq} 2> {log}'

## A rule to clean fastq files with fastp SE
## [J. WANG] fastp uses up to 16 threads
rule fastp_SE:
    input:
        fastq_0="fastq/{sample}_0.fastq.gz",
    output:
        fastq_clean='fastq_clean/{sample}_0.fastq.gz',
        html_report='fastp_reports/{sample}_fastp_report.html',
        json_report='fastp_reports/{sample}_fastp_report.json'
    log:
        "logs/fastp/{sample}_fastp.html.log"
    params:
        queue = "mediumq",
        fastp = config["APP_FASTP"],
        adapters = config["FASTP_ADAPTERS"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        '{params.fastp} --thread {threads} --dont_overwrite -i {input.fastq_0} -o {output.fastq_clean} --compression 9 --adapter_fasta {params.adapters} --trim_poly_g --trim_poly_x --length_required 25 --overrepresentation_analysis  --html {output.html_report} --json {output.json_report} 2> {log}'
        
## A rule to clean fastq files with fastp PE
## [J. WANG] add -t {threads} in shell
##           raise threads 1 to 16
##           raise mem_mb 5000 to 51200
##           change longq to mediumq
rule fastp_PE:
    input:
        fastq_1="fastq/{sample}_1.fastq.gz",
        fastq_2="fastq/{sample}_2.fastq.gz"
    output:
        fastq_clean_1='fastq_clean/{sample}_1.fastq.gz',
        fastq_clean_2='fastq_clean/{sample}_2.fastq.gz',
        html_report='fastp_reports/{sample}_fastp_report.html',
        json_report='fastp_reports/{sample}_fastp_report.json'
    log:
        "logs/fastp/{sample}_fastp.html.log"
    params:
        queue = "mediumq",
        fastp = config["APP_FASTP"],
        adapters = config["FASTP_ADAPTERS"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        '{params.fastp} --thread {threads} --dont_overwrite -i {input.fastq_1} -o {output.fastq_clean_1} -I {input.fastq_2} -O {output.fastq_clean_2} --compression 9 --adapter_fasta {params.adapters} --trim_poly_g --trim_poly_x --length_required 25 --overrepresentation_analysis  --html {output.html_report} --json {output.json_report} 2> {log}'

## A rule to generate fastq quality control for cleaned fastq
## [J. WANG] add -t {threads} in shell
##           raise threads 1 to 16
##           raise mem_mb 2000 to 51200
##           change mediumq to shortq
rule fastqc_clean:
    input:
        fastq='fastq_clean/{fastq_sample}.fastq.gz'
    output:
        'fastq_QC_clean/{fastq_sample}_fastqc.html',
        'fastq_QC_clean/{fastq_sample}_fastqc.zip'
    log:
        "logs/fastq_QC_clean/{fastq_sample}_fastqc.html.log"
    threads : 16
    resources:
        mem_mb = 51200
    params:
        queue = "shortq",
        fastqc = config["APP_FASTQC"],
        adapters = config["FASTQC_ADAPTERS"]
    shell:
        '{params.fastqc} -t {threads} -a {params.adapters} -o fastq_QC_clean/ {input} 2> {log}'

        
## A rule to map single-end RNA sample using STAR
rule star_single:
    input:
        reads1 = "fastq_clean/{sample}_0.fastq.gz"
    output:
        bam = "bam/{sample}.bam"
    params:
        queue = "mediumq",
        star  = config["APP_STAR"],
        index = config["STAR_INDEX"],
        tmp_directory_for_index = "TMP_STAR/tmp_genome_for_star_{sample}/",
        genome_fasta = config["GENOME_FASTA"],
        star_sjdbOverhang = config["STAR_SJDBOVERHANG"],
        star_gtf = config["STAR_GTF"]
    log:
        "logs/bam/{sample}.bam.log"
    threads: 16
    resources:
        mem_mb = 51200
    shell:
        "mkdir -p {params.tmp_directory_for_index};"
        "{params.star}"
        " --readFilesCommand zcat"
        " --outStd Log"
        " --runThreadN {threads}"
        " --genomeDir {params.index}"
        " --readFilesIn {input.reads1}"
        " --sjdbOverhang {params.star_sjdbOverhang}"
        " --sjdbGTFfile {params.star_gtf}"
        " --twopassMode Basic"
        " --outSAMtype BAM SortedByCoordinate"
        " --outSAMunmapped Within"
        " --outSAMstrandField intronMotif"
        " --outSAMmultNmax 5"
        " --outSJfilterReads Unique"
        " --quantMode GeneCounts"
        " --outFileNamePrefix {params.tmp_directory_for_index}{wildcards.sample}_TMP_; "
        " mv {params.tmp_directory_for_index}{wildcards.sample}_TMP_Aligned.sortedByCoord.out.bam {output.bam} 2>> {log};"

## A rule to map paired-end RNA sample using STAR
## [J. WANG] ATTENTION :
##      if threads number is more than 12, then you need to add the command "ulimit -n 10000" 
##      at the begin of the shell script. 
rule star_paired:
    input:
        reads1 = "fastq_clean/{sample}_1.fastq.gz",
        reads2 = "fastq_clean/{sample}_2.fastq.gz"
    output:
        bam = "bam/{sample}.bam"
    params:
        queue = "mediumq",
        star  = config["APP_STAR"],
        index = config["STAR_INDEX"],
        tmp_directory_for_index = "TMP_STAR/tmp_genome_for_star_{sample}/",
        genome_fasta = config["GENOME_FASTA"],
        star_sjdbOverhang = config["STAR_SJDBOVERHANG"],
        star_gtf = config["STAR_GTF"],
        start_outSAMmultNmax = config["START_OUTSAMMULTNMAX"]
    log:
        "logs/bam/{sample}.bam.log"
    threads: 32
    resources:
        mem_mb = 102400
    shell:
        "ulimit -n 10000 2>> {log}; "
        "mkdir -p {params.tmp_directory_for_index} 2>> {log} ; "
        "{params.star}"
        " --readFilesCommand zcat"
        " --outStd Log"
        " --runThreadN {threads}"
        " --genomeDir {params.index}"
        " --readFilesIn {input.reads1} {input.reads2}"
        " --sjdbOverhang {params.star_sjdbOverhang}"
        " --sjdbGTFfile {params.star_gtf}"
        " --twopassMode Basic"
        " --outSAMtype BAM SortedByCoordinate"
        " --outSAMunmapped Within"
        " --outSJfilterReads Unique"
        " --quantMode GeneCounts"
        " --outSAMstrandField intronMotif"
        " --outSAMmultNmax 5"
        " --outFileNamePrefix {params.tmp_directory_for_index}{wildcards.sample}_TMP_ 2>> {log} ; "
        " mv {params.tmp_directory_for_index}{wildcards.sample}_TMP_Aligned.sortedByCoord.out.bam {output.bam} 2>> {log};"
      
## A rule to generate bam index with samtools
rule indexbam:
    input:
        bam = "bam/{sample}.bam"
    output:
        bai = "bam/{sample}.bam.bai"
    threads : 1
    resources:
        mem_mb = 2000
    params:
        queue = "shortq"
    log:
        "logs/bam/{sample}.bam.bai.log"
    shell:
        "samtools index {input} 2> {log}"

## A rule to check mapping metrics with samtools flagstat
rule samtools_flagstat:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai"
    output:
        "mapping_QC/{sample}_flagstat.txt"
    log:
        "logs/mapping_QC/{sample}.flagstat.log"
    threads : 1
    params:
        queue = "shortq"
    resources:
        mem_mb = 2000
    shell:
        "samtools flagstat {input.bam} > {output} 2> {log}"

        
if config["SPECIES"] == "Humain":
    # A rule to check mapping coverage  with rseqc
    rule rseqc_geneBody_coverage:
        input:
            bam = "bam/{sample}.bam",
            bai = "bam/{sample}.bam.bai",
            rseqc_ref = config["RSEQC_REF"]
        output:
            o1 = "rseqc_geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf"
            # o2 = "rseqc_geneBody_coverage/{sample}.geneBodyCoverage.r",
            # o3 = "rseqc_geneBody_coverage/{sample}.geneBodyCoverage.txt"
        log:
            "logs/rseqc_geneBody_coverage/{sample}.log"
        threads : 1
        params:
            queue = "mediumq",
            genbd = config["APP_GENE_BD_CVR"]
        resources:
            mem_mb = 5120
        shell:
            "cd rseqc_geneBody_coverage/ ;  {params.genbd} -r {input.rseqc_ref} -l 500 -i ../{input.bam} -o {wildcards.sample} >2 ../{log}"
    
    #A rule to check transcript integrity number with rseqc
    rule rseqc_tin:
        input:
            bam = "bam/{sample}.bam",
            bai = "bam/{sample}.bam.bai",
            rseqc_ref = config["RSEQC_REF"]
        output:
            o1 = "rseqc_tin/{sample}.tin.xls",
            o2 = "rseqc_tin/{sample}.summary.txt"
        log:
            "logs/rseqc_tin/{sample}.log"
        threads : 1
        params:
            queue = "mediumq",
            tin   = config["APP_TIN"]
        resources:
            mem_mb = 5000
        shell:
            "cd rseqc_tin/ ; {params.tin} -r {input.rseqc_ref} -c 30 -i ../{input.bam} 2> ../{log}"
        
    # A rule to check transcript integrity number with rseqc
    rule rseqc_readDuplication:
        input:
            bam = "bam/{sample}.bam",
            bai = "bam/{sample}.bam.bai",
        output:
            o1 ="rseqc_read_duplication/{sample}.DupRate_plot.pdf",
            o2 ="rseqc_read_duplication/{sample}.pos.DupRate.xls",
            o3 ="rseqc_read_duplication/{sample}.seq.DupRate.xls",
            o4 ="rseqc_read_duplication/{sample}.DupRate_plot.r"
        log:
            "logs/rseqc_readDuplication/{sample}.log"
        threads : 1
        params:
            queue = "mediumq",
            read_dup = config["APP_READ_DUP"]
        resources:
            mem_mb = 51200
        shell:
            "cd rseqc_read_duplication/; {params.read_dup} -i ../{input.bam} -o {wildcards.sample} 2> ../{log}"

## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_genesID:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        gene_id_count = "HTSeq_geneID_count/{sample}_geneID_count.table"
    log:
        "logs/HTSeq_geneID_count/{sample}_geneID_count.table"
    threads : 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20000
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i gene_id"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.gene_id_count} 2> {log}"

## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_transcriptID:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        transcript_id_count = "HTSeq_transcriptID_count/{sample}_transcriptID_count.table"
    log:
        "logs/HTSeq_transcriptID_count/{sample}_transcriptID_count.table"
    threads: 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i transcript_id"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.transcript_id_count} 2> {log}"
    
## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_genesNAME:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        gene_name_count = "HTSeq_geneNAME_count/{sample}_geneNAME_count.table"
    log:
        "logs/HTSeq_geneNAME_count/{sample}_geneNAME_count.table"
    threads : 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20000
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i gene_name"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.gene_name_count} 2> {log}"

## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_transcriptNAME:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        transcript_name_count = "HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.table"
    log:
        "logs/HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.table"
    threads : 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness  = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20000
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i transcript_name"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.transcript_name_count} 2> {log}"

# A rule to quantify reads per annotation and produce RPKM, with QTLtools
# rule quantification_with_QTLtools:
    # input:
        # bam = "bam/{sample}.bam",
        # bai = "bam/{sample}.bam.bai",
        # qtltool_gtf = config["QTLTOOL_GTF"]
    # output:
        # exon_count = "QTLtools/{sample}.exon.count.bed",
        # exon_rpkm = "QTLtools/{sample}.exon.rpkm.bed",
        # gene_count = "QTLtools/{sample}.gene.count.bed",
        # gene_rpkm = "QTLtools/{sample}.gene.rpkm.bed"
    # params:
        # out_prefix = "QTLtools/{sample}"
    # log:
        # "logs/QTLtools/{sample}.flagstat.log"
    # threads : 1
    # resources:
        # mem_mb = 20000
    # shell:
        # "QTLtools_1.0_CentOS6.8_x86_64 quan" 
        # " --bam {input.bam}"
        # " --gtf {input.qtltool_gtf}"
        # " --filter-mapping-quality 255"
        # " --filter-mismatch 5"
        # " --rpkm"
        # " --out {params.out_prefix} 2> {log}"
