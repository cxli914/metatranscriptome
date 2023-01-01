#! /usr/bin/env nextflow

/*
=====================================================================================
                              metatranscriptome.nf                                      
=====================================================================================
    A standardized reference-based meta-transcriptomics data processing pipeline.                                     
=====================================================================================
*/

//

params.help = false
if (params.help) {
    HelpMessage()
    exit 0
}

def HelpMessage() {
    log.info """
    Please run the command: [ nextflow run metatranscriptome.nf --help ] to print the help message.
    
    Usage:
    
    nextflow run metatranscriptome.nf --stdin1 <reads 1> --stdin2 <reads 2> --rrna_db <rrna database> --ref <reference database> --gtf <GTF file> --kraken2_db <kraken2 databse> <Options> <Functions>
        
    ====================
    Required arguments: 
    ====================

    --stdin1         Primary reads input, fastq file. For example: --stdin1 /path/to/1.fastq.gz.

    --stdin2         For paired reads in two files, fastq file. For example: --stdin2 /path/to/2.fastq.gz.

    --rrna_db        rrna database for the rrna removal in sequence reads. Note that this parameter needs the filename prefix of rrna database (minus trailing .X.bt2). For example, --rrna_db /path/to/rrnadb_index/rrna_db. See bowtie2 -x for details.

    --ref            Specify a reference sequence for reads mapping, fasta file. For example, --ref /path/to/reference.fa.

    --gtf            Specify a GTF file as the input for featureCounts. For example, --gtf /path/to/example.gtf

    --kraken2_db     kraken2 database for kraken2. Note that this parameter needs a folder that contains kraken2 database. For example, --kraken2_db /path/to/kraken2_db. See kraken2 -db for details. 
    
    ====================
    Options arguments:
    ====================
    
    --sample        Specify the prefix of output files. Default: the prefix of the name of stdin1 (exclude suffix).

    --outdir        Specify an output folder. Default: the results folder in the current path where the command is executed.

    --thread        Set the thread. Default: 40.

    --cache         The mode to store the process results to a local cache, see https://www.nextflow.io/docs/latest/basic.html for details. Default: "deep".
    
    --publish_mode  The publication mode of output files, see https://www.nextflow.io/docs/latest/basic.html for details. Default: "copy".

    --help          Print this help message.

    ====================
    Functions arguments:
    ====================
        
    ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    These parameters are built-in functions of Nextflow, they can generate some visual graphics, which or show the total time consumption of the pipeline, or show the time consumption, memory occupation, cpu usage of each process. 
    Interested can add these parameters to observe relative information. See https://www.nextflow.io/docs/latest/basic.html for details.
    ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    -with-timeline  It renders a timeline.html file that records the time, memory consumption of different processes. 

    -with-report    It generates a report.html file that records the single core CPU Usage, execution time, memory occupation and Disk read write information of different processes.

    -with-trace     It creates an execution tracing file that contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used.

    -with-dag       It outputs the pipeline execution DAG. It creates a file named dag.dot containing a textual representation of the pipeline execution graph in the DOT format.
        
    -resume         It means only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

    """.stripIndent()
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    SET UP CONFIGURATION VARIABLES                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////


params.sample = false
params.stdin1 = null
params.stdin2 = null
params.outdir = "./results"


params.rrna_db = null
params.ref = null 
params.gtf = null
params.kraken2_db = null

params.thread = 40
params.cache = "deep"
params.publish_mode = "copy"

//When the command line lacks the corresponding parameter, give the necessary warning
def Warnings() {
    if (params.stdin1 == null) {
        println "\nError:\n     The standard input file 1 is not found, please specify it by adding the --stdin1 parameter in the command-line!\n"
        exit 0
    }
    if (params.stdin2 == null) {
        println "\nError:\n     The standard input file 2 is not found, please specify it by adding the --stdin2 parameter in the command-line!\n"
        exit 0
    }
    if (params.rrna_db == null) {
        println "\nError:\n     The rrna database is not found, please specify it by adding the --rrna_db parameter in the command-line!\n"
        exit 0
    }
    if (params.ref == null) {
        println "\nError:\n     The reference database is not found, please specify it by adding the --ref parameter in the command-line!\n"
        exit 0
    }
    if (params.gtf == null) {
        println "\nError:\n     The GTF file is not found, please specify it by adding the --gtf parameter in the command-line!\n"
        exit 0
    }
    if (params.kraken2_db == null) {
        println "\nError:\n     The database for kraken2 is not found, please specify it by adding the --kraken2_db parameter in the command-line!\n"
        exit 0
    }
}
Warnings()

def Setsample() {
    if (params.sample == false) {
        infile = file("$params.stdin1")
        sample = infile.simpleName
    }
    else if (params.sample) {
        sample = params.sample
    }
}
Setsample()

def Makedirs() {
    outdir = file("$params.outdir")
    if (outdir.exists() == false) {
        outdir.mkdir()
    }

    qc_outdir = file("${params.outdir}/01_trim_qc")
    qc_outdir.mkdir()

    fastqc_outdir = file("${params.outdir}/02_fastqc_qe")
    fastqc_outdir.mkdir()
    
    rrnaRemove_outdir = file("${params.outdir}/03_bowtie2_rrnaRemove")
    rrnaRemove_outdir.mkdir()
    
    readsMapping_outdir = file("${params.outdir}/04_bbmap_map")
    readsMapping_outdir.mkdir()
    
    featureCounts_outdir = file("${params.outdir}/05_featureCounts_tpm")
    featureCounts_outdir.mkdir()

    kraken2_outdir = file("${params.outdir}/06_kraken2")
    kraken2_outdir.mkdir()
}
Makedirs()

Channel.fromPath(params.stdin1, checkIfExists: true)
       .ifEmpty{exit 1, "The $params.stdin1 file is empty!"}
       .view{"$workflow.start - INFO - Load standard input file 1: $it"}
       .set{stdin1}

Channel.fromPath(params.stdin2, checkIfExists: true)
       .ifEmpty{exit 1, "The $params.stdin2 file is empty!"}
       .view{"$workflow.start - INFO - Load standard input file 2: $it"}
       .set{stdin2}


// quality control for raw squencing data. Tools: Trimmomatic

process qualityControl {

    cache params.cache
    publishDir "$qc_outdir", mode: params.publish_mode

    input:
    path r1 from stdin1
    path r2 from stdin2

    output:
    path "*"
    path "${sample}_01_trim_1.fastq.gz" into trim_fq1
    path "${sample}_01_trim_2.fastq.gz" into trim_fq2

    script:
    """
    c1=${sample}_01_trim_1.fastq.gz
    c2=${sample}_01_trim_2.fastq.gz
    t1=${sample}_01_trim_unpaired_1.fastq.gz
    t2=${sample}_01_trim_unpaired_2.fastq.gz

    java -jar /lustre/home/acct-ioozy/ioozy-user3/FQL/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $r1 $r2 \${c1} \${t1} \${c2} \${t2} ILLUMINACLIP:/lustre/home/acct-ioozy/ioozy-user3/FQL/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 HEADCROP:10
    """
}

// copy quality-controlled fastq file into multiple channels
trim_fq1.into{trim_fq1_01; trim_fq1_02}
trim_fq2.into{trim_fq2_01; trim_fq2_02}

// quality evaluation of quality-controlled fastq file
// tools: fastqc

process qualityEvaluate {

    cache params.cache
    publishDir "$fastqc_outdir", mode: params.publish_mode

    input:
    path c1 from trim_fq1_01
    path c2 from trim_fq2_01

    output:
    path "*"

    script:
    """
    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/FastQC/fastqc -f fastq $c1 -o ./ -t $params.thread 
    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/FastQC/fastqc -f fastq $c2 -o ./ -t $params.thread
    """
}

// rrna removal of quality-controlled meta-transcriptomics fastq file
// tools: bowtie2

process rrnaRemove {

    cache params.cache
    publishDir "$rrnaRemove_outdir", mode: params.publish_mode

    input:
    path c1 from trim_fq1_02
    path c2 from trim_fq2_02

    output:
    path "*"
    path "${sample}_filter_1.fastq.gz" into filter_fq1
    path "${sample}_filter_2.fastq.gz" into filter_fq2

    script:
    """
    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 -x $params.rrna_db -1 $c1 -2 $c2 --no-head --no-unal --no-mixed -p $params.thread --un-conc-gz ${sample}_filter.fastq.gz --al-conc-gz ${sample}_rRNA.fastq.gz -S ${sample}_rRNA.sam

    rm ${sample}_rRNA.sam

    mv ${sample}_filter.fastq.1.gz ${sample}_filter_1.fastq.gz
    mv ${sample}_filter.fastq.2.gz ${sample}_filter_2.fastq.gz
    """
}

// copy quality-controlled and rrna-removed fastq file into multiple channels
filter_fq1.into{filter_fq1_01; filter_fq1_02}
filter_fq2.into{filter_fq2_01; filter_fq2_02}

// reads mapping
// tools: script - bbmap.sh & samtools

process readsMapping {

    cache params.cache
    publishDir "$readsMapping_outdir", mode: params.publish_mode

    input:
    path f1 from filter_fq1_01
    path f2 from filter_fq2_01

    output:
    path "*"
    path "${sample}_04_bbmap.sorted.bam" into mapped_bam

    script:
    """
    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/bbmap/bbmap.sh in=$f1 in2=$f2 ref=$params.ref nodisk threads=$params.thread covstats=${sample}_depth_bbmap.depth out=${sample}_04_bbmap_temp.sam

    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/samtools-1.11/samtools view -bS -S ${sample}_04_bbmap_temp.sam -o ${sample}_04_bbmap.bam
    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/samtools-1.11/samtools sort -@ $params.thread -O bam ${sample}_04_bbmap.bam -o ${sample}_04_bbmap.sorted.bam

    rm ${sample}_04_bbmap_temp.sam
    """
}

// calculate and obtain gene count matrix
// tools: featureCounts

process featureCounts {
    
    cache params.cache
    publishDir "$featureCounts_outdir", mode: params.publish_mode

    input:
    path bam from mapped_bam

    output:
    path "*"

    script:
    """
    /lustre/home/acct-ioozy/ioozy-user3/FQL/software/subread-2.0.1-Linux-x86_64/bin/featureCounts -a $params.gtf -p $bam -t transcript -g gene_id -T $params.thread -o ${sample}_bp_trim.txt
    python /lustre/home/acct-ioozy/ioozy-user3/FQL/SMTD/best_pipeline/06_tpm.py ${sample}_bp_trim.txt ${sample}_bp_trim_tpm.txt
    """
}

process kraken2 {

    cache params.cache
    publishDir "$kraken2_outdir", mode: params.publish_mode

    input:
    path f1 from filter_fq1_02
    path f2 from filter_fq2_02
    
    output:
    path "*"

    script:
    """
    mkdir trim

    /lustre/home/acct-ioozy/ioozy-user3/anaconda3/envs/kraken2/bin/kraken2 --db $params.kraken2_db --output ${sample}_06_kraken2_nr.out --report ${sample}_06_kraken2_nr.report --paired --gzip-compressed $f1 $f2 --threads $params.thread

    for list in `cat /lustre/home/acct-ioozy/ioozy-user3/FQL/SMTD/D_kraken2/sh/bac_20220922.txt`
    do
        type=`echo \$list | awk -F ":" '{print \$2}'`
        name=`echo \$list | awk -F ":" '{print \$1}'`
        cat ${sample}_06_kraken2_nr.report | awk -F "\\t" '{if(\$4=="'\$type'") print \$0}'| grep \$name | awk -F "\\t" '{print \$1}'>> ./trim/${sample}_summary_ratio.txt
        cat ${sample}_06_kraken2_nr.report | awk -F "\\t" '{if(\$4=="'\$type'") print \$0}'| grep \$name | awk -F "\\t" '{print \$2}'>> ./trim/${sample}_summary_count.txt
    done
    """
}
