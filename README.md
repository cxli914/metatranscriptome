# An Instruction on the Analysis of Example Meta-transcriptomics Datasets
## Tools required
```shell
Trim
Python
Bowtie2 
featureCounts
```
## Data preparation
Please execute the following command in your terminal to clone the metatranscriptome repository to your own machine. 
```shell
git clone https://github.com/cxli914/metatranscriptome.git
```
Then, download the example meta-transcriptomics data and the corresponding parameter files such as rRNA database, reference database. Please store them in the `data` folder.

## Data analysis
Execute the following command in your terminal to start the analysis of meta-transcriptomics data:
```shell
nextflow run metatranscriptome.nf --stdin1 /path/to/1.fastq.gz --stdin2 /path/to/2.fastq.gz --rrna_db /path/to/rrnadb_index/rrna_db --ref /path/to/reference.fa --gtf /path/to/example.gtf --kraken2_db /path/to/kraken2_db
```
**Note:** This step will take about **two hours**. The output results will be stored in the folder named `results` by default. Please refer to the **Help Message** section or execute `nextflow run metatranscriptome.nf --help` in the terminal to view the detailed information of parameter passing.

# Help Message
This help message can also be obtained by executing the following command in the terminal：
```shell
nextflow run metatranscriptome.nf --help
```
## Command line for meta-transcriptomics data analysis:
```
nextflow run metatranscriptome.nf --stdin1 <reads 1> --stdin2 <reads 2> --rrna_db <rrna database> --ref <reference database> --gtf <GTF file> --kraken2_db <kraken2 databse> <Options> <Functions> ```
## Library-based mode: 
```
## Parameters descriptions
### Required arguments
|parameters|descriptions|
|---|---|
|--stdin1|Primary reads input, fastq file. For example: --stdin1 /path/to/1.fastq.gz.|
|--stdin2|For paired reads in two files, fastq file. For example: --stdin2 /path/to/2.fastq.gz.|
|--rrna_db|rrna database for the rrna removal in sequence reads. Note that this parameter needs the filename prefix of rrna database (minus trailing .X.bt2). For example, --rrna_db /path/to/rrnadb_index/rrna_db. See bowtie2 -x for details.|
|--ref|Specify a reference sequence for reads mapping, fasta file. For example, --ref /path/to/reference.fa.|
|--gtf|Specify a GTF file as the input for featureCounts. For example, --gtf /path/to/example.gtf|
|--kraken2_db|kraken2 database for kraken2. Note that this parameter needs a folder that contains kraken2 database. For example, --kraken2_db /path/to/kraken2_db. See kraken2 -db for details.|

### Options arguments
|parameters|descriptions|
|---|---|
|--sample|Specify the prefix of output files. Default: the prefix of the name of stdin1 (exclude suffix).|
|--outdir|Specify an output folder. Default: the results folder in the current path where the command is executed.|
|--thread| Set the thread. Default: 40.|
|--cache| The mode to store the process results to a local cache, see [nextflow document](https://www.nextflow.io/docs/latest/basic.html) for details. Default: "deep".|
|--publish_mode|The publication mode of output files, see [nextflow document](https://www.nextflow.io/docs/latest/basic.html) for details. Default: "copy".|
|--help|Print this help message.|

### Functions arguments
These parameters are built-in functions of Nextflow, they can generate some visual graphics, which or show the total time consumption of the pipeline, or show the time consumption, memory occupation, cpu usage of each process. Interested can add these parameters to observe relative information.
|parameters|descriptions|
|---|---|
|-with-timeline|It renders a timeline.html file that records the time, memory consumption of different processes.|
|-with-report|It generates a report.html file that records the single core CPU Usage, execution time, memory occupation and Disk read write information of different processes.|
|-with-trace|It creates an execution tracing file that contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used.|
|-with-dag|It outputs the pipeline execution DAG. It creates a file named dag.dot containing a textual representation of the pipeline execution graph in the DOT format.|
|-resume|It means only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead. Also, the pipeline can be restarted by add the parameter when any disconnection of the network or server occurs.|
