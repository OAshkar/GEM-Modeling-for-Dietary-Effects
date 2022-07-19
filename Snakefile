################################
# Adapted from this script:
# https://github.com/bartongroup/MG_GlycoTreg/blob/bfbf5282af38ae4190ff3ea55269d1546f30e3ad/Snakefile
################################

import re
from snakemake.io import *
import snakemake
import glob

# Genome URLs: ensembl
genomeURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa.gz"
gtfURL = "http://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf"
cdnaURL = "http://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa"
cdsURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz"


# Directories
fastqDir = "fastq/"
genomeDir = "genome/"
starIndexDir = "starindex/"
salmonIndexDir = "salmonindex/"

# Genome files
genomeFile = genomeDir + "Mus_musculus.GRCm38.dna_rm.primary_assembly.fa"
gtfFile = genomeDir + "Mus_musculus.GRCm38.93.gtf"
cdsFile = genomeDir + "Mus_musculus.GRCm38.cds.all.fa"
transcriptomeFile = genomeDir + "Mus_musculus.GRCm38.cdna.all.fa"
genomeIndex = genomeFile + ".fai"
genomeSizeFile = genomeFile + ".txt"


# Files to test if an index was created
starIndexTestFile = starIndexDir + "chrName.txt"
salmonIndexTestFile = salmonIndexDir + "indexing.log"

DIR = "/home/omar/Downloads/Data/ObOb/rawdata/"
(SAMPLES, LANES, READS) = snakemake.io.glob_wildcards(DIR+ "fastq/{sample}_{lane}_{read}_001.fastq.gz")

# SAVE computation and time
SAMPLES = list(dict.fromkeys(SAMPLES))
LANES = list(dict.fromkeys(LANES))
READS = list(dict.fromkeys(READS))


# Standardize filenames
SAMPLESsum  = list(set([re.sub("(ob_ob|WT)_(ND|HFD|HDF)_\d_(\D\D).*", "\\1_\\2_\\3",x) for x in SAMPLES]))
SAMPLESsum = list(set([re.sub("HDF", "HFD" ,x ) for x in SAMPLESsum]))
#print("FASTQ Number %d" % len(SAMPLES))
#print("Conditions Number %d" % len(SAMPLESsum))

# Lists for rule all
QCS1 = expand("qc/qc1/{sample}_{lane}_{read}_001_fastqc.html", sample=SAMPLES, lane = LANES, read = READS)
QCS2 = expand("qc/qc2/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=READS)
QCS3 = expand("qc/qc3/{sample}_{pair}_val_{num}_fastqc.html", sample=SAMPLES, pair=READS, num = [1,2])
MULTIQC1 = ["multiqc/multiqc1.html"]
MULTIQC2 = ["multiqc/multiqc2.html"]
trim = [expand(fastqDir + "trimmed/{sample}.R1_val_1.fq.gz", sample = SAMPLES)]
MULTIQC3 = ["multiqc/multiqc3.html"]

QUANTS = expand("salmon/{sample}/quant.sf", sample=SAMPLES)
COUNTS = ["counts/TPMaverage.csv", "counts/Count.Rdata"]
medianCounts = "counts/medianCounts.csv"
DGE_re = "paper/figure/Anno.8Genes.png"
mouse_char = ["paper/tex/MiceCharactersticsTable.tex", "paper/tex/summaryChar.tex"]

modelNames = snakemake.io.glob_wildcards("models/fastc_median_plus/model_{modelNames}.mat") #FIXME I need to set a real names
modelNames = modelNames.modelNames
modelNamesUnctrl = [x for x in modelNames if "WT_ND" not in x]
rule final: # requires all outputs as input
     input:
         QCS1,
         MULTIQC1,
         QCS2,
         MULTIQC2,
         trim,
         # MULTIQC3,
         QUANTS,
         COUNTS,
         medianCounts,
         DGE_re,
         mouse_char,
         expand("models/fastc_median/model_{X}.mat", X = modelNames),
         samples = expand("samplesMedianRes2/model_{X}.csv", X = modelNames),
         tidysamples = "samplesAllTidy_ACHR2.Rdata",
         xgb_tune = "tmp/xgb_res.RDS" ,
         xgb_final = "tmp/final_fit_xgb.RDS",
         rmta =  expand("rMTAResults/{X}/rmtares.csv", X = modelNamesUnctrl)



## fastqc1
####################################################################
# Quality control
rule fastqc1:
    input:
        R1 = DIR+"fastq/{sample}_R1_001.fastq.gz",
        R2 = DIR+"fastq/{sample}_R2_001.fastq.gz"
    output:
        #QCS1
        R1 = "qc/qc1/{sample}_R1_001_fastqc.html",
        R2 = "qc/qc1/{sample}_R2_001_fastqc.html"
    threads: 2
    shell:
        "fastqc -o qc/qc1 --threads {threads} -f fastq {input.R1} {input.R2}"

rule multiqc1:
    input: QCS1
    output: MULTIQC1
    shell:
        """
        #mkdir -p multiqc/multiqc1
        multiqc -f --filename multiqc1 --outdir multiqc qc/qc1 -i "All samples before trimming or merging"
        """

####################################################################
# Merge lanes
rule combine_lanes:
    # Add file one by one to the output target
    output:
        r1 = "fastq/{sample}_R1.fastq.gz",
        r2 = "fastq/{sample}_R2.fastq.gz"
    shell:
        """
        cat {DIR}/fastq/{wildcards.sample}_L*_R1_001.fastq.gz > {output.r1}
        cat {DIR}/fastq/{wildcards.sample}_L*_R2_001.fastq.gz > {output.r2}
        """
####################################################################
# Quality control
rule fastqc2:
    input:
        R1 = "fastq/{sample}_R1.fastq.gz",
        R2 = "fastq/{sample}_R2.fastq.gz"
    output:
        "qc/qc2/{sample}_R1_fastqc.html",
        "qc/qc2/{sample}_R2_fastqc.html"
    threads: 2
    shell:
        "fastqc -o qc/qc2 --threads {threads} -f fastq {input.R1} {input.R2}"


rule multiqc2:
    input: QCS2
    output: MULTIQC2
    shell:
        """
        #mkdir -p multiqc/multiqc2
        multiqc -f --filename multiqc2 --outdir multiqc qc/qc2 -i "merged data before trimming"
        """

####################################################################
# trimming

rule trim_galore_pe:
    input:
        fastqDir + "{sample}_R1.fastq.gz",
        fastqDir + "{sample}_R2.fastq.gz"
    output:
        fastqDir + "trimmed/{sample}.R1_val_1.fq.gz",
        fastqDir + "trimmed/{sample}.R1.fastq.gz_trimming_report.txt",
        fastqDir + "trimmed/{sample}.R2_val_2.fq.gz",
        fastqDir + "trimmed/{sample}.R2.fastq.gz_trimming_report.txt",
    log:
        "logs/trim_galore/{sample}.log",
    shell:
         """
         ~/NoApps/TrimGalore-0.6.6/trim_galore {input} -o "fastq/trimmed" --paired --fastqc --fastqc_args "-o qc/qc3"
         """

# rule trim:
#     input:
#         ["fastq/raw/{sample}.{lane}.R1.fastq.gz",
#             "fastq/raw/{sample}.{lane}.R2.fastq.gz"]
#     output:
#         fastq1=temp("fastq/trimmed/{sample}.{lane}.R1.fastq.gz"),
#         fastq2=temp("fastq/trimmed/{sample}.{lane}.R2.fastq.gz"),
#         qc="qc/cutadapt/{sample}.{lane}.txt"
    # params:
    #     extra='--fastqc --fastqc_args "--outdir qc/qc3',
    #
    # wrapper:
    #     "0.80.3/bio/trim_galore/pe"
# rule multiqc3:
#     input: QCS3
#     output: MULTIQC3
#     shell:
#         """
#         #mkdir -p multiqc
#         multiqc -f --multiqc3 report --outdir multiqc qc/qc3
#        """


####################################################################
# Load genome files

rule load_db_files:
    output: genomeFile, gtfFile, transcriptomeFile, cdsFile
    shell:
        """
        wget {genomeURL} -O - | gunzip -c > {genomeFile}
        wget {gtfURL} -O - | gunzip -c > {gtfFile}
        wget {cdnaURL} -O - | gunzip -c > {transcriptomeFile}
        wget {cdsURL} -O - | gunzip -c > {cdsFile}
        """

####################################################################
# Salmon

rule salmon_index:
    input: transcriptomeFile
    output: salmonIndexTestFile
    log: "logs/salmon_index.log"
    shell:
        "salmon index -t {input} -i {salmonIndexDir} &> {log}"


rule salmon_quant:
    input:
        R1 = fastqDir + "trimmed/{sample}.R1_val_1.fq.gz",
        R2 = fastqDir + "trimmed/{sample}.R2_val_2.fq.gz",
        #R1 = fastqDir + "{sample}_R1.fastq.gz",
        #R2 = fastqDir + "{sample}_R2.fastq.gz",
        testfile = salmonIndexTestFile
    output: "salmon/{sample}/quant.sf"
    params:
        prefix = "salmon/{sample}"
    threads: 12
    log: "logs/salmon_{sample}_quant.log"
    shell:
        """
        salmon quant \
        --index {salmonIndexDir} \
        --libType A \
        --numBootstraps 100 \
        --threads {threads} \
        -1 {input.R1} -2 {input.R2} \
        --output {params.prefix} &> {log}
        """

rule extract_TPM:
    input: QUANTS
    output: COUNTS, medianCounts
    shell:
        """
        Rscript tx2median.R --genome {gtfFile} --directory salmon --out {COUNTS}
        """

rule DGE:
    input: COUNTS
    output: DGE_re
    shell:
        """
        Rscript DGE.R
        """
rule MouseChar:
    output: mouse_char
    shell:
        """
        Rscript char.R
        """

rule fastc_median:
    input:
        medianCounts
    # log:
    #     "logs/fastc.log"
    output:
        expand("models/fastc_median/model_{X}.mat", X = modelNames)
    shell:
        """
        /home/omar/NoApps/MATLAB/R2019b/bin/matlab -nosplash -nodesktop -r "run('rFastcromicsMedian.m');exit;"
        """

rule sampling: # this step will also write new models with fixed constarints
    input:
        "models/fastc_median/{X}.mat"
    output:
        models = "models/fastc_median2/model_{X}.mat",
        samples = "samplesMedianRes2/model_{X}.csv"
    shell:
        """
        python3 simulations2.py "models/fastc_median/{wildcards.X}"
        """

rule tidy_samples:
    input:
        rules.final.input.samples
    output:
        rules.final.input.tidysamples
    shell:
        """
        Rscript sim_tidying.R
        """

# rule dflux_cluster:

# rule dflux_wilcox:
#     input:
#          rules.sampling.output.samples


rule dflux_xgb_tuning:
   threads: 10
   resources:
       mem_mb=25000
   input:
       rules.final.input.tidysamples
   output:
       rules.final.input.xgb_tune
   shell:
       """
       Rscript dflux_xgb_prep.R
       """

rule dflux_xgb_final:
    threads: 10
    resources:
        mem_mb=30000
    input:
        rules.dflux_xgb_tuning.output
    output:
        rules.final.input.xgb_final
    shell:
        """
        Rscript dflux_xgb_final_fit.R
        """

rule rMTA:
    threads: 10 # give 6 cores for every job . It already takes almost 100% for CPU for single run!try to ignore complete also
    resources:
        mem_mb=25000
    log:
        "rMTAResults/{X}/rmta.log"
    input:
        "models/fastc_median2/model_{X}.mat"
    output:
        "rMTAResults/{X}/rmtares.csv"
    shell:
        """
        /home/omar/NoApps/MATLAB/R2019b/bin/matlab -nosplash -nodesktop -r "rmta('{wildcards.X}');exit"
        """
