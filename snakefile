#this snakemake pipeline is used to analyse the performance of isONform,
#especially on bigger datasets (>1k reads)

#to activate snakemake environment: conda activate snake
#run via: snakemake --cores 1
#original position 100kSIRV6/0
#SAMPLES = list(range(0,25)
#SAMPLES =["0","1","2"]
#wildcard_constraints: id ="\d+",
IDS, = glob_wildcards("spoa{id}merged.fa")
#max_id = max(IDS)

rule all:
    input:
        #expand("alignment_{sample}_spoa.sam",sample=IDS)
        #expand("analysis{sample}spoa.csv",sample=IDS),expand("nroccs_{sample}.txt",sample=IDS)

#rule all:
#    input:
#        expand("alignment_0_spoa.sam",sample=samples)"""

rule align_spoa:
    input: "cluster0Original_Transcriptome.fa","spoa{sample}merged.fa"
    output: "alignment_{sample}_spoa.sam"
    shell: "minimap2 -ax map-ont {input} > {output}"

rule count_spoa_instances:
    input: db = "cluster0Original_Transcriptome.fa",
            map ="mapping{sample}.txt",
            align="alignment_{sample}_spoa.sam"
    output:"analysis{sample}spoa.csv"
    shell:"python get_error_rates_mapping.py {input.db} {input.map} {input.align} {output}"

rule align_batch:
    input: "cluster0Original_Transcriptome.fa","{sample}_batchfile.fa"
    output: "alignment{sample}batch.sam"
    shell: "minimap2 -ax map-ont  {input} > {output}"

rule error_rates_batch:
    input:"cluster0Original_Transcriptome.fa","alignment{sample}batch.sam"
    output:"analysis{sample}_batch.csv"
    shell:"python get_error_rates_original.py {input} {output}"

rule count_id_abundances_batch:
    input:"analysis{sample}_batch.csv"
    output: "nroccs_{sample}.txt"
    shell:"./count_SIRV_id_abundances.sh  {input} {output}"

