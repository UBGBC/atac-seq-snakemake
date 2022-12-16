#!/usr/bin/py
import pandas as pd
import os

# trim -> align -> map -> index -> remove mito -> index -> shift reads -> peak call -> bigWig -> filter peaks -> motif analysis

configfile: "config.json"
localrules: all, mkdir

wildcard_constraints:
    sample_id = '\w+'

df = pd.read_csv(config["meta_file"], sep='\t', header=0, index_col=0)
sample_ids = list(df.index)

def get_pair_gz(sample_id):
    dir = config["raw_fastq_gz_dir"]
    return tuple(os.path.join(dir, df.loc[str(sample_id), x]) for x in ('ForwardFastqGZ', 'ReverseFastqGZ'))

def get_forward_primer(sample_id):
    return df.loc[sample_id]["Adapter_1"]

def get_reverse_primer(sample_id):
    return df.loc[sample_id]["Adapter_2"]

###

rule all:
    input:expand("{dir}/{sample_id}.sort.shifted.bam.bai", dir=config["dir_names"]["shifted_dir"],sample_id=sample_ids)
    run:
        for sample in sample_ids:
            print("Wrapping up pipeline")

rule mkdir:
    output: touch(config["file_names"]["mkdir_done"])
    params: dirs = list(config["dir_names"].values())
    resources: time_min=10, mem_mb=2000, cpus=1
    shell: "mkdir -p {params.dirs}"

rule trim:
    input: 
        rules.mkdir.output,
        all_read1 = lambda wildcards: get_pair_gz(wildcards.sample_id)[0],
        all_read2 = lambda wildcards: get_pair_gz(wildcards.sample_id)[1]
    resources: time_min=360, mem_mb=2000, cpus=1
    output: 
        trimmed_read1 = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R1.fastq.gz",
        trimmed_read2 = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R2.fastq.gz",
        trimmed_stats = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.stats"
    version: config["tool_version"]["cutadapt"]
    params:
        adapter1=lambda wildcards: get_forward_primer(wildcards.sample_id),
        adapter2=lambda wildcards: get_reverse_primer(wildcards.sample_id)
    shell: "cutadapt -m 15 -a {params.adapter1} -A {params.adapter2} -n 2 -o {output.trimmed_read1} -p {output.trimmed_read2} {input.all_read1} {input.all_read2} >{output.trimmed_stats}"

rule map:
    input:
        p1 = rules.trim.output.trimmed_read1,
        p2 = rules.trim.output.trimmed_read2
    resources: time_min=360, mem_mb=20000, cpus=6
    output:
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam",
        stats = config["dir_names"]["mapped_dir"] + "/{sample_id}.stats",
    params:
        threads = config["params"]["bowtie2"]["threads"],
        map_all = config["params"]["bowtie2"]["all"],
        reference = config["params"]["bowtie2"]["bowtie2_reference"]
    shell:
        """
        bowtie2 --threads {params.threads} -1 {input.p1} -2 {input.p2} -x {params.reference} 2> {output.stats} | samtools sort -T {output.mapped_bam_file}.tmp -O bam -o {output.mapped_bam_file}
        """

rule index:
    input:
        mapped_bam_file = rules.map.output.mapped_bam_file
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        indexed_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam.bai"
    shell:
        """
        samtools index {input.mapped_bam_file} {output.indexed_bam_file}
        """

rule rmmito:
    input:
        rules.index.output.indexed_bam_file,
        mapped_bam_file = rules.map.output.mapped_bam_file
    resources: time_min=360, mem_mb=2000, cpus=6
    output:
        nomito_mapped_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.nomito.bam"
    params:
        reference = config["params"]["samtools"]["samtoolsview_ref"]
    shell:
        """
        samtools view -h {input.mapped_bam_file} -L {params.reference} -b | samtools sort -T {output.nomito_mapped_bam_file}.tmp -O bam -o {output.nomito_mapped_bam_file}
        """

rule remove_dup:
    input:
        nomito_mapped_bam_file = rules.rmmito.output.nomito_mapped_bam_file
    resources: time_min=360, mem_mb=12000, cpus=6
    output:
        marked_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.marked.bam",
        dedup_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.dedup.bam",
        dedup_stats = config["dir_names"]["filtered_dir"] + "/{sample_id}.dedup.metrics.txt"
    params:
        picard_path = config["params"]["picard"]["picard_path"]
    shell:
        """
        java -jar {params.picard_path}/picard.jar MarkDuplicates INPUT={input.nomito_mapped_bam_file} OUTPUT={output.marked_bam_file} METRICS_FILE={output.dedup_stats};
        samtools view -bh -F 1024 {output.marked_bam_file} > {output.dedup_bam_file};
        """

rule index2:
    input:
        dedup_bam_file = rules.remove_dup.output.dedup_bam_file
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        indexed_dedup_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.dedup.bam.bai"
    shell:
        """
        samtools index {input.dedup_bam_file} {output.indexed_dedup_file}
        """

rule shift:
    input:
        rules.index2.output.indexed_dedup_file,
        dedup_bam_file = rules.remove_dup.output.dedup_bam_file
    resources: time_min=360, mem_mb=2000, cpus=6
    output:
        shifted_bam_file = config["dir_names"]["shifted_dir"] + "/{sample_id}.shifted.bam",
        metrics = config["dir_names"]["shifted_dir"] + "/{sample_id}.shiftmetrics.txt"
    shell:
        """
        alignmentSieve -b {input.dedup_bam_file} --ATACshift --minMappingQuality 10 --filterMetrics {output.metrics} -o {output.shifted_bam_file}
        """

rule sort_shifted:
    input:
        shifted_bam_file = rules.shift.output.shifted_bam_file
    resources: time_min=360, mem_mb=2000, cpus=6
    output:
        sorted_shifted_bam_file = config["dir_names"]["shifted_dir"] + "/{sample_id}.sort.shifted.bam"
    shell:
        """
        samtools sort {input.shifted_bam_file} -O bam -o {output.sorted_shifted_bam_file}
        """

rule index_shifted:
    input:
        sorted_shifted_bam_file = rules.sort_shifted.output.sorted_shifted_bam_file
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        indexed_shifted_file = config["dir_names"]["shifted_dir"] + "/{sample_id}.sort.shifted.bam.bai"
    shell:
        """
        samtools index {input.sorted_shifted_bam_file} {output.indexed_shifted_file}
        """
