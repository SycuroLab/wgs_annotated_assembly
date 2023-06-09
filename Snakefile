# ***************************************
# * Snakefile for wgs_annotated_assembly pipeline *
# ***************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
import os

os.environ["GTDBTK_DATA_PATH"] = "/bulk/IMCshared_bulk/shared/dbs/gtdbtk-1.5.0/db"

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# The forward read number for output files.
forward_read_num = config["forward_read_suffix"].split(".",1)[0]

# The reverse read number for output files.
reverse_read_num = config["reverse_read_suffix"].split(".",1)[0]


# **** Rules ****

rule all:
    input:
        os.path.join(config["output_dir"],"multiqc","multiqc_report_raw.html"),
        expand(os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"), sample=SAMPLES),
        os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html"),
        expand(config["output_dir"]+"/genomes/{sample}/assembly/spades/scaffolds.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/assembly/{sample}_genome.fa",sample=SAMPLES),        
        expand(config["output_dir"]+"/genomes/{sample}/quast/transposed_report.tsv", sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/prokka/{sample}.fna",sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/prokka/{sample}.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/prokka/{sample}.faa",sample=SAMPLES),
	expand(config["output_dir"]+"/genomes/{sample}/extracted_sequences/cpn60_metadata.csv",sample=SAMPLES),
##        expand(config["output_dir"]+"/{sample}/metaerg/data/all.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/checkm/checkm.tsv",sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/gtdbtk/gtdbtk.bac120.summary.tsv",sample=SAMPLES),
        expand(config["output_dir"]+"/genomes/{sample}/eggnog_mapper/{sample}.emapper.annotations",sample=SAMPLES),
#        os.path.join(config["output_dir"],"metadata_files","all_merged_assembly_analysis_metadata.tsv")

rule merge_reads:
    input:
        r11 = config["input_dir"]+"{sample}_L001"+config["forward_read_suffix"],
        r12 = config["input_dir"]+"{sample}_L001"+config["reverse_read_suffix"],
        r21 = config["input_dir"]+"{sample}_L002"+config["forward_read_suffix"],
        r22 = config["input_dir"]+"{sample}_L002"+config["reverse_read_suffix"]
    output:
        o1=config["output_dir"] + "/merged_data/{sample}_R1.fastq",
        o2=config["output_dir"] + "/merged_data/{sample}_R2.fastq"
    params:
        r11 = config["input_dir"]+"{sample}_L001_R1.fastq",
        r12 = config["input_dir"]+"{sample}_L001_R2.fastq",
        r21 = config["input_dir"]+"{sample}_L002_R1.fastq",
        r22 = config["input_dir"]+"{sample}_L002_R2.fastq"
    shell:
#        "gunzip  {input.r11};"
#        "gunzip  {input.r12};"
#        "gunzip  {input.r21};"
#        "gunzip {input.r22};"
        "cat {params.r11} {params.r21} > {output.o1};"
        "cat {params.r12} {params.r22} > {output.o2};"

rule fastqc_raw:
    input:
        r1 = config["output_dir"] + "/merged_data/{sample}_R1.fastq",
        r2 = config["output_dir"] + "/merged_data/{sample}_R2.fastq"
    output:
        r1 = os.path.join(config["output_dir"],"fastqc_raw","{sample}"+forward_read_num+"_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"fastqc_raw","{sample}"+reverse_read_num+"_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"fastqc_raw")
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_raw:
    input:
        r1 = expand(os.path.join(config["output_dir"],"fastqc_raw","{sample}"+forward_read_num+"_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"fastqc_raw","{sample}"+reverse_read_num+"_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_raw.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"fastqc_raw/"),
        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -c utils/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_raw.html"

rule prinseq:
    input:
        r1 = config["output_dir"] + "/merged_data/{sample}_R1.fastq",
        r2 = config["output_dir"] + "/merged_data/{sample}_R2.fastq"
    params:
        prefix = os.path.join(config["output_dir"],"prinseq","{sample}_filtered")
    output:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    conda: "utils/envs/prinseq_env.yaml"
    shell:
            "perl utils/scripts/prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} "
            "-trim_left {config[trimleft]} -trim_right {config[trimright]} "
            "-out_good {params.prefix} -out_bad null -lc_method {config[lc_method]} -lc_threshold {config[lc_threshold]} "
            "-derep {config[derep]} -trim_qual_type {config[trim_qual_type]} -trim_qual_window "
            "{config[trim_qual_window]} -trim_qual_step {config[trim_qual_step]} "
            "-trim_qual_rule {config[trim_qual_rule]} -trim_qual_left {config[trim_qual_left]} "
            "-trim_qual_right {config[trim_qual_right]} -min_len {config[minlength]} "
            "-ns_max_n {config[maxn]}"

rule fastqc_prinseq_filt:
    input:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/")
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_prinseq_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -c utils/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_prinseq_filtered.html"

rule spades_assembly:
    input:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        spades_scaffolds_assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","spades/scaffolds.fasta")
    params:
        memory_in_gb = config["memory_in_gb"],
        threads = config["assembler_threads"],
        sample_assembly_dir = os.path.join(config["output_dir"],"genomes","{sample}","assembly","spades")
    conda: "utils/envs/spades_env.yaml"
#    conda: "spades_env"
    shell:
        "spades.py -t {params.threads} -m {params.memory_in_gb} -o {params.sample_assembly_dir} -1 {input.r1} -2 {input.r2}"

rule filter_sequences_by_length:
    input:
        spades_scaffolds_assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","spades/scaffolds.fasta")
    output:
         renamed_genome_assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa")
    params:
        min_scaffold_length = config["min_scaffold_length"]
    conda: "utils/envs/biopython_env.yaml"
    shell:
        "python utils/scripts/filter_sequences_by_length.py -i {input.spades_scaffolds_assembly_file} -l {params.min_scaffold_length} -o {output.renamed_genome_assembly_file}"
    
rule quast:
    input:
        assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa")
    output:
        quast_transposed_report_file = os.path.join(config["output_dir"],"genomes","{sample}","quast","transposed_report.tsv")
    params:
        quast_dir = os.path.join(config["output_dir"],"genomes","{sample}","quast"),
        threads = config["quast_threads"]
    conda: "utils/envs/quast_env.yaml"
    shell:
       "quast.py --output-dir {params.quast_dir} --threads {params.threads} {input.assembly_file}"

rule prokka:
    input:
        assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa")
    output:
        prokka_fna_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.gff"),
        prokka_faa_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.faa")
    params:
        prokka_dir = os.path.join(config["output_dir"],"genomes","{sample}","prokka"),
        threads = config["prokka_threads"],
	prefix = "{sample}"
    conda: "utils/envs/prokka_env.yaml"
    shell:
       "prokka --outdir {params.prokka_dir} --prefix {params.prefix} {input.assembly_file} --cpus {params.threads} --rfam 1 --force"

rule extract_marker_sequences:
    input:
        prokka_fna_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.gff")
    output:
        extracted_marker_seqs_csv_file = os.path.join(config["output_dir"],"genomes","{sample}","extracted_sequences","cpn60_metadata.csv"),
    params:
        extracted_sequences_dir = os.path.join(config["output_dir"],"genomes","{sample}","extracted_sequences"),
    conda: "utils/envs/biopython_env.yaml"
    shell:
        "python utils/scripts/extract_marker_sequences.py --fasta_infile {input.prokka_fna_file} --gff_infile {input.prokka_gff_file} --output_dir {params.extracted_sequences_dir}"

#rule metaerg:
#    input:
#        assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa")
#    output:
##        metaerg_fna_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","metaerg","{sample}_genome.fna"),
#        metaerg_gff_file = os.path.join(config["output_dir"],"genomes","{sample}","metaerg","data","all.gff")
#    params:
#        metaerg_dir = os.path.join(config["output_dir"],"genomes","{sample}","metaerg"),
#        metaerg_database_path = config["metaerg_database_path"],
#        locustag = "{sample}",
#        threads = config["metaerg_threads"]
#    shell:
#       "singularity run -H $HOME -B {params.metaerg_database_path}:/NGStools/metaerg/db -B /work:/work -B /bulk:/bulk /global/software/singularity/images/software/metaerg2.sif /NGStools/metaerg/bin/metaerg.pl --mincontiglen 200 --gcode 11 --gtype meta --minorflen 180 --cpus {params.threads} --evalue 1e-05 --identity 20 --coverage 70 --locustag {params.locustag} --force --outdir {params.metaerg_dir} {input.assembly_file}"

rule checkm:
    input:
        assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa")
    output:
        checkm_table_file = os.path.join(config["output_dir"],"genomes","{sample}","checkm","checkm.tsv")
    params:
        checkm_database = config["checkm_database_path"],
        checkm_dir = os.path.join(config["output_dir"],"genomes","{sample}","checkm"),
	threads = config["checkm_threads"]
    conda: "utils/envs/checkm_env.yaml"
    shell:
       "checkm data setRoot {params.checkm_database}; "
       "mkdir -p {params.checkm_dir}; "
       "filename=$(basename {input.assembly_file}); "
       "cp {input.assembly_file} {params.checkm_dir}/$filename; "
       "checkm lineage_wf -t {params.threads} -x fa --tab_table --file {output.checkm_table_file} {params.checkm_dir} {params.checkm_dir}; "

rule gtdbtk:
    input:
        assembly_file = os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa")
    output:
        gtdbtk_file = os.path.join(config["output_dir"],"genomes","{sample}","gtdbtk","gtdbtk.bac120.summary.tsv")
    params:
       gtdbtk_data_path = config["gtdbtk_database_path"],
       gtdbtk_dir = os.path.join(config["output_dir"],"genomes","{sample}","gtdbtk"),
       threads = config["gtdbtk_threads"]
    conda: "utils/envs/gtdbtk_env.yaml"
    shell:
       "GTDBTK_DATA_PATH=\"{params.gtdbtk_data_path}\"; "
       "filename=$(basename {input.assembly_file}); "
       "cp {input.assembly_file} {params.gtdbtk_dir}/$filename; "
       "gtdbtk classify_wf --genome_dir {params.gtdbtk_dir} --extension \"fa\" --cpus {params.threads} --out_dir {params.gtdbtk_dir}; "


rule eggnog_mapper:
    input:
       prokka_faa_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.faa"), 
    output:
       eggnog_mapper_file = os.path.join(config["output_dir"],"genomes","{sample}","eggnog_mapper","{sample}.emapper.annotations"),
    params:
       eggnog_mapper_db = config["eggnog_mapper_db"],
       eggnog_mapper_output_file_prefix = os.path.join(config["output_dir"],"genomes","{sample}","eggnog_mapper","{sample}"),
       threads = config["eggnog_mapper_threads"]
    conda: "utils/envs/eggnog_mapper_env.yaml"
    shell:
       "python /bulk/IMCshared_bulk/shared/shared_software/eggnog-mapper/emapper.py -i {input.prokka_faa_file} --itype proteins --cpu {params.threads} --data_dir {params.eggnog_mapper_db} --output {params.eggnog_mapper_output_file_prefix}; "

#rule dbcan:
#    input:
#       prokka_faa_file = os.path.join(config["output_dir"],"genomes","{sample}","prokka","{sample}.faa"),
#    output:
#       dbcan_overview_file = os.path.join(config["output_dir"],"genomes","{sample}","eggnog_mapper","{sample}.emapper.annotations"),
#    params:
#       eggnog_mapper_db = config["eggnog_mapper_db"],
#       eggnog_mapper_output_file_prefix = os.path.join(config["output_dir"],"genomes","{sample}","eggnog_mapper","{sample}"),
#       threads = config["eggnog_mapper_threads"]
#    conda: "utils/envs/eggnog_mapper_env.yaml"
#    shell:
#       "python /bulk/IMCshared_bulk/shared/shared_software/eggnog-mapper/emapper.py -i {input.prokka_faa_file} --itype proteins --cpu {params.threads} --data_dir {params.eggnog_mapper_db} --output {params.eggnog_mapper_output_file_prefix}; "

#rule merge_assembly_analysis_data:
#    input:
#        expand(os.path.join(config["output_dir"],"genomes","{sample}","assembly","{sample}_genome.fa"),sample=SAMPLES)
#    output:
#        merged_metadata_file = os.path.join(config["output_dir"],"metadata_files","all_merged_assembly_analysis_metadata.tsv"),

#    params:
#       input_dir = os.path.join(config["output_dir"],"genomes"),
#       output_dir = os.path.join(config["output_dir"],"metadata_files")
#    shell:
#       "python utils/scripts/merge_assembly_analysis_data.py --input_dir {params.input_dir} --output_dir {params.output_dir}"


