# ***************************************
# * Snakefile for genome_assembly_analysis pipeline *
# ***************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
import os

os.environ["GTDBTK_DATA_PATH"] = "/bulk/IMCshared_bulk/shared/dbs/gtdbtk-1.5.0/db"

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:



        os.path.join(config["output_dir"],"multiqc","multiqc_report_raw.html"),
#   expand(config["input_dir"]+"/assembly_analysis/{sample}.fa", sample=SAMPLES),
        expand(os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"), sample=SAMPLES)
#        os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html"),
        expand(config["output_dir"]+"/{sample}/assembly/spades/scaffolds.fasta",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/spades/long_scaffolds.fasta",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/spades/mapped_spades_assembly_reads.sam",sample=SAMPLES),
 #       expand(config["output_dir"]+"/{sample}/assembly/spades/unmapped_spades_assembly_reads_R1.fastq",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/spades/unmapped_spades_assembly_reads_R2.fastq",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/megahit/fixed_headers_final.contigs.fa",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/final_assembly.fasta",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/{sample}_metagenome.fasta",sample=SAMPLES),
        
#        expand(config["output_dir"]+"/assembly_analysis/{sample}/quast/transposed_report.tsv", sample=SAMPLES),
#        expand(config["output_dir"]+"/assembly_analysis/{sample}/prokka/{sample}.fna",sample=SAMPLES),
#        expand(config["output_dir"]+"/assembly_analysis/{sample}/prokka/{sample}.gff",sample=SAMPLES),
#    expand(config["output_dir"]+"/assembly_analysis/{sample}/extracted_sequences/cpn60_metadata.csv",sample=SAMPLES),
##        expand(config["output_dir"]+"/assembly_analysis/{sample}/metaerg/data/all.gff",sample=SAMPLES),
 #       expand(config["output_dir"]+"/assembly_analysis/{sample}/checkm/checkm.tsv",sample=SAMPLES),
 #       expand(config["output_dir"]+"/assembly_analysis/{sample}/gtdbtk/gtdbtk.bac120.summary.tsv",sample=SAMPLES)


rule fastqc_raw:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
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

#rule prinseq:
#    input:
#        r1 = os.path.join(config["output_dir"],"cutadapt","{sample}_r1_trimmed.fastq") if config["run_cutadapt"] else os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
 #       r2 = os.path.join(config["output_dir"],"cutadapt","{sample}_r2_trimmed.fastq") if config["run_cutadapt"] else os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
#    params:
#        prefix = os.path.join(config["output_dir"],"prinseq","{sample}_filtered")
#    output:
#        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
#        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
#    conda: "utils/envs/prinseq_env.yaml"
#    shell:
#            "perl utils/scripts/prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} "
#            "-trim_left {config[trimleft]} -trim_right {config[trimright]} "
#            "-out_good {params.prefix} -out_bad null -lc_method {config[lc_method]} -lc_threshold {config[lc_threshold]} "
#            "-derep 1 -trim_qual_type {config[trim_qual_type]} -trim_qual_window "
#            "{config[trim_qual_window]} -trim_qual_step {config[trim_qual_step]} "
#            "-trim_qual_rule {config[trim_qual_rule]} -trim_qual_left {config[trim_qual_left]} "
#            "-trim_qual_right {config[trim_qual_right]} -min_len {config[minlength]} "
#            "-ns_max_n {config[maxn]}"

#rule fastqc_prinseq_filt:
#    input:
#        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
#        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
#    output:
#        r1 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"),
#        r2 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html")
#    params:
#        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/")
#    conda: "utils/envs/fastqc_env.yaml"
#    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

#rule multiqc_prinseq_filt:
#    input:
#        r1 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"), sample=SAMPLES),
#        r2 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html"), sample=SAMPLES)
#    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html")
#    params:
#        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/"),
#        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
#    conda: "utils/envs/multiqc_env.yaml"
#    shell: "multiqc -c utils/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_prinseq_filtered.html"

#rule spades_assembly:
#    input:
#        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
#        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
#    output:
#        spades_scaffolds_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/scaffolds.fasta")
#    params:
#        memory_in_gb = config["memory_in_gb"],
#        threads = config["assembler_threads"],
#        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly","spades")
#    conda: "utils/envs/spades_env.yaml"
#    shell:
#        "spades.py -t {params.threads} -m {params.memory_in_gb} -o {params.sample_assembly_dir} -1 {input.fastq_read1} -2 {input.fastq_read2}"

#rule filter_long_scaffolds:
#    input:
#        spades_scaffolds_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/scaffolds.fasta")
#    output:
#        long_scaffolds_spades_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/long_scaffolds.fasta")
#    params:
#        min_sequence_length = config["min_sequence_length"]
#    conda: "utils/envs/biopython_env.yaml"
#    shell:
#        "python utils/scripts/filter_sequences_by_length.py -i {input.spades_scaffolds_assembly_file} -l {params.min_sequence_length} -o {output.long_scaffolds_spades_file}"
        
#rule map_reads_to_assembly:
#    input:
#        long_scaffolds_spades_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/long_scaffolds.fasta"),
#        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
#        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
#    output:
#        spades_scaffolds_sam_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/mapped_spades_assembly_reads.sam")
#    params:
#        threads = config["assembler_threads"],
#        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly")
#    conda: "utils/envs/bwa_env.yaml"
#    shell:
#        "bwa index {input.long_scaffolds_spades_file}; "
#        "bwa mem -t {params.threads} {input.long_scaffolds_spades_file} {input.fastq_read1} {input.fastq_read2} > {output.spades_scaffolds_sam_file}; "

#rule recover_unmapped_spades_reads:
#    input:
#        spades_scaffolds_sam_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/mapped_spades_assembly_reads.sam")
#    output:
#        unmapped_spades_scaffolds_read1_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/unmapped_spades_assembly_reads_R1.fastq"),
#        unmapped_spades_scaffolds_read2_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/unmapped_spades_assembly_reads_R2.fastq")
#    params:
#        threads = config["assembler_threads"],
#        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly"),
#        spades_scaffolds_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/mapped_spades_assembly_reads.bam"),
#        unmapped_spades_scaffolds_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/unmapped_spades_assembly_reads.bam")
#    conda: "utils/envs/samtools_env.yaml"
#    shell:
#        "samtools view -S -b {input.spades_scaffolds_sam_file} --threads {params.threads} -o {params.spades_scaffolds_bam_file}; "
#        "samtools view -b -f 4 {params.spades_scaffolds_bam_file} --threads {params.threads} -o {params.unmapped_spades_scaffolds_bam_file}; "
#        "samtools fastq {params.unmapped_spades_scaffolds_bam_file} --threads {params.threads} -1 {output.unmapped_spades_scaffolds_read1_fastq_file} -2 {output.unmapped_spades_scaffolds_read2_fastq_file}; "

#rule megahit_unmapped_spades_reads:
#    input:
#        unmapped_spades_scaffolds_read1_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/unmapped_spades_assembly_reads_R1.fastq"),
#        unmapped_spades_scaffolds_read2_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/unmapped_spades_assembly_reads_R2.fastq")
#    output:
#        fixed_megahit_final_contigs_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/fixed_headers_final.contigs.fa")
#    params:
#        megahit_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/final.contigs.fa"),
#        threads = config["assembler_threads"],
#        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly"),
#        memory_in_gb = config["memory_in_gb"]
#    conda: "utils/envs/megahit_env.yaml"
#    shell:
#        "megahit -1 {input.unmapped_spades_scaffolds_read1_fastq_file} -2 {input.unmapped_spades_scaffolds_read2_fastq_file} -t {params.threads} -m {params.memory_in_gb}000000000 -f -o {params.sample_assembly_dir}/megahit; "
#        "sed 's/ /_/g' {params.megahit_assembly_file} > {output.fixed_megahit_final_contigs_file} "

#rule combine_and_sort_assembly:
#    input:
#        long_scaffolds_spades_file = os.path.join(config["output_dir"],"{sample}","assembly","spades/long_scaffolds.fasta"),
#        fixed_megahit_final_contigs_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/fixed_headers_final.contigs.fa")
#    output:
        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta")
#    params:
#        threads = config["assembler_threads"],
#        combined_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","combined_file.fasta"),
#        long_scaffolds_megahit_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/fixed_headers_long.contigs.fa"),
#        min_sequence_length = config["min_sequence_length"],
#        memory_in_gb = config["memory_in_gb"]

#    conda: "utils/envs/biopython_env.yaml"
#    shell:
#        "python utils/scripts/filter_sequences_by_length.py -i {input.fixed_megahit_final_contigs_file} -l {params.min_sequence_length} -o {params.long_scaffolds_megahit_file}; "
#        "cat {input.long_scaffolds_spades_file} {params.long_scaffolds_megahit_file} > {params.combined_assembly_file}; "
#        "python utils/scripts/sort_assembly_by_length.py -i {params.combined_assembly_file} -o {output.metagenome_final_assembly_file}; "
        
#rule rename_final_assembly_file:
#    input:
#        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta")
    output:
#        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
#    conda: "utils/envs/biopython_env.yaml"
#    shell:
#        "python utils/scripts/fix_fasta_header_length.py -i {input.metagenome_final_assembly_file} -o {output.renamed_metagenome_assembly_file}"

#rule quast:
#    input:
#        assembly_file = os.path.join(config["input_dir"], "{sample}.fa")
#    output:
#        quast_transposed_report_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","quast","transposed_report.tsv")
#    params:
#        quast_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","quast"),
#        threads = config["quast_threads"]
#    conda: "utils/envs/quast_env.yaml"
#    shell:
#       "quast.py --output-dir {params.quast_dir} --threads {params.threads} {input.assembly_file}"

#rule prokka:
#    input:
#        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
#    output:
#        prokka_fna_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.fna"),
#        prokka_gff_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.gff")
#    params:
#        prokka_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka"),
#        threads = config["prokka_threads"],
#	prefix = "{sample}"
#    conda: "utils/envs/prokka_env.yaml"
#    shell:
#       "prokka --metagenome --outdir {params.prokka_dir} --prefix {params.prefix} {input.assembly_file} --cpus {params.threads} --rfam 1 --force"

#rule extract_marker_sequences:
#    input:
#        prokka_fna_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.fna"),
#        prokka_gff_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.gff")
#    output:
#        extracted_marker_seqs_csv_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","extracted_sequences","cpn60_metadata.csv"),
#    params:
#        extracted_sequences_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","extracted_sequences"),
#    conda: "utils/envs/biopython_env.yaml"
#    shell:
#        "python utils/scripts/extract_marker_sequences.py --fasta_infile {input.prokka_fna_file} --gff_infile {input.prokka_gff_file} --output_dir {params.extracted_sequences_dir}"

#rule metaerg:
#    input:
#        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
#    output:
##        metaerg_fna_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","assembly","metaerg","{sample}_metagenome.fna"),
#        metaerg_gff_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","metaerg","data","all.gff")
#    params:
#        metaerg_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","metaerg"),
#        metaerg_database_path = config["metaerg_database_path"],
#        locustag = "{sample}",
#        threads = config["metaerg_threads"]
#    shell:
#       "singularity run -H $HOME -B {params.metaerg_database_path}:/NGStools/metaerg/db -B /work:/work -B /bulk:/bulk /global/software/singularity/images/software/metaerg2.sif /NGStools/metaerg/bin/metaerg.pl --mincontiglen 200 --gcode 11 --gtype meta --minorflen 180 --cpus {params.threads} --evalue 1e-05 --identity 20 --coverage 70 --locustag {params.locustag} --force --outdir {params.metaerg_dir} {input.assembly_file}"

#rule checkm:
#    input:
#        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
#    output:
#        checkm_table_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","checkm","checkm.tsv")
#    params:
#        checkm_database = config["checkm_database_path"],
#        checkm_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","checkm"),
#	threads = config["checkm_threads"]
#    conda: "utils/envs/checkm_env.yaml"
#    shell:
#       "checkm data setRoot {params.checkm_database}; "
#       "mkdir -p {params.checkm_dir}; "
#       "filename=$(basename {input.assembly_file}); "
#       "cp {input.assembly_file} {params.checkm_dir}/$filename; "
#       "checkm lineage_wf -t {params.threads} -x fa --tab_table --file {output.checkm_table_file} {params.checkm_dir} {params.checkm_dir}; "

#rule gtdbtk:
#    input:
#        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
#    output:
#        gtdbtk_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","gtdbtk","gtdbtk.bac120.summary.tsv")
#    params:
#       gtdbtk_data_path = config["gtdbtk_database_path"],
#       gtdbtk_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","gtdbtk"),
#       threads = config["gtdbtk_threads"]
#    conda: "utils/envs/gtdbtk_env.yaml"
#    shell:
#       "GTDBTK_DATA_PATH=\"{params.gtdbtk_data_path}\"; "
#       "filename=$(basename {input.assembly_file}); "
#       "cp {input.assembly_file} {params.gtdbtk_dir}/$filename; "
#       "gtdbtk classify_wf --genome_dir {params.gtdbtk_dir} --extension \"fa\" --cpus {params.threads} --out_dir {params.gtdbtk_dir}; "