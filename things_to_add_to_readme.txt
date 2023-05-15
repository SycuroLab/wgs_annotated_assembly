# To install multiqc
Place this into bash_profile or bash_rc. Or in the rule itself. Probably just bash_profile.
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

 1043  gunzip -c SRR7659115_pass_1.fastq.gz | sed 's;\(^@SRR7659115\.[0-9]*\)\.[0-9].*;\1/1;g' > SRR7659115_pass_1.fastq
 1044  gunzip -c SRR7659115_pass_1.fastq.gz | sed 's;\(^@SRR7659115\.[0-9]*\)\.[0-9].*;\1/2;g' > SRR7659115_pass_2.fastq

