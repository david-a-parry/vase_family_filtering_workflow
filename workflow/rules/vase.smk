rule tabix_variants:
    input:
        config['vcf'],
    output:
        config['vcf'] + '.tbi',
    log:
        "logs/tabix/original.log",
    conda:
        "../envs/vase.yaml"
    shell:
        "tabix -fp vcf {input} 2> {log}"


checkpoint write_contigs:
    input:
        config['vcf'],
        config['vcf'] + '.tbi',
    output:
        "resources/contigs.txt",
    log:
        "logs/tabix/write_contigs.log",
    conda:
        "../envs/vase.yaml"
    shell:
        "tabix -l {input[0]} > {output} 2> {log}"


rule tabix_annotated_variants:
    input:
        "results/vase_annotated/vase_anno.vcf.gz",
    output:
        "results/vase_annotated/vase_anno.vcf.gz.tbi",
    log:
        "logs/tabix/annotated.log",
    conda:
        "../envs/vase.yaml"
    shell:
        "tabix -fp vcf {input} 2> {log}"


rule vase_annotate:
    input:
        vcf=config['vcf'],
        tbi=config['vcf'] + '.tbi',
        contigs="resources/contigs.txt",
    output:
        temp("results/vase_annotated/vase_anno.{contig}.vcf.gz")
    params:
        get_vase_annot_params
    conda:
        "../envs/vase.yaml"
    log:
        "logs/vase/vase_annot.{contig}.log"
    resources:
        mem_mb=4096
    shell:
        "bcftools view -O u {input.vcf} {wildcards.contig} | vase {params} -i - -o {output} 2> {log}"


rule vase_filter:
    input:
        vcf="results/vase_annotated/vase_anno.vcf.gz",
    output:
        vcf="results/vase_filtered/vase_filtered.{seg}.vcf.gz",
    params:
        get_vase_filter_params
    conda:
        "../envs/vase.yaml"
    log:
        "logs/vase/vase_filter_{seg}.log"
    shell:
        "vase {params} -i {input.vcf} -o {output.vcf} 2> {log}"


rule merge_annotated_variants:
    input:
        vcfs=lambda w: expand(
            "results/vase_annotated/vase_anno.{contig}.vcf.gz",
            contig=get_contigs()
        ),
    output:
        vcf="results/vase_annotated/vase_anno.vcf.gz",
    log:
        "logs/picard/merge_vase_annotated.log",
    wrapper:
        "0.84.0/bio/picard/mergevcfs"


if len(get_segregation_modes()) > 1:
  rule merge_filtered_variants:
      input:
          vcfs=lambda w: expand(
              "results/vase_filtered/vase_filtered.{seg}.vcf.gz",
              seg=get_segregation_modes(),
          ),
      output:
          vcf="results/vase_filtered/vase_filtered.all.vcf.gz",
      log:
          "logs/picard/merge_vase_filtered.log",
      wrapper:
          "0.84.0/bio/picard/mergevcfs"
else:
  rule rename_single_output:
      input:
          vcfs=lambda w: expand(
              "results/vase_filtered/vase_filtered.{seg}.vcf.gz",
              seg=get_segregation_modes(),
          ),
      output:
          "results/vase_filtered/vase_filtered.all.vcf.gz"
      log:
          "logs/picard/merge_vase_filtered.log",
      shell:
        "mv {input} {output} 2> {log}"


rule create_variant_report:
    input:
        vcf="results/vase_filtered/vase_filtered.all.vcf.gz",
    output:
        "results/report/vase_filtered.report.{ext}"
    conda:
        "../envs/vase.yaml"
    log:
        "logs/reporter/{ext}_report.log",
    params:
        get_reporter_args()
    shell:
        "vase_reporter {params} -o {wildcards.ext} {input} {output} 2> {log}"

