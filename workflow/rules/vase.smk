import glob


def get_vase_annot_params(wildcards):
    args = "-g {} -d {}".format(config['gnomad'],
                                get_variation_vcf())
    if config.get("cadd_dir"):
        args += " --cadd_dir {}".format(config.get("cadd_dir"))
    if config.get("cadd_files"):
        args += " --cadd_files {}".format(config.get("cadd_files"))
    if config.get("splice_ai_files"):
        args += " --splice_ai_vcfs {}".format(config.get("splice_ai_files"))
    if config.get("splice_ai_dir"):
        svcfs = glob.glob("{}/*.vcf.gz".format(config.get("splice_ai_dir")))
        args += " --splice_ai_vcfs {}".format(" ".join(svcfs))
    if config.get("cadd_files") or config.get("cadd_dir"):
        args += " --missing_cadd_scores results/vase_annotated/missing_cadd_scores.{}.vcf.gz".format(wildcards.contig)
    if config.get("splice_ai_files") or config.get("splice_ai_dir"):
        args += " --missing_splice_ai_scores results/vase_annotated/missing_spliceai_scores.{}.vcf.gz".format(wildcards.contig)
    if config.get("extra"):
        args += " {}".format(config.get("extra"))
    return args


def get_vase_filter_params(wildcards):
    args = '--ped {}'.format(config['ped'])
    if wildcards.seg == 'recessive':
        args += ' --recessive'
    else:
        args += ' --de_novo'
    args += '--freq {}'.format(config['filtering']['freq'][wildcards.seg])
    args += '--csq {}'.format(config['filtering']['csq'])
    arg2flag = {'het_vaf': '--het_ab',
                'hom_vaf': '--hom_ab',
                'con_het_vaf': '-con_het_ab',
                'con_hom_vaf': '-con_hom_ab',
                'max_con_ref_vaf': '-con_ref_ab',
                'gq': '--gq',
                'dp': '--dp'}
    for k, flag in arg2flag.items():
        if k in config['filtering']:
            args += ' {} {}'.format(flag, config['filtering'][k])
    return args


rule tabix_annotated_variants:
    input:
        config['vcf'],
    output:
        "results/annotated/all.vcf.gz.tbi",
        config['vcf'] + '.tbi',
    log:
        "logs/tabix/variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "0.84.0/bio/tabix"


rule vase_annotate:
    input:
        vcf=config['vcf'],
        tbi=config['vcf'] + '.tbi',
    output:
        temp("results/vase_annotated/all.{contig}.vcf.gz")
    params:
        get_vase_annot_params
    conda:
        "../envs/vase.yaml"
    log:
        "logs/vase/vase_annot.{contig}.log"
    resources:
        mem_mb=4096
    shell:
        "bcftools view -O u {input.vcf} {wildcards.contig} | vase {params} -i - -o {output}"


rule vase_filter:
    input:
        vcf="results/vase_annotated/all.vcf.gz",
    output:
        "results/vase_filtered/all.{seg}.vcf.gz",
    params:
        get_vase_annot_params
    conda:
        "../envs/vase.yaml"
    log:
        "logs/gatk/vase_annot_{seg}.log"
    shell:
        "vase {params} -i {input.vcf} -o {output.vcf}"


rule merge_annotated_variants:
    input:
        vcfs=lambda w: expand(
            "results/vase_annotated/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/vase_annotated/all.vcf.gz",
    log:
        "logs/picard/merge-vase_annotated.log",
    wrapper:
        "0.84.0/bio/picard/mergevcfs"


rule merge_filtered_variants:
    input:
        vcfs=lambda w: expand(
            "results/vase_filtered/all.{seg}.vcf.gz",
            seg=['recessive', 'de_novo', 'dominant']
        ),
    output:
        vcf="results/vase_filtered/all.vcf.gz",
    log:
        "logs/picard/merge-vase_fitlered.log",
    wrapper:
        "0.84.0/bio/picard/mergevcfs"
