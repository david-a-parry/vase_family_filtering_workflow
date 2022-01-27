import pandas as pd
import glob
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")



###### Config file and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")


def get_contigs():
    with checkpoints.write_contigs.get().output[0].open() as contigs:
        return pd.read_csv(contigs, usecols=[0], dtype=str, header=None,
                           squeeze=True)


def get_segregation_modes():
    if 'filtering' in config:
        return [x for x in ['recessive', 'dominant', 'de_novo'] if x in
                config['filtering']]
    return None


def get_vase_annot_params(wildcards):
    args = "--log_progress"
    if 'gnomad' in config:
        args += " -g {}".format(config['gnomad'])
    if 'dbsnp' in config:
        args += " -d {}".format(config['dbsnp'])
    if config.get("cadd_files") or config.get("cadd_dir"):
        args += " --missing_cadd_scores results/vase_annotated/missing_cadd_scores.{}.vcf.gz".format(wildcards.contig)
        if config.get("cadd_dir"):
            args += " --cadd_dir {}".format(config.get("cadd_dir"))
        if config.get("cadd_files"):
            args += " --cadd_files {}".format(config.get("cadd_files"))
    if config.get("splice_ai_files") or config.get("splice_ai_dir"):
        args += " --missing_splice_ai_scores results/vase_annotated/missing_splice_ai_scores.{}.vcf.gz".format(wildcards.contig)
        if config.get("splice_ai_files"):
            args += " --splice_ai_vcfs {}".format(config.get("splice_ai_files"))
        if config.get("splice_ai_dir"):
            svcfs = glob.glob("{}/*.vcf.gz".format(config.get("splice_ai_dir")))
            args += " --splice_ai_vcfs {}".format(" ".join(svcfs))
    if 'annotating' in config:
        if config['annotating'].get("extra"):
            args += " {}".format(config['annotating'].get("extra"))
    return args


def get_vase_filter_params(wildcards):
    mode = wildcards.seg
    args = '--log_progress --ped {}'.format(config['ped'])
    if mode == 'recessive':
        args += ' --recessive'
    elif mode == 'dominant':
        args += ' --dominant'
    else:
        args += ' --de_novo'
    args += ' --freq {}'.format(config['filtering'][mode]['freq'])
    args += ' --csq {}'.format(config['filtering'][mode]['csq'])
    arg2flag = {'het_vaf': '--het_ab',
                'hom_vaf': '--hom_ab',
                'con_het_vaf': '-con_het_ab',
                'con_hom_vaf': '-con_hom_ab',
                'max_con_ref_vaf': '-con_ref_ab',
                'gq': '--gq',
                'dp': '--dp'}
    for k, flag in arg2flag.items():
        if k in config['filtering'][mode]:
            args += ' {} {}'.format(flag, config['filtering'][mode][k])
    return args


def get_reporter_args():
    args = "-p {} {}".format(config['ped'], config['report']['extra'])
    if config.get("splice_ai_files") or config.get("splice_ai_dir"):
        args += " -i SpliceAI"
    return args
