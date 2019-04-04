#!/usr/bin/env python3

import pandas
import pathlib2
import multiprocessing


#############
# FUNCTIONS #
#############

def write_config_file(fastq, threads, config_string, config_file):
    '''
    Accept fastq file, threads config string and output location and write
    config
    '''
    my_fastq = str(pathlib2.Path(fastq).resolve())
    my_conf = config_string.format(my_fastq, threads)
    with open(config_file, 'wt') as f:
        f.write(my_conf)
    return True


def read_goi_exons(goi_exons_file):
    csv = pandas.read_csv(goi_exons_file)
    region_list = ['--region {0}:{1}-{2}'.format(
        row[1][0], row[1][1], row[1][2])
        for row in csv.iterrows()]
    return(' '.join(region_list))


###########
# GLOBALS #
###########

run_info_file = 'data/SraRunInfo.csv'
csd_locus_fasta = 'data/csd.fasta'
honeybee_ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
bbduk_ref = '/sequencing_artifacts.fa.gz'
bbduk_adaptors = '/adapters.fa'
meraculous_config_file = 'src/meraculous_config.txt'
meraculous_threads = 50

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
mer_container = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
spades_container = 'shub://TomHarrop/singularity-containers:spades_3.12.0'
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
vcflib_container = 'shub://TomHarrop/singularity-containers:vcflib_1.0.0-rc2'
bioc_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.7'

########
# MAIN #
########

# get a list of all names
run_info = pandas.read_csv(run_info_file)
all_samples = sorted(set(run_info['SampleName']))

# TEMPORARY # subset samples to get the pipeline to run quickly
# all_samples = [x for x in all_samples if x.startswith('REUd')]

# read the meraculous config
with open(meraculous_config_file, 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())


#########
# RULES #
#########

rule target:
    input:
        # expand(('output/050_meraculous/'
        #         '{sample}/'
        #         'meraculous_final_results/final.scaffolds.fa'),
        #        sample=all_samples),
        # expand('output/040_norm/{sample}_kmer_plot.pdf',
        #        sample=all_samples),
        # expand('output/040_norm/{sample}_kmer_plot.pdf',
        #        sample=all_samples),
        # expand('output/070_bwa/{sample}_marked.bam',
        #        sample=all_samples),
        'output/080_freebayes/variants_filtered.vcf'

# run freebayes
rule filter_vcf:
    input:
        'output/080_freebayes/variants.vcf'
    output:
        'output/080_freebayes/variants_filtered.vcf'
    params:
        filter = "QUAL > 100"
    log:
        'output/000_logs/080_freebayes/freebayes_filter.log'
    singularity:
        vcflib_container
    shell:
        'vcffilter -f "{params.filter}" {input} > {output} 2> {log}'


rule freebayes:
    input:
        bam = expand('output/070_bwa/{sample}_marked.bam',
                     sample=all_samples),
        bai = expand('output/070_bwa/{sample}_marked.bam.bai',
                     sample=all_samples),
        fa = honeybee_ref,
        goi_exons = 'output/080_freebayes/goi_exons.csv'
    output:
        vcf = 'output/080_freebayes/variants.vcf'
    params:
        ploidy = '1',
        region = lambda wildcards, input: read_goi_exons(input.goi_exons)
    log:
        'output/000_logs/080_freebayes/freebayes.log'
    benchmark:
        'output/000_benchmarks/080_freebayes/freebayes.tsv'
    singularity:
        freebayes_container
    shell:
        'freebayes '
        '{params.region} '
        '--ploidy {params.ploidy} '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'

rule extract_goi_exons:
    input:
        gff = ('data/GCF_003254395.2_Amel_HAv3.'
               '1_genomic.cds-exon-mrna-gene.gff3'),
    output:
        goi_exons = 'output/080_freebayes/goi_exons.csv'
    params:
        goi = 'XM_026439432.1'  # CSD transcript ID
    log:
        'output/000_logs/080_freebayes/extract_goi_exons.log'
    singularity:
        bioc_container
    script:
        'src/extract_goi_exons.R'

# mark duplicates
rule index_bamfile:
    input:
        'output/070_bwa/{sample}_marked.bam'
    output:
        'output/070_bwa/{sample}_marked.bam.bai'
    log:
        'output/000_logs/070_bwa/{sample}_index.log',
    threads:
        2
    singularity:
        samtools_container
    shell:
        'samtools index -@ {threads} {input} 2> {log}'

rule markdup:
    input:
        'output/070_bwa/{sample}.sam'
    output:
        sorted = temp('output/070_bwa/{sample}_sorted.bam'),
        marked = 'output/070_bwa/{sample}_marked.bam'
    threads:
        5
    log:
        s = 'output/000_logs/070_bwa/{sample}_sort.log',
        m = 'output/000_logs/070_bwa/{sample}_markdup.log'
    benchmark:
        'output/000_benchmarks/070_bwa/{sample}_markdup.tsv'
    singularity:
        samtools_container
    shell:
        'samtools fixmate '
        '-m '
        '-O BAM '
        '-@ {threads} '
        '{input} '
        '- '
        '2> {log.s} '
        '| '
        'samtools sort '
        '-o {output.sorted} '
        '-O BAM '
        '-l 0 '
        '-@ {threads} '
        '- '
        '2>> {log.s} '
        '; '
        'samtools markdup '
        '-@ {threads} '
        '-s '
        '{output.sorted} '
        '{output.marked} '
        '2> {log.m}'

# map each individual to reference contig
rule bwa:
    input:
        fq = 'output/010_trim-decon/{sample}.fq.gz',
        index = expand('output/070_bwa/honeybee_ref.fasta.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/070_bwa/{sample}.sam'
    params:
        prefix = 'output/070_bwa/honeybee_ref.fasta',
        rg = '\'@RG\\tID:{sample}\\tSM:{sample}\''
    threads:
        multiprocessing.cpu_count()
    log:
        'output/000_logs/070_bwa/{sample}.log'
    benchmark:
        'output/000_benchmarks/070_bwa/{sample}.tsv'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '-p '
        '-R {params.rg} '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2> {log}'


rule index:
    input:
        honeybee_ref
    output:
        expand('output/070_bwa/honeybee_ref.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/070_bwa/honeybee_ref.fasta'
    threads:
        1
    log:
        'output/000_logs/070_bwa/index.log'
    benchmark:
        'output/000_benchmarks/070_bwa/index.tsv'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input} '
        '2> {log}'

# run a meraculous assembly for each indiv
rule spades:
    input:
        fq = 'output/030_bbsplit-csd/{sample}.fq.gz'
    output:
        'output/060_spades/{sample}/scaffolds.fasta'
    params:
        outdir = 'output/060_spades/{sample}'
    threads:
        meraculous_threads
    priority:
        10
    log:
        'output/000_logs/060_spades/{sample}.log'
    benchmark:
        'output/000_benchmarks/060_spades/{sample}.tsv'
    singularity:
        spades_container
    shell:
        'spades.py '
        '-t {threads} '
        # '--careful '
        '--pe1-12 {input.fq} '
        '-o {params.outdir} '
        '&> {log}'


# run a meraculous assembly for each indiv
rule meraculous:
    input:
        fq = 'output/040_norm/{sample}.fq.gz',
        config = 'output/050_meraculous/{sample}/config.txt'
    output:
        contigs = ('output/050_meraculous/'
                   '{sample}/'
                   'meraculous_final_results/final.scaffolds.fa')
    params:
        outdir = 'output/050_meraculous/{sample}'
    threads:
        meraculous_threads
    priority:
        10
    log:
        'output/000_logs/050_meraculous/{sample}.log'
    benchmark:
        'output/000_benchmarks/050_meraculous/{sample}.tsv'
    singularity:
        mer_container
    shell:
        'run_meraculous.sh '
        '-dir {params.outdir} '
        '-config {input.config} '
        '-cleanup_level 2 '
        '&> {log}'

rule meraculous_config:
    input:
        fq = 'output/040_norm/{sample}.fq.gz'
    output:
        config = 'output/050_meraculous/{sample}/config.txt'
    params:
        threads = meraculous_threads
    priority:
        10
    run:
        write_config_file(
            input.fq,
            params.threads,
            meraculous_config_string,
            output.config)

# normalise coverage to 10x and kick out reads < k
rule plot_kmer_coverage:
    input:
        hist_before = 'output/040_norm/{sample}_hist.txt',
        hist_after = 'output/040_norm/{sample}_hist_out.txt',
        peaks = 'output/040_norm/{sample}_peaks.txt'
    output:
        plot = 'output/040_norm/{sample}_kmer_plot.pdf'
    threads:
        1
    priority:
        2
    log:
        log = 'output/logs/040_norm/{sample}_plot-kmer-coverage.log'
    singularity:
        r_container
    script:
        'src/plot_kmer_coverage.R'

rule bbnorm:
    input:
        fq = 'output/030_bbsplit-csd/{sample}.fq.gz'
    output:
        fq_reformat = temp('output/040_norm/{sample}_reformat.fq'),
        fq_norm = 'output/040_norm/{sample}.fq.gz',
        fq_toss = 'output/040_norm/{sample}_toss.fq.gz',
        hist = 'output/040_norm/{sample}_hist.txt',
        hist_out = 'output/040_norm/{sample}_hist_out.txt',
        peaks = 'output/040_norm/{sample}_peaks.txt'
    log:
        reformat = 'output/000_logs/040_norm/{sample}_reformat.log',
        norm = 'output/000_logs/040_norm/{sample}_bbnorm.log'
    benchmark:
        'output/000_benchmarks/040_norm/{sample}.tsv'
    params:
        target = 7                     # 7x 31-mer ~ 10 read depth
    threads:
        2
    priority:
        2
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'threads={threads} '
        'in={input.fq} '
        'int=t '
        'out={output.fq_reformat} '
        'minlength=71 requirebothbad=f '
        '2> {log.reformat} '
        '; '
        'bbnorm.sh '
        'in={output.fq_reformat} '
        'threads={threads} '
        '-Xmx50g '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target={params.target} '
        'min=3 '                       # 3x 31-mer ~ 3 read depth
        'peaks={output.peaks} '
        '2> {log.norm} '

# attempt to k-mer match reads with bbsplit
rule bbsplit_csd:
    input:
        fq = 'output/010_trim-decon/{sample}.fq.gz',
        ref = csd_locus_fasta
    output:
        fq = 'output/030_bbsplit-csd/{sample}.fq.gz',
        stats = 'output/031_bbsplit-csd_stats/{sample}.txt'
    log:
        'output/000_logs/030_bbsplit-csd/{sample}.log'
    benchmark:
        'output/000_benchmarks/030_bbsplit-csd/{sample}.tsv'
    threads:
        2
    priority:
        1
    singularity:
        bbduk_container
    shell:
        'bbsplit.sh '
        'ref_csd={input.ref} '
        'in={input.fq} '
        'interleaved=t '
        'out_csd={output.fq} '
        'refstats={output.stats} '
        'outu=/dev/null '
        'zl=9 '
        '-Xmx50g '
        '2> {log} '

# 01 trim and decontaminate reads
rule trim_decon:
    input:
        r1 = 'data/reads/{sample}_1.fastq.gz',
        r2 = 'data/reads/{sample}_2.fastq.gz'
    output:
        fq = 'output/010_trim-decon/{sample}.fq.gz',
        f_stats = 'output/011_trim-decon_stats/{sample}_filter-stats.txt',
        t_stats = 'output/011_trim-decon_stats/{sample}_trim-stats.txt'
    log:
        filter = 'output/000_logs/010_trim-decon/{sample}_filter.log',  
        trim = 'output/000_logs/010_trim-decon/{sample}_trim.log',
        repair = 'output/000_logs/010_trim-decon/{sample}_repair.log'
    benchmark:
        'output/000_benchmarks/010_trim-decon/{sample}.tsv'
    params:
        filter = bbduk_ref,
        trim = bbduk_adaptors
    threads:
        2
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'repair.sh '
        'in=stdin.fastq '
        'zl=9 '
        '-Xmx50g '
        'out={output.fq} '
        '2> {log.repair} '
