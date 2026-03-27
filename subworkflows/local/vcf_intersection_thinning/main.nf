/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_ISEC                                  } from '../../../modules/nf-core/bcftools/isec'
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP_ISEC          } from '../../../modules/nf-core/bcftools/bgzip/main'
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP_EXCLUDE       } from '../../../modules/nf-core/bcftools/bgzip/main'
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP_THIN          } from '../../../modules/nf-core/bcftools/bgzip/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_ISEC          } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_EXCLUDE       } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_THIN          } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_ISEC          } from '../../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_EXCLUDE       } from '../../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_THIN          } from '../../../modules/nf-core/bcftools/stats/main'
include { MAKE_SCAFFOLD_BED as MAKE_SCAFFOLD_EXCLUDE_BED } from '../../../modules/local/make_scaffold_bed/'
include { MAKE_SCAFFOLD_BED as MAKE_SCAFFOLD_INCLUDE_BED } from '../../../modules/local/make_scaffold_bed/'
include { VCFTOOLS as VCFTOOLS_EXCLUDE                   } from '../../../modules/nf-core/vcftools/'
include { VCFTOOLS as VCFTOOLS_THIN                      } from '../../../modules/nf-core/vcftools/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow VCF_INTERSECTION_THINNING {
    take:
    vcf_tool1 // channel: [ meta, vcf]
    tbi_tool1 // channel: [ meta, tbi]
    vcf_tool2 // channel: [ meta, vcf]
    tbi_tool2 // channel: [ meta, tbi]
    intervals // channel: [ meta, bed, number_of_intervals]
    fasta     // channel: [ meta, fasta]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    vcf_tool1_prepared = vcf_tool1
        .map { meta, vcf ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), vcf)
        }
    vcf_tool2_prepared = vcf_tool2
        .map { meta, vcf ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), vcf)
        }
    tbi_tool1_prepared = tbi_tool1
        .map { meta, tbi ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), tbi)
        }
    tbi_tool2_prepared = tbi_tool2
        .map { meta, tbi ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), tbi)
        }

    ch_vcfs = vcf_tool1_prepared
        .mix(vcf_tool2_prepared)
        .groupTuple()
        .map { meta, vcfs -> tuple(meta, vcfs) }

    ch_tbis = tbi_tool1_prepared
        .mix(tbi_tool2_prepared)
        .groupTuple()
        .map { meta, tbis -> tuple(meta, tbis) }

    ch_isec_input = ch_vcfs
        .join(ch_tbis)
        .map { meta, vcfs, tbis -> tuple(meta, vcfs, tbis, [], [], []) }

    // Run BCFTOOLS_ISEC
    BCFTOOLS_ISEC(ch_isec_input)

    // Collect intersection output
    def has_include = params.include_scaffolds
    def has_exclude = params.exclude_scaffolds
    def need_scaffold_filter = has_include || has_exclude
    intersection = BCFTOOLS_ISEC.out.results
        .map { meta, dir ->
            def file_common = file("${dir}/0002.vcf")
            tuple(meta, file_common)
        }
        .branch { _tuple ->
            filter:      need_scaffold_filter
            passthrough: !need_scaffold_filter
        }

    // bgzip + index the isec VCF (plain .vcf → .vcf.gz) for BCFTOOLS_STATS
    ch_isec_vcf = intersection.filter
        .mix(intersection.passthrough)
        .map { meta, vcf -> tuple(meta + [id: meta.id + "_isec"], vcf) }

    BCFTOOLS_BGZIP_ISEC(ch_isec_vcf)
    versions = versions.mix(BCFTOOLS_BGZIP_ISEC.out.versions)

    BCFTOOLS_INDEX_ISEC(BCFTOOLS_BGZIP_ISEC.out.output)
    versions = versions.mix(BCFTOOLS_INDEX_ISEC.out.versions)

    ch_isec_vcf_tbi = BCFTOOLS_BGZIP_ISEC.out.output
        .join(BCFTOOLS_INDEX_ISEC.out.tbi)

    // Run BCFTOOLS_STATS on isec output
    BCFTOOLS_STATS_ISEC(
        ch_isec_vcf_tbi,
        [[id: 'no_regions'], []],
        [[id: 'no_targets'], []],
        [[id: 'no_samples'], []],
        [[id: 'no_exons'],   []],
        fasta
    )
    versions      = versions.mix(BCFTOOLS_STATS_ISEC.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS_ISEC.out.stats.map { _meta, stats -> stats })

    def include_scaffolds = normalize_scaffold_param(params.include_scaffolds)
    def exclude_scaffolds = normalize_scaffold_param(params.exclude_scaffolds)
    include_ch = include_scaffolds ? channel.value(include_scaffolds) : channel.empty()
    exclude_ch = exclude_scaffolds ? channel.value(exclude_scaffolds) : channel.empty()

    // Make scaffold BED for inclusion
    bed_include = intervals
        .map { meta, bed_file, _number_of_intervals ->
            tuple(meta + [id: "include_scaffolds"], bed_file)
        }
    MAKE_SCAFFOLD_INCLUDE_BED(
        include_ch,
        bed_include
    )

    // Make scaffold BED for exclusion
    bed_exclude = intervals
        .map { meta, bed_file, _number_of_intervals ->
            tuple(meta + [id: "exclude_scaffolds"], bed_file)
        }
    MAKE_SCAFFOLD_EXCLUDE_BED(
        exclude_ch,
        bed_exclude
    )

    bed = MAKE_SCAFFOLD_INCLUDE_BED.out.bed
        .mix(MAKE_SCAFFOLD_EXCLUDE_BED.out.bed)
        .map { _meta, bed_file -> bed_file }.collect()

    vcftools_exclude_input = intersection.filter
        .map { meta, vcf_file ->
            tuple(meta + [id: meta.id + "_scaffolds_excluded"], vcf_file)
        }

    VCFTOOLS_EXCLUDE(
        vcftools_exclude_input,
        bed,
        [] // diff_variant_file: unused
    )
    versions = versions.mix(VCFTOOLS_EXCLUDE.out.versions)

    // bgzip + index the scaffolds-EXCLUDE VCF for BCFTOOLS_STATS
    BCFTOOLS_BGZIP_EXCLUDE(VCFTOOLS_EXCLUDE.out.vcf)
    versions = versions.mix(BCFTOOLS_BGZIP_EXCLUDE.out.versions)

    BCFTOOLS_INDEX_EXCLUDE(BCFTOOLS_BGZIP_EXCLUDE.out.output)
    versions = versions.mix(BCFTOOLS_INDEX_EXCLUDE.out.versions)

    ch_EXCLUDE_vcf_tbi = BCFTOOLS_BGZIP_EXCLUDE.out.output
        .join(BCFTOOLS_INDEX_EXCLUDE.out.tbi)

    // Run BCFTOOLS_STATS on scaffolds-EXCLUDE output
    BCFTOOLS_STATS_EXCLUDE(
        ch_EXCLUDE_vcf_tbi,
        [[id: 'no_regions'], []],
        [[id: 'no_targets'], []],
        [[id: 'no_samples'], []],
        [[id: 'no_exons'],   []],
        fasta
    )
    versions      = versions.mix(BCFTOOLS_STATS_EXCLUDE.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS_EXCLUDE.out.stats.map { _meta, stats -> stats })

    vcf_cleaned = VCFTOOLS_EXCLUDE.out.vcf
        .mix(intersection.passthrough)

    vcftools_thin_input = vcf_cleaned
        .map { meta, vcf_file ->
            tuple(meta + [id: meta.id + "_thinned"], vcf_file)
        }

    VCFTOOLS_THIN(
        vcftools_thin_input,
        [], // bed: unused
        []  // diff_variant_file: unused
    )
    versions = versions.mix(VCFTOOLS_THIN.out.versions)

    // bgzip + index the THIN VCF for BCFTOOLS_STATS
    BCFTOOLS_BGZIP_THIN(VCFTOOLS_THIN.out.vcf)
    versions = versions.mix(BCFTOOLS_BGZIP_THIN.out.versions)

    BCFTOOLS_INDEX_THIN(BCFTOOLS_BGZIP_THIN.out.output)
    versions = versions.mix(BCFTOOLS_INDEX_THIN.out.versions)

    ch_THIN_vcf_tbi = BCFTOOLS_BGZIP_THIN.out.output
        .join(BCFTOOLS_INDEX_THIN.out.tbi)

    // Run BCFTOOLS_STATS on THIN output
    BCFTOOLS_STATS_THIN(
        ch_THIN_vcf_tbi,
        [[id: 'no_regions'], []],
        [[id: 'no_targets'], []],
        [[id: 'no_samples'], []],
        [[id: 'no_exons'],   []],
        fasta
    )
    versions      = versions.mix(BCFTOOLS_STATS_THIN.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS_THIN.out.stats.map { _meta, stats -> stats })

    emit:
    intersection = VCFTOOLS_THIN.out.vcf
    versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def normalize_scaffold_param(param) {
    if (!param)
        return []

    // CASE 1: Array of strings
    if (param instanceof List)
        return param*.trim()

    // CASE 2: Single string
    if (param instanceof String) {
        // CASE 2a: a single string
        def f = file(param)
        if (f.exists()) {
            return f.readLines()
                .findAll { tuple -> tuple && !tuple.startsWith("#") }
                .collect { tuple -> tuple.trim().tokenize()[0] }
        }

        // CASE 2b: a single string
        return [param.trim()]
    }

    return []
}
