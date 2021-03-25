include: "rules/common.smk"


##### Target rules #####


rule all:
    input:
        "annotated/all.vcf.gz",
        "qc/multiqc.html",
        "plots/depths.png",
        "plots/allele-freqs.png",


##### Modules #####


include: "rules/ref.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
include: "rules/visualization.smk"
