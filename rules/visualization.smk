rule visualization:
    input:
        calls="annotated/all.vcf.gz"
    output:

    params:
        ref = config["ref"]["species"],
        outdir = directory("visualization/")
    log:
        "logs/visualization.log"
    shell:
        "Rscripts scripts/mutationalParrerns.r --input {input.calls} --ref {params.ref} --outdir {params.outdir}"