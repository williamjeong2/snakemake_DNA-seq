rule visualization:
    input:
        calls="annotated/all.vcf.gz"
    output:
        done = "visualization/idone"
    params:
        ref = config["ref"]["species"],
        outdir = directory("visualization/")
    log:
        "logs/visualization.log"
    shell:
        "rm {output.done} && Rscripts scripts/mutationalParrerns.r --input {input.calls} --ref {params.ref} --outdir {params.outdir} && touch {output.done}"