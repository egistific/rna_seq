rule salmon_index:
    input:
        config["data"]["ref"]["transcriptome"]
    output:
        directory(f"{OUTDIR}/salmon_index")
    params:
        kmer=config["params"]["kmer"]
    threads:
        config["params"]["threads"]["salmon"]
    conda:
        "../../envs/salmon.yaml"
    log:
        "logs/salmon_index.log"
    shell:
        "salmon index -t {input} -i {output} -k {params.kmer}  -p {threads} > {log} 2>&1"


rule salmon_quant:
    input:
        r1=f"{OUTDIR}/trimmed/{{sample}}_r1.trimmed.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_r2.trimmed.fastq.gz",
        index=f"{OUTDIR}/salmon_index"
    output:
        quant=f"{OUTDIR}/quant/{{sample}}/quant.sf"
    params:
        outdir=f"{OUTDIR}/quant/{{sample}}",
        libtype=config["params"]["salmon_libtype"]
    threads:
        config["params"]["threads"]["salmon"]
    conda:
        "../../envs/salmon.yaml"
    log:
        "logs/{sample}/salmon_quant.log"
    shell:
        """
        salmon quant -i {input.index} -l {params.libtype} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} --validateMappings \
            -o {params.outdir} \
            --gcBias > {log} 2>&1
        """