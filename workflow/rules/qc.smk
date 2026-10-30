rule fastqc:
    input:
        unpack(fastqs)
    output:
        html1=f"{OUTDIR}/qc/fastqc/{{sample}}_r1_fastqc.html",
        html2=f"{OUTDIR}/qc/fastqc/{{sample}}_r2_fastqc.html"
    params:
        outdir=f"{OUTDIR}/qc/fastqc"
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/{sample}/fastqc.log"
    shell:
        """
        mkdir -p {params.outdir}
        # symlink the inputs
        ln -sf "$(readlink -f {input.r1})" {params.outdir}/{wildcards.sample}_r1.fastq.gz
        ln -sf "$(readlink -f {input.r2})" {params.outdir}/{wildcards.sample}_r2.fastq.gz
        fastqc --outdir {params.outdir} {params.outdir}/{wildcards.sample}_r1.fastq.gz {params.outdir}/{wildcards.sample}_r2.fastq.gz > {log} 2>&1
        """


rule fastp:
    input:
        unpack(fastqs)
    output:
        # r1 and r2 are temporary files that are deleted after they have been used to save space 
        r1=temp(f"{OUTDIR}/trimmed/{{sample}}_r1.trimmed.fastq.gz"),
        r2=temp(f"{OUTDIR}/trimmed/{{sample}}_r2.trimmed.fastq.gz"),
        html=f"{OUTDIR}/qc/fastp/{{sample}}.fastp.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}.fastp.json"
    threads: config["params"]["threads"]["fastp"]
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/{sample}/fastp.log"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json} --thread {threads} > {log} 2>&1
        """


rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastqc/{{sample}}_r1_fastqc.html", sample=SAMPLES),
        expand(f"{OUTDIR}/qc/fastqc/{{sample}}_r2_fastqc.html", sample=SAMPLES),
        expand(f"{OUTDIR}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
        expand(f"{OUTDIR}/quant/{{sample}}/quant.sf", sample=SAMPLES)
    output:
        html=f"{OUTDIR}/multiqc/multiqc_report.html"
    params:
        searchdir=f"{OUTDIR}",
        outdir=f"{OUTDIR}/multiqc"
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc {params.searchdir} --outdir {params.outdir} --force > {log} 2>&1
        """
