rule deseq2:
    input:
        quant=expand(f"{OUTDIR}/quant/{{sample}}/quant.sf", sample=SAMPLES),
        tx2gene=config["data"]["ref"]["tx2gene"],
        samples=config["data"]["samples"]
    output:
        results=f"{OUTDIR}/deseq2/results.tsv",
        ma_plot=f"{OUTDIR}/deseq2/ma_plot.png",
        pca_plot=f"{OUTDIR}/deseq2/pca_plot.png",
        significant_genes=f"{OUTDIR}/deseq2/significant_genes.tsv",
        validation=f"{OUTDIR}/deseq2/validation.txt"
    params:
        quantdir=f"{OUTDIR}/quant",
        formula=config["study_design"]["formula"],
        condition_col=config["study_design"]["condition_col"],
        reference_level=config["study_design"]["reference_level"],
        fdr=config["params"]["fdr_threshold"],
        min_count=config["params"]["min_count"],
        lfc_threshold=config["params"]["lfc_threshold"],
        known_de_genes=",".join(config["params"]["known_de_genes"])
    conda:
        "../../envs/deseq2.yaml"
    log:
        "logs/deseq2.log"
    shell:
        """
        Rscript workflow/scripts/run_deseq2.R \
            --quantdir {params.quantdir} \
            --samples {input.samples} \
            --tx2gene {input.tx2gene} \
            --formula '{params.formula}' \
            --condition-col {params.condition_col} \
            --reference-level {params.reference_level} \
            --fdr {params.fdr} \
            --min-count {params.min_count} \
            --lfc-threshold {params.lfc_threshold} \
            --known-de-genes {params.known_de_genes} \
            --results {output.results} \
            --ma-plot {output.ma_plot} \
            --pca-plot {output.pca_plot} \
            --significant-genes {output.significant_genes} \
            --validation {output.validation} > {log} 2>&1
        """


rule ora:
    input:
        results=f"{OUTDIR}/deseq2/results.tsv"
    output:
        results=f"{OUTDIR}/ora/ora_results.tsv"
    params:
        outdir=f"{OUTDIR}/ora",
        fdr=config["params"]["fdr_threshold"]
    conda:
        "../../envs/deseq2.yaml"
    log:
        "logs/ora.log"
    shell:
        """
        Rscript workflow/scripts/run_ora.R \
            --results {input.results} \
            --fdr {params.fdr} \
            --outdir {params.outdir} > {log} 2>&1
        """

