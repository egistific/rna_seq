rule reformat_transcriptome:
    input:
        raw=config["data"]["ref"]["raw"]
    output:
        transcriptome=config["data"]["ref"]["transcriptome"],
        tx2gene=config["data"]["ref"]["tx2gene"]
    log:
        "logs/reformat_transcriptome.log"
    shell:
        """
        python workflow/scripts/reformat_transcriptome.py \
            --raw {input.raw} --transcriptome {output.transcriptome} --tx2gene {output.tx2gene} > {log} 2>&1
        """