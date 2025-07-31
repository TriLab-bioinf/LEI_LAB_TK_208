# Run multiqc report

rule multiqc:
    input:
        i1 = lambda wildcards: expand("results/1-trim/{s}.P.R{r}.fastq.gz", s = samples, r = [1,2]) if (config["trimming"] == True) else []
    output: "results/3-multiqc/multiqc_report.html"
    resources:
        partition = "quick",
        mem_mb = 4000
    benchmark:
        "benchmarks/3-multiqc/multiqc.tsv"
    log:
        logfile = "logs/3-multiqc/multiqc.log"
    shell:
        """
        module load multiqc

        multiqc -f -d -o results/3-multiqc results/ > {log.logfile} 2>&1
        """