# configfile: "config.yaml"

rule firststep:
    input:
        fasta="../data/conserved/{subtype}_CR.fasta",
        csv="../data/variable/{subtype}_VR.csv"
    output:
        fasta="../data/output/conserved/{subtype}_filtered.fasta",
        csv="../data/output/variable/{subtype}_filtered.csv",
    shell:
        "python {params.scripts}/filter.py {input.fasta} {input.csv} {output.fasta} {output.csv}"