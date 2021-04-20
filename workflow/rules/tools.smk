rule ceq8000_to_tsv:
    conda: "workflow/envs/tools.yml"
    input:
        "resources/ceq8000/{folder}/{sample}.txt"
    output:
        "resources/ceq8000tsv/{folder}/{sample}.tsv"
    shell:
        "python scripts/tools/ceq8000_to_tsv.py {input} {output}"
