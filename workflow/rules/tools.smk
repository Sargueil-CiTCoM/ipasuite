rule ceq8000_to_tsv:
    conda: "workflow/envs/tools.yml"
    input:
        "resources/ceq8000/{folder}/{samplerx}.txt"
        "resources/ceq8000/{folder}/{samplebg}.txt"
    output:
        "resources/fluorescence/{folder}/{sample}.tsv"
    shell:
        "python scripts/tools/ceq8000_to_tsv.py {input} {output}"

rule generate_project_qushape:
    conda: "workflow/envs/tools.yml"
    input:
        rx = "resources/fluorescence/{folder}/{samplerx}.tsv"
        bg = "resources/fluorescence/{folder}/{samplebg}.tsv"
        refseq = "resources/sequences/{refseq}.fa"
        refproj "resources/qushape/{refproj}.qushape"
    params:
        refseq = lambda input, wildcards: '--refseq {input.refseq}' if wildcards.refseq 
        refproj = lambda input, wildcards, '--qushape_refproj {input.refproj}' if wildcards.refproj
    output:
        "resources/qushape/{folder}/{sample}.qushape"
    shell:
        "python scripts/tools/qushape_proj_generator.py {input.rx} {input.bg} {params.refseq} {params.refproj} --output {output}"


