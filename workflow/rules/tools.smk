rule import_raw_probe_data:
    input: get_raw_probe_input 
    output: expand("resources/{fluodir}/{rna_id}" + config["format"]["condition"] + "_{replicate}.{fluoext}.tsv", fluodir = config["folders"][config["rawdata"]["type"]], fluoext=config["rawdata"]["type"], allow_missing=True)
    shell:
        "cp {input} {output}"

rule import_raw_control_data:
    input: get_raw_control_input 
    output:
        expand("resources/{fluodir}/{rna_id}" + config["format"]["control_condition"] + "_{replicate}.{fluoext}.tsv", fluodir = config["folders"][config["rawdata"]["type"]],fluoext=config["rawdata"]["type"], allow_missing=True)
    shell:
        "cp {input} {output}"

rule fluo_ceq8000:
    conda: "workflow/envs/tools.yml"
    input: expand("ressources/{folder}/{sample}.fluo-ceq8000.tsv", folder=config["folders"]["fluo-ceq8000"], allow_missing=True)
    output: expand("resources/{folder}/{sample}.tsv", folder=config["folders"]["fluo-ce"], allow_missing=True)
    shell:
        "python scripts/tools/ceq8000_to_tsv.py {input} {output}"

rule generate_project_qushape:
    conda: "workflow/envs/tools.yml"
    input:
        rx = expand("resources/{folder}/{rna_id}" + config["condition"] + "_{replicate}.fluo-ce.tsv", folder = config["folders"]["fluo-ce"], allow_missing=True),
        bg = expand("resources/{folder}/{rna_id}" + config["control_condition"] + "_{replicate}.fluo-ce.tsv", folder = config["folders"]["fluo-ce"], control = config["rowdata"]["control"], allow_missing=True),
        refseq = get_qushape_refseq, 
        refproj = get_qushape_refproj
    params:
        refseq = lambda input, wildcards: '--refseq {input.refseq}' if wildcards.refseq,
        refproj = lambda input, wildcards: '--qushape_refproj {input.refproj}' if wildcards.refproj
    output: expand("results/{folder}/{rna_id}" + config["format"]["condition"] + "_{replicate}.qushape.tsv", folder = config["folders"]["qushape"] , allow_missing=True)
    shell:
        "python scripts/tools/qushape_proj_generator.py {input.rx} {input.bg} {params.refseq} {params.refproj} --output {output}"
