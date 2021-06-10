CONDITION = config["format"]["condition"]
CTRL_CONDITION = config["format"]["control_condition"]
RAW_DATA_TYPE = config["rawdata"]["type"]
CONTROL = config["rawdata"]["control"]
FOLDERS = config["folders"]
MESSAGE = config["format"]["message"]
TOOLS = "workflow/scripts/tools/"


rule import_raw_probe_data:
    input: get_raw_probe_input
    output: construct_path(RAW_DATA_TYPE, results_dir = False)
    shell:
        "cp {input} {output}"

rule import_raw_control_data:
    input: get_raw_control_input 
    output: construct_path(RAW_DATA_TYPE, control = True, results_dir = False)
    shell:
        "cp {input} {output}"

rule fluo_ceq8000:
    conda: "../scripts/tools/conda-env.yml"
    input: construct_path(step="fluo-ceq8000", results_dir=False) 
    output: construct_path(step="fluo-ce", results_dir=False) 
    shell:
        "python " + TOOLS + "ceq8000_to_tsv.py {input} {output}"

rule generate_project_qushape:
    conda: "../scripts/tools/conda-env.yml"
    input:
        rx = construct_path("fluo-ce", results_dir = False),
        bg = construct_path("fluo-ce", control = True, results_dir = False),
        refseq = get_qushape_refseq, 
        refproj = get_qushape_refproj
    params:
        refseq=lambda wildcards, input: expand('--refseq {refseq}', refseq=input.refseq) if not input.refseq is [] else "",
        refproj=lambda wildcards, input: expand('--refproj {refproj}', refproj=input.refproj) if not input.refproj is [] else "",
        ddNTP= "--ddNTP " + config["qushape"]["ddNTP"] if "ddNTP" in config["qushape"] else ""
        #TODO channels
        #channels=  "--channels" + config["qushape"]["ddNTP"] if "ddNTP" in config["qushape"] else ""
    #output: protected(construct_path("qushape", ext=".qushape"))
    output: construct_path("qushape", ext=".qushape")
    shell:
        "python " + TOOLS + "qushape_proj_generator.py {input.rx} {input.bg} {params.refseq} {params.refproj} {params.ddNTP} --output {output}"

rule extract_reactivity:
    conda:  "../scripts/tools/conda-env.yml"
    input: construct_path("qushape", ext=".qushape")
    output: construct_path("reactivity") 
    message: "Extracting reactivity from QuShape for" + MESSAGE + " - replicate {wildcards.replicate}"
    shell:
        "python " + TOOLS + "qushape_extract_reactivity.py {input} --output {output}"

rule normalize_reactivity:
    conda:  "../scripts/tools/conda-env.yml"
    input: construct_path("reactivity")
    output: construct_path("normreact")
    message: "Normalizing reactivity for" + MESSAGE + " - replicate {wildcards.replicate}"
    shell:
        "python " + TOOLS + "normalize_reactivity.py {input} {output}"

rule aggregate_reactivity:
    conda:  "../scripts/tools/conda-env.yml"
    input: lambda wildcards: expand(construct_path("normreact"), replicate=get_replicates(wildcards, qushape_analysed = True), allow_missing=True)
    output: 
        full= construct_path("aggreact", replicate = False), 
        compact = construct_path("aggreact-ipanemap", replicate=False, ext=".txt")
    message: "Aggregating normalized reactivity for " + MESSAGE
    shell:
        "python "+ TOOLS + "aggregate_reactivity.py {input} {output.full} --ipanemap_output {output.compact}"

#rule ipanemap:
#    conda: "../envs/ipanemap.yml"
#    input: construct_path("aggreact-ipanemap", replicate = False)
#    output: "python workflow/scripts/IPANEMAP/IPANEMAP.py"
#    
