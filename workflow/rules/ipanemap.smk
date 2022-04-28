import copy
import yaml
VARNA = "workflow/scripts/VARNA/build/jar/VARNAcmd.jar"

def generate_conditions_config(pool_id, config):
    conditions = {}
    for pool in config["ipanemap"]["pools"]:
        if pool["id"] == pool_id:
            for cond in pool["conditions"]:
                react_file = expand("_" + config["format"]["condition"], **cond)
                conditions['pool_id'] = {"reactivity_file": react_file}
            break
    return conditions

def generate_ipanemap_config_file(configfile_path, pool_id):
    
    ipanemap_config = copy.deepcopy(config["ipanemap"]["config"])
    conditions = generate_conditions_config(pool_id, config)

    config["conditions"] = conditions

    with open(configfile_path, "w") as file:
        yaml.dump(ipanemap_config, file, default_flow_style=False)
        return configfile_path
    return None
    

# configfile_path: str,
# input_softdir: str,
# input_rnaseq: str,
# output_dir: str,
# input_conditions: [str],
# input_harddir: str,
# log_file: str,
# config):




rule configure_ipanemap:
    output:
        cfg=expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}.yaml",
            folder=config["folders"]["ipanemap-config"],
            allow_missing=True,
        ),
    log:
        "logs/ipanemap-config-{rna_id}_pool_{pool_id}.log",
    run:
        generate_ipanemap_config_file(configfile_path=output.cfg[0])

checkpoint ipanemap:
    conda:
        "../envs/ipanemap.yml"
    threads: 8
    input:
        config=expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}.yaml",
            folder=config["folders"]["ipanemap-config"],
            allow_missing=True,
        ),
        files=get_ipanemap_inputs,
    output:
        directory(
            expand(
                f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}",
                folder=config["folders"]["ipanemap-out"],
                allow_missing=True,
            )
        ),
    log:
        "logs/ipanemap-out-{rna_id}_pool_{pool_id}.log",
    shell:
        f"ipanemap -f {{input.config}} | tee {{log}}"


rule structure:
    conda:
        "../envs/ipanemap.yml"
    threads: 8
    input:
        expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}/"
            f"{{rna_id}}_pool_{{pool_id}}_optimal_{idx}.dbn",
            folder=config["folders"]["ipanemap-out"],
            allow_missing=True,
        ),
    output:
        expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx, \d+}}.dbn",
            folder=config["folders"]["structure"],
            allow_missing=True,
        ),
    log:
        "logs/ipanemap-{rna_id}_pool_{pool_id}_{idx}.log",
    shell:
        "cp {input} {output} &> {log}"


rule varna_color_by_condition:
    conda:
        "../envs/ipanemap.yml"
    input:
        struct=expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx}}.dbn",
            folder=config["folders"]["structure"],
            allow_missing=True,
        ),
        aggreact=construct_path("aggreact-ipanemap", show_replicate = False,
                ext=".txt")
    params:
        colorstyle= f"-colorMapStyle '{config['varna']['colormapstyle']}'",
    output:
        varna=expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx, \d+}}_cond_{{conditions}}.varna",
            folder=config["folders"]["varna"],
            conditions=CONDITION,
            allow_missing=True,
        ),
        svg=report(expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx, \d+}}_cond_{{conditions}}.svg",
            folder=config["folders"]["varna"],
            conditions=CONDITION,
            allow_missing=True,
        ), category="6.-Secondary structure", subcategory="{rna_id} - {pool_id}"),
    log:
        f"logs/varna-{{rna_id}}_pool_{{pool_id}}_{{idx}}_{CONDITION}.log",
    shell:
        f"java -jar {VARNA} -i {{input.struct}} -o {{output.varna}}" 
        f" {{params.colorstyle}} -colorMap {{input.aggreact}}"
        f" -title '{{wildcards.probe}} - {{wildcards.pool_id}}_"
        f"{{wildcards.idx}} - {{wildcards.rna_id}}' "
        f"&> {{log}};"
        f"java -jar {VARNA} -i {{input.struct}} -o {{output.svg}}" 
        f" {{params.colorstyle}} -colorMap {{input.aggreact}}"
        f" -title '{{wildcards.probe}} - {{wildcards.pool_id}}_"
        f"{{wildcards.idx}} - {{wildcards.rna_id}}' "
        f"&> {{log}};"


# Not working while Varna can't save several rna inside one session file
rule varna_pool_concat:
    conda:
        "../envs/ipanemap.yml"
    input:
        get_varna_pool_concat_inputs
    output:
        varna=expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx, \d+}}.varna",
            folder=config["folders"]["varna"],
            allow_missing=True,
        ),
    params: 
        inputs=lambda wildcards, input: "".join(f" -i {ipt}" for ipt in input)
    log:
        f"logs/varna-{{rna_id}}_pool_{{pool_id}}_{{idx}}.log",
    shell:
        f"java -jar {VARNA} {{params.inputs}} -o {{output.varna}} &> {{log}};"



#run varna_all_conditions:
