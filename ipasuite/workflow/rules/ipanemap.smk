import copy
import yaml
import json

def generate_conditions_config(pool_id, config):
    conditions = {}
    rna_id = None
    for pool in config["ipanemap"]["pools"]:
        if pool["id"] == pool_id:
            rna_id = pool["rna_id"]
            if "conditions" in pool:
                for cond in pool["conditions"]:
                    cond_name = expand(config["format"]["condition"], **cond)[0]
                    react_file = expand(construct_path("aggreact-ipanemap", show_replicate=False,
                    ext=".shape"), rna_id=rna_id, **cond)[0]
                    conditions[cond_name] = {"reactivity_file": react_file}
            if "external_conditions" in pool:
                for cond in pool["external_conditions"]:
                    cond_name = cond["name"]
                    react_file = cond["path"]
                    conditions[cond_name] = {"reactivity_file": react_file}
            if "alignements" in pool:
                print(pool["alignements"])
                for name, file in pool["alignements"].items():
                    conditions[name] = {"alignement_file": file}
            break
    return conditions, rna_id

def generate_ipanemap_config_file(configfile_path, pool_id, output_dir):

    ipanemap_config = json.loads(json.dumps(config["ipanemap"]["config"]))
    conditions, rna_id = generate_conditions_config(pool_id, config)
    ipanemap_config["conditions"] = conditions
    ipanemap_config["output_dir"] = output_dir
    ipanemap_config["sequence_file"] = config["sequences"][rna_id]

    ipanemap_config["tmp_dir"] = f"{output_dir}/tmp"
    ipanemap_config["format"] = {
            'dbn_file_pattern':
            f"{output_dir}/{rna_id}_pool_{pool_id}_optimal_{{idx}}.dbn",
            'dbn_centroid_file_pattern':
            f"{output_dir}/{rna_id}_pool_{pool_id}_centroid_{{idx}}.dbn",
            }




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
        "results/logs/ipanemap-config-{rna_id}_pool_{pool_id}.log",
    run:
        generate_ipanemap_config_file(
            configfile_path=output.cfg[0],
            pool_id=wildcards.pool_id,
            output_dir=expand(
                f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}",
                folder=config["folders"]["ipanemap-out"],
                rna_id=wildcards.rna_id,
                pool_id=wildcards.pool_id,
                allow_missing=True,
            )[0])

checkpoint ipanemap:
    #conda:
    #    "../envs/ipanemap.yml"
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
        "results/logs/ipanemap-out-{rna_id}_pool_{pool_id}.log",
    shell:
        f"ipanemap -f {{input.config}} --log {{log}}"

rule structure:
    """
    Copy Ipanemap's best and second best predicted cluster centroid to
    the ipanemap output folder
    """
    #conda:
    #    "../envs/ipanemap.yml"
    threads: 8
    input:
        dir = expand(
                f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}",
                folder=config["folders"]["ipanemap-out"],
                allow_missing=True,
            ),
    params:
        input_folder=config["folders"]["ipanemap-out"],

    output:
        expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{k}}.dbn",
            folder=config["folders"]["structure"],
            k = [1,2],
            allow_missing=True,),

    log:
        log = "results/logs/ipanemap-{rna_id}_pool_{pool_id}.log",

    run:
        import os
        import re
        import shutil

        def gen_centroid_name(k,tag="_centroid"):
            return f"{wildcards.rna_id}_pool_{wildcards.pool_id}{tag}_{k}.dbn"

        def kbest_centroids(dir, k):
            import glob
            centroids = list(glob.glob(os.path.join(dir,gen_centroid_name("*"))))
            pc_pairs=list()
            for i in range(1,len(centroids)+1):
                centroid = gen_centroid_name(i)
                with open(os.path.join(dir,centroid)) as fh:
                    probability = float(fh.readline().strip().split(' ')[-1])
                    pc_pairs.append((probability,i))
            pc_pairs = sorted(pc_pairs,key=lambda x:-x[0])
            print(pc_pairs)
            return [ c for (p,c) in pc_pairs[:k] ]

        source_dir = str(input.dir[0])
        target_dir = os.path.join(RESULTS_DIR,config["folders"]["structure"])
        for i,j in enumerate(kbest_centroids(source_dir,2)):
            shutil.copy(os.path.join(source_dir,gen_centroid_name(j)),os.path.join(target_dir,gen_centroid_name(i+1,tag="")))

rule varna_color_by_condition:
    #conda:
    #    "../envs/ipanemap.yml"
    input:
        struct_optimal = expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_1.dbn",
            folder=config["folders"]["structure"],
            allow_missing=True,
        ),
        struct_2 = expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_2.dbn",
            folder=config["folders"]["structure"],
            allow_missing=True,
        ),
        aggreact=construct_path("aggreact-ipanemap", show_replicate = False,
                ext=".shape"),
    params:
        colorstyle= f"-colorMapStyle '{config['varna']['colormapstyle']}'",
    output:
        varna_optimal = expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_1_cond_{{conditions}}.varna",
            folder=config["folders"]["varna"],
            conditions=CONDITION,
            allow_missing=True,
        ),
        varna_2 = expand(
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_2_cond_{{conditions}}.varna",
            folder=config["folders"]["varna"],
            conditions=CONDITION,
            allow_missing=True,
        ),
        svg_optimal = report(expand(
            f"{RESULTS_DIR}/figures/{{folder}}/{{rna_id}}_pool_{{pool_id}}_1_cond_{{conditions}}.svg",
            folder=config["folders"]["varna"],
            conditions=CONDITION,
            allow_missing=True,
        ), category="6.-Secondary structure", subcategory="{rna_id} - {pool_id}"),
        svg_2 = report(expand(
            f"{RESULTS_DIR}/figures/{{folder}}/{{rna_id}}_pool_{{pool_id}}_2_cond_{{conditions}}.svg",
            folder=config["folders"]["varna"],
            conditions=CONDITION,
            allow_missing=True,
        ), category="6.-Secondary structure", subcategory="{rna_id} - {pool_id}"),
    log:
        log1 = f"results/logs/varna-{{rna_id}}_pool_{{pool_id}}_1_{CONDITION}.log",
        log2 = f"results/logs/varna-{{rna_id}}_pool_{{pool_id}}_2_{CONDITION}.log",
    shell:
        f"varna -i {{input.struct_optimal}} -o {{output.varna_optimal}}"
        f" {{params.colorstyle}} -colorMap {{input.aggreact}}"
        f" -title '{{wildcards.probe}} - {{wildcards.pool_id}} - {{wildcards.rna_id}}_optimal' &> {{log.log1}};\n"
        f"varna -i {{input.struct_optimal}} -o {{output.svg_optimal}}"
        f" {{params.colorstyle}} -colorMap {{input.aggreact}}"
        f" -title '{{wildcards.probe}} - {{wildcards.pool_id}} - {{wildcards.rna_id}}_optimal' &> {{log.log1}};\n"
        f"varna -i {{input.struct_2}} -o {{output.varna_2}}"
        f" {{params.colorstyle}} -colorMap {{input.aggreact}}"
        f" -title '{{wildcards.probe}} - {{wildcards.pool_id}} - {{wildcards.rna_id}}_centrioid' &> {{log.log2}};\n"
        f"varna -i {{input.struct_2}} -o {{output.svg_2}}"
        f" {{params.colorstyle}} -colorMap {{input.aggreact}}"
        f" -title '{{wildcards.probe}} - {{wildcards.pool_id}} - {{wildcards.rna_id}}_centrioid' &> {{log.log2}};"

# Not working while Varna can't save several rna inside one session file
rule varna_pool_concat:
    #conda:
    #    "../envs/ipanemap.yml"
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
        "results/logs/varna-{rna_id}_pool_{pool_id}_{idx}.log",
    shell:
        "varna {params.inputs} -o {output.varna} &> {log};"



#run varna_all_conditions:
