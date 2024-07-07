if "refactor_addpositions" in config and config["refactor_addpositions"]:

    rule add_pos_ceq8000:
        input: construct_path(step="fluo-ceq8000", results_dir=False)
        output: construct_path(step="fluo-ceq8000", results_dir=False, force_split_seq=True)
        log: construct_path("fluo-ceq8000", ext=".log", log_dir=True, force_split_seq=True)
        message: f"Adding rt_position to ceq8000 filename"
        shell: "mv -n {input} {output} &> {log}"

    rule add_pos_fluo_ce:
        input: construct_path(step="fluo-ce", results_dir=False)
        output: construct_path(step="fluo-ce", results_dir=False, force_split_seq=True)
        log: construct_path("fluo-ce", ext=".log", log_dir=True, force_split_seq=True)
        message: f"Adding rt_position to fluo-ce filename"
        shell: "mv -n {input} {output} &> {log}"

    rule add_pos_qushape:
        input: construct_path("qushape", ext=".qushape")
        output: construct_path("qushape", ext=".qushape", force_split_seq=True)
        log:construct_path("qushape", ext=".log", log_dir=True, force_split_seq=True)
        message: f"Adding rt_position to qushape filename"
        shell: "mv -n {input} {output} &> {log}"

if "refactor_rename" in config and config["refactor_rename"]:
    rule rename_ceq8000:
        #input: lambda wildcards, output: expand(construct_path(step="fluo-ceq8000", results_dir=False, previous=True), **get_additionnal_wildcards(wildcards), allow_missing=True)
        #input:  construct_path("fluo-ceq8000", merged_conditions=True, previous=True )
        output: construct_path("fluo-ceq8000",results_dir=False )
        log: construct_path("fluo-ceq8000", ext=".log", log_dir=True)
        params:
            config_file = f"{workflow.configfiles[0]}",
            full_condition_format = lambda wildcards: construct_full_condition_path(wildcards, step="fluo-ceq8000", results_dir=False), 
            full_prev_condition_format = lambda wildcards: construct_full_condition_path(wildcards, previous=True, step="fluo-ceq8000", results_dir=False), 
            rna_id = lambda wildcards: f"--rna_id {wildcards.rna_id}",
            replicate = lambda wildcards: f"--replicate {wildcards.replicate}",
            #conditions = lambda wildcards: f"--conditions {wildcards.conditions}"
            #config_file = f"--config_file={workflow.configfiles[0]}"
        message: f"Rename ceq8000 file"
        #shell: "shapece_renamer --new-file-path {output} {params} &> {log}"
        shell: "shapece_renamer {output} {params} &> {log}"

    rule rename_fluo_ce:
        #input: construct_path(step="fluo-ce", results_dir=False, previous=True)
        output: construct_path(step="fluo-ce", results_dir=False, )
        log: construct_path("fluo-ce", ext=".log", log_dir=True, )
        params:
            config_file = f"{workflow.configfiles[0]}",

            full_condition_format = lambda wildcards: construct_full_condition_path(wildcards, step="fluo-ce", results_dir=False), 
            full_prev_condition_format = lambda wildcards: construct_full_condition_path(wildcards, previous=True, step="fluo-ce", results_dir=False),
            rna_id = lambda wildcards: f"--rna_id {wildcards.rna_id}",
            replicate = lambda wildcards: f"--replicate {wildcards.replicate}",
            #conditions = lambda wildcards: f"--conditions {wildcards.conditions}"
            #config_file = f"--config_file={workflow.configfiles[0]}"
        message: f"Rename fluo-ce file"
        #shell: "shapece_renamer --new-file-path {output} {params} &> {log}"
        shell: "shapece_renamer {output} {params} &> {log}"

    rule rename_qushape:
        #input: lambda wildcards: expand(construct_path("qushape", ext=".qushape", previous=True),**get_additionnal_wildcards(wildcards), allow_missing=True)
        output: construct_path("qushape", ext=".qushape")
        params:
            config_file = f"{workflow.configfiles[0]}",
            full_condition_format = lambda wildcards: construct_full_condition_path(wildcards, step="qushape", ext=".qushape" ), 
            full_prev_condition_format = lambda wildcards: construct_full_condition_path(wildcards, previous=True, step="qushape", ext=".qushape"), 
            rna_id = lambda wildcards: f"--rna_id {wildcards.rna_id}",
            replicate = lambda wildcards: f"--replicate {wildcards.replicate}",
            #conditions = lambda wildcards: f"--conditions {wildcards.conditions}"
            #config_file = f"--config_file={workflow.configfiles[0]}"
        log:construct_path("qushape", ext=".log", log_dir=True)
        message: f"Rename qushape filename"
        shell: "shapece_renamer --new-file-path {output} {params} &> {log}"
