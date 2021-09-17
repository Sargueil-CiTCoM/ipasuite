if config["allow_auto_import"]:
    rule import_raw_probe_data:
        input: get_raw_probe_input
        output: protected(construct_path(RAW_DATA_TYPE, results_dir = False))
        log: construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True) 
        message: "Importing raw data from external source: " + MESSAGE + " replicate {wildcards.replicate}"
        shell:
            "cp {input} {output} &> {log}"
    
    rule import_raw_control_data:
        input: get_raw_control_input 
        output: protected(construct_path(RAW_DATA_TYPE, control = True, results_dir = False))
        log: construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, control= True) 
        message: "Importing raw data from external source: " + MESSAGE + " replicate {wildcards.replicate}"
        shell:
            "cp {input} {output} &> {log}"
    rule import_external_qushape:
        input: get_external_qushape
        output: protected(construct_path("qushape", ext=".qushape", results_dir = True))
        log: construct_path('qushape', ext=".log", log_dir=True) 
        message: "Importing from external source: " + MESSAGE + " replicate {wildcards.replicate}"
        shell:
            "cp {input} {output} &> {log}"

    ruleorder: import_external_qushape > generate_project_qushape

rule importraw:
    input: get_all_raw_outputs()

rule importqushape:
    input: ancient(get_all_qushape_outputs())