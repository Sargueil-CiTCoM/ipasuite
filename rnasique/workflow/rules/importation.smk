if config["allow_auto_import"] and not ('refactor_rename' in config and
        config['refactor_rename']) and not ('refactor_addpositions' in config and
        config['refactor_addpositions']):

    rule import_raw_probe_data:
        input:
            get_raw_probe_input,
        output:
            protected(construct_path(RAW_DATA_TYPE, split_seq=True)),
        log:
            construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, split_seq=True),
        message:
            (
                "Importing raw data from external source: "
                + MESSAGE
                + " replicate {wildcards.replicate}"
            )
        shell:
            "cp '{input}' '{output}' &> {log}"

    rule import_raw_control_data:
        input:
            get_raw_control_input,
        output:
            protected(construct_path(RAW_DATA_TYPE, control=True, split_seq=True)),
        log:
            construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, control=True, split_seq=True),
        message:
            (
                "Importing raw data from external source: "
                + MESSAGE
                + " replicate {wildcards.replicate}"
            )
        shell:
            "cp '{input}' '{output}' &> {log}"

    rule import_external_qushape:
        input:
            get_external_qushape,
        output:
            construct_path("qushape", ext=".qushape", split_seq=True),
        log:
            construct_path("qushape", ext=".log", log_dir=True, split_seq=True),
        message:
            (
                "Importing from external source: "
                + MESSAGE
                + " replicate {wildcards.replicate}"
            )
        shell:
            "cp '{input}' '{output}' &> {log}"

    ruleorder: import_external_qushape > generate_project_qushape


rule importraw:
    input:
        get_all_raw_outputs(),


rule importqushape:
    input:
        ancient(get_all_qushape_outputs()),
    
