if config["allow_auto_import"] and not ('refactor_rename' in config and
        config['refactor_rename']) and not ('refactor_addpositions' in config and
        config['refactor_addpositions']):
    if "fluo-ce" in RAW_DATA_TYPE:
        rule import_raw_control_data:
            input:
                get_raw_control_input,
            output:
                protected(construct_path(RAW_DATA_TYPE, control=True, split_seq=False)),
            log:
                construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, control=True, split_seq=False),
            message:
                (
                    "Importing raw data from external source: "
                    + MESSAGE
                    + " replicate {wildcards.replicate}"
                )
            shell:
                "cp '{input}' '{output}' &> {log}"

        rule import_raw_probe_data:
            input:
                get_raw_probe_input,
            output:
                protected(construct_path(RAW_DATA_TYPE, split_seq=False)),
            log:
                construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, split_seq=False),
            message:
                (
                    "Importing raw data from external source: "
                    + MESSAGE
                    + " replicate {wildcards.replicate}"
                )
            shell:
                " cp '{input}' '{output}' &> {log}"


    if "fluo-fsa" in RAW_DATA_TYPE:
        rule import_raw_control_data:
            input:
                get_raw_control_input,
            output:
                protected(construct_path(RAW_DATA_TYPE, control=True, split_seq=False, ext=".fsa")),
            log:
                construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, control=True, split_seq=False),
            message:
                (
                    "Importing raw data from external source: "
                    + MESSAGE
                    + " replicate {wildcards.replicate}"
                )
            shell:
                "cp '{input}' '{output}' &> {log}"

        rule import_raw_probe_data:
            input:
                get_raw_probe_input,
            output:
                protected(construct_path(RAW_DATA_TYPE, split_seq=False, ext=".fsa")),
            log:
                construct_path(RAW_DATA_TYPE, ext=".log", log_dir=True, split_seq=False),
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
            get_external_qushape
        output:
            construct_path("qushape", ext=".qushapey", split_seq=True),
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


rule import_external_map:
    input:
        get_external_map
    output:
        full= construct_path("aggreact", show_replicate = False),
        shape_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".shape"),
        shape_IP_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".ip.shape"),
        map_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".map"),
        #relation_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".csv"),
        #plot =report(construct_path("aggreact", ext=".aggreact.svg",
        #    show_replicate=False, figure=True),
        #    category="4-Aggregated reactivity", subcategory=CONDITION) if config["aggregate"]["plot"] else [],
        #fullplot = report(construct_path("aggreact", ext=".aggreact.full.svg",
        #    show_replicate=False, figure=True),
        #    category="4-Aggregated reactivity", subcategory=CONDITION) if config["aggregate"]["plot"] else [],
    log:
        construct_path("aggreact-ipanemap", ext=".log", log_dir=True, split_seq=True, show_replicate=False),
    message:
        (
            "Importing from external map file: "
            + MESSAGE
        )
    shell:
        "touch {output.full} &> {log}; "
        "cp '{input}' '{output.map_file}' &>> {log}; "
        "cut -f1,2 '{input}' >'{output.shape_file}' 2>> {log}; "
        "cp '{output.shape_file}' '{output.shape_IP_file}' &>> {log}"


rule importraw:
    input:
        get_all_raw_outputs(),


rule importqushape:
    input:
        ancient(get_all_qushape_outputs()),
