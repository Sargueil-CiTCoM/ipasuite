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


def aggreact_tsv_from_map(map_file,aggreact_tsv_file):
    import pandas as pd
    import os

    input = pd.read_csv(map_file, names=['reactivity','error','nucleotide'], sep='\t')

    map_file_base = os.path.basename(map_file)

    result = list()
    for index,row in input.iterrows():
        df = pd.DataFrame({
            'seqNum': index,
            'sequence': row['nucleotide'],
            map_file_base: row['reactivity'],
            'mean': row['reactivity'],
            'stdev': 0,
            'sem': row['error'],
            'mad': row['error'],
            'used_values': 1,
            'desc': 'accepted' if row['reactivity']>-1 else 'undetermined'
        }, index=[index])
        result.append(df)
    result_df =pd.concat(result)
    result_df.to_csv(aggreact_tsv_file, index=False, sep='\t')

rule import_external_map:
    input:
        get_external_map
    output:
        map_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".map"),
        shape_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".shape"),
        shape_IP_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".ip.shape"),
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
        "cp '{input}' '{output.map_file}' &> {log}; "
        "cut -f1,2 '{input}' >'{output.shape_file}' 2>> {log}; "
        "cp '{output.shape_file}' '{output.shape_IP_file}' &>> {log}"

rule generate_aggreact_tsv_from_map:
    input:
        external = get_external_map,
        internal = construct_path("aggreact-ipanemap", show_replicate=False, ext=".map"),
    output:
        construct_path("aggreact", show_replicate = False),
    run:
        aggreact_tsv_from_map(input.internal[0],output[0])


rule importraw:
    input:
        get_all_raw_outputs(),


rule importqushape:
    input:
        ancient(get_all_qushape_outputs()),
