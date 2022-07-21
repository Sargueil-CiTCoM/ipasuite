

if RAW_DATA_TYPE == "fluo-ceq8000":
    rule fluo_ceq8000:
        input: construct_path(step="fluo-ceq8000")
        output: protected(construct_path(step="fluo-ce"))
        log: construct_path('fluo-ce', ext=".log", log_dir=True)
        message: f"Converting ceq8000 data for qushape: {MESSAGE} replicate"
                 f" {{wildcards.replicate}}"
        shell:
            f"ceq8000_to_tsv {{input}} {{output}} &> {{log}}"

if config["qushape"]["use_subsequence"]:
    rule split_fasta:
        input:
            ancient(get_refseq)
        output:
            f"{RESULTS_DIR}/{config['folders']['subseq']}/{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.fasta"
    
        message:
            f"Fragmenting fasta file {{input}} "
            f"from {{wildcards.rt_end_pos}} to {{wildcards.rt_begin_pos}}"
    
        log:
            f"logs/{config['folders']['subseq']}/{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.log"
        shell:
            f"split_fasta {{input}} {{output}} --begin "
            f"{{wildcards.rt_end_pos}}"
            f" --end {{wildcards.rt_begin_pos}}"

rule generate_project_qushape:
    input:
        rx = ancient(construct_path("fluo-ce")),
        bg = ancient(construct_path("fluo-ce", control = True )),
        refseq = ancient(lambda wildcards: get_subseq(wildcards, split_seq=True)),
        refproj = ancient(get_qushape_refproj)
    message: f"Generate QuShape project for {MESSAGE}"
             f"- replicate {{wildcards.replicate}}"
    params:
        refseq=lambda wildcards, input: f"--refseq={input.refseq}",
        refproj=lambda wildcards, input: expand('--refproj={refproj}', refproj=input.refproj)[0] if len(input.refproj) > 0 else "",
        ddNTP=lambda wildcards: f"--ddNTP={get_ddntp_qushape(wildcards)}",
        channels= construct_dict_param(config["qushape"], "channels"),
        overwrite="--overwrite=untreated"
    log: construct_path('qushape', ext=".log", log_dir=True, split_seq=True)
    output: construct_path("qushape", ext=".qushape", split_seq=True)
    shell:
        f"qushape_proj_generator {{input.rx}} {{input.bg}}"
        f" {{params}} --output={{output}} &> {{log}}"

rule extract_reactivity:
    input: 
        qushape = construct_path("qushape", ext=".qushape", split_seq=True),
        refseq = ancient(lambda wildcards: get_subseq(wildcards, split_seq=True))
    output:
        react=construct_path("reactivity", split_seq=True),
        plot=report(construct_path("reactivity", ext=".reactivity.svg",  figure=True, split_seq=True), category="3.1-Reactivity", subcategory=CONDITION),
        #,protect = protected(construct_path("qushape", ext=".qushape"))
    params:
        rna_file=lambda wildcards, input: f"--rna_file {input.refseq}" if config["qushape"]["check_integrity"] else "",
        plot_title=lambda wildcards: f"--plot_title='Reactivity of {MESSAGE.format(wildcards=wildcards)} - replicate {wildcards.replicate}'"
    message: f"Extracting reactivity from QuShape for {MESSAGE}"
             f"- replicate {{wildcards.replicate}}"
    log: construct_path('reactivity', ext=".log", log_dir=True, split_seq=True)

    shell:
        f"qushape_extract_reactivity {{input.qushape}} {{params}}"
        f" --output={{output.react}} --plot={{output.plot}} &> {{log}}"

rule normalize_reactivity:
    input: construct_path("reactivity", split_seq=True)
    output:
        nreact=construct_path("normreact", split_seq=True),
        plot=report(construct_path("normreact", ext=".normreact.svg",
            figure=True, split_seq=True) ,
                category="3.2-Normalized reactivity", subcategory=CONDITION) if
        config["normalization"]["plot"] else [],
    message: f"Normalizing reactivity for {MESSAGE}"
             f" - replicate {{wildcards.replicate}}"
    log: construct_path('normreact', ext=".log", log_dir=True, split_seq=True)
    params:
        react_nuc = get_reactive_nucleotides,
        st_perc = construct_param(CNORM, "stop_percentile"),
        low_norm_reac_thres = construct_param(CNORM, "low_norm_reactivity_threshold"),
        norm_methods = construct_list_param(CNORM, "norm_methods"),
        snorm_out_perc= construct_param(CNORM, "simple_outlier_percentile"),
        snorm_term_avg_perc= construct_param(CNORM, "simple_norm_term_avg_percentile"),
        plot =  lambda wildcards, output: f"--plot={output.plot}" if config["normalization"]["plot"] else "",
        plot_title=lambda wildcards: f"--plot_title='Normalized Reactivity of {MESSAGE.format(wildcards=wildcards)} - replicate {wildcards.replicate}'"
    shell:
        f"normalize_reactivity {{params}} {{input}}"
        f" --output={{output.nreact}}  &> {{log}}"

if config["qushape"]["use_subsequence"]:
    rule align_reactivity_to_ref:
        input: unpack(get_align_reactivity_inputs)
        output: construct_path("alignnormreact")
        params:
            rt_end_pos = get_align_begin,
            #rna_end = get_align_end

        log: construct_path('alignnormreact', ext=".log", log_dir=True)
        shell:
            f"shift_reactivity {{input.norm}} {{input.refseq}}"
            f" {{output}} --begin {{params.rt_end_pos}} &> {{log}}"
            #f" --end {{params.rna_end}} &> {{log}}"
     


rule aggregate_reactivity:
    input:
        norm= lambda wildcards: expand(construct_path(aggregate_input_type()),
                replicate=get_replicate_list(wildcards), allow_missing=True),
        refseq = lambda wildcards: get_refseq(wildcards)
    output:
        full= construct_path("aggreact", show_replicate = False),
        shape_file = construct_path("aggreact-ipanemap", show_replicate=False,
                ext=".shape"),
        map_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".map"),
        plot =report(construct_path("aggreact", ext=".aggreact.svg",
            show_replicate=False, figure=True),
            category="4-Aggregated reactivity", subcategory=CONDITION) if config["aggregate"]["plot"] else [],
        fullplot = report(construct_path("aggreact", ext=".aggreact.full.svg",
            show_replicate=False, figure=True),
            category="4-Aggregated reactivity", subcategory=CONDITION) if config["aggregate"]["plot"] else [],

    #message: f"Aggregating normalized reactivity for {MESSAGE}"
    log: construct_path('aggreact', ext=".log", log_dir=True, show_replicate=False)
    params:
        norm_method= construct_normcol(),
        minndp = construct_param(config["aggregate"], "min_ndata_perc"),
        mindndp = construct_param(config["aggregate"], "min_nsubdata_perc"),
        maxmp = construct_param(config["aggregate"], "max_mean_perc"),
        mind = construct_param(config["aggregate"], "min_dispersion"),
        plot = lambda wildcards, output: f"--plot={output.plot} --fullplot={output.fullplot}" if
        config["aggregate"]["plot"] else "",
        plot_title=lambda wildcards: f"--plot_title='Aggregated reactivity of {MESSAGE.format(wildcards=wildcards)}'"
        #refseq = lambda wildcards, input: expand('--refseq={refseq}', refseq=input.refseq)[0] if len(input.refseq) > 0 else ""
    shell:
        f"aggregate_reactivity {{input.norm}}"
        f" --output={{output.full}} {{params}}"
        f" --shape_output={{output.shape_file}}"
        f" --map_output={{output.map_file}}"
        f" --err_on_dup={config['aggregate']['err_on_dup']} &> {{log}}"



rule footprint:
    input:
        get_footprint_inputs,
    output:
        tsv=f"{RESULTS_DIR}/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}.tsv",
        plot=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}.svg", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
    log: "logs/footprint_{{rna_id}}_footprint_{foot_id}.log"
    params:
        ttest_pvalue_thres = construct_param(config["footprint"]["config"],
                "ttest_pvalue_thres"),
        deviation_type = construct_param(config["footprint"]["config"],
                "deviation_type"),
        diff_thres = construct_param(config["footprint"]["config"],
                "diff_thres"),
        ratio_thres = construct_param(config["footprint"]["config"],
                "ratio_thres"),
        cond1_name = lambda wildcards: f"--cond1_name='{get_footprint_condition_name(wildcards, 1)}'",
        cond2_name = lambda wildcards: f"--cond2_name='{get_footprint_condition_name(wildcards, 2)}'",
        plot_title = lambda wildcards: f"--plot_title='Compared reactivity between {get_footprint_condition_name(wildcards, 1)} and {get_footprint_condition_name(wildcards, 2)}'"
    shell:
        f"footprint {{input}}"
        f" --output={{output.tsv}} {{params}}"
        f" --plot={{output.plot}} --plot_format=svg "
#rule ipanemap:
#    conda: "../envs/ipanemap.yml"
#    input: construct_path("aggreact-ipanemap", replicate = False)
#    output: "python workflow/scripts/IPANEMAP/IPANEMAP.py"
#
