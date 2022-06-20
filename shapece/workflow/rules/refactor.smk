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

