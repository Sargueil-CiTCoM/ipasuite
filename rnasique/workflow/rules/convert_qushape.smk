rule qushape_yamlify:
    conda: "../envs/qushape.yml"
    input: ancient(construct_path("qushape", ext=".qushape", split_seq=True))
    output: construct_path("qushape", ext=".qushapey", split_seq=True)
    log: construct_path('qushape', ext="_qushapey.log", log_dir=True, split_seq=True)
    priority: 50
    message: f"Converting format to Yaml Qushape file for {MESSAGE}"
             f"- replicate {{wildcards.replicate}}"
    shell:
        f"qushape_yamlify {{input}} {{output}} &> {{log}}"
