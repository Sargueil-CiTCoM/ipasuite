def generate_ipanemap_config_from_snakeconf(
        configfile_path: str,
        input_softdir: str,
        input_rnaseq: str,
        output_dir: str,
        input_harddir: str,
        log_file: str,
        config):

    generate_ipanemap_config(
        configfile_path = configfile_path, 
        input_softdir = input_softdir, 
        input_rnaseq = input_rnaseq, 
        output_dir = output_dir,
        input_conditions = config["ipanemap"]["conditions"], 
        input_harddir = input_harddir, 
        log_file = log_file,
        sampling= config["ipanemap"]["sampling"]["enable"],
        sampling_nstructures = config["ipanemap"]["sampling"]["nstructures"],
        sampling_temperature = config["ipanemap"]["sampling"]["temperature"], 
        sampling_m = config["ipanemap"]["sampling"]["m"],
        sampling_b = config["ipanemap"]["sampling"]["b"],
        clustering_max_diameter_thres = config["ipanemap"]["clustering"]["max_diameter_threshold"],
        clustering_max_avg_diameter_thres = config["ipanemap"]["clustering"]["max_avg_diameter_threshold"],
        pareto_perc = config["ipanemap"]["pareto"]["percentage"],
        pareto_zcutoff = config["ipanemap"]["pareto"]["zcutoff"],
        draw_models = config["ipanemap"]["visualization"]["draw_models"],
        draw_centroids = config["ipanemap"]["visualization"]["draw_centroids"],
        show_probing = config["ipanemap"]["visualization"]["show_probing"]) 

def generate_ipanemap_config(
        configfile_path: str
        input_softdir: str,
        input_rnaseq: str,
        output_dir: str,
        input_conditions: [str],
        input_harddir: str = None,
        log_file: str = "ipanemap.log" 
        sampling: bool = True,
        sampling_nstructures: int = 1000,
        sampling_temperature: float = 37,
        sampling_m: float = 1.3,
        sampling_b: float = -0.4,
        clustering_max_diameter_thres: int = 7,
        clustering_max_avg_diameter_thres: float = 7,
        pareto_perc: float = 20,
        pareto_zcutoff: float = 0.05,
        draw_models: bool = True,
        draw_centroids: bool = True,
        show_probing: bool = True
        ):
    """generate_ipanemap_config.
    :param input_softdir:
    :type input_softdir: str
    :param output_dir:
    :type output_dir: str
    :param input_conditions:
    :type input_conditions: [str]
    :param input_harddir:
    :type input_harddir: str
    :param log_file:
    :type log_file: str
    :param sampling:
    :type sampling: bool
    :param sampling_nstructures:
    :type sampling_nstructures: int
    :param sampling_temperature:
    :type sampling_temperature: float
    :param sampling_m:
    :type sampling_m: float
    :param sampling_b:
    :type sampling_b: float
    :param clustering_max_diameter_thres:
    :type clustering_max_diameter_thres: int
    :param clustering_max_avg_diameter_thres:
    :type clustering_max_avg_diameter_thres: float
    :param pareto_perc:
    :type pareto_perc: float
    :param pareto_zcutoff:
    :type pareto_zcutoff: float
    :param draw_models:
    :type draw_models: bool
    :param draw_centroids:
    :type draw_centroids: bool
    :param show_probing:
    :type show_probing: bool
    """
    ipanemap_conf = configparser.ConfigParser()
    ipanemap_conf['Input'] = {}
    ipanemap_conf['Input']['HardConstraintsDir'] = input_harddir
    ipanemap_conf['Input']['SoftConstraintsDir'] = input_softdir
    ipanemap_conf['Input']['Conditions'] = input_conditions 
    ipanemap_conf['Input']['RNA'] = input_rnaseq

    ipanemap_conf['Paths'] = {}
    ipanemap_conf['Paths']['WorkingDir'] = output_dir
    ipanemap_conf['Paths']['LogFile'] = log_file

    ipanemap_conf['Sampling'] = {}
    ipanemap_conf['Sampling']['DoSampling'] = sampling 
    ipanemap_conf['Sampling']['NumStructures'] = sampling_nstructures
    ipanemap_conf['Sampling']['Temperature'] = sampling_temperature
    ipanemap_conf['Sampling']['m'] = sampling_m 
    ipanemap_conf['Sampling']['b'] = sampling_b

    ipanemap_conf['Clustering'] = {}
    ipanemap_conf['Clustering']['MaxDiameterThreshold'] = clustering_max_diameter_thres
    ipanemap_conf['Clustering']['MaxAverageDiameterThreshold'] = clustering_max_avg_diameter_thres

    ipanemap_conf['Pareto'] = {}
    ipanemap_conf['Pareto']['Percent'] = pareto_perc
    ipanemap_conf['Pareto']['ZCutoff'] = pareto_zcutoff 

    ipanemap_conf['Visualization'] = {}
    ipanemap_conf['Visualization']['DrawModels'] = draw_models 
    ipanemap_conf['Visualization']['DrawCentroids'] = draw_centroids 
    ipanemap_conf['Visualization']['ShowProbing'] = show_probing
    
    
    with open(configfile_path, 'w') as file:
        ipanemap_conf.write(file)
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
    output: construct_path("ipanemap-config"), replicate = False)
    run:
        generate_ipanemap_config_from_snakeconf(
                configfile_path = output,
                input_softdir = config["folders"]["aggreact-ipanemap"]
                input_rnaseq =
                output_dir = temp(directory("output"))
                input_harddir = None,
                log_file = None,
                config)
        
        

rule ipanemap:
    conda: "../envs/ipanemap.yml"
    input: 
        config= construct_path("ipanemap-config"), replicate = False)
    output: "python workflow/scripts/IPANEMAP/IPANEMAP.py --config {input.config}"
