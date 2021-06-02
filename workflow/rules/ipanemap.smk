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
    config = configparser.ConfigParser()
    config['Input'] = {}
    config['Input']['HardConstraintsDir'] = input_harddir
    config['Input']['SoftConstraintsDir'] = input_softdir
    config['Input']['Conditions'] = input_conditions 
    config['Input']['RNA'] = input_rnaseq

    config['Paths'] = {}
    config['Paths']['WorkingDir'] = output_dir
    config['Paths']['LogFile'] = log_file

    config['Sampling'] = {}
    config['Sampling']['DoSampling'] = sampling 
    config['Sampling']['NumStructures'] = sampling_nstructures
    config['Sampling']['Temperature'] = sampling_temperature
    config['Sampling']['m'] = sampling_m 
    config['Sampling']['b'] = sampling_b

    config['Clustering'] = {}
    config['Clustering']['MaxDiameterThreshold'] = clustering_max_diameter_thres
    config['Clustering']['MaxAverageDiameterThreshold'] = clustering_max_avg_diameter_thres

    config['Pareto'] = {}
    config['Pareto']['Percent'] = pareto_perc
    config['Pareto']['ZCutoff'] = pareto_zcutoff 

    config['Visualization'] = {}
    config['Visualization']['DrawModels'] = draw_models 
    config['Visualization']['DrawCentroids'] = draw_centroids 
    config['Visualization']['ShowProbing'] = show_probing
    
    
    with open(configfile_path, 'w') as file:
        config.write(file)
        return configfile_path
    return None
