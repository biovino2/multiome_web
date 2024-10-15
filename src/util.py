"""Utility functions.

Ben Iovino  09/30/24    CZ-Biohub
"""


def get_datasets() -> 'dict[str, str]':
    """Returns a dictionary of timepoints and their corresponding datasets.

    Returns:
        dict_timepoints: Dictionary of timepoints and datasets.
    """

    timepoints = {'10 hours post fertilization': 'TDR126',
                '12 hours post fertilization': 'TDR127',
                '14 hours post fertilization': 'TDR128',
                '16 hours post fertilization': 'TDR118',
                '19 hours post fertilization': 'TDR125',
                '24 hours post fertilization': 'TDR124'}
    
    return timepoints


def get_timepoints() -> 'dict[str, str]':
    """Returns a dictionary of datasets and their corresponding timepoints.

    Returns:
        dict_datasets: Dictionary of datasets and timepoints.
    """

    datasets = {'TDR126': '10 hours post fertilization',
                'TDR127': '12 hours post fertilization',
                'TDR128': '14 hours post fertilization',
                'TDR118': '16 hours post fertilization',
                'TDR125': '19 hours post fertilization',
                'TDR124': '24 hours post fertilization'}
    
    return datasets


def get_timepoints_abbr() -> 'dict[str, str]':
    """Returns a dictionary of datasets and their corresponding timepoints.

    Returns:
        dict_datasets: Dictionary of datasets and timepoints.
    """

    datasets = {'TDR126': '10hpf',
                'TDR127': '12hpf',
                'TDR128': '14hpf',
                'TDR118': '16hpf',
                'TDR125': '19hpf',
                'TDR124': '24hpf'}
    
    return datasets


def define_color_dict() -> 'dict[str:str]':
    """Return a dictionary of colors for each cell type.

    Returns:
        dict[str, str]: The dictionary of colors.
    """

    cell_type_color_dict = {
        'NMPs': '#8dd3c7',
        'PSM': '#008080',
        'differentiating_neurons': '#bebada',
        'endocrine_pancreas': '#fb8072',
        'endoderm': '#80b1d3',
        'enteric_neurons': '#fdb462',
        'epidermis': '#b3de69',
        'fast_muscle': '#df4b9b',
        'floor_plate': '#d9d9d9',
        'hatching_gland': '#bc80bd',
        'heart_myocardium': '#ccebc5',
        'hemangioblasts': '#ffed6f',
        'hematopoietic_vasculature': '#e41a1c',
        'hindbrain': '#377eb8',
        'lateral_plate_mesoderm': '#4daf4a',
        'midbrain_hindbrain_boundary': '#984ea3',
        'muscle': '#ff7f00',
        'neural': '#e6ab02',
        'neural_crest': '#a65628',
        'neural_floor_plate': '#66a61e',
        'neural_optic': '#999999',
        'neural_posterior': '#393b7f',
        'neural_telencephalon': '#fdcdac',
        'neurons': '#cbd5e8',
        'notochord': '#f4cae4',
        'optic_cup': '#c0c000',
        'pharyngeal_arches': '#fff2ae',
        'primordial_germ_cells': '#f1e2cc',
        'pronephros': '#cccccc',
        'somites': '#1b9e77',
        'spinal_cord': '#d95f02',
        'tail_bud': '#7570b3'
    }

    return cell_type_color_dict
