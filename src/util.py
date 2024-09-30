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
