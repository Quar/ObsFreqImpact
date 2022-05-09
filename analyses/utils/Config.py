from os import path
from glob import glob
from utils.Constants import *
from collections import namedtuple

OUTPUT_ROOT = 'path/to/output/root'

def get_root():
    return path.dirname(path.dirname(
            path.abspath(__file__)
    ))


def get_rpath(relative_path_to_root):
    return path.join(get_root(), relative_path_to_root)


def get_opath(file_name):
    return path.join(get_rpath("fig"), file_name)


def get_contact_records(dataset_name):
    basedir = f'{OUTPUT_ROOT}/pp_data/rssi_-80/contact_matrix'
    file_name = f'pp_{dataset_name}.csv'

    result_file = path.join(basedir, file_name)

    return result_file


def get_infection_pair_files(data_collection, downsample_method, dataset_name, disease_name, obs_interval=None):
    basedir_template = '{OUTPUT_ROOT}/{downsample_method}/rssi_-80/infection_pairs/{data_collection}'
    data_collections = {
        'origin': 'data_origin',
        'repeat4times': 'data_repeat4times',
    }
    downsample_methods = {
        'snapshot': 'pp_newpara_snapshot',
        'upperbound': 'pp_newpara_upperbound'
    }

    basedir = basedir_template.format(
        downsample_method=downsample_methods[downsample_method],
        data_collection=data_collections[data_collection],
    )

    expt_idx = DISEASE_NAMES.index(disease_name)

    result_folder = path.join(basedir, f'pp_{dataset_name}_{data_collection}')

    if obs_interval is not None:
        glob_str = f"*e{EXPT_EX[expt_idx]}_samp{obs_interval}_*.csv"
    else:
        glob_str = f"*e{EXPT_EX[expt_idx]}_*.csv"
    result_files = glob(
        path.join(
            path.abspath(result_folder),
            glob_str
        )
    )

    return result_files


def get_individual_infection_counts(dataset_name, expt_idx, obs_interval, data_collection, downsample_method):

    data_collections = {
        'origin': 'data_origin',
        'repeat4times': 'data_repeat4times'
    }

    downsample_methods = {
        'snapshot': 'pp_newpara_snapshot',
        'upperbound': 'pp_newpara_upperbound'
    }

    assert(data_collection in data_collections)
    assert(downsample_method in downsample_methods)

    expt_id = expt_idx + 1  # in results, expt_id starts from 1, expt_idx starts from 0

    csv_filename = f'res_infection_counts_pp_{dataset_name}_repeat4timese{expt_id}_samp{obs_interval}.csv'

    basedir = '{OUTPUT_ROOT}/{}/rssi_-80/res_infection_counts/{}/'.format(
        downsample_methods[downsample_method],
        data_collections[data_collection],
    )

    result_file = path.join(
        basedir, csv_filename
    )
    return result_file