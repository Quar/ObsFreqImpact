import pandas as pd
from utils import *


def node_pair_incidences_count(
        result_file,
        population_size,
        data_type=INFECTIONS_PAIR_DATATYPE,
        node_pair_col_indicies=INFECTION_NODE_PAIR_COL_INDICES):
    """

    :param result_file: result file contains successful
    :param population_size: population size corresponding to result_file
    :param data_type: data_type of result_file
    :param node_pair_col_indices: tuple of col-idx (from, to) in result_file
    :return: 1D np array of successful infection counts corresponding to canonical node-pairs
    """
    df = pd.read_csv(result_file, dtype=data_type)
    groups_keys = [
            (x, y)
            for x in range(0, population_size)
            # for y in range(x, population_size)
            for y in range(0, population_size)
            if x != y
    ]

    def get_node_pair(record):
        node_pair = [record[i] for i in node_pair_col_indicies]
        return tuple(sorted(node_pair))
    # get_node_pair = partial(get_node_pair_with_col_indicies, node_pair_col_indicies=node_pair_col_indicies)

    # avoid error when apply to matrix with 0 row.
    if df.shape[0] > 0:
        df['node_pair'] = df.apply(get_node_pair, axis=1)
    else:
        df['node_pair'] = []

    return np.array(df.groupby("node_pair").size().reindex(groups_keys, fill_value=0))