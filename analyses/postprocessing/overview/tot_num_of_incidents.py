import logging
import pandas as pd
from natsort import natsorted
import multiprocessing as mp
from functools import partial
from utils import *


logger = logging.getLogger(__name__)



def get_distance_matrix(dataset_name, expt_id, population_size):
    logger.info(
        "Processing {} experiment {} pop_size {} ...".format(
            dataset_name,
            expt_id,
            population_size
        )
    )
    result_files = get_infection_pair_files(dataset_name, expt_id)
    dist_mat = calc_pairwise_distance_matrix(result_files, population_size)
    return (dataset_name, expt_id), dist_mat


def plot_matrix(distance_matrix,
                axish,
                cmap=CMAP):
    ploth = axish.matshow(distance_matrix, cmap=cmap)
    # ploth.set_xticklabels(obsFreq)
    # ploth.set_yticklabels(obsFreq)
    axish.axis('off')
    return ploth


def plot_matrices(ds_name=DS_NAME,
                  disease_names=DISEASE_NAMES,
                  expt_idx=range(11),
                  pop_sizes=POP_SIZE):
    figh = plt.figure(figsize=(21, 8))
    outer_grid = gridspec.GridSpec(
            1, 5,
            wspace=0.2, hspace=0.1
    )

    def get_inner_grid(idx_in_outer_grid):
        return gridspec.GridSpecFromSubplotSpec(
            4, 3,
            subplot_spec=outer_grid[idx_in_outer_grid],
            wspace=0.1, hspace=0.1
        )   # an array of inner_grid coords

    for d in range(len(ds_name)):
        for e in expt_idx:
            population_size = pop_sizes[DS_NAME[d]]
            axh = plt.Subplot(figh, get_inner_grid(d)[e])
            _, dist_mat = get_distance_matrix(ds_name[d], e, population_size)
            subploth = plot_matrix(dist_mat, axh)
            axh.set_title(disease_names[e])
            figh.add_subplot(axh)

    # plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
    figh.savefig(get_opath('distance_matrices.pdf'))

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    plot_matrices()