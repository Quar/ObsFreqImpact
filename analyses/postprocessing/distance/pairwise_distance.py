from sklearn.metrics.pairwise import pairwise_distances
import logging
import matplotlib.pyplot as plt
from natsort import natsorted
import multiprocessing as mp
from functools import partial
from utils import *
from scipy.spatial.distance import minkowski
from sklearn.preprocessing import normalize
from data_reader import node_pair_incidences_count
import pickle


logger = logging.getLogger(__name__)

sns.set_style({
    'axes.grid' : False,
    'axes.edgecolor': '1.0',
})


def calc_pairwise_distance_matrix(result_files, population_size):
    result_files = natsorted(result_files)
    logging.debug(
        "calc_pairwise_distance_matrix(m ): result_files = {}".format(result_files)
    )
    pool = mp.Pool()
    f = partial(node_pair_incidences_count, population_size=population_size)
    tup = pool.map(f, result_files)
    pool.close()
    pool.join()
    vector_array = np.row_stack(tup)

    def weighted_minkowski(baseline, downsampled):
        # assert population_size > 0
        return minkowski(baseline, downsampled, p=1, w=normalize([baseline], norm='l1').ravel()) / population_size

    dist_mat = pairwise_distances(vector_array, n_jobs=-1, metric=weighted_minkowski)
    return dist_mat


def get_distance_matrix(data_collection, downsample_method, dataset_name, disease_name, population_size):
    logger.info(
        "Processing {} experiment {} pop_size {} ...".format(
            dataset_name,
            DISEASE_NAMES.index(disease_name),
            POP_SIZE[dataset_name]
        )
    )
    result_files = get_infection_pair_files(data_collection, downsample_method, dataset_name, disease_name)
    dist_mat = calc_pairwise_distance_matrix(result_files, population_size)
    return (dataset_name, disease_name), dist_mat


def plot_matrix(distance_matrix,
                axish,
                cmap=CMAP,
                vmax=None,
                show_value=True):
    distance_matrix_to_degree = distance_matrix
    ploth = axish.matshow(
        distance_matrix_to_degree,
        cmap=cmap,
        vmin=0, vmax=vmax or distance_matrix_to_degree.max(),
    )
    # ploth.set_xticklabels(obsFreq)
    # ploth.set_yticklabels(obsFreq)
    axish.axis('off')

    if show_value:
        nrow, ncol = distance_matrix.shape
        logger.info(
            "distance matrix shape: ({},{})".format(
                nrow, ncol
            )
        )
        for i in range(nrow):
            for j in range(ncol):
                # lbl = "{:.0f}".format(np.sqrt(distance_matrix[i, j]))
                lbl = "{:.0f}".format(distance_matrix_to_degree[i, j])  # cosine distance to deg
                print(lbl, end=" ")
                axish.text(i, j, lbl,
                           ha="center", va="center", color="black",
                           fontsize=6)
        print("")
    return ploth


def plot_matrices_with_colorbar(
        data_collection,
        downsample_method,
        ds_names=DS_NAME,
        disease_names=DISEASE_NAMES,
        cached_pickle_file=None,
        save_to_pickle_file=None,
        output_figure_path='distance_matrices.pdf'):

    if cached_pickle_file is None:
        distance_matricies = {}
        vmax = 0

        for ds_name in ds_names:
            distance_matricies[ds_name] = {}
            for disease_name in disease_names:
                population_size = POP_SIZE[ds_name]
                _, dist_mat = get_distance_matrix(data_collection, downsample_method, ds_name, disease_name, population_size)
                distance_matricies[ds_name][disease_name] = dist_mat
                vmax = max(vmax, dist_mat.max())

        if save_to_pickle_file is not None:
            with open(save_to_pickle_file, 'wb') as pickle_file:
                pickle.dump(
                    dict(
                        distance_matricies=distance_matricies,
                        vmax=vmax,
                    ),
                    pickle_file,
                )
    else:
        with open(cached_pickle_file, 'rb') as pickle_file:
            cached_data = pickle.load(pickle_file)
            distance_matricies = cached_data['distance_matricies']
            vmax = cached_data['vmax']

    fig = plt.figure(figsize=(24, 8))
    fig_gridspec = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[1, .05, .1])

    fig_gridspec_mat, fig_gridspec_title, fig_gridspec_colorbar = fig_gridspec

    subfig_mat = fig.add_subfigure(fig_gridspec_mat)
    subfig_title = fig.add_subfigure(fig_gridspec_title)
    subfig_colorbar = fig.add_subfigure(fig_gridspec_colorbar)

    gridspec_dataasets = subfig_mat.add_gridspec(
        1, len(ds_names),
        wspace=0.2,
        hspace=0.1,
        bottom=.05,
    )

    gridspec_title = subfig_title.add_gridspec(
        1, len(ds_names),
        wspace=0.2,
        # top=0.01,
        # hspace=0.1,
        # bottom=.5,
    )

    unified_vmax = 8 # unify vmax for both snapshot and upperbound
    # unified_cmap = sns.light_palette((260, 100, 51), input="husl", as_cmap=True)
    # unified_cmap = sns.color_palette('Spectral_r', as_cmap=True)
    unified_cmap = sns.color_palette('RdYlGn_r', as_cmap=True)

    for dataset_name, gs_title, gs in zip(ds_names, list(gridspec_title), list(gridspec_dataasets)):
        dataset_ax = subfig_title.add_subplot(gs_title)
        dataset_ax.axis('off')
        # dataset_ax.set_title(dataset_name.upper(), fontsize=16, fontname='Helvetica Neue')
        dataset_ax.set_title(dataset_name.upper(), fontsize=16, fontname='Arial')
        subgridspec_diseases = gs.subgridspec(4, 3, wspace=0.1, hspace=0.1)
        for disease_name, gs_sub in zip(disease_names, list(subgridspec_diseases)):
            ax = subfig_mat.add_subplot(gs_sub)
            dist_mat = distance_matricies[dataset_name][disease_name]
            # subploth = plot_matrix(dist_mat, ax, vmax=vmax)
            subploth = plot_matrix(dist_mat, ax, vmax=unified_vmax, cmap=unified_cmap, show_value=False)  # unify vmax for both snapshot and upperbound
            ax.set_title(
                # PRINTABLE_DISEASE_NAMES[disease_name],
                disease_name,
                y=1,
                # fontname='Ubuntu Condensed',
                # fontname='Ubuntu',
                fontname='Arial',
                fontsize=16,
            )

    # norm = plt.Normalize(0, vmax)
    norm = plt.Normalize(0, unified_vmax)

    mappable = plt.cm.ScalarMappable(norm, unified_cmap)

    ax_colorbar = subfig_colorbar.add_subplot()
    # ax_colorbar = subfig_colorbar.add_axes([0.2, 0.4, 0.6, 0.5])

    plt.colorbar(
        mappable,
        cax=ax_colorbar,
        use_gridspec=True,
        orientation='horizontal',
        # shrink=0.6,
    )
    # ax_colorbar.tick_params(direction='in', length=20, bottom=False, top=True, labeltop=True, labelbottom=False, labelsize=12)
    ax_colorbar.tick_params(direction='in', length=20, labelsize=15)
    # ax_colorbar.set_title('Weighted Minkowski Distance', fontsize=14, fontname='Ubuntu')
    ax_colorbar.set_title('Weighted Minkowski Distance', fontsize=15, fontname='Arial')

    # plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
    fig.tight_layout()

    # subfig_colorbar.subplots_adjust(left=0.09, bottom=.4, right=0.1, top=0.41)
    # subfig_colorbar.subplots_adjust(bottom=0.49)
    subfig_colorbar.subplots_adjust(bottom=0.49, right=.8, left=.2)  # the parameter names are confusing, just a hack
    fig.savefig(output_figure_path)
    return fig


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    normalized = True
    output_filename_template = 'l1_wminkowski_distance_matricies_{data_collection}_{downsample_method}.pdf'

    selected_disease_names = [
        "flu", "SARS", "fifth", "pertussis",
        "chickenpox", "MERS", "diphtheria", "measles",
        "covid19", "covid19alpha", "covid19beta", "covid19delta",
    ]


    for data_collection in ['repeat4times']:
        for downsample_method in ('snapshot', 'upperbound'):
            logger.info('Processing DataCollection={}, DownsamplingMethod={}'.format(data_collection, downsample_method))
            experiment_config = {
                'data_collection': data_collection,
                'downsample_method': downsample_method,
                'disease_names': selected_disease_names,
            }
            output_figure_name = output_filename_template.format(**experiment_config)
            output_cache_name = output_figure_name.replace(".pdf", ".pickle")
            figh = plot_matrices_with_colorbar(
                **experiment_config,
                cached_pickle_file=get_opath(output_cache_name),
                # save_to_pickle_file=get_opath(output_cache_name),
                output_figure_path=get_opath(output_figure_name),
            )
            figh.show()
    exit(0)