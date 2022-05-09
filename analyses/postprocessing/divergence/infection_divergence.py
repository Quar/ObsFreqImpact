import pandas as pd
import matplotlib.pyplot as plt
from functools import partial
from scipy.stats import entropy
import multiprocessing as mp
import logging
from itertools import product as catesian_product
from itertools import repeat
from utils import *
import pickle
import matplotlib.patches as mpatches

logger = logging.getLogger(__name__)

sns.set_style({
    'axes.grid' : False,
    'axes.edgecolor': '1.0',

})


def individual_infection_prob(dataset_name, expt_idx, obs_interval, data_collection, downsample_method, smooth=True):

    def get_denominator(pop_size, replication=REPLICATION):
        return pop_size * replication

    if obs_interval == 0:
        return pd.DataFrame(
            {'infection_prob': 1.0 / POP_SIZE[dataset_name]},
            index=range(POP_SIZE[dataset_name])
        )

    df = pd.read_csv(
        get_individual_infection_counts(dataset_name, expt_idx, obs_interval, data_collection, downsample_method),
        sep=' ',
        header=None,
        names=['infection_count', 'node_id'],
        skipinitialspace=True
    )

    if smooth:
        df = df.set_index(['node_id']).reindex(range(POP_SIZE[dataset_name]), fill_value=0)
        df['infection_prob'] = (df['infection_count'] + 1) / get_denominator(POP_SIZE[dataset_name])
        assert (all(df.infection_prob > 0))
    else:
        df = df.set_index(['node_id']).reindex(range(POP_SIZE[dataset_name]), fill_value=0)
        df['infection_prob'] = df['infection_count'] / get_denominator(POP_SIZE[dataset_name])

    return df['infection_prob']


def kl_divergence(pk_obs_interval, qk_obs_interval,
                  dataset_name, expt_idx,
                  data_collection, downsample_method,
                  base=2, smooth=True):
    logger.info(
        'calc kl-divergence: dataset = {}, disease = {}, obs_interval = ({}, {}) ...'.format(
            dataset_name,
            DISEASE_NAMES[expt_idx],
            pk_obs_interval,
            qk_obs_interval
        )
    )
    f = partial(
        individual_infection_prob,
        dataset_name=dataset_name,
        expt_idx=expt_idx,
        data_collection=data_collection,
        downsample_method=downsample_method,
        smooth=smooth,
    )
    pk = f(obs_interval=pk_obs_interval).values.flatten()
    qk = f(obs_interval=qk_obs_interval).values.flatten()

    try:
        res = entropy(pk, qk, base)
    except Exception as err:
        logger.debug((
            f'\ndataset={dataset_name}'
            f', disease={DISEASE_NAMES[expt_idx]}'
            f', obs_interval=({pk_obs_interval}'
            f',{qk_obs_interval})'
            f'\n{err}\npk.values:{type(pk)} = {pk}'
            f'\n----------'
            f'\nqk.values:{type(qk)} = {qk}'
            f'\n==========='
            f'\ntry trunc to float32'
            f'\n{entropy(pk.astype("float32"), qk.astype("float32"), base)}'
            f'\ntry sync to float64'
            f'\n{entropy(pk.astype("float64"), qk.astype("float64"), base)}'
        ))
        res = entropy(pk.astype("float64"), qk.astype("float64"), base)

    res = entropy(pk.astype("float32"), qk.astype("float32"), base)

    logger.debug(f'KL@({dataset_name}, {expt_idx}, {pk_obs_interval}, {qk_obs_interval}) = {res}')

    return res


def kl_divergence_variation(dataset_name, expt_idx, data_collection, downsample_method, ref_unif=False, unif_pk=False):
    """
    ref_unif=True, unif_pk=False:  D_KL( obs || unif )
    ref_unif=True, unif_pk=True:  D_KL( unif || obs )
    ref_unif=False             :  D_KL( obs{10,30,..360} || obs{5} )

    :param dataset_name:
    :param expt_idx:
    :param data_collection:
    :param downsample_method:
    :param ref_unif: whether refers to the uniform distribution as the theoretical distribution
    :param unif_pk: whether use uniform distribution as "observed" and "obs" as theoretical, this will reverse the D_KL.
    :return:
    """
    def get_f_kld():
        pf_kld = partial(
            kl_divergence,
            dataset_name=dataset_name,
            data_collection=data_collection,
            downsample_method=downsample_method,
            expt_idx=expt_idx,
        )

        if ref_unif and not unif_pk:
            return lambda obs_i: pf_kld(pk_obs_interval=obs_i, qk_obs_interval=0)
        elif ref_unif and unif_pk:
            return lambda obs_i: pf_kld(pk_obs_interval=0, qk_obs_interval=obs_i)
        elif not ref_unif:
            return lambda obs_i: pf_kld(pk_obs_interval=OBS_INTERVAL[0], qk_obs_interval=obs_i)

        raise('Unable to determine kld function')

    f_kld = get_f_kld()
    klds = [f_kld(obs_i) for obs_i in OBS_INTERVAL[1:]]

    return pd.DataFrame(
        {
            'dataset': dataset_name,
            'disease': DISEASE_NAMES[expt_idx],
            'downsample_method': downsample_method,
            'data_collection': data_collection,
            'obs_interval': OBS_INTERVAL[1:],
            'kld': klds,
        }
    )


def calc_df_kld(ds_name_expt_idx_pair, data_collection, downsample_method, ref_unif=False, unif_pk=False):
    pool = mp.Pool()
    pf = partial(
        kl_divergence_variation,
        data_collection=data_collection,
        downsample_method=downsample_method,
        ref_unif=ref_unif,
        unif_pk=unif_pk,
    )
    kld_set = pool.starmap(pf, ds_name_expt_idx_pair)
    pool.close()
    pool.join()
    return pd.concat(kld_set, ignore_index=True)


def plot_kld_set(kld_set):
    x_coord = repeat(OBS_INTERVAL[1:], len(ds_name_expt_idx_pair))
    # cmap = plt.cm.hsv(np.linspace(0, 1, len(EXPT_IDX)))
    colors = [
        disease_col
        for _ in DS_NAME
        for disease_col in COLOR12 #iter(cmap)
    ]
    shapes = [
        '{}'.format(dataset)
        for dataset in ['o', '^', 's', 'd', '>']
        for _ in DISEASE_NAMES
    ]
    labels = [
        '{}'.format(disease)
        for _ in DS_NAME
        for disease in DISEASE_NAMES
    ]

    fig = plt.figure(figsize=(14, 8))
    # sns.set_palette("husl")
    phs = []
    for x, y, m, c, l in zip(x_coord, kld_set, shapes, colors, labels):
        ph = plt.plot(x, y, marker=m, color=c, label=l)
        phs.append(ph)

    legend_handles = [
        mpatches.Patch(color=colors, label=label)
        for colors, label in zip(COLOR12, DISEASE_NAMES)
    ]

    plt.legend(handles=legend_handles)
    plt.xlabel('Observation Interval (Unit: minutes)')
    plt.ylabel('KL-divergence on Pairwise Contacts')
    plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)

    return plt


def dump_df_kld(df_kld, downsample_method, data_collection):
    opath = f'../../ranal/kld/data/kld_{downsample_method}_{data_collection}.pickle'
    logger.info(f'write to {opath}')
    with open(opath, "bw") as f:
        pickle.dump(df_kld, f)


def save_df_kld(df_kld, downsample_method, data_collection):
    opath = f'../../ranal/kld/data/kld_{downsample_method}_{data_collection}.csv'
    df_kld.to_csv(opath, index=False)



def save_fig(plt, downsample_method, data_collection, formats=['pdf', 'png']):
    base_dir = '../../ranal/kld/data'
    output_path_template = f'{base_dir}/kl_infection_prob_{downsample_method}_{data_collection}.{{fmt}}'

    for fmt in formats:
        opath = output_path_template.format(fmt=fmt)
        logger.info(f'write to {opath}')
        plt.savefig(opath)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    selected_disease_names = [
        "flu", "SARS", "fifth", "pertussis",
        "chickenpox", "MERS", "diphtheria", "measles",
        "covid19", "covid19alpha", "covid19beta", "covid19delta",
    ]

    selected_experiment_ids = [ DISEASE_NAMES.index(disease) for disease in selected_disease_names]

    ds_name_expt_idx_pair = list(catesian_product(DS_NAME, selected_experiment_ids))

    for data_collection in ['repeat4times']:
        for downsample_method in DOWNSAMPLE_METHODS:
            df_kld = calc_df_kld(
                ds_name_expt_idx_pair,
                data_collection=data_collection,
                downsample_method=downsample_method,
                ref_unif=False, unif_pk=True,
            )
            print(df_kld.sample(n=5))
            save_df_kld(df_kld, downsample_method, data_collection)
            dump_df_kld(df_kld, downsample_method, data_collection)
            #plt = plot_kld_set(df_kld)
            #save_fig(plt, downsample_method, data_collection)
