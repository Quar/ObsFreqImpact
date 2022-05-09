import numpy as np
import seaborn as sns

OBS_INTERVAL = [5, 10, 30, 60, 90, 180, 360]

EXPT_IDX = range(0, 15)

EXPT_EX = [
    *range(1, 16),   # for new run
]

REPLICATION = 30

DS_NAME = ['shed1', 'shed2', 'shed7', 'shed8', 'shed9']

POP_SIZE = dict(zip(DS_NAME, [39, 32, 61, 74, 78]))

EXPT_PERIOD_IN_DAYS = np.array([28, 28, 28, 31, 38])

TOT_NUM_SLOTS = EXPT_PERIOD_IN_DAYS * 24 * 60 / 5

DISEASE_NAMES = [
    "flu", "norovius", "SARS", "fifth", "pertussis",
    "chickenpox", "MERS", "H7N9", "diphtheria",
    "meningococcal", "measles", "covid19", "covid19alpha", "covid19beta", "covid19delta"
]   # 12 in total, excluded duplication

# "e1": "flu",
# "e2": "norovirus",
# "e3": "sars",
# "e4": "fifth",
# "e5": "pertussis",
# "e6": "chickenpox",
# "e7": "mers",
# "e8": "h7n9",
# "e9": "diphtheria",
# "e10": "meningococcal",
# "e11": "measles",
# "e12": "covid19"

DATA_COLLECTIONS = 'origin', 'repeat4times'

DOWNSAMPLE_METHODS = 'snapshot', 'upperbound'


INFECTIONS_PAIR_DATATYPE = {
    'timestampInMinutes': np.float64,
    'target': np.int32,
    'source': np.int32,
    'infectionResult': np.int32
}

INFECTION_NODE_PAIR_COL_INDICES = (1, 2)

RANDOM_SEEDS = [(seed ^ 0x5DEECE66D) & ((1 << 48) - 1) for seed in range(1, 31)]


CMAP = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
# CMAP = sns.light_palette("#3498db", as_cmap=True)

COLOR12 = ['#e6194b', '#3cb44b', '#ffe119', '#0082c8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#d2f53c', '#fabebe', '#008080', "#000000"]


PRINTABLE_DISEASE_NAMES = {
    "flu": "Flu",
    "norovius": "Norovirus",
    "SARS": "SARS",
    'fifth': "Fifth",
    "pertussis": "Pertussis",
    "chickenpox": "Chickenpox",
    "MERS": "MERS",
    "H7N9": "H7N9",
    "diphtheria": "Diphtheria",
    "meningococcal": "Meningococcal",
    "measles": "Measles",
    "covid19": "COVID19",
    "covid19alpha": "COVID19-Alpha",
    "covid19beta": "COVID19-Beta",
    "covid19delta": "COVID19-Delta"
}