from __future__ import print_function
import keras
from keras.models import Sequential, Model, load_model
from keras import backend as K

import tensorflow as tf

import isolearn.keras as iso

import sys
import os
import pandas as pd

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from aparent.predictor import *


data = sys.argv[1]
model = sys.argv[2]
outpath = sys.argv[3]
prefix = sys.argv[4]

aparent_model = load_model(model)
aparent_encoder = get_aparent_encoder(lib_bias=4)


df = pd.read_csv(data, delimiter = "\t")


data = []
for index, row in df.iterrows():
    id = row["er_id"]
    seq = row["x"]
    peak_ixs, polya_profile = find_polya_peaks(
        aparent_model,
        aparent_encoder,
        seq,
        sequence_stride=5,
        conv_smoothing=True,
        peak_min_height=0.01,
        peak_min_distance=50,
        peak_prominence=(0.01, None)
    )
    peaks = ','.join(str(e) for e in peak_ixs)

    peak_iso_scores = score_polya_peaks(
        aparent_model,
        aparent_encoder,
        seq,
        peak_ixs,
        sequence_stride=1,
        strided_agg_mode='max',
        iso_scoring_mode='both',
        score_unit='log'
    )
    pas_score = ','.join(str(e) for e in peak_iso_scores)

    data.append([id,peaks,pas_score])

res_df = pd.DataFrame(data, columns = ["er_id", "pas_location", "pas_score"])

res_df.to_csv(outpath + "/" + prefix + "_aparent.txt", sep="\t", index=False, header=True)
