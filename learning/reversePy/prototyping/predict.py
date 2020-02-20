import os
import re
import json

from pepms2 import PeptideMS2Options, PeptideMS2Predictor
from pepms2.modeling import build_model_from_weights
from pepms2.utils import save_data_json, load_peptides_csv

from peprt import PeptideRTPredictor
from peprt.models import max_sequence_length

import pandas as pd

ms2options = PeptideMS2Options.default()

c2Predictor = PeptideMS2Predictor(ms2options)
c3Predictor = c2Predictor
rtPredictor = PeptideRTPredictor()

c2Predictor.model = build_model_from_weights(options=ms2options, weights_path='c2Model.hdf5')
c3Predictor.model = build_model_from_weights(options=ms2options, weights_path='c3Model.hdf5')
rtPredictor.load_weights('rtMode.hdf5')


f = open(pepList.peptide.csv, "r")
for line in f:
    peptide = line.split()[1]
    c2Pred = c2Predictor.predict(peptide, None)
    c3Pred = c3Predictor.predict(peptide, None)
    rtPred = rtPredictor.predict(peptide, None)
