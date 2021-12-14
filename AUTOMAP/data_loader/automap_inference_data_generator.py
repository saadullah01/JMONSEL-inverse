import numpy as np
import tensorflow as tf
import mat73

import sys
import os

from scipy.io import loadmat

class InferenceDataGenerator:
    def __init__(self, config):
        self.config = config

        inference_in_file = os.path.join(self.config.data_dir,self.config.inference_input)
        inference_out_file = os.path.join(self.config.data_dir,self.config.inference_target_output)

        print('*** LOADING INFERENCE INPUT DATA ***')
        inference_in_dict = loadmat(inference_in_file)
        
        print('*** LOADING INFERENCE OUTPUT DATA ***')
        inference_out_dict = loadmat(inference_out_file)

        inference_in_key = list(inference_in_dict.keys())[3]
        inference_out_key = list(inference_out_dict.keys())[3]
        
        self.input = inference_in_dict[inference_in_key]
        self.output = inference_out_dict[inference_out_key]
        
        self.input = np.transpose(inference_in_dict[inference_in_key])
        self.output = np.transpose(inference_out_dict[inference_out_key])

        self.len = self.input.shape[0]

    def next_batch(self, ind_start, batch_size):
        idx = np.arange(ind_start,ind_start+batch_size)
        yield self.input[idx], self.output[idx]
