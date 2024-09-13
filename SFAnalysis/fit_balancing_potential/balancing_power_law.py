import sys
import os
import pandas as pd
sys.path.append('../code/')
import sfanalysis as sf
import matplotlib.pyplot as plt
from visualisations import *

deg_dir = 'degseqs/'

deg_df = sf.organize_degree_sequences(deg_dir)

analysis_df = sf.analyze_degree_sequences(deg_dir,deg_df)
hyps_df = sf.categorize_networks(analysis_df)