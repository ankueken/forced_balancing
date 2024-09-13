import sys
import pandas as pd
sys.path.append('../code/')
import sfanalysis as sf

# location of gml files to analyze
gml_dir = 'gmls/'
# location to write degree sequences
deg_dir = 'degseqs/'
# make catalog of gmls and write degree sequence files
# each row of deg_df is a degree sequence file
deg_df = sf.write_degree_sequences(gml_dir, deg_dir)
# analyze all degree sequences (this will take a while for many or large data sets)
analysis_df = sf.analyze_degree_sequences(deg_dir, deg_df)
# categorize networks (by unique gml file) into scale-free categories
hyps_df = sf.categorize_networks(analysis_df)