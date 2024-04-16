# webapp_rnadiff

A streamlit application using Sequana to explore differential analysis results (using sequana.rnadiff.RNADiffResults object).

## Pickling an object
from sequana.rnadiff import RNADiffResults
import pickle

rd = RNADiffResults("rnadiff_folder", condition="my_condition")
with open("test.pkl", "wb") as f:
	pickle.dump(rd, f)
