# Modelling protein translocation

## Environment setup

1. Make sure Python is installed (preferably a version > 3.9)
2. Create a virtual environment like this: ```python3 -m venv preferred-venv-location```
3. Activate the environment you just created by running ```source path-to-your-venv/bin/activate```
4. Install all necessary dependencies by running: ```pip3 install -r requirements.txt```

### Parameter exploration

1. The parameter values for the model can be set in the ```parameters.py``` file. Whether all are used, depends on the rate functions defined in ```model.py```.
2. Define your rate functions, including multipliers (Linear, Exponential, Gaussian-like, etc.)
3. Define the base nuclear import and export ranges over which you would like to simulate
4. Run ```model.py``` (output plots will be written to ```../Model_vs_Reference/{ref_trace_file}/rates/{averages_file}/```)

The model will simulate all possible base rate and multiplier combinations, generate a plot for it, and report
the plot number that has the highest similarity score.
