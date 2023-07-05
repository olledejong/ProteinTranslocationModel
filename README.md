# Modelling protein translocation

## Environment setup

1. Make sure Python is installed (preferably a version > 3.9)
2. Create a virtual environment like this: ```python3 -m venv preferred-venv-location```
3. Activate the environment you just created by running ```source path-to-your-venv/bin/activate```
4. Install all necessary dependencies by running: ```pip3 install -r requirements.txt```

## Running simulations

1. The parameter values for the model can be set in the ```parameters.py``` file. Whether all are used, depends on the rate functions defined in ```model.py```.
2. Define your rate functions (Linear, Exponential, Gaussian-like, etc.)
3. Run ```model.py``` (output plots will be locally written to ```./output/volume_file_used/```)
