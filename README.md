# CAMB-MultifluidDE
## Authors: Vivian Miranda, João Victor Rebouças
## Description:
This is a CAMB modification that makes it better to add multiple components to dark energy.
## Installation:
Pull this repository and install in the same way as the original CAMB
## Usage:
### Python Interface:
All dark energy models from the original CAMB still work. To use a model with multiple fluids, select `dark_energy_model = 'MultiFluidDE'`. `num_of_components` selects the number of components you want to use (so far only accepts 2 components). To select the specific models, use `models = [<model1>, <model2>, <model3>, <model4>]` (needs to have 4 components). The models are:
	- First component, Late DE models: 1 - constant w; 2 - w0wa
	- Second component, Early DE models: 1 - fluid EDE; 2 - scalar field
Each model has parameters you need to input in the `set_params` function:
 - w0 and wa for late DE
 - zc, fde_zc, theta_i and wn for fluid EDE or scalar field EDE
 - use_zc, initial_phi, m, f, V0 for scalar field EDE