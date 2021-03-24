# streaming_pca
Python and C implementation for MCU-compatible pca streaming methods.

The repository is composed of:
- micro_test.ipynb: jupyter-notebook crafted in order to automatically communicate with MCU. This can be used to test the MCU and assess the performance of the PCA streaming methods.
- micro_src: this folder contains the C code to be loaded on MCU in order to work with micro_test.ipynb. Further instructions on the usage of the code can be found in the readme inside the folder.
- micro_pylib: this folder contains custom python libraries used by the jupyter notebooks
