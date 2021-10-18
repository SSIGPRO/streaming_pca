## Streaming Algorithms for Subspace Analysis: an Edge- and Hardware-oriented review

## Paper:
Marchioni, A., Prono, L., Mangia, M., Pareschi, F., Rovatti, R., & Setti, G. (2021). [Streaming Algorithms for Subspace Analysis: an Edge-and Hardware-oriented review](https://www.techrxiv.org/articles/preprint/Streaming_Algorithms_for_Subspace_Analysis_an_Edge-_and_Hardware-oriented_review/14694420).

## Abstract
Subspace analysis is a basic tool for coping with high-dimensional data and is becoming a fundamental step in early processing of many signals elaboration tasks. Nowadays trend  of collecting huge quantities of usually very redundant data by means of decentralized systems suggests these techniques be deployed as close as possible to the data sources. Regrettably, despite its conceptual simplicity, subspace analysis is ultimately equivalent to eigenspace computation and brings along non-negligible computational and memory requirements. To make this fit into typical systems operating at the edge, specialized streaming algorithms have been recently devised that we here classify and review giving them a coherent description, highlighting features and analogies, and easing comparisons. Implementation of these methods is also tested on a common edge digital hardware platform to estimate not only abstract functional  and  complexity  features, but also more practical running times and memory footprints on which compliance with real-world applications hinges.

## The repository
This repository contains a MCU-compatible C implementation of pca streaming methods.

The repository is composed of:
- micro_test.ipynb: jupyter-notebook crafted in order to automatically communicate with MCU. This can be used to test the MCU and assess the performance of the PCA streaming methods.
- micro_src: this folder contains the C code to be loaded on MCU in order to work with micro_test.ipynb. Further instructions on the usage of the code can be found in the readme inside the folder.
- micro_pylib: this folder contains custom python libraries used by the jupyter notebooks
