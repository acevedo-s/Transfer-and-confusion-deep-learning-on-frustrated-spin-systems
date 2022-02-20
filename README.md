# Transfer and confusion deep learning on frustrated spin systems

All the codes uploaded in this repository are available for use and distribute. In either case, please cite our work using  

@article{CORTE2021110702,

title = {Exploring neural network training strategies to determine phase transitions in frustrated magnetic models},

journal = {Computational Materials Science},

volume = {198},

pages = {110702},

year = {2021},

issn = {0927-0256},

doi = {https://doi.org/10.1016/j.commatsci.2021.110702},

url = {https://www.sciencedirect.com/science/article/pii/S0927025621004298},

author = {I. Corte and S. Acevedo and M. Arlego and C.A. Lamas},

keywords = {Frustrated magnetism, Machine learning, Ising model, Honeycomb lattice, Square lattice, Neural networks}

}




Codes written by Santiago Acevedo (acevedo-s in Github): 

-- Honey-IsingJ2-AF-ML-T.c runs several independent Monte carlo simulations of the J1-J2 Ising model on the honeycomb lattice. It calculates the main thermodynamics and exports both thermodynamics and raw configurations. 

-- Honeycomb_Dense.ipynb  Dense feed forward neural network that makes the classification from Fig.2 in our work.

-- triangular_dense_correlations.ipynb Dense feed forward neural network that makes the classification from Fig.3 in our work, using both spin configurations and spin correlations.

-- Model_transfer_Honeycomb.ipynb CNN that makes the classification from Fig. 5 in our work. 

-- Model_lattice_transfer_Neel.ipynb CNN that makes the classification from Fig. 8 in our work. 
