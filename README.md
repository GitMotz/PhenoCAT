# PhenoCAT (**Pheno**typic **C**lustering using **A**I **T**echnology)

Propulsion Academy Final Project DSWD-2019-05

Authors: Doris Berchtold, Simon Käppeli


### Goals:

Classification of ~4 million unlabeled images of human nuclei (derived from image-based siRNA screens, Berchtold et al., Molecular Cell 2019) using deep learning approaches

Identify gene perturbations that result in a specific nuclear phenotype (nucleolar caps: diamond-ring-like constellation of nucleoli and Cajal bodies) 



### Strategy:

unsupervised DL approaches:

- ClusterGAN (https://github.com/sudiptodip15/ClusterGAN)
- Autoencoders (Variational AE, DeepAE, Convolutional AE)
- Clustering of transfer values

supervised DL approaches:
- label subset of images: start with few labels=classes, then refine. 
- fine tune pre-trained CNNs (VGG16, ResNet, Inception v3) to classify unlabeled images


### Data set description

The original data set contains 3.7 mio images of individual cell nuclei from a biological screening experiment. Red and green fluorescence show specific structures, where nucleoli are red and appear roughly circular, and Cajal bodies are small, bright green and appear as dots.

A randomly sampled subset ("small data set") contains ~38K images and will be used to explore methods

Image sizes: 256x256 pixels, RGB


### Data preparation

see README file in folder 'image preparation'
