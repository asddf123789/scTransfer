# scTransfer
## CSE546 Fall 2018
### Xingfan Huang xh11@cs.washington.edu

#### 1. Introduction
In this project we explore the extent to which we can learn a model on existing annotated datasets, apply to a new dataset and transfer annotations to the new dataset. We show that we can detect novel celltypes in new datasets and transfer cell type annotations cell by cell to new datasets.

#### 2. Files

We provide the following files:

- pancreas.R: Comparing our unsupervised learning method and that provided by Stuart and Butler et al., 2018 on human pancreas dataset.
- plotconfusion.R: Plot confusion matrix given transferred labels.	
- preprocess_SeuratVer3.R: Preprocess mouse tissue data. 
- supervised_example.py: Transfer labels with supervised learning, the example given is on mouse kidney dataset with KNN classfication. 
- novelty_detect.py: Test novelty detection method on mouse kidney dataset.