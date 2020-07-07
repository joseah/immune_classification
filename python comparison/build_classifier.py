# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 11:37:16 2019

@author: Lieke
"""

from sklearn import svm
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind
from newick import *
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

@ignore_warnings(category=ConvergenceWarning)
def train_tree(data, labels, tree, classifier = 'svm_occ', dimred = True):
    '''
    Train the hierarchical classifier. 
    
    Parameters
    ----------
    data: training data (cells x genes)
    labels: labels of the training data
    tree: classification tree (build for the training data using newick.py)
    classifier: which classifier to use ('svm' or 'svm_occ')
    dimred: if dimensionality reduction should be applied
    
    Return
    ------
    tree: trained classification tree
    '''
    
    numgenes = np.shape(data)[1]
    
    if(dimred == True):
        pca = PCA(n_components = 0.9)
        pca.fit(data)
        tree[0].set_pca(pca, None) #Save PCA transformation to the root node, so we can apply it to a test set
        data = pca.transform(data)
        data = pd.DataFrame(data)
    
    #recursively train the classifiers for each node in the tree
    for n in tree[0].descendants:
        group = train_node(data, labels, n, classifier, dimred, numgenes)

    return tree


def train_node(data, labels, n, classifier, dimred, numgenes):
    '''
    Train a linear of one-class SVM for a node.
    
    Parameters
    ----------
    data: training data
    labels: labels of the training data
    n: node to train the classifier for
    classifier: which classifier to use
    dimred: dimensionality reduction
    numgenes: number of genes in the training data
    
    Return
    ------
    group: vector which indicates the positive samples of a node
    
    '''
    
    group = np.zeros(len(labels), dtype = int)

    if n.is_leaf:
        group[np.where(labels == n.name)[0]] = 1
    else:
        if n.name != None:
            group[np.where(labels == n.name)[0]] = 1
        for j in n.descendants:
            group_new = train_node(data, labels, j, classifier, dimred, numgenes)
            group[np.where(group_new == 1)[0]] = 1
            
    if(dimred):
        data = find_pcs(data, labels, group, n, numgenes)
    
    if(classifier == 'svm'):
        train_svm(data, labels, group, n)
    else:
        train_occ(data, labels, group, n)
        
    return group


def find_pcs(data, labels, group, n, numgenes):
    '''
    Do a ttest between the positive and negative samples to find explaining pcs
    '''
    
    group = find_negativesamples(labels, group, n)
    
    # positive samples
    this_class = np.where(group == 1)[0]
    this_data = data.iloc[this_class]
    
    # negative samples
    other_class = np.where(group == 2)[0]
    other_data = data.iloc[other_class]
    
    statistic, pvalue = ttest_ind(this_data, other_data, equal_var = False)

    explaining_pcs = np.where(pvalue*numgenes <= 0.05)[0]  
    
    ## If there are no explaining PCs, just pick the 5 top ones
    if len(explaining_pcs) == 0:
        explaining_pcs = np.argsort(pvalue)[:5]
    
    data = data.iloc[:,explaining_pcs]
    
    # Save the explaining pcs in the tree
    n.set_pca(None, explaining_pcs)
    
    return data


def train_svm(data, labels, group, n):
    '''
    Train a linear svm and attach to the node
    
    Parameters:
    -----------
    data: training data
    group: indicating which cells of the training data belong to that node
    n: node
    '''
    
    # group == 1 --> positive samples
    # group == 2 --> negative samples
    group = find_negativesamples(labels, group, n) 
    idx_svm = np.where((group == 1) | (group == 2))[0]
    data_svm = data.iloc[idx_svm]
    group_svm = group[idx_svm]
    
    clf = svm.LinearSVC(random_state=1).fit(data_svm, group_svm)
    n.set_classifier(clf) #save classifier to the node
    
        

def train_occ(data, labels, group, n):
    '''
    Train a one-class-classifier SVM and attach to the node
    
    Parameters:
    ----------
    data: training data
    group: indicating which cells of the training data belong to that node
    n: node
    '''
    
    data_group = data.iloc[np.where(group == 1)[0]]
    
    clf = svm.OneClassSVM(gamma = 'scale', nu = 0.05).fit(data_group)
    n.set_classifier(clf) #save classifier to the node
    
    return 

def find_negativesamples(labels, group, n):
    
    a = n.ancestor
    
    # Find the 'sister' nodes
    for i in a.descendants:
        if(i.name != n.name):
            for j in i.walk():
                group[np.where(labels == j.name)[0]] = 2
                
    # If we find no sisters, we compare with the other samples of its ancestor
    if(len(np.where(group == 2)[0])) == 0:
        group[np.where(labels == a.name)[0]] = 2

    return group









