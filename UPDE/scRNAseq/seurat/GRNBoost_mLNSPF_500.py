#!/usr/bin/env python

#Run as a python script
#cd <arboreto repo>/resources/dream5/net1
#python run_grnboost2

#Author: Joern Pezoldt
#Date: 20.08.2018

#Function:
# 1) Perform GRNBoost



import os
import pandas as pd


from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

#Set directory
wd = '/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/scenic/mLN_pilot/int'
#Set directories for example data
net1_ex_path = wd + '/1.1_exprMatrix_filtered_t.txt'
net1_tf_path = wd + '/1.2_inputTFs.txt'

#Load data
ex_matrix = pd.read_csv(net1_ex_path, sep='\t')

#shape of matrix
ex_matrix.shape

#head of matrix
ex_matrix.head()

#load TF list from file
tf_names = load_tf_names(net1_tf_path)
#Quick inspection
tf_names[:5]
len(tf_names)

#Set computational local environment
# n ... number of nodes / 2
from distributed import LocalCluster, Client
local_cluster = LocalCluster(n_workers=12, threads_per_worker=8)
custom_client = Client(local_cluster)
custom_client

#Start Job
%%time
network = grnboost2(expression_data=ex_matrix,
                    tf_names=tf_names,
                    client_or_address=custom_client)

#QC job
network.head()
len(network)

#Save output
network.to_csv('ex_02_network.tsv', sep='\t', header=False, index=False)

#close client
custom_client.close()
local_cluster.close()





					




