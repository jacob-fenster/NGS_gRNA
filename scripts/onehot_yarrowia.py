# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:36:59 2022

@author: andre
"""
import tensorflow
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.wrappers.scikit_learn import KerasRegressor,KerasClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor,MLPClassifier
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np
from sklearn import tree
from sklearn.ensemble import RandomForestRegressor
from sklearn import decomposition
import matplotlib.pyplot as plt


def kmer1(guide_list):
    # convert base position to individual features
    kmer1_A = []
    kmer1_T = []
    kmer1_G = []
    kmer1_C = []
    vectors = {}
    vectors['A'] = [1,0,0,0]
    vectors['T'] = [0,1,0,0]
    vectors['G'] = [0,0,1,0]
    vectors['C'] = [0,0,0,1] 
    for guide in guide_list:
        onehot = np.zeros((len(guide),4))
        p = 0
        for base in guide:
            onehot[p] = vectors[base]
            p += 1
        onehot_transp = onehot.T
        kmer1_A.append(onehot_transp[0])
        kmer1_T.append(onehot_transp[1])
        kmer1_G.append(onehot_transp[2])
        kmer1_C.append(onehot_transp[3])
    return kmer1_A, kmer1_T, kmer1_G, kmer1_C

def kmer2(guide_list):
    # Convert position of twomers to features
    
    twomer_list = []
    l = np.zeros(19)
    for guide in guide_list:
        twomers = {'AA':l,'AT':l,'AG':l,'AC':l,'TA':l,'TT':l,'TG':l,'TC':l,'GA':l,'GT':l,'GG':l,'GC':l,'CA':l,'CT':l,'CG':l,'CC':l}
        for idx in range(19):
            mer = twomers[guide[idx:(idx+2)]].copy()
            mer[idx] = 1
            twomers[guide[idx:(idx+2)]] = mer
        twomer_list.append(twomers)
    return twomer_list

def kmer3(guide_list):
    # Convert position of 3-mers to features
    
    threemer_list = []   
    l = np.zeros(18)        #Only 18 3-mer positions now.
    
    # Generate 3mer list
    bases = ['A','T','G','C']
    threemers = []
    for a in bases:
        for b in bases:
            for c in bases:
                threemers.append(a + b + c)
    
    for guide in guide_list:
        # Generate initial dictionary
        merdict = {}   
        for p in threemers:
            merdict[p] = l
        # Change 0 to 1 based on corresponding threemer
        for idx in range(18):
            mer = merdict[guide[idx:(idx+3)]].copy()
            mer[idx] = 1
            merdict[guide[idx:(idx+3)]] = mer
        threemer_list.append(merdict)
    # The naming is opposite sorry
    return threemer_list, threemers

def kmer4(guide_list):
    fourmers = []
    l = np.zeros(17)
    
    bases = ['A','T','G','C']
    fourmer_list = []
    for a in bases:
        for b in bases:
            for c in bases:
                for d in bases:
                    fourmer_list.append(a + b + c + d)
    
    for guide in guide_list:
        merdict = {}
        for p in fourmer_list:
            merdict[p] = l
        for idx in range(17):
            mer = merdict[guide[idx:(idx+4)]].copy()
            mer[idx] = 1
            merdict[guide[idx:(idx+4)]] = mer
        fourmers.append(merdict)
    return fourmers, fourmer_list
            


def simple_model():
    # assemble the structure
    model = Sequential()
    model.add(Dense(10, input_dim=80, kernel_initializer='normal', activation='relu'))
    model.add(Dense(1, kernel_initializer='normal'))
    # compile the model
    model.compile(loss='mean_squared_error', optimizer='adam')
    return model

two_list = ['AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC']

#guides = pd.read_csv('depleted_guides.csv')
#guide_list = guides['gRNA_sequence'].to_list()
#guide_ids = guides['Guide No. '].to_list()
#guide_list = np.array(guide_list)

guides = pd.read_csv('depleted_putida_sample.csv')
guide_list = guides['gRNA_sequence'].to_list()
guide_ids = guides['gRNA'].to_list()    #or 'gRNA'
ones = kmer1(guide_list)
data = pd.DataFrame()
data['guide'] = guide_ids
data['sequence'] = guide_list

twomers = kmer2(guide_list)
threemers,three_list = kmer3(guide_list)
fourmers,four_list = kmer4(guide_list)
# Swap these when using sample data or yarrowia data
#data['deltadelta'] = guides['delta1'] - guides['delta2']
data['deltaLFC'] = guides['deltaLFC']


# Convert the single base vectors to unique columns
bases = ['A','T','G','C']
cnames = []
for base in range(4):
    basename = bases[base]
    for pos in range(20):
        As = []
        for i in range(len(guides)):
            As.append(ones[base][i][pos])
        cname = basename + str(pos + 1)
        cnames.append(cname)
        data[cname] = As
data_nn = data.copy()       
# Convert the twomer arrays into individual column features
for value in range(16):        
    for x in range(19):
        unique_column = []
        for i in range(len(guides)):
            mer = two_list[value]
            unique_column.append(twomers[i][mer][x])
        cname = mer + str(x+1)
        cnames.append(cname)
        data[cname] = unique_column

# Same thing, but for threemers
for value in range(64):
    for x in range(18):
        unique_column = []
        for i in range(len(guides)):
            mer = three_list[value]
            unique_column.append(threemers[i][mer][x])
        cname = mer + str(x+1)
        cnames.append(cname)
        data[cname] = unique_column
        
# one more time
for value in range(256):
    for x in range(17):
        col = []
        for i in range(len(guides)):
            mer = four_list[value]
            col.append(fourmers[i][mer][x])
        cname = mer + str(x+1)
        cnames.append(cname)
        data[cname] = col
            
   



X = data_nn.loc[:, ~data_nn.columns.isin(['guide','sequence','deltaLFC','deltadelta'])]
## Swap these too
#Y = data['deltadelta']
Y = data_nn['deltaLFC']

X_train, X_test, Y_train, Y_test = train_test_split(X, Y,
                                                    test_size=0.20,
                                                    random_state=2022)

estimator = KerasRegressor(build_fn=simple_model, epochs=500, batch_size=25,verbose=1) # what does verbose=0 do?
history = estimator.fit(X_train, Y_train, validation_split=0.20)
#print(estimator.model.get_weights())

### Use the neural net to predict new toxic sequences?

trainerror = []
testerror = []
"""
trees = np.arange(2,6,1)
model = tree.DecisionTreeRegressor()

# loop over depth of tree
for t in trees:
    model = RandomForestRegressor(max_depth=t,oob_score=True,n_estimators=100)
    model.fit(X_train,Y_train)
    trainerror.append(model.score(X_train,Y_train))
    testerror.append(model.oob_score_)

plt.figure(figsize=(8,4))
plt.subplot(121)
plt.plot(trees,trainerror,marker='o',label='train error')
plt.plot(trees,testerror,marker="s",label='test error')
plt.legend()
plt.xlabel('Max tree depth')
"""
#model = RandomForestRegressor(max_depth=3,oob_score=True,n_estimators=5000)
#model.fit(X_train,Y_train)

#importance = model.feature_importances_
#results = pd.DataFrame(cnames)
#results['importance'] = importance




#data_pca = X_train.copy()
#pca = decomposition.PCA(copy=False,n_components=4)
#results = pca.fit_transform(data_pca)
#print("PCA explained variance ratio=\n",pca.explained_variance_ratio_)
#print("PCA singular values=\n",pca.singular_values_)
#print("PCA components=\n",pca.components_)




