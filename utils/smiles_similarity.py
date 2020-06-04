# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import sdf
import urllib
import os
import numpy as np
import pandas as pd
import sys
from sys import getsizeof
import scipy
from scipy import stats

import csv
#import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SDMolSupplier

from rdkit.Chem.Fingerprints import FingerprintMols


#Load Data
path1 = input("Path for list 1: ")
path2= input("Path for list 2: ")

smis1=pd.read_csv(path1)
smis1=smis1["smiles"]
smis2=pd.read_csv(path2)
smis2=smis2["smiles"]
l1=len(smis1)
l2=len(smis2)
l=l1*l2
ll=round(l/20)

#Get molecules from smiles
bad1=[]
molecules1=[]
for i,smi in enumerate(smis1):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print('smile with number:',i,'in list 1 could not be converted to molecule')
        bad1.append(i)
        continue
    molecules1.append(m)
    
bad2=[]
molecules2=[]
for i,smi in enumerate(smis2):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print('smile with number:',i,'in list 2 could not be converted to molecule')
        bad2.append(i)
        continue
    molecules2.append(m)

indices_A = [i for i, item in enumerate(smis1) if item in set(smis2)]
indices_B = [i for i, item in enumerate(smis2) if item in set(smis1)]
#smis1=[]
#smis2=[]    


#Final output matrix
similarity=np.zeros(shape=(l1,l2),dtype=np.float32)

from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem

print('Begining fingerprint calculation...wait')
fps_topol1=[FingerprintMols.FingerprintMol(x) for x in molecules1]
fps_maccs1=[MACCSkeys.GenMACCSKeys(x) for x in molecules1]
fps_pairs1 = [Pairs.GetAtomPairFingerprint(x) for x in molecules1]
fps_tts1 = [Torsions.GetTopologicalTorsionFingerprintAsIntVect(x) for x in molecules1]
fps_ecfp4_1 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024) for x in molecules1]
fps_ecfp6_1 = [AllChem.GetMorganFingerprintAsBitVect(x,3,nBits=1024) for x in molecules1]
fps_fcfp4_1 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024,useFeatures=True) for x in molecules1]
fps_fcfp6_1 = [AllChem.GetMorganFingerprintAsBitVect(x,3,nBits=1024,useFeatures=True) for x in molecules1]
print('Begining fingerprint calculation...50%')
fps_topol2=[FingerprintMols.FingerprintMol(x) for x in molecules2]
fps_maccs2=[MACCSkeys.GenMACCSKeys(x) for x in molecules2]
fps_pairs2 = [Pairs.GetAtomPairFingerprint(x) for x in molecules2]
fps_tts2 = [Torsions.GetTopologicalTorsionFingerprintAsIntVect(x) for x in molecules2]
fps_ecfp4_2 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024) for x in molecules2]
fps_ecfp6_2 = [AllChem.GetMorganFingerprintAsBitVect(x,3,nBits=1024) for x in molecules2]
fps_fcfp4_2 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024,useFeatures=True) for x in molecules2]
fps_fcfp6_2 = [AllChem.GetMorganFingerprintAsBitVect(x,3,nBits=1024,useFeatures=True) for x in molecules2]
print('Begining fingerprint calculation...done\n')

for j in bad1:
    fps_topol1.insert(j,1)
    fps_maccs1.insert(j,1)
    fps_pairs1.insert(j,1)
    fps_tts1.insert(j,1)
    fps_ecfp4_1.insert(j,1)
    fps_ecfp6_1.insert(j,1)
    fps_fcfp4_1.insert(j,1)
    fps_fcfp6_1.insert(j,1)

for j in bad2:
    fps_topol2.insert(j,1)
    fps_maccs2.insert(j,1)
    fps_pairs2.insert(j,1)
    fps_tts2.insert(j,1)
    fps_ecfp4_2.insert(j,1)
    fps_ecfp6_2.insert(j,1)
    fps_fcfp4_2.insert(j,1)
    fps_fcfp6_2.insert(j,1)
    
print('Begining of fingerprints similarity calculation\n')



k=0
for i in range(l1):
    for j in range(l2):
        if not((i in bad1) or (j in bad2)):
            similarities_topol=((DataStructs.FingerprintSimilarity(fps_topol1[i],fps_topol2[j]))>=0.75)+0
            similarities_maccs=((DataStructs.FingerprintSimilarity(fps_maccs1[i],fps_maccs2[j]))>=0.85)+0
            similarities_pairs=((DataStructs.DiceSimilarity(fps_pairs1[i],fps_pairs2[j]))>=0.7)+0
            similarities_tts=((DataStructs.DiceSimilarity(fps_tts1[i],fps_tts2[j]))>=0.7)+0
            similarities_ecfp4=((DataStructs.FingerprintSimilarity(fps_ecfp4_1[i],fps_ecfp4_2[j]))>=0.65)+0 
            similarities_ecfp6=((DataStructs.FingerprintSimilarity(fps_ecfp6_1[i],fps_ecfp6_2[j]))>=0.6)+0
            similarities_fcfp4=((DataStructs.FingerprintSimilarity(fps_fcfp4_1[i],fps_fcfp4_2[j]))>=0.65)+0 
            similarities_fcfp6=((DataStructs.FingerprintSimilarity(fps_fcfp6_1[i],fps_fcfp6_2[j]))>=0.6)+0
            similarity[i][j]=(((0.5*(similarities_ecfp4+similarities_ecfp6)+0.5*(similarities_fcfp4+similarities_fcfp6)+0.5*(similarities_tts+similarities_pairs)+similarities_maccs+similarities_topol)/5)>=0.5)+0
        k+=1
    if k%ll==0:
        print('running:',(k/l)*100,'%')
        #for other similarity metrics use for example DataStructs.FingerprintSimilarity(fps[0],fps[1], metric=DataStructs.DiceSimilarity)

similarity[indices_A,indices_B]=2
similarity[bad1,:]=666
similarity[:,bad2]=666

print('End of fingerprints similarity calculation')
bad1=[]
bad2=[]
indices_A=[]
indices_B=[]
#Write in csv
print('Saving results in csv...')
df_similarity = pd.DataFrame(similarity)
similarity=[]
df_similarity.to_csv('similarity.csv')
print('END')

