## Libraries loading
import numpy as np
from nanocompore.SampCompDB import SampCompDB
from sklearn.preprocessing import StandardScaler
from collections import *
from nanocompore.common import *
import sys
import ast
from tqdm import tqdm

## Nanocompore database loading
dbdir = '/path/to/nanocompore/results/folder'
reference_transcriptome='/path/to/the/reference/transcriptome/fasta'

## TARDBP
tx='ENSMUST00000084125::4:148612381-148627019'

## db database loading
db = SampCompDB(dbdir+"/out_SampComp.db", fasta_fn=reference_transcriptome)

## Selection of significant positions
db.list_significant_positions(tx,'GMM_logit_pvalue',thr=0.05)

## Positions of interest for TARDBP
# Transcriptomic coordinates
txPosList=[2387, 2471, 2596, 2616, 2788]
# Genomic coordinates
genePosList=[148617470, 148617386, 148617261, 148617241, 148617069]

## For each site of interest save a set of files with: Condition, Class, Cluster probabilities, Dwell time and Current intensity
for i in [0, 1, 2, 3, 4]:
    txPos=txPosList[i]
    genePos=genePosList[i]
    print(txPos)
    print(genePos)
    data = db[tx][txPos]['data']
    # Condition labels defined according to the samples.txt file
    condition_labels = tuple(["WT", "KD"])
    # Sample labels defined according to the samples.txt file
    sample_labels = ['WT_1', 'WT_2', 'WT_3', 'KD_1', 'KD_2', 'KD_3']
    global_intensity = np.concatenate(([v['intensity'] for v in data[condition_labels[0]].values()]+[v['intensity'] for v in data[condition_labels[1]].values()]), axis=None)
    global_dwell = np.concatenate(([v['dwell'] for v in data[condition_labels[0]].values()]+[v['dwell'] for v in data[condition_labels[1]].values()]), axis=None)
    global_dwell_linear = global_dwell
    global_dwell = np.log10(global_dwell)
    # Parameters
    X = StandardScaler().fit_transform([(i, d) for i,d in zip(global_intensity, global_dwell)])    
    # Classifications
    Y = [ k for k,v in data[condition_labels[0]].items() for _ in v['intensity'] ] + [ k for k,v in data[condition_labels[1]].items() for _ in v['intensity'] ]
    # Model prediction - clusters probabilities
    model=db[tx][txPos]['txComp']['GMM_model']['model']
    y_pred_prob = model.predict_proba(X)
    y_pred = model.predict(X)
    # Saving - keep the names coherent with those required in singleMolecule.R
    np.savetxt('singleMoleculeCond_'+str(txPos)+'_'+str(genePos)+'.csv', Y, delimiter=",", fmt='%s')
    np.savetxt('singleMoleculeClass_'+str(txPos)+'_'+str(genePos)+'.csv', y_pred, delimiter=",")
    np.savetxt('singleMoleculeProba_'+str(txPos)+'_'+str(genePos)+'.csv', y_pred_prob, delimiter=",")
    np.savetxt('singleDwellTimes_'+str(txPos)+'_'+str(genePos)+'.csv', global_dwell_linear, delimiter=",")
    np.savetxt('singleIntensityTimes_'+str(txPos)+'_'+str(genePos)+'.csv', global_intensity, delimiter=",")