# -*- coding: utf-8 -*-
# Another alignment code in a difference approach
# By Alireza Dantism
# Created on Tuesday, July 21, 2022
# Version 1

from Bio.Blast import NCBIWWW
import numpy as np

class ScoreParams:

    def __init__(self,gap,match,mismatch):
        self.gap = gap
        self.match = match
        self.mismatch = mismatch
        
    def getMatrix(self,sizeX,sizeY):
        matrix = []
        for i in range(len(sizeY)+1):
            subMatrix = []
            for j in range(len(sizeX)+1):
                subMatrix.append(0)
            matrix.append(subMatrix)
        return matrix
        
    def getTraceBackMatrix(self,sizeX,sizeY):
        matrix = []
        for i in range(len(sizeY)+1):
            subMatrix = []
            for j in range(len(sizeX)+1):
                subMatrix.append(0)
            matrix.append(subMatrix)
        for j in range(1,len(sizeX)+1):
            matrix[0][j] = 'left'
        for i in range(1,len(sizeY)+1):
            matrix[i][0] = 'up'
        matrix[0][0] = 'done'
        return matrix


def load_primers(file_name):
    primers = []
    fileReader = open(file_name)
    while( True == True ):
        line = fileReader.readline()
        if (not line):
            break
        line = line.replace("\n","")
        primers.append(line)
    return primers
    
def read_sequence_file(fileName):
    fileReader = open(fileName)
    line = fileReader.readline()
    line = line.replace("\n","")
    sequence = ""
    while( True == True ):
        line = fileReader.readline()
        if (not line):
            break
        sequence+= line.replace("\n","")
    return sequence



def align_primer_sequence(sequence,primer):
    match_score = 10
    mismatch_score = -5

    positions = []
    identities = []
    scores = []

    mismatch = 0
    primer_len = len(primer)
    sequence_len = len(sequence)

    for position,_ in enumerate(sequence[0:sequence_len-primer_len]):
        mismatch = 0
        score = 0
        window = sequence[position:position+primer_len]

        for primer_base,seq_base in zip(primer,window):
            if primer_base ==seq_base:
                score+=match_score
            else:
                mismatch+=1
                score+=mismatch_score
        identity =1 - (mismatch/primer_len)
        score = primer_len*match_score + mismatch*mismatch_score

        if  identity > 0.8:
            positions.append(position)
            identities.append(identity)
            scores.append(score)
            
    return np.array(positions),np.array(identities),np.array(scores)


sequence = read_sequence_file("data/sequence1.fasta")
primers = load_primers("data/primers.csv")

print("\n \n Hi Dear user, welcome to the program \n ")

for primer in primers:
    
    print ("\n -- Result for primer {} ".format(primer))
    positions,identities,scores = align_primer_sequence(sequence, primer)
    sortindex = np.argsort(scores)

    for pos,iden,score in zip(positions[sortindex],identities[sortindex],scores[sortindex]):
        print("{} : {} , {}".format(pos,iden,score))
        
    
