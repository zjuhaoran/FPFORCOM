#FF算法 
import numpy as np
import pandas as pd 
import scipy.fftpack as fftpack

from matplotlib import pyplot as plt

from aaindex import aaindex1

def get_numerical_sequence(sequence,dictionary):
    return  [dictionary[aminoacid] for aminoacid in sequence]

def fft_plt_y(sequence,draw=False):

    fft_y=fftpack.fft(sequence) 
    abs_y=np.abs(fft_y)         
    normalization_y= np.true_divide(abs_y, max(abs_y))  

    x = np.linspace(0,1,len(sequence))         
    normalization_x=x   

    half_x = normalization_x[range(int(len(sequence)/2))]                                  
    normalization_half_y = normalization_y[range(int(len(sequence)/2))]      
    if draw:
        plt.figure()
        plt.plot(half_x,normalization_half_y,'b')
        plt.title('Normalized Spectrum',fontsize=16,color='blue')
        plt.show()

    return normalization_half_y 


def encoding(seq,aaindexvalue,N=None):
    if N == None:
        sequence = get_numerical_sequence(seq,aaindexvalue)
        sequence = list(np.subtract(np.array(sequence), np.mean(sequence, axis=0))) # subtract mean to avoid 0 effects
    else:
        sequence.extend([0]*(N-len(sequence)))   
    return fft_plt_y(sequence,draw=False)


def get_aai_encoding(mutdict,N=None,indices=None):
       
        #validate AAI indices are present in the input parameter, if not raise error
        if (indices == None or indices == ""):
            raise ValueError('AAI indices input parameter cannot be None or empty.')

        #check input indices is of correct type (str/list), if not raise type error
        if (not isinstance(indices, str) and (not isinstance(indices, list))):
            raise TypeError("Input indices parameter must be a string or list, got {}.".format(type(indices)))

        #cast index string to list, split multiple indices using comma
        if (isinstance(indices, str)):
            if (',' in indices):
                indices = indices.split(',')  #split on ',' just in case multiple indices passed in as str
            else:
                indices = [indices]
        #get seqs from mut_dict
        mut_seq = [val for val in mutdict.values()] 

        #create zeros numpy array to store encoded sequence output
        encoded_aai_ = np.zeros((len(mutdict), len(mut_seq[0])*len(indices)))

        
        #if multiple indices used then calculate AAI index encoding for each and
        #then concatenate after each calculation

        for index in range(0, len(indices)):
            
            #get values from aaindex record using its accession number
            encoded_aai = aaindex1[indices[index]].values

            #initialise temp arrays to store encoded sequences
            temp_seq_vals = []
            temp_all_seqs = []

            #iterate through each protein sequence and amino acid, getting the AAI index encoding value
            for protein in mut_seq:
                temp_seq_vals=encoding(protein,encoded_aai,N)
                #append encoding and reset temp array
                temp_all_seqs.append(temp_seq_vals)
            
            #convert list of lists into array
            temp_all_seqs = np.array(temp_all_seqs, dtype="float32")

            #in first iteration through indices (index=0) set encoded_aai_ to zero-initialised 
            #numpy array, else concatenate to the array in previous iteration
            if (index == 0):
                encoded_aai_ = temp_all_seqs
            else:
                encoded_aai_ = np.concatenate((encoded_aai_, temp_all_seqs), axis=1)
    
        return encoded_aai_