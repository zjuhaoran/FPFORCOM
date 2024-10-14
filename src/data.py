
import itertools as it 
import re

def getxulie(file):
    seq =''
    with open(file,'r') as f: 
        lines = f.readlines()
        for line in lines:
            line =line.strip()
            if '>' not in line:
                seq = seq+line
    return seq

def getKey(dic,value):
    if value not in dic.values():
        return None
    result=set()
    for key in dic:
        if dic[key] == value:
            result.add(key)
    result_s = list(result)
    result_s.sort(key = lambda x:int(re.findall(r"\d+",x)[0]))
    return result_s

def getseqs(wt,mutants,q=1):
    mutant_seqs ={}
    wt_list = list(wt)
    if q ==1 : #Generate the mutation sequence for each mutant in the list of mutants
        for mutant in mutants:
            mutant = mutant.strip()
            n = re.findall("\d+",mutant)
            if mutant[1].isdigit():
                mutantsite = re.findall(r'[A-Z]',mutant)
                seqs = wt_list[:]
                for i in range(len(n)): 
                    idx = int(n[i])-1
                    if seqs[idx] == mutantsite[2*i]:
                        seqs[idx]= mutantsite[2*i+1]
                    else:
                        print("Mutant Error(AA index not match ):",mutant)          
            else:
                seqs = wt_list  
            mutant_seqs[mutant] = ''.join(seqs)
    if q==2:  #Generate the co-occurring mutation sequences in the list of mutants.
        mutant_seqs ={}
        seqs = list(wt)  
        for mutant in mutants:
            mutant = mutant.strip()
            n = re.findall("\d+",mutant) 
            mutantsite = re.findall(r'[A-Z]',mutant)
            
            for i in range(len(n)): 
                    idx = int(n[i])-1
                    if seqs[idx] == mutantsite[2*i]:
                        seqs[idx]= mutantsite[2*i+1]
                    else:
                        print("Mutant Error(AA index not match ):",mutant)         

        mutant_seqs['/'.join(mutants)] = ''.join(seqs)

    return mutant_seqs


def getmutimutants(wt,predict_m):
    s = list(it.product((0,1),repeat=len(predict_m)))
    multi_m ={}
    for nub in s:
        if sum(nub)<=1:
            continue
        dictionary = dict(zip(predict_m, nub))  
        mutants = getKey(dictionary,1)
        nubm = set(re.findall("\d+",''.join(mutants))  ) ## Avoid the simultaneous occurrence of mutations at the same position when combining
        if len(nubm) != len(mutants):
            continue
        multi_m.update(getseqs(wt,mutants,2))
    return multi_m


def getsinglemutants(wt_seq,length):
    # Specify the mutation sites to be predicted
    # To achieve saturated single-point mutations of a specified length sequence, encode all possible single-point mutation sequences in a sequence of aimed length
    AAscan_m={}
    AAlist=['A','R','D','C','Q','E','H','I','G','N','L','K','M','F','P','S','T','W','Y','V']
    for mutAA in AAlist:
        for i in range(length):
            seq = list(wt_seq)
            name = seq[i]+str(i+1)+mutAA
            seq[i]=mutAA
            seqs =''.join(seq)
            AAscan_m[name]=seqs
    return AAscan_m