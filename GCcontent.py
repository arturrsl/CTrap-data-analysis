#!/usr/bin/env python
# coding: utf-8

# In[67]:


import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


from Bio.SeqUtils import GC
from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord


# ## Load the SnapGene sequence and plot GC content 

# In[ ]:


filepath = '/Users/Artur/Desktop/lambda.dna'
dictionary = snapgene_file_to_dict(filepath)
seqrecord = snapgene_file_to_seqrecord(filepath)


# In[42]:


DNAseq = seqrecord[:] 
print("This sequence contains " + str(len(DNAseq.seq)) + ' base pairs')
print(DNAseq.seq)


# In[84]:


binsize = 600
GCcontent = []

GC(str(DNAseq.seq[0:5]))
   
for i in range(0,int(len(DNAseq.seq)/binsize)):
    #a = GC(str(DNAseq.seq[(i*binsize):((i+1)*binsize)]))
    GCcontent = np.append(GCcontent,[GC(str(DNAseq[(i*binsize):((i+1)*binsize)]))])

plt.plot(GCcontent, 'g'); 


# In[ ]:




