#!/usr/bin/env python
# coding: utf-8

# Artur Kaczmarczyk
# August 2019

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio import SeqIO

from Bio.SeqUtils import GC
from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord


# ## Load the SnapGene sequence and plot GC content 

# In[2]:


#filepath = '/Users/Artur/Desktop/lambda.dna'     # if Snap Gene file is going to be loaded
#dictionary = snapgene_file_to_dict(filepath)
#seq_record = snapgene_file_to_seqrecord(filepath)


# loading from the genbank database
#Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="J02459.1") as handle:
    seq_record = SeqIO.read(handle, "gb") #using "gb" as an alias for "genbank"
#print("%s with %i features" % (seq_record.id, len(seq_record.features)))
#print(seq_record)



# In[19]:


DNA = seq_record[:] 
print("This sequence contains " + str(len(DNA.seq)) + ' base pairs')
print(DNA.seq)


# In[28]:


binsize = 1000
GCcontent = []

#GC(str(DNAseq.seq[0:5]))
   
for i in range(0,int(len(DNAseq.seq)-250)):
    #a = GC(str(DNAseq.seq[(i*binsize):((i+1)*binsize)]))
    if i>(len(DNAseq.seq)-binsize):
        GCcontent = np.append(GCcontent,[GC(str(DNA.seq[i:]))])
    else:
        GCcontent = np.append(GCcontent,[GC(str(DNA.seq[i:(i+binsize)]))])

        
ATcontent = (100-GCcontent)*0.01

fig, ax = plt.subplots(figsize=(20,5))

ax.tick_params(direction='out', length=6, width=2, labelsize= 15,
               grid_color='r', grid_alpha=0.5)
ax.plot(ATcontent, linewidth=2.0,color='g')
ax.set_xlabel("position (bp)", fontsize=15)
ax.set_ylabel("frequency", fontsize=15)
np.savetxt('AT_content_binsize_' + (str(binsize)) + 'nt.csv', ATcontent, delimiter=",")
plt.savefig('ATcontent.png', dpi=300, bbox_inches='tight')


# In[ ]:




