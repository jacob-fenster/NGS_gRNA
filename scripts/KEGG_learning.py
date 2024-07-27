from bioservices.kegg import KEGG
import pandas as pd
import os, glob
k = KEGG()
k.organism = 'ppu'
# get all pathway Ids from P. putida KT2440 ('ppu')
k.pathwayIds

# find gene ID and description 
k.find('ppu', 'mreB')

# get a lot of information about an entry
k.get('ppu:PP_0933')
breakpoint()
# get pathway by gene, returns None is no pathway information 
k.get_pathway_by_gene('pcaI', 'ppu') #ppu00362
# get information about the pathway 
data = k.get('ppu00362')
dict_data = k.parse(data)