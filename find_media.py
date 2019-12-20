#!/bin/python3
from cobra import Model,Reaction,Metabolite
import cobra.io

filename = "Examples/iAF1260.xml"
model = cobra.io.read_sbml_model(filename)
reactions = model.reactions

for r in reactions:
    mets = r.metabolites
    one_side = True
    for met in mets:
        one_side = one_side&(mets[met]<0)
    if one_side:
        print(r.name)
        print(r)
