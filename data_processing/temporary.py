import pickle

results=pickle.load(open('../Results/rna_data_format.p','rb'))

dav,dst,dref=results

print len(dav)

