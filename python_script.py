import numpy as np


maf = np.linspace(0,2,41)
float_formatter = lambda x: "%.2f" % x
maf  =np.array([ float_formatter(i) for i in maf])


max_likelihood=list()


likelihood_file = "loglikeli_maf_{0}_ldsc.npy"

for a in maf:
    file_load = likelihood_file.format(a)
    likelihood_array = np.load(file_load)
    max_likelihood.append(likelihood_array[::-1][0])

maf[np.argmax(max_likelihood)]

