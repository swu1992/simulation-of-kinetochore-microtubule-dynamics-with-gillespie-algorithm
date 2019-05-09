import numpy as np
import matplotlib.pyplot as plt

nDOF = 6
exp_data = [1, 1.95, 0.98, 0.25, 0.7, 0.2]
sim_data = np.zeros(nDOF)
sim_data2 = np.zeros(nDOF)
for i in range(nDOF):
    fn_data = './data/summary_%d.dat' % i
    data = np.loadtxt(fn_data)
    sim_data[i] = np.mean(data[:,2])
    sim_data2[i] = sim_data[i]/sim_data[0]
 
ind1 = np.arange(nDOF)-0.1
ind2 = np.arange(nDOF)+0.1
p1 = plt.bar(ind1, exp_data, width=0.2,color='r')
p2 = plt.bar(ind2, sim_data2, width=0.2, color='g')

plt.xlabel('RNAi treatment')
plt.ylabel('Relative flux rate')
plt.xticks(np.arange(nDOF), ('Control', '67A', '59C', 'Mast', 'EB1', 'Msps'))
plt.yticks(np.arange(0, 2.5, 0.5))
plt.legend((p1[0], p2[0]), ('Experiment', 'Simulation'))

plt.show()



    

