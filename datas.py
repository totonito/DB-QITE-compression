import matplotlib.pyplot as plt   
from matplotlib import pyplot
# ###NON DEGEN
# t = [0.10,0.10,0.15,0.15, 0.2, 0.3 , 0.3 , 0.4]
# rep = [1,1,1,1, 1, 1, 1 , 1]                 #size of the ansatz
# num_qubits = [9,8,7, 6 ,5, 4, 3 , 2]
# num_circuits = [200,200,100,80, 50, 40, 40, 40]     # 40 default, increased if no convergeance, decreased if it takes too much time
# initial_depth= [345217,85674,21017,5108, 1126, 255, 35,10]
# compressed_depth = [12,11,10 ,9 ,8, 7, 6 , 5]      # if goal is only to create ground state from a FIXED initial state 1, 0, 0 ...

###NON DEGEN
#qubits variables
t = [0.10]
rep = [1,1,1,1, 1, 1, 1 , 1]                 #size of the ansatz
num_qubits = [2,3,4,5,6,7,8,9]
num_circuits = [200]     # 40 default, increased if no convergeance, decreased if it takes too much time
initial_depth= [10,35,255,1134,5113,20959,85674,345217]
compressed_depth = [5,6,7,8,9,10,11,12]      # if goal is only to create ground state from a FIXED initial state 1, 0, 0 ...
#T_variable
t_variable = [0.6,0.5,0.4,0.3,0.2,0.1]
depth_t = [1132,1124,1132,1139,1134,1134]
#num_circuits variable
num_circuitsvar = [25,50, 75, 100, 125, 150, 175, 200]
depth_num = [1128 ,1135 ,1133 ,1133 ,1138 ,1140 ,1134,1134]

#############################""

#Depth of compressed circuit vs k
k = [10, 25, 50, 75, 100, 125, 150]
compressed_depth_k = [24, 20, 24, 24, 24, 20, 8]

plt.figure()
plt.plot(k,compressed_depth_k,'.',markersize = 12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.xlabel('DB-QITE iterations',fontsize = 14)
plt.ylabel('compressed depth',fontsize = 14)


plt.figure()
plt.plot(num_qubits,initial_depth,'+',label ='depth of the original circuit')
plt.plot(num_qubits,compressed_depth,'+',label ='depth of the compressed circuit')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()

pyplot.yscale('log')
plt.xlabel('Number of qubits',fontsize = 14)
plt.ylabel('depth',fontsize = 14)
plt.legend()

plt.figure()
plt.title('compressed depth vs number of qubit')
plt.plot(num_qubits,compressed_depth,'+')
plt.xlabel('Number of qubits')
plt.ylabel('compressed depth')

plt.figure()
plt.title('t')
plt.plot(t_variable,depth_t,'+')
plt.xlabel('t')
plt.ylabel('depth')

plt.figure()
plt.title('DBQITE iterations')
plt.plot(num_circuitsvar,depth_num,'+')
plt.xlabel('t')
plt.ylabel('depth')



plt.show()

## how much compression can we given no convergeance 