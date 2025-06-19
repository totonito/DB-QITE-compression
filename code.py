import numpy as np
import matplotlib.pyplot as plt   
from math import sqrt,pi
from qiskit import *
from scipy.linalg import expm,svd
import time
from qiskit.circuit.library import UnitaryGate
from qiskit_aer import QasmSimulator, Aer, AerSimulator
from scipy.sparse.linalg import eigsh
from qiskit_addon_aqc_tensor import generate_ansatz_from_circuit
from qiskit.circuit.library import EfficientSU2




def closest_unitary(U):
    # Compute the SVD
    X, _, Yh = svd(U)
    # Reconstruct the closest unitary
    return X @ Yh

#########################################################
start_time = time.time()  # Record the start time for execution time

Z = np.array([(1,0), (0,-1)])
X = np.array([(0,1), (1,0)])


#### Initial conditions #####
situation = "Ising_nondegenerate"
### ising_nondegen : for t = 0.4, stuck at E1, t = 1 stuck at E3 Straight up MAGIC, 

t = 0.1
num_circuits = 150
n_qubits = 5

stopping_point = 0.0001
rep = 1

def operator(pauli, i, N):
    left = np.identity(2**i)
    right = np.identity(2 ** (N - i - 1))
    mat = np.kron(np.kron(left, pauli), right)
    return mat

match situation:
    case "easy":
        initial_vector = 1/sqrt(2)*np.array([1,1]).reshape(-1,1)
        initial_state = initial_vector@initial_vector.T
        H = np.array([(1,0), (0,-1)])
        Emin = -1

    case "Ising_nondegenerate":
        J = 1       # Interaction strength
        h = 1       # Magnetic field strength

        ### ALL QUBITS INITIALIZED TO |+X>
        #initial_vector = (1/pow(2,n_qubits/2) * np.ones(pow(2,n_qubits))).reshape(-1, 1)
        ### ALL QUBITS INITIALIZED TO |0>
        initial_vector = np.zeros(pow(2,n_qubits)).reshape(-1,1)
        initial_vector[0,0] = 1 
        initial_state = initial_vector@initial_vector.T
        H = 0
        for i in range (1 , n_qubits):
            H += J * operator(Z, 0, n_qubits) @ operator(Z, i, n_qubits)
            H += h * operator(X, i, n_qubits)
        eigenvalues = eigsh(H, k=64, which="SA", return_eigenvectors=False)

    case "Ising_degenerate":    
        J = 1       # Interaction strength
        h = 1       # Magnetic field strength
        ### ALL QUBITS INITIALIZED TO |+X>
        #initial_vector = 1/pow(2,n_qubits/2) * np.ones(pow(2,n_qubits)).reshape(-1, 1)
        ### ALL QUBITS INITIALIZED TO |0>
        initial_vector = np.zeros(pow(2,n_qubits)).reshape(-1,1)
        initial_vector[0,0] = 1
        initial_state = initial_vector@initial_vector.T
        H = 0
        for i in range (0 , n_qubits):
            j = (i+1)%n_qubits
            H += J * operator(Z, i, n_qubits) @ operator(Z, j, n_qubits)
            H += h * operator(X, i, n_qubits)
        eigenvalues = eigsh(H, k=64, which="SA", return_eigenvectors=False)
        
    case _:
        print("TYPO error")
        raise SystemExit(1)


##### EVO_MATRIX####
H_evo = expm(1j * t * H)

###### REFLECTION MATRIX ######
Ref = expm(1j*t*initial_state)

###### COMPUTATION OF THE STATES #####

def DBQITE(iteration): 
    gate = np.identity(2**n_qubits)
    for i in range (1,iteration+1):
        if i%5 == 0:
            gate = closest_unitary(H_evo @ gate @ Ref  @ gate.conj().T @ H_evo.conj() @ gate)
        else:
            gate = H_evo @ gate @ Ref  @ gate.conj().T @ H_evo.conj() @ gate
    return gate
        


def plot_Energy_theory():
    Energy_theory = []
    for i in range(0,num_circuits+1):
        tampon = np.trace(H @ DBQITE(i)@ initial_state @DBQITE(i).conjugate().T).real
        #print("Iteration =",i ,"Energy =",tampon)
        Energy_theory.append(tampon)
    plt.plot(Energy_theory, marker='o')
    plt.xlabel('Iterations',fontsize=14)
    plt.ylabel('Energy',fontsize=14)
    for eigenvalue in eigenvalues:
        plt.axhline(y=eigenvalue, color='r', linestyle='--')#label=f'Energy ={eigenvalue}') #uncomment to see the values
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid()
    plt.show()
plot_Energy_theory()


gates = []
Energy_circuit = []
for i in range (num_circuits):
    if (i == num_circuits-1): # honteux, je sais
        qc = QuantumCircuit(n_qubits)
        qc.append(UnitaryGate(DBQITE(i)),range(n_qubits))
        backend = Aer.get_backend('statevector_simulator')

        #qc = transpile(qc,backend, optimization_level=3)
        qc = transpile(qc, basis_gates=['cx', 'rz', 'ry', 'rx'], optimization_level=3)
        result = backend.run(qc).result()
        # Get the statevector
        statevector = result.get_statevector(qc)
        # Compute the expectation value of the Hamiltonian
        expectation_value = np.real(statevector.expectation_value(H))
        Energy_circuit.append(expectation_value)
        gates.append(qc.depth())
    

#print('Energy for the initial circuit:',expectation_value)
#plt.plot(gates)
#plt.show()


def plot_Energies():
    plt.plot(Energy_circuit, marker='o')
    plt.xlabel('Iterations')
    plt.ylabel('Energy')
    for eigenvalue in eigenvalues:
        plt.axhline(y=eigenvalue, color='r', linestyle='--')#label=f'Energy ={eigenvalue}') #uncomment to see the values
    plt.legend()
    plt.grid()
    plt.show()
#plot_Energies()







circuit = EfficientSU2(n_qubits, reps=rep) #TO BE INCREASED IF IT DOES NOT CONVERGE
circuit = transpile(circuit, basis_gates=['cx', 'rz', 'ry', 'rx'], optimization_level=3)
param = list(circuit.parameters)
theta = np.random.uniform(0, 2 * np.pi, len(param)) #random initialisation of parameters



print(f"Initial circuit: depth {qc.depth()}")
print(f"Ansatz circuit: depth {circuit.depth()}, with {len(param)} parameters")



from qiskit_addon_aqc_tensor.simulation import tensornetwork_from_circuit
from scipy.optimize import OptimizeResult, minimize
from qiskit_addon_aqc_tensor.objective import MaximizeStateFidelity

#### Optimization of the parameters:
simulator_settings = AerSimulator(
    method="matrix_product_state",
    matrix_product_state_max_bond_dimension=100,
)


aqc_target_mps = tensornetwork_from_circuit(qc, simulator_settings)
objective = MaximizeStateFidelity(aqc_target_mps, circuit, simulator_settings)



def callback(intermediate_result: OptimizeResult):
    print(f"Intermediate result: Fidelity {1 - intermediate_result.fun:.8}")
    if intermediate_result.fun < stopping_point:
        # Good enough for now
        raise StopIteration


result = minimize(
    objective.loss_function,
    theta,
    method="L-BFGS-B",
    jac=True,
    options={"maxiter": 100},
    callback=callback,
)
if result.status not in (
    0,
    1,
    99,
):  # 0 => success; 1 => max iterations reached; 99 => early termination via StopIteration
    raise RuntimeError(f"Optimization failed: {result.message} (status={result.status})")

print(f"Done after {result.nit} iterations.")
aqc_final_parameters = result.x

final_circuit = circuit.assign_parameters(aqc_final_parameters)




end_time = time.time()  # Record the end time
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.6f} seconds")





## VERIFICATION that the circuit still creates the ground state



final_circuit = transpile(final_circuit, basis_gates=['cx', 'rz', 'ry', 'rx'], optimization_level=3)
#final_circuit = transpile(final_circuit,backend, optimization_level=3)
print('initial depth = ', qc.depth())
print('final depth = ', final_circuit.depth())
final_circuit.draw("mpl", fold=-1)
plt.show()



result = backend.run(final_circuit).result()
statevector = result.get_statevector(final_circuit)
expectation_value = np.real(statevector.expectation_value(H))
print('Energy for the optimized circuit :',expectation_value)
print('Energy for the initial circuit :',Energy_circuit[-1])

