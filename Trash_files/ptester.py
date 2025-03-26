  
import numpy as np
import sympy as sp
from sympy import log, diff
from scipy import integrate
from scipy.special import j0

def hard_core_approach(sigmai, rhos, flag):
    
    variables = [[sp.symbols(f"n{j}_{i} ") for j in range(6)] for i in range(len(sigmai))]
    etas = [ variables[i][3] for i in range(len(sigmai))]
    

    fac1 = 1 - sum(etas)
    fac2 = sum(etas[i] for i in range(len(sigmai)) if flag[i] == 1)
    
    
    phi0 = fac1 * sp.log(1 - fac2) + fac2
    
    
    diff_1 = [sp.diff(phi0, etas[i]) for i in range(len(sigmai))]
    diff_2 = [[sp.diff(phi0, etas[i], etas[j]) for j in range(len(sigmai))] for i in range(len(sigmai))]
    diff_3 = [[[sp.diff(phi0, etas[i], etas[j], etas[k]) for k in range(len(sigmai))] for j in range(len(sigmai))] for i in range(len(sigmai))]

  
    
    

    phi1 = sum(variables[i][0] * diff_1[i] for i in range(len(sigmai)))
    phi2 = sum((variables[i][1] * variables[j][2] - variables[i][4] * variables[j][5]) * diff_2[i][j] for i in range(len(sigmai)) for j in range(len(sigmai)))
    phi3 = sum((1/(8*np.pi)) * (variables[i][2]*variables[j][2]*variables[k][2]/3.0  - variables[i][2] *variables[j][5] * variables[k][5] + (3/2) * ( variables[i][5] *variables[j][5] * ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[k][2]/2 ) - ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[i][2]/2) * ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[j][2]/2) * ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[k][2]/2 ) + 2 * ( ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[i][2]/2 )/2 *((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[j][2]/2 )/2 *((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[k][2]/2 ))    )  ) * diff_3[i][j][k] for i in range(len(sigmai)) for j in range(len(sigmai)) for k in range(len(sigmai)))
    
    total_phi = phi1+phi2+phi3
    
    
    functions = [[sp.diff(total_phi, variables[i][j]) for j in range(6)] for i in range(len(sigmai))]
    functions_func = [[sp.lambdify(variables, functions[i][j], 'numpy') for j in range(6)] for i in range(len(sigmai))]

    return functions



sigmai = [1.0, 1.0, 1.0]

flag = [0, 0, 1]

rhos = [2.4, 2.4, 0.001]

functions =  hard_core_approach (sigmai, rhos, flag)

print(functions)
