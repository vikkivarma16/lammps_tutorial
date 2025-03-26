  
import numpy as np
import sympy as sp
from sympy import log, diff, lambdify
from scipy import integrate
from scipy.special import j0





def free_energy_mean_field( epsilonij, sigmaij, interaction_type_ij, rhos):
    '''
    k_values = np.linspace(0.00001, 5, 100)
    
    
    new_variables= [[sp.symbols(f"pot_{i}_{j}") for j in range(len(sigmaij))] for i in range(len(sigmaij))]
    
    densities = [sp.symbols(f"rho_{i}") for i in range(len(epsilonij))]


    print(densities)
    fext =  sum(0.5 * new_variables[i][j] * densities[i] * densities[j] for j in range(len(sigmaij)) for i in range(len(sigmaij)))
            
   
   
    t_variables  = []
    
    for j in range(len(epsilonij)):
        t_variables.append(densities[j])     
    for i in range(1):
        for j in range(len(epsilonij)):
            t_variables.append (new_variables[i][j])
            
   
    
    parf = [diff(fext, densities[i]) for i in range(len(epsilonij))]
    
    '''
    xt, yt, zt = sp.symbols("xt yt zt") 
    
    
    fxt = xt**2 + yt + zt
    
    func_fxt  =  sp.lambdify((xt, yt, zt), fxt, 'numpy')
    
    
    result = func_fxt(2, 3, 4)  # xt=2, yt=3, zt=4
    print("Function output:", result)

    
    print("debugging part",func_fxt)
    
    '''
    print(t_variables)
    
    print(parf[0])

    parf_func = sp.lambdify(t_variables, parf[0], 'numpy') 
    print(parf_func)
    
    '''
    parf_func = 0
    return parf_func
        
        
epsilonij = [[0, 0, 1], [0, 0, 1], [0, 0, 1]]
sigmaij =  [[0, 0, 1], [0, 0, 1], [1, 2, 3]] 
interaction_type_ij =  [['gs', 'gs', 'gs'], ['gs', 'gs', 'gs'], ['gs', 'gs', 'gs']] 
rhos = [1, 2, 2]
        
pdphi_mf = free_energy_mean_field( epsilonij, sigmaij, interaction_type_ij, rhos)
    
print (pdphi_mf)




'''

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

#print(functions)
'''
