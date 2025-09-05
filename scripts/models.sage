# SPDX-License-Identifier: MIT
# This file contains material from:
#   six-worlds-anemoi (commit <SHA>), MIT License — see third_party/six-worlds-anemoi/LICENCE
#   anemoi-hash/anemoi-hash (commit <SHA>), MIT License — see third_party/anemoi-hash/LICENSE.md
# Local modifications by <Your Name>, <Year>: <one-line summary>.

def model_F_CICO(P, n_rounds, final_ll=True, ordering=None, debug=False):
    """
    Model F from the Anemoi paper. Without final linear layer. 
    CICO variables (first l inputs, first l outputs) inserted directly.
    
    Inputs:
    P         Anemoi permutation
    n_rounds  Number of rounds
    ordering  1: x1 > y0 > x2 > y1 > ... > x{N-1} > y{N-2} > y{N-1} > y{N} (custom, good) 
              2: x1 > x2 > ... > x{N-1} > y0 > y1 > ... > y{N}             (original, ok) 
              3: x1 < x2 < ... < x{N-1} < y0 < y1 < ... < y{N}             (original reversed, bad)
    debug     Debug flag for intermediate output
    """
    P.n_rounds = n_rounds
    
    # If no ordering is specified, use custom (good) one.
    if ordering is None:
        ordering = 1
    
    if debug:
        print("\n" + "="*100)
        print(P)
    
    #------------------------------------------------------------------------------------------
    # Variables (all variables, including CICO variables and unused output variables)
    #------------------------------------------------------------------------------------------
    
    # 2l input variables per round, 2l additional variables after last round
    num_vars = 2*P.n_cols*(P.n_rounds+1)
    
    x_vars = []
    y_vars = []
    all_vars = []
    
    for r in range(0, P.n_rounds+1):
        x_vars.append([f"X{r:02d}{i:02d}" for i in range(0, P.n_cols)])
        y_vars.append([f"Y{r:02d}{i:02d}" for i in range(0, P.n_cols)])
        all_vars += x_vars[-1]
        all_vars += y_vars[-1]
    
    assert(len(all_vars) == num_vars)
    
    # Polynomial ring
    poly_ring = PolynomialRing(base_ring=P.F, names=all_vars, order='degrevlex')
    
    if debug:
        print("-"*100)
        print(poly_ring)
    
    # Dictionary to use for equations
    dict_vars = {"X" : [], "Y" : []}
    for r in range(0, P.n_rounds+1):
        dict_vars["X"].append([])
        dict_vars["Y"].append([])
        for i in range(0, P.n_cols):
            dict_vars["X"][r].append(poly_ring(f"X{r:02d}{i:02d}"))
            dict_vars["Y"][r].append(poly_ring(f"Y{r:02d}{i:02d}"))
    
    if debug:
        for var_name in dict_vars:
            print(dict_vars[var_name])
            
    #------------------------------------------------------------------------------------------
    # Equations
    #------------------------------------------------------------------------------------------
    
    # Equations that all the intermediate values must satisfy. Implicitely relies on open Flystel.
    all_eqs = []
    
    # 2l equations per round, in last round only equations for l rightmost variables
    num_eqs = 2*P.n_cols*(P.n_rounds)
    
    # Helper functions
    E = lambda x : (x)**P.alpha  # alpha fixed to what is defined by permutation P
    Q = lambda x, c : P.beta * (x)**P.QUAD + c  # either squared or cubed, depending on characteristic (or custom)
    Qd = lambda x : Q(x, P.delta) 
    Qg = lambda x : Q(x, P.F.zero()) # gamma is zero
    
    for r in range(0, P.n_rounds):
        # Round inputs
        x, y = dict_vars["X"][r], dict_vars["Y"][r]
        
        # Round constant addition
        x = [x[i] + P.C[r][i] for i in range(0, P.n_cols)]
        y = [y[i] + P.D[r][i] for i in range(0, P.n_cols)]
        
        # Linear Layer (outputs = inputs to open flystel)
        x, y = P.linear_layer(x, y)
        
        # Outputs of open flystel (= round outputs)
        u = dict_vars["X"][r+1]
        v = dict_vars["Y"][r+1]
        
        # Last full round and final linear layer
        if (r == P.n_rounds - 1) and final_ll:
            # Invert last linear layer
            u, v = P.inv_linear_layer(u, v)

        # Exploit open flystel verification property: H(x,y) = u,v  <=>  V(v,y) = (x,u)
        for i in range(0, P.n_cols):
            all_eqs.append(Qg(y[i]) + E(y[i]-v[i]) - x[i])
            all_eqs.append(Qd(v[i]) + E(y[i]-v[i]) - u[i])
    
    # Final linear layer application
    #if final_ll:
    #    x, y = dict_vars["X"][P.n_rounds], dict_vars["Y"][P.n_rounds]
    #    x, y = P.linear_layer(x, y)
    #    u, v = dict_vars["X"][P.n_rounds+1], dict_vars["Y"][P.n_rounds+1]
    #    for i in range(0, P.n_cols):
    #        all_eqs.append(x[i]-u[i])
    #        all_eqs.append(y[i]-v[i])
    
    assert(len(all_eqs) == num_eqs)
    
    if debug:
        print("-"*100)
        print("Number of equations: ", num_eqs)
        print("Number of variables: ", num_vars)
        print("-"*100)
        #for eqn in all_eqs:
        #    print(f"deg = {eqn.degree()}, ", eqn)
        print(f"Degrees: {[f.degree() for f in all_eqs]}")
    
    #------------------------------------------------------------------------------------------
    # CICO constraints according to anemoi paper (l inputs x, l outputs x)
    #------------------------------------------------------------------------------------------
    
    cico_vars = {cico_var : P.F.zero() for cico_var in dict_vars["X"][0] + dict_vars["X"][-1]}
    
    if debug:
        print("-"*100)
        print("CICO variables: ", cico_vars)
    
    F_CICO = [f.specialization(cico_vars) for f in all_eqs]
    
    #------------------------------------------------------------------------------------------
    # Variable ordering
    #------------------------------------------------------------------------------------------
    # REMARK: ONLY CHECKED FOR l=1
    #print(f"Apply ordering {ordering}")
    
    # ORDERING 1: custom (good)
    if (P.n_cols == 1) and ordering == 1:
        # x1 > y0 > x2 > y1 > ... > x{N-1} > y{N-2} > y{N-1} > y{N}
        V_CICO = flatten([[x,y] for x,y in zip(flatten(dict_vars['X'][1:-1]),flatten(dict_vars['Y'][:-2]))]) 
        V_CICO += dict_vars['Y'][-2] 
        V_CICO += dict_vars['Y'][-1]
    # ORDERING 2: original (ok)
    elif (P.n_cols == 1) and ordering == 2:
        # x1 > x2 > ... > x{N-1} > y0 > y1 > ... > y{N}
        V_CICO = flatten(dict_vars['X'][1:-1]) + flatten(dict_vars['Y'])
    # ORDERING 3: original reversed (bad)
    elif (P.n_cols == 1) and ordering == 3:
        # x1 < x2 < ... < x{N-1} < y0 < y1 < ... < y{N}
        V_CICO = flatten(dict_vars['X'][1:-1]) + flatten(dict_vars['Y'])
        V_CICO.reverse()
    else:
        print(f"Ordering ({ordering}) not implemented for l={P.n_cols}")
        exit(-1)
    
    # Reinterpret polynomials wrt new variable ordering
    poly_ring = PolynomialRing(base_ring=P.F, names=V_CICO, order='degrevlex')
    F_CICO = [poly_ring(f) for f in F_CICO]
    #V_CICO = [poly_ring(v) for v in V_CICO] # never used
    dict_vars['X'][0] = [0]*P.n_cols
    dict_vars['X'][-1] = [0]*P.n_cols
    
    if debug:
        print("-"*100)
        print("Number of equations: ", len(F_CICO))
        print("Number of variables: ", len(V_CICO))
        print("-"*100)
        print(poly_ring)
    
    return F_CICO, dict_vars


#=============================================================================================
def model_P_CICO(P, n_rounds, final_ll=True, ordering=1, debug=False):
    """
    Model P from the Anemoi paper. Without or without final linear layer. 
    CICO variables (first l inputs, first l outputs) inserted directly.
    
    Inputs:
    P         Anemoi permutation
    n_rounds  Number of rounds
    ordering  1: y0 > sr > ... > s1 (custom, good) 
              2: y0 > s1 > ... > sr (original, bad) 
              3: sr > ... > s1 > y0 (custom, bad)
    debug     Debug flag for intermediate output
    """   
    P.n_rounds = n_rounds
    
    # If no ordering is specified, use custom (good) one.
    if ordering is None:
        ordering = 1
    
    if debug:
        print("\n" + "="*100)
        print(P)
    
    #------------------------------------------------------------------------------------------
    # Variables
    #------------------------------------------------------------------------------------------
    # 2l input variables, l additional variables per round (one for each flystel)
    num_vars = 2*P.n_cols + P.n_cols*P.n_rounds
    
    all_vars = []
    
    for r in range(0, P.n_rounds+1):
        if r == 0:
            all_vars += [f"X{r:02d}{i:02d}" for i in range(0, P.n_cols)]
            all_vars += [f"Y{r:02d}{i:02d}" for i in range(0, P.n_cols)]
        else:
            all_vars += [f"S{r:02d}{i:02d}" for i in range(0, P.n_cols)]
    
    assert(len(all_vars) == num_vars)
    
    # Polynomial ring
    poly_ring = PolynomialRing(base_ring=P.F, names=all_vars, order='degrevlex')
    
    if debug:
        print("-"*100)
        print(poly_ring)
    
    # Dictonary to use for equations
    dict_vars = {"X" : {}, "Y" : {}, "S" : {}}
    for r in range(0, P.n_rounds+1):
        if r == 0:
            dict_vars["X"][r] = [poly_ring(f"X{r:02d}{i:02d}") for i in range(0, P.n_cols)]
            dict_vars["Y"][r] = [poly_ring(f"Y{r:02d}{i:02d}") for i in range(0, P.n_cols)]
        else:
            dict_vars["S"][r] =  [poly_ring(f"S{r:02d}{i:02d}") for i in range(0, P.n_cols)]
    if debug:
        for var_name in dict_vars:
            print(var_name, dict_vars[var_name])
            
    #------------------------------------------------------------------------------------------
    # Equations
    #------------------------------------------------------------------------------------------
    # Equations that all the intermediate values must satisfy. Implicitely relies on open Flystel.
    all_eqs = []
    
    # 2l equations per round
    num_eqs = P.n_cols*P.n_rounds
    
    # Helper functions
    E = lambda x : (x)**P.alpha  # alpha fixed to what is defined by permutation P
    Q = lambda x, c : P.beta * (x)**P.QUAD + c  # either squared or cubed, depending on characteristic (or custom)
    Qd = lambda x : Q(x, P.delta) 
    Qg = lambda x : Q(x, P.F.zero()) # gamma is zero
    
    x, y = dict_vars["X"][0], dict_vars["Y"][0]
    
    for r in range(1, P.n_rounds+1):
        # Inputs to open flystel (= round inputs + round constand addition + linear layer)
        x = [x[i] + P.C[r-1][i] for i in range(0, P.n_cols)]
        y = [y[i] + P.D[r-1][i] for i in range(0, P.n_cols)]
        x, y = P.linear_layer(x, y)
        
        # Each state gets modelled via state variable si
        s = dict_vars["S"][r]
        all_eqs += [E(s[i]) - x[i] + Qg(y[i]) for i in range(0, P.n_cols)]
        
        # Outputs of open flystel
        v = [y[i] - s[i] for i in range(0, P.n_cols)]
        x = [(x[i] - Qg(y[i])) + Qd(v[i]) for i in range(0, P.n_cols)]
        y = v
    
    assert(len(all_eqs) == num_eqs)
    
    if debug:
        print("-"*100)
        print("Number of equations: ", num_eqs)
        print("Number of variables: ", num_vars)
        print("-"*100)
        #for eqn in all_eqs:
        #    print(f"deg = {eqn.degree()}, ", eqn)
        print(f"Degrees: {[f.degree() for f in all_eqs]}")
    
    #------------------------------------------------------------------------------------------
    # CICO constraints according to anemoi paper (l inputs x, l outputs x)
    #------------------------------------------------------------------------------------------
    
    # constained input
    cico_vars = {cico_var : P.F.zero() for cico_var in dict_vars["X"][0]}
    
    # constrained output
    if final_ll:
        x, y = P.linear_layer(x, y)
    all_eqs += x
    
    F_CICO = [f.specialization(cico_vars) for f in all_eqs]
    
    #------------------------------------------------------------------------------------------
    # Fast variable ordering
    #------------------------------------------------------------------------------------------
    # REMARK: ONLY CHECKED FOR l=1
    #print(f"Apply ordering {ordering}")
    
    # ORDERING 1: custom (good)
    if (P.n_cols == 1) and ordering == 1:
        # y0 > sr > ... > s1
        V_CICO = flatten([dict_vars['S'][r] for r in range(1,P.n_rounds + 1)])
        V_CICO = dict_vars['Y'][0] + V_CICO[::-1]
    # ORDERING 2: original (bad)
    elif (P.n_cols == 1) and ordering == 2:
        # y0 > s1 > ... > sr
        V_CICO = dict_vars['Y'][0] + flatten([dict_vars['S'][r] for r in range(1,P.n_rounds + 1)])
    # ORDERING 3: custom (bad)
    elif (P.n_cols == 1) and ordering == 3:
        # sr > ... > s1 > y0
        V_CICO = flatten([dict_vars['S'][r] for r in range(1,P.n_rounds + 1)])
        V_CICO = V_CICO[::-1] + dict_vars['Y'][0]
    else:
        print(f"Ordering ({ordering}) not implemented for l={P.n_cols}")
        exit(-1)
        
    # Reinterpret polynomials wrt new variable ordering
    poly_ring = PolynomialRing(base_ring=P.F, names=V_CICO, order='degrevlex')
    F_CICO = [poly_ring(f) for f in F_CICO]
    #V_CICO = [poly_ring(v) for v in V_CICO] # never used
    dict_vars['X'][0] = [0]*P.n_cols
    
    if debug:
        print("-"*100)
        print("Number of equations: ", len(F_CICO))
        print("Number of variables: ", len(V_CICO))
        print("-"*100)
        print(poly_ring)
    
    return F_CICO, dict_vars