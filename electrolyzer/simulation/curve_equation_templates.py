
def p_to_i_cubic_constant_temp(P_T,p1,p2,p3,p4,p5,p6): 
    #calculates i-v curve coefficients given the stack power and stack temp
    pwr,tempc=P_T
    i_stack=p1*(pwr**3) + p2*(pwr**2) +  (p3*pwr) + (p4*pwr**(1/2)) + p5
    return i_stack 

def p_to_i_squared_variable_temp(P_T,a,b,c,d,e,f): 
    """
    Given a power input (kW), temperature (C), and set of coefficients, returns
    current (A).  Coefficients can be determined using non-linear least squares
    fit (see `Stack.create_polarization`).
    """
    pwr,tempc=P_T
    i = (a * (pwr**2)) + (b * tempc**2) + (c * pwr * tempc) + (d * pwr) + (e * tempc) + f
    return i