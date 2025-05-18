"""
pyengtoolkit - concrete_shear 
---
Calculation of shear resistance for concrete sections with different formulation and national standard
Using the "rendered" versions it is possible to show the output of calculation on Jupyter.
"""
from handcalcs.decorator import handcalc
from math import sqrt

def concrete_EC2_NotReinf(f_ck: float, a_s: float, b_w: float, d: float, gamma_c: float, sigma_cp: float = 0):
    """
    Calculation of shear resistance for concrete section without shear reinforcement, according to EC2 §
    """
    C_Rdc = 0.18 / gamma_c
    if sqrt(1+(200/d)) < 2: 
        k = sqrt(1+(200/d))
    elif sqrt(1+(200/d)) >= 2: 
        k = 2
    if a_s / (b_w * d) <= 0.02:
        rho_l = a_s / (b_w * d)
    elif a_s / (b_w * d) > 0.02:
        rho_l = 0.02   
    k_1 = 0.15
    v_min = 0.035 * k**(3/2) * sqrt(f_ck)
    V_Rdc = max( (C_Rdc * k * (100* rho_l * f_ck)**(1/3) + k_1 * sigma_cp) * b_w * d , (v_min + k_1 * sigma_cp) * b_w * d)
    return V_Rdc

def concrete_NBR_laje(f_ck: float, a_s: float, b_w: float, d: float, gamma_c: float, sigma_cp: float = 0):
    """
    Calculation of shear resistance for slabs according to NBR 6118/2023 § 19.4.1
    """
    f_ctm = 0.3*f_ck**(2/3)
    f_ctk_inf = 0.7 * f_ctm
    f_ctd = f_ctk_inf / gamma_c
    tau_Rd = 0.25 * f_ctd
    if a_s / (b_w * d) <= 0.02:
        rho_1 = a_s / (b_w * d)
    elif a_s / (b_w * d) > 0.02:
        rho_1 = 0.02
    k = 1
    V_Rd1 = (tau_Rd * k * (1.2 + 40 * rho_1) + 0.15 * sigma_cp) * b_w * d
    return V_Rd1



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rendered versions for handcalcs output

@handcalc(override= "", precision= 2, left= "", right= "", jupyter_display=True)
def concrete_EC2_NotReinf_render(f_ck: float, a_s: float, b_w: float, d: float, gamma_c: float, sigma_cp: float = 0):
    """
    Calculation of shear resistance for concrete section without shear reinforcement, according to EC2 §
    """
    C_Rdc = 0.18 / gamma_c
    if sqrt(1+(200/d)) < 2: 
        k = sqrt(1+(200/d))
    elif sqrt(1+(200/d)) >= 2: 
        k = 2
    if a_s / (b_w * d) <= 0.02:
        rho_l = a_s / (b_w * d)
    elif a_s / (b_w * d) > 0.02:
        rho_l = 0.02   
    k_1 = 0.15
    v_min = 0.035 * k**(3/2) * sqrt(f_ck)
    V_Rdc = max( (C_Rdc * k * (100* rho_l * f_ck)**(1/3) + k_1 * sigma_cp) * b_w * d , (v_min + k_1 * sigma_cp) * b_w * d)
    return V_Rdc

@handcalc(override= "", precision= 2, left= "", right= "", jupyter_display=True)
def concrete_NBR_laje_render(f_ck: float, a_s: float, b_w: float, d: float, gamma_c: float, sigma_cp: float = 0):
    """
    Calculation of shear resistance for slabs according to NBR 6118/2023 § 19.4.1
    """
    f_ctm = 0.3*f_ck**(2/3)
    f_ctk_inf = 0.7 * f_ctm
    f_ctd = f_ctk_inf / gamma_c
    tau_Rd = 0.25 * f_ctd
    if a_s / (b_w * d) <= 0.02:
        rho_1 = a_s / (b_w * d)
    elif a_s / (b_w * d) > 0.02:
        rho_1 = 0.02
    k = 1
    V_Rd1 = (tau_Rd * k * (1.2 + 40 * rho_1) + 0.15 * sigma_cp) * b_w * d
    return V_Rd1
