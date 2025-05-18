"""
pyengtoolkit - steel_beam
---
Calculation of resistance and stability for steel beams, based on Eurocode 3.
It is possible to show inside Jupyter Notebook the Latex output of calculations with the optional parameter "render" in every method of SteelBeam class.
"""

from handcalcs.decorator import handcalc
from math import pi
from math import sqrt
import os
import json

import forallpeople 
forallpeople.environment('structural', top_level=True)
from forallpeople import Physical

# Import the database of steel profiles
dirname = os.path.dirname(__file__)
database = json.load(open(os.path.join(dirname, 'steel_profilesdb_EC.json'), 'r'))

# Create lists of profile type and list                    
profile_type = []
for profile_t in database:
    profile_type.append(profile_t)
profile_list = []
for profile_t in profile_type:
    for prof in database[profile_t]:
        profile_list.append(prof) 


class SteelBeam:
    def __init__(self,
                 length,
                 elastic_modulus,
                 f_yk,
                 profile, 
                 section_area=0,
                 section_area_shear_y=0,
                 section_area_shear_z=0,
                 section_inertia_y=0,
                 section_inertia_z=0,
                 section_w_pl_y=0,
                 section_w_pl_z=0):
        """
        Class that represents the steel beam object.
        
        Parameters
        ------------
        length: Length of the beam;
        elastic_modulus: Elastic modulus of the steel;
        f_yk: yielding tension for the steel;
        profile: The name of the section considered. It is possible to assign a 'User defined' section or choose a profile contained in 'steel_profilesdb_EC' database.
            In case the variable profile is set to 'User defined' the following parameters can be defined: 
            section_area = gross area of the steel section
            section_area_shear_y = area for shear calculation along the y-axis
            section_area_shear_z = area for shear calculation along the z-axis
            section_inertia_y = inertia of the steel section around the y-axis (principal)
            section_inertia_z = inertia of the steel section around the z-axis (secondary)
            section_w_pl_y = plastic modulus of resistance around the y-axis
            section_w_pl_z = plastic modulus of resistance around the z-axis

        Methods
        ------------    
        The SteelBeam object has the following methods:
        
        --- RESISTANCE
            - SteelBeam.normal_force_tension: TO BE CONTINUED
        
        --- STABILITY
            - SteelBeam.: TO BE CONTINUED
        """
        self.length = length
        self.elastic_modulus = elastic_modulus
        self.f_yk = f_yk

        self.profile = profile
        if self.profile == 'User defined':
            self.section_area = section_area
            self.section_area_shear_y = section_area_shear_y
            self.section_area_shear_z = section_area_shear_z
            self.section_inertia_y = section_inertia_y
            self.section_inertia_z = section_inertia_z
            self.section_w_pl_y = section_w_pl_y
            self.section_w_pl_z = section_w_pl_z
        elif self.profile in profile_list:
            for value_type in database:
                for prof in database[value_type]:
                    if self.profile == prof:
                        # Check if the variables are unit-aware
                        if isinstance(self.length, Physical): 
                            self.section_area = float(database[value_type][prof]['A']) * mm**2
                            self.section_area_shear_y = float(database[value_type][prof]['Avy']) * mm**2
                            self.section_area_shear_z = float(database[value_type][prof]['Avz']) * mm**2
                            self.section_inertia_y = float(database[value_type][prof]['Iy']) * mm**4
                            self.section_inertia_z = float(database[value_type][prof]['Iz']) * mm**4
                            self.section_w_pl_y = float(database[value_type][prof]['Wpl_y']) * mm**3
                            self.section_w_pl_z = float(database[value_type][prof]['Wpl_z']) * mm**3
                        else:
                            self.section_area = float(database[value_type][prof]['A']) 
                            self.section_area_shear_y = float(database[value_type][prof]['Avy']) 
                            self.section_area_shear_z = float(database[value_type][prof]['Avz'])
                            self.section_inertia_y = float(database[value_type][prof]['Iy'])
                            self.section_inertia_z = float(database[value_type][prof]['Iz'])
                            self.section_w_pl_y = float(database[value_type][prof]['Wpl_y'])
                            self.section_w_pl_z = float(database[value_type][prof]['Wpl_z'])
        else:
            raise ValueError(f"""The profile is not present in the current database. Please use 'User defined'!""")

    def __repr__(self):
        return f"""SteelBeam(Profile='{self.profile}' | L={self.length} | E={self.elastic_modulus} | f_yk={self.f_yk}
        A={self.section_area} | Iy={self.section_inertia_y} | Iz={self.section_inertia_z}
        Wpl_y={self.section_w_pl_y} | Wpl_z={self.section_w_pl_z})"""   
    
# ----- RESISTANCE
    
    def normal_force_tension(self, a_net: float, render = False, prefix: str = '', precision: int = 3):
        """
        ..........
        """
        gamma_m0 = 1.00
        gamma_m2 = 1.25

        if render==False:
            n_pl = self.section_area * self.f_yk / gamma_m0
            n_u = 0.9 * a_net * self.f_yk / gamma_m2
            normal_force_tension = min(n_pl, n_u)
            return normal_force_tension
            
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A, f_yk, A_net):
                    """
                    """
                    N_plRd = (A * f_yk / gamma_m0).prefix(prefix)
                    N_uRd = (0.9 * A_net * f_yk / gamma_m2).prefix(prefix)
                    N_tRd = min(N_plRd, N_uRd)
                return render_instance(self.section_area, self.f_yk, a_net)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A, f_yk, A_net):
                    """
                    """
                    N_plRd = A * f_yk / gamma_m0
                    N_uRd = 0.9 * A_net * f_yk / gamma_m2
                    N_Rd = min(N_plRd, N_uRd)
                return render_instance(self.section_area, self.f_yk, a_net)

    def normal_force_compression(self, render = False, prefix: str = '', precision: int = 3):
        gamma_m0 = 1.00

        if render==False:
            normal_force_compression = self.section_area * self.f_yk / gamma_m0
            return normal_force_compression
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A, f_yk):
                    """
                    """
                    N_cRd = (A * f_yk / gamma_m0).prefix(prefix)
                return render_instance(self.section_area, self.f_yk)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A, f_yk):
                    """
                    """
                    N_cRd = A * f_yk / gamma_m0
                return render_instance(self.section_area, self.f_yk)  
    
    def bending_moment_y(self, render = False, prefix: str = '', precision: int = 3):
        gamma_m0 = 1.0
        if render==False:
            bending_moment_y = self.section_w_pl_y * self.f_yk / gamma_m0
            return bending_moment_y
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(W_ply, f_yk):
                    """
                    """
                    M_cRdy = (W_ply * f_yk / gamma_m0).prefix(prefix)
                return render_instance(self.section_w_pl_y, self.f_yk)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(W_ply, f_yk):
                    """
                    """
                    M_cRdy = W_ply * f_yk / gamma_m0
                return render_instance(self.section_w_pl_y, self.f_yk)
    
    def bending_moment_z(self, render = False, prefix: str = '', precision: int = 3):
        gamma_m0 = 1.0
        if render==False:
            bending_moment_z = self.section_w_pl_z * self.f_yk / gamma_m0
            return bending_moment_z
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(W_plz, f_yk):
                    """
                    """
                    M_cRdz = (W_plz * f_yk / gamma_m0).prefix(prefix)
                return render_instance(self.section_w_pl_z, self.f_yk)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(W_plz, f_yk):
                    """
                    """
                    M_cRdz = W_plz * f_yk / gamma_m0
                return render_instance(self.section_w_pl_z, self.f_yk)    

### Bending and axial force
    ## Uniaxial bending moment

    ## Biaxial bending moment


    def shear_y(self, render = False, prefix: str = '', precision: int = 3):
        gamma_m0 = 1.0
        if render==False:
            shear_y = self.section_area_shear_y * (self.f_yk / sqrt(3)) / gamma_m0
            return shear_y
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A_vy, f_yk):
                    """
                    """
                    V_cRdy = (A_vy * (f_yk / sqrt(3)) / gamma_m0).prefix(prefix)
                return render_instance(self.section_area_shear_y, self.f_yk)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A_vy, f_yk):
                    """
                    """
                    V_cRdy = A_vy * (f_yk / sqrt(3)) / gamma_m0
                return render_instance(self.section_area_shear_y, self.f_yk)            

    def shear_z(self, render = False, prefix: str = '', precision: int = 3):
        gamma_m0 = 1.0
        if render==False:
            shear_z = self.section_area_shear_z * (self.f_yk / sqrt(3)) / gamma_m0
            return shear_z
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A_vz, f_yk):
                    """
                    """
                    V_cRdz = (A_vz * (f_yk / sqrt(3)) / gamma_m0).prefix(prefix)
                return render_instance(self.section_area_shear_y, self.f_yk)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A_vz, f_yk):
                    """
                    """
                    V_cRdz = A_vz * (f_yk / sqrt(3)) / gamma_m0
                return render_instance(self.section_area_shear_y, self.f_yk)            

### Bending and shear

### Torsion

# ----- STABILITY

    def eulerian_critic_load(self, render = False, prefix: str = '', precision: int = 3):
        """
        """
        if render==False:    
            n_cr_y = (pi**2 * self.elastic_modulus * self.section_inertia_y) / (self.length)**2
            n_cr_z = (pi**2 * self.elastic_modulus * self.section_inertia_z) / (self.length)**2
            n_cr = min(n_cr_y, n_cr_z)
            return n_cr
        elif render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(E, I_y, I_z, L_0):
                    """
                    """
                    N_cry = ((pi**2 * E * I_y) / L_0**2).prefix(prefix)
                    N_crz = ((pi**2 * E * I_z) / L_0**2).prefix(prefix)
                    N_cr = (min(N_cry, N_crz)).prefix(prefix)
                return render_instance(self.elastic_modulus, self.section_inertia_y, self.section_inertia_z, self.length)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(E, I_y, I_z, L_0):
                    """
                    """
                    N_cry = ((pi**2 * E * I_y) / L_0**2)
                    N_crz = ((pi**2 * E * I_z) / L_0**2)
                    N_cr = min(N_cry, N_crz)
                return render_instance(self.elastic_modulus, self.section_inertia_y, self.section_inertia_z, self.length)            

    def slenderness_relative(self):
        """"
        """
        lambda_relative = sqrt(self.section_area * self.f_yk / self.eulerian_critic_load())
        return lambda_relative

    def normal_force_buckling(self, buckling_curve: str, render = False, prefix: str = '', precision: int = 3):
        """
        """
        imperfection_factors = {
        'a0': 0.13,
        'a': 0.21,
        'b': 0.34,
        'c': 0.49,
        'd': 0.76
        }
        alpha = imperfection_factors.get(buckling_curve)
        gamma_m_1 = 1.1
        lambd = self.slenderness_relative()
        if render==False:
            psi = 0.5*(1+alpha*(lambd-0.2)+lambd**2)
            chi = 1/ (psi+sqrt(psi**2 - lambd**2))
            N_bRd = chi * self.section_area * self.f_yk / gamma_m_1
            return N_bRd
        if render==True:
            try:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A, f_yk, alpha, lamb):
                    """
                    """
                    Phi = 0.5*(1+alpha*(lamb-0.2)+lamb**2)
                    chi = 1/ (Phi+sqrt(Phi**2 - lamb**2))
                    N_bRd = (chi * A * f_yk / gamma_m_1).prefix(prefix)
                return render_instance(self.section_area, self.f_yk, alpha, lambd)
            except:
                @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
                def render_instance(A, f_yk, alpha, lamb):
                    """
                    """
                    Phi = 0.5*(1+alpha*(lamb-0.2)+lamb**2)
                    chi = 1/ (Phi+sqrt(Phi**2 - lamb**2))
                    N_bRd = chi * A * f_yk / gamma_m_1
                return render_instance(self.section_area, self.f_yk, alpha, lambd)
        


    # def compression_critic_load():
    #     """
    #     """

    #     chi =
    #     csi =
    #     phi =
    #     return
    # to be continued