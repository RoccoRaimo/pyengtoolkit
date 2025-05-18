"""
pyengtoolkit - stay_cables
---
Calculation of stay cables characteristics.
Using the "rendered" versions it is possible to show the output of calculation on Jupyter.
"""
from handcalcs.decorator import handcalc
from math import sqrt
from math import pi

class Stay:
    def __init__(self,
                 x_d, y_d, z_d, 
                 x_p, y_p, z_p,
                 section_area,
                 elastic_modulus,
                 f_yk,
                 gamma):
        """
        Class Object that represents a stay cable.
        The cable's characteristics are defined by the following parameters:

        Parameters
        -----------------------------------------
        x_d, y_d, z_d: numbers that represents the 3 coordinates of the deck point;
        x_p, y_p, z_p: numbers that represents the 3 coordinates of the pylon point;
        section_area: area of the cross section of stay cable;
        elastic_modulus: elastic modulus of steel for stay cable;
        f_yk: yielding stress for stay cable;
        gamma: density for the material of the stay cable.

        The Stay object has the following methods:

        - Stay.length_0: length of the stay cable in undeformed condition (property)
        - Stay.deformation_elastic():
        - Stay.deformation_geometric():
        - Stay.elastic_modulus_dischinger():

        """
        self.deck_coord = [x_d, y_d, z_d]
        self.pylon_coord = [x_p, y_p, z_p]
        self.section_area: float = section_area
        self.elastic_modulus: elastic_modulus = 195000 #MPa
        self.f_yk: float = f_yk
        self.gamma: gamma = 78.5 #kN/m3

    def __repr__(self):
        return f'Stay(A={self.section_area}, E={self.elastic_modulus}, f_yk={self.f_yk})'    

    @property
    def length_0(self):
        """
        Calculus of length in underformed condition.
        """
        length_0 = sqrt(
            (self.deck_coord[0]-self.pylon_coord[0])**2 +
            (self.deck_coord[1]-self.pylon_coord[1])**2 +
            (self.deck_coord[2]-self.pylon_coord[2])**2
        )
        return length_0
    
    def deformation_elastic(self, n_pre, n_post):
        """
        Calculate the elastic deformation of the cable for a chosen delta of normal force.

        Parameters
        ------------------------
        n_pre: value of normal force before the stress
        n_post: value of normal force after the stress
        """
        delta_n = n_post - n_pre
        deformation = delta_n * self.length_0 / (self.elastic_modulus * self.section_area)
        return deformation

    def deformation_geometric(self, deform_pre: list[list[float]], deform_post: list[list[float]]):
        """
        Calculate the geometric deformation of the cable for specified vectors of deformation.

        Parameters
        ------------------------
        deform_pre: list[list[float]] - [[dx_d, dy_d, dz_d],[dx_p, dy_p, dz_p]]  | List of deformations for deck and pylon before the stress of stay cable
        deform_post: list[list[float]] - [[dx_d, dy_d, dz_d],[dx_p, dy_p, dz_p]]  | List of deformations for deck and pylon after the stress of stay cable
        """
        coord_pre_d = [
        self.deck_coord[0] + deform_pre[0][0],
        self.deck_coord[1] + deform_pre[0][1],
        self.deck_coord[2] + deform_pre[0][2]
        ]

        coord_pre_p = [
        self.pylon_coord[0] + deform_pre[1][0],
        self.pylon_coord[1] + deform_pre[1][1],
        self.pylon_coord[2] + deform_pre[1][2]
        ]

        # Length before the stress of the cable
        length_1 = sqrt(
            (coord_pre_d[0] - coord_pre_p[0])**2 +
            (coord_pre_d[1] - coord_pre_p[1])**2 +
            (coord_pre_d[2] - coord_pre_p[2])**2
        )

        coord_post_d = [
        self.deck_coord[0] + deform_post[0][0],
        self.deck_coord[1] + deform_post[0][1],
        self.deck_coord[2] + deform_post[0][2]
        ]

        coord_post_p = [
        self.pylon_coord[0] + deform_post[1][0],
        self.pylon_coord[1] + deform_post[1][1],
        self.pylon_coord[2] + deform_post[1][2]
        ]

        # Length after the stress of the cable
        length_2 = sqrt(
            (coord_post_d[0] - coord_post_p[0])**2 +
            (coord_post_d[1] - coord_post_p[1])**2 +
            (coord_post_d[2] - coord_post_p[2])**2
        )

        deformation = length_2 - length_1

        return deformation


    # def elastic_modulus_dischinger(self, L_H: float, gamma: float, sigma_ed: float):
    #     """
    #     """
    #     E_eq = self.elastic_modulus / (1+ (gamma * L_H)**2 * self.elastic_modulus / (12 * sigma_ed**3))
    #     return E_eq

    def elastic_modulus_dischinger(self, L_H: float, gamma: float, sigma_Ed: float, 
                                   prefix: str = '', precision: int = 2):
        """
        """
        try:
            @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
            def render_instance(E, L_H, gamma, sigma_Ed):
                """
                """
                E_Dischinger = (E / (1+ (gamma * L_H)**2 * E / (12 * sigma_Ed**3))).prefix(prefix)
            return render_instance(self.elastic_modulus, L_H, self.gamma, sigma_Ed)
        except:
            @handcalc(override= "", precision= precision, left= "", right= "", jupyter_display=True)
            def render_instance(E, L_H, gamma, sigma_Ed):
                """
                """
                E_Dischinger = (E / (1+ (gamma * L_H)**2 * E / (12 * sigma_Ed**3)))
            return render_instance(self.elastic_modulus, L_H, self.gamma, sigma_Ed)
            