"""
pyengtoolkit - PT_losses
---
Calculation of prestressed losses due to short and long term condition, based on prescriptions of Eurocode 2.
"""

import math
import scipy.interpolate   

class Tendon:
    def __init__(self,
                 f_pu,
                 elastic_modulus,
                 area_singlestrand,
                 strands_number,
                friction_coefficient,
                unin_angular_change,
                anchor_set
                 ):
        """
        A class that represents the tendon object.

        TO BE CONTINUED
        """
        self.f_pu = f_pu
        self.elastic_modulus = elastic_modulus
        self.area_singlestrand = area_singlestrand
        self.strands_number = strands_number
        self.friction_coefficient = friction_coefficient
        self.unin_angular_change = unin_angular_change
        self.anchor_set = anchor_set
        self.area_total = area_singlestrand * strands_number
        self.vertical_profile = []
        self.horizontal_profile = []

    def __repr__(self):
        return f"""Tendon(f_pu = {self.f_pu} | E = {self.elastic_modulus} 
        Ap_i = {self.area_singlestrand} | n_s = {self.strands_number} | Ap = {self.area_total}
        mu = {self.friction_coefficient} | k = {self.unin_angular_change} | Ds = {self.anchor_set}"""

    def set_vertical_profile(self, stations: list[float], elevations: list[float]):
        """
        """
        if len(stations) != len(elevations):
            raise ValueError(f"The number of stations and elevations must be the same!")
        
        self.vertical_profile = list(zip(stations, elevations))
        self.vertical_profile.sort(key=lambda x: x[0])
        
        # Verify ascending order
        for i in range(1, len(self.vertical_profile)):
            if self.vertical_profile[i][0] <= self.vertical_profile[i-1][0]:
                raise ValueError("The stations in vertical profile must be in ascending order")
    
    def set_horizontal_profile(self, stations: list[float], offsets: list[float]):
        """
        """
        if len(stations) != len(offsets):
            raise ValueError(f"The number of stations and offsets must be the same!")

        self.horizontal_profile = list(zip(stations, offsets))
        self.horizontal_profile.sort(key=lambda x: x[0])
        
        # Verify ascending order
        for i in range(1, len(self.horizontal_profile)):
            if self.horizontal_profile[i][0] <= self.horizontal_profile[i-1][0]:
                raise ValueError("The stations in horizontal profile must be in ascending order")

    # Combine vertical and horizontal profiles
    def tendon_profile(self):
        """
        """
        all_stations = set([s for s, _ in self.vertical_profile] + [s for s, _ in self.horizontal_profile])
        all_stations = sorted(all_stations)

        # Create interpolation functions using scipy
        if self.vertical_profile:
            v_stations, v_elevs = zip(*self.vertical_profile)
            v_interp = scipy.interpolate.interp1d(v_stations, v_elevs, kind='linear', fill_value='extrapolate')
        
        if self.horizontal_profile:
            h_stations, h_offsets = zip(*self.horizontal_profile)
            h_interp = scipy.interpolate.interp1d(h_stations, h_offsets, kind='linear', fill_value='extrapolate')
        
        # Interpolate for each station
        self.combined_profile = []
        for station in all_stations:
            elev = float(v_interp(station)) if self.vertical_profile else 0
            offset = float(h_interp(station)) if self.horizontal_profile else 0
            self.combined_profile.append((station, offset, elev))
        
        return self.combined_profile

    def angular_changes(self):
        """
        Cumulative angular change for the tendon profile (combine vertical and horizontal angular change)
        """

        # Generic parabola equation
        # y = ax^2 + bx + c
        
        # Tangent segments
        # a = 0 | b = (y2-y1) / (x2-x1) | c = -bx1 + y1

        # Parabola segments
        # TBC

        self.angular_change_v = []
        a_v = []
        b_v = []
        c_v = []
        self.angular_change_h = []
        a_h = []
        b_h = []
        c_h = []

        for count, point in enumerate(self.combined_profile):
            a_v.append(0)
            b_v.append((self.combined_profile[count][2] - self.combined_profile[count-1][2]) / (self.combined_profile[count][0] - self.combined_profile[count-1][0]))
            c_v.append(- b_v[count] * self.combined_profile[count][0] + self.combined_profile[count][2])
            if count in [0, 1]:
                angular_v = 0
            else:
                angular_v = abs(math.atan(b_v[count]) - math.atan(b_v[count - 1]))
            self.angular_change_v.append(angular_v)

            a_h.append(0)
            b_h.append((self.combined_profile[count][1] - self.combined_profile[count-1][1]) / (self.combined_profile[count][0] - self.combined_profile[count-1][0]))
            c_h.append(- b_h[count] * self.combined_profile[count][0] + self.combined_profile[count][1])
            if count in [0, 1]:
                angular_h = 0
            else:
                angular_h = abs(math.atan(b_h[count]) - math.atan(b_h[count - 1]))
            self.angular_change_h.append(angular_h)

        self.angular_change = []
        for count, angle in enumerate(self.angular_change_v):
            if count in [0, 1]:
                self.angular_change.append(0)
            else:
                self.angular_change.append(math.acos(math.cos(self.angular_change_v[count])*math.cos(self.angular_change_h[count])) + self.angular_change[count-1])

        return self.angular_change

    def calculate_losses_short(self, 
                               stress_ratio: float,
                                jacking_stress_nr: int = 2,
                                stress_order: list[int] = [1, 2],
                                stress_ratio_order: list[float] = [1, 1]):
        """
        """
        self.segment_length = []
        for count, point in enumerate(self.combined_profile):
            if count ==0:
                self.segment_length.append(0)
            else:
                self.segment_length.append(self.combined_profile[count][0] - self.combined_profile[count-1][0])
        
        # jacking_stress = stress_ratio * self.f_pu
        # friction_losses = []
        # for count, length in enumerate(self.segment_length):
        #     friction_loss = jacking_stress * (1-math.exp(-self.friction_coefficient * (self.angular_change[count] + self.unin_angular_change * self.combined_profile[count][0])))
        #     friction_losses.append(friction_loss)

        # Define the function for determination of friction losses
        def friction_losses(combined_profile: list[list[float]],
                            segment_length: list[float],
                            jacking_stress: float,
                            friction_coefficient: float,
                            angular_change: float,
                            unin_angular_change: float):
                friction_losses = []
                for count, length in enumerate(segment_length):
                    friction_loss = jacking_stress * (1-math.exp(-friction_coefficient * (angular_change[count] + unin_angular_change * combined_profile[count][0])))
                    friction_losses.append(friction_loss)
                return friction_losses
        
        for count in [0, jacking_stress_nr]:
            jacking_stress_i = self.f_pu * 
            friction_losses_i = friction_losses(self.combined_profile, self.segment_length, )


        # Apply the function according to stress steps, order and ratio
        #  
            
        # if stress_pull == 'Both-End':
        #     if stress_order == '1-2':
        #         pass
        # return friction_losses

# Short term losses (Initial Losses)

    ## Tendon profiles
        ### Vertical profile
        ### Horizontal profile
        ### Combined profile

    ## Loss due to the intermediate elastic shortening of concrete (not x dependent)

    ## Loss due to friction (x dependent)

    ## Loss due to wedge draw-in of anchorage devices (not x dependent)


#