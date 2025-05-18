"""
pyengtoolkit - effective_section_composite
---
Calculation of effective width of composite (steel and concrete) sections, based on prescriptions of Eurocode 4.
"""

from dataclasses import dataclass
import pandas as pd
import numpy as np
from matplotlib.figure import Figure

def b_e_i(l_e: float, b_i: float):
    """
    A function that takes the equivalent length le and b_i and returns the equivalent width of the concrete flange on each side of the web, b_e_i.
    """
    b_e_i = min(l_e/8, b_i)
    return b_e_i

def beta_i(l_e: float, b_e_i: float):
    """
    A function that takes l_e and b_e_i and returns the i-th values of beta.
    """
    beta_i = min(0.55 + 0.025* l_e/b_e_i, 1)
    return beta_i 

def b_eff(b_0: float, beta_1: float, b_e_1: float, beta_2: float, b_e_2: float):
    """
    A function that takes b_0, beta_i and b_e_i to calculate the total effective width b_eff
    """
    b_eff = b_0 + beta_1 * b_e_1 + beta_2 * b_e_2
    return b_eff

@dataclass
class Bridge():
    """
    A class that represent a bridge geometric configuration.

    ---------------------------------------------
    Parameters
    n_spans: a int representing the number of spans
    spans_length: a list representing the length of spans
    spans_type: a list representing the type of spans
            span_type: 
                1 for "Simply supported - End"
                2 for "Simply supported - Central"
                3 for "Cantilever"
    """
    n_spans: int
    spans_length: list
    spans_type: list
    b_0: list
    b_1: list
    b_2: list

    @property
    def spans_coord(self):
        """
        A function that defines the coordinates of every span for the defined Bridge class.
        """
        spans_coord = [0]
        for count, length in enumerate(self.spans_length, start=1):
            span_coord_previous = spans_coord[count-1]
            delta_length = length         
            span_coord_i = span_coord_previous + delta_length
            spans_coord.append(span_coord_i)
        return spans_coord

    @property
    def poi_coord(self):
        """
        A function that defines the coordinates of points of interest (end-support, mid-span, internal supports, cantilevers, etc.)
        """
        if self.n_spans == 1:
            poi_coord = [0]
            poi_coord.append(self.spans_length[0]/2)
            poi_coord.append(self.spans_length[0])
        if self.n_spans == 2:
            poi_coord = [0]
            poi_coord.append(self.spans_length[0]/2)
            poi_coord.append(self.spans_length[0])
            poi_coord.append(self.spans_length[0] + self.spans_length[1]/2)
            poi_coord.append(self.spans_length[0] + self.spans_length[1])
        else:
            #Start span
            poi_coord = [0]
            poi_coord.append(self.spans_length[0]/2)
            poi_coord.append(self.spans_length[0])
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                poi_coord.append(self.spans_length[count-1] + self.spans_length[count]/2)
                poi_coord.append(self.spans_length[count-1] + self.spans_length[count])
            #End span
            poi_coord.append(poi_coord[-1] + self.spans_length[-1]/2)
            poi_coord.append(poi_coord[-1] + self.spans_length[-1]/2)
        return poi_coord

    @property
    def key_type(self):
        """
        A function that defines the key type for every point of interest.
        """
        if self.n_spans == 1:
            key_type = [1, 1, 1]
        if self.n_spans == 2:
            if self.spans_type[0] == 1 and self.spans_type[1] == 1:
                key_type = [1, 1, 2, 1, 1]
            elif self.spans_type[0] == 3 and self.spans_type[1] == 1:
                key_type = [4, 4, 4, 1, 1]
            elif self.spans_type[0] == 1 and self.spans_type[1] == 3:
                key_type = [1, 1, 4, 4, 4]        
        else:
            #Start span
            if self.spans_type[0] == 1:
                key_type = [1, 1, 2]
            elif self.spans_type[0] == 3:
                key_type = [4, 4, 4]
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                if span_type != 2:
                    raise ValueError(f"The span type for central span must be 2, not {span_type}") 
                else:
                    key_type.append(3)
                    key_type.append(2)
            #End span
            if self.spans_type[-1] == 1:
                key_type.append(1)
                key_type.append(1)
            elif self.spans_type[-1] == 3:
                key_type[-1] = 4
                key_type.append(4)
                key_type.append(4)
        return key_type

    @property
    def spans_lengths_poi(self):
        """
        A function that defines the spans lengths for every point of interest.
        """
        if self.n_spans == 1:
            spans_lengths_poi = [self.spans_length[0], self.spans_length[0], self.spans_length[0]]
        if self.n_spans == 2:
            spans_lengths_poi = [self.spans_length[0], self.spans_length[0], self.spans_length[1], self.spans_length[1], self.spans_length[1]]
        else:
            #Start span
            spans_lengths_poi = [self.spans_length[0], self.spans_length[0]]
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                spans_lengths_poi.append(self.spans_length[count])
                spans_lengths_poi.append(self.spans_length[count])
                spans_lengths_poi.append(self.spans_length[count])
            #End span
            spans_lengths_poi.append(self.spans_length[-1])
            spans_lengths_poi.append(self.spans_length[-1])
        return spans_lengths_poi

    @property
    def l_e(self):
        """
        A function that calculates the equivalent length for every point of interest previously defined.
        """
        l_e = []
        for count, key_type in enumerate(self.key_type):
            if key_type == 1:
                l_e_i = 0.85 * self.spans_lengths_poi[count]
            elif key_type == 2:
                l_e_i = 0.25 * (self.spans_lengths_poi[count-1] + self.spans_lengths_poi[count])
            elif key_type == 3:
                l_e_i = 0.7 * self.spans_lengths_poi[count]
            elif key_type == 4:
                l_e_i = 2 * self.spans_lengths_poi[count]    
            l_e.append(l_e_i)
        return l_e

    @property
    def b_0_poi(self):
        """
        """
        if self.n_spans == 1:
            b_0_poi = [self.b_0[0], self.b_0[0], self.b_0[0]]
        if self.n_spans == 2:
            b_0_poi = [self.b_0[0], self.b_0[0], self.b_0[1], self.b_0[1], self.b_0[1]]
        else:
            #Start span
            b_0_poi = [self.b_0[0], self.b_0[0]]
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                b_0_poi.append(self.b_0[count])
                b_0_poi.append(self.b_0[count])
                b_0_poi.append(self.b_0[count])
            #End span
            b_0_poi.append(self.b_0[-1])
            b_0_poi.append(self.b_0[-1])
        return b_0_poi

    @property
    def b_1_poi(self):
        """
        """
        if self.n_spans == 1:
            b_1_poi = [self.b_1[0], self.b_1[0], self.b_1[0]]
        if self.n_spans == 2:
            b_1_poi = [self.b_1[0], self.b_1[0], self.b_1[1], self.b_1[1], self.b_1[1]]
        else:
            #Start span
            b_1_poi = [self.b_1[0], self.b_1[0]]
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                b_1_poi.append(self.b_1[count])
                b_1_poi.append(self.b_1[count])
                b_1_poi.append(self.b_1[count])
            #End span
            b_1_poi.append(self.b_1[-1])
            b_1_poi.append(self.b_1[-1])
        return b_1_poi

    @property
    def b_2_poi(self):
        """
        """
        if self.n_spans == 1:
            b_2_poi = [self.b_2[0], self.b_2[0], self.b_2[0]]
        if self.n_spans == 2:
            b_2_poi = [self.b_2[0], self.b_2[0], self.b_2[1], self.b_2[1], self.b_2[1]]
        else:
            #Start span
            b_2_poi = [self.b_2[0], self.b_2[0]]
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                b_2_poi.append(self.b_2[count])
                b_2_poi.append(self.b_2[count])
                b_2_poi.append(self.b_2[count])
            #End span
            b_2_poi.append(self.b_2[-1])
            b_2_poi.append(self.b_2[-1])
        return b_2_poi

    @property
    def b_eff(self):
        """
        """
        b_e_1 = []
        for count, l_e in enumerate(self.l_e):
            b_e_1_i = b_e_i(l_e, self.b_1_poi[count])
            b_e_1.append(b_e_1_i)
        b_e_2 = []
        for count, l_e in enumerate(self.l_e):
            b_e_2_i = b_e_i(l_e, self.b_2_poi[count])
            b_e_2.append(b_e_2_i)
        beta_1 = []
        beta_2 = []
        b_eff_tot = []
        for count, value in enumerate(self.b_0_poi):
            beta_1_i = 1 
            beta_1.append(beta_1_i)
            beta_2_i = 1
            beta_2.append(beta_2_i)

        if self.key_type[0] == 1:
            beta_1[0] = beta_i(self.l_e[0], b_e_1[0])
            beta_2[0] = beta_i(self.l_e[0], b_e_2[0])
        elif self.key_type[-1] == 1:
            beta_1[-1] = beta_i(self.l_e[-1], b_e_1[-1])
            beta_2[-1] = beta_i(self.l_e[-1], b_e_2[-1])

        for count, value in enumerate(self.b_0_poi):
            b_eff_i = b_eff(self.b_0_poi[count], beta_1[count], b_e_1[count], beta_2[count], b_e_2[count])
            b_eff_tot.append(b_eff_i)
        return [b_e_1, b_e_2, beta_1, beta_2, b_eff_tot]

    # def b_eff_table(self):
    #     """
    #     A function that returns a table that contains all the informations related to effective width of section,
    #     along all the segments of the bridge.
    #     """
    #     #First the coordinates of spans, length_0 and total coordinates are defined
    #     x_span = self.spans_coord
    #     x_le = self.l_e
    #     x_tot = np.unique(x_span + x_le).tolist()

    #     #Initial and final coordinates for every segment of spans coord
    #     x_span_in = []
    #     x_span_fin = []
    #     for count, coord in enumerate(x_span[0:(len(x_span)-1)]):
    #         x_span_in.append(x_span[count])
    #         x_span_fin.append(x_span[count+1])

    #     #Initial and final coordinates for every segment of l_eff segments
    #     x_le_in = []
    #     x_le_fin = []
    #     if self.n_spans == 1:
    #         x_le_in = x_span_in
    #         x_le_fin = x_span_fin
    #     else:
    #         for count, coord in enumerate(x_le[0:(len(x_le)-1)]):
    #             x_le_in.append(x_le[count])
    #             x_le_fin.append(x_le[count+1])

    #     #Initial and final coordinates for every segment of total segments
    #     x_tot_in = []
    #     x_tot_fin = []
    #     for count, coord in enumerate(x_tot[0:(len(x_tot)-1)]):
    #         x_tot_in.append(x_tot[count])
    #         x_tot_fin.append(x_tot[count+1])
        
    #     #Number of total segments
    #     le_tot_segments = []
    #     for count in range(1,len(x_tot)):
    #         le_tot_segments.append(count)

    #     #Values of l_eff
    #     le = []
    #     for count, x in enumerate(x_le_in):
    #         le_i = x_le_fin[count] - x_le_in[count]
    #         le.append(le_i)

    #     #Values of length eff for every segment of x_tot
    #     le_tot = []
    #     #Values of le in x_tot
    #     for value1 in x_tot_fin:
    #         for count, value2 in enumerate(x_le_fin):
    #             if value1 <= value2:
    #                 le_tot.append(le[count])
    #                 break
        
    #     #Values of b_0 for every segment of x_tot
    #     b_0_tot = []
    #     #Values of l_eff in x_tot
    #     for value1 in x_tot_fin:
    #         for count, value2 in enumerate(x_span_fin):
    #             if value1 <= value2:
    #                 b_0_tot.append(self.b_0[count])
    #                 break

    #     b_1_tot = []
    #     b_2_tot = []
    #     #Values of b_1 in x_tot
    #     for value1 in x_tot_fin:
    #         for count, value2 in enumerate(x_span_fin):
    #             if value1 <= value2:
    #                 b_1_tot.append(self.b_1[count])
    #                 b_2_tot.append(self.b_2[count])
    #                 break

    #     b_eff_1 = []
    #     b_eff_2 = []
    #     for count, length in enumerate(le_tot):
    #         b_eff_1_i = b_e_i(length,b_1_tot[count])
    #         b_eff_2_i = b_e_i(length,b_2_tot[count])
    #         b_eff_1.append(b_eff_1_i)
    #         b_eff_2.append(b_eff_2_i)
        
    #     b_eff_list = []
    #     for count, b in enumerate(b_eff_1):
    #         b_eff_list.append(b_eff(b_0_tot[count], b_eff_1[count], b_eff_2[count]))

    #     #Dictionary for the DataFrame is defined
    #     dict = {"Segment": le_tot_segments,
    #             "x_in [m]": x_tot_in,
    #             "x_fin [m]": x_tot_fin,
    #             "l_0 [m]": le_tot,
    #             "b_1 [m]": b_1_tot,
    #             "b_2 [m]": b_2_tot,
    #             "b_eff_1 [m]": b_eff_1,
    #             "b_eff_2 [m]": b_eff_2,
    #             "b_0 [m]": b_0_tot,
    #             "b_eff [m]": b_eff_list
    #             }

    #     df = pd.DataFrame(data=dict)
    #     df_indexed = df.set_index(keys="Segment")
    #     return df_indexed #l0_tot_segments, x_l0_in, x_l0_fin, l0, l0_tot, b_1_tot, b_eff_1, b_eff_list

    # def plot_b_eff(self):
    #     """
    #     A function that plot the values of b_eff throughout the axis of the bridge.
    #     """
    #     #First the data are defined
    #     df = self.b_eff_table()
    #     x_in = df["x_in [m]"].values.tolist()
    #     x_fin = df["x_fin [m]"].values.tolist()
    #     l0 = df["l_0 [m]"].values.tolist()
    #     b_eff = df["b_eff [m]"].values.tolist()
    
    #     x_axis = []
    #     y_axis = []
    #     for count, value in enumerate(x_in):
    #         x_axis.append(x_in[count])
    #         x_axis.append(x_fin[count])
    #         y_axis.append(b_eff[count])
    #         y_axis.append(b_eff[count])
    #     fig = Figure()
    #     ax = fig.gca()
    #     ax.step(x=x_axis, y=y_axis, color='firebrick')
    #     ax.set_xlabel("x [m]")
    #     ax.set_ylabel("Effective width b_eff [m]")
    #     ax.vlines(self.spans_coord, 0, max(y_axis), colors='lightgrey',linestyles='dashed')
    #     ax.grid(color='0.8', linestyle='--', linewidth=0.5)
    #     return fig