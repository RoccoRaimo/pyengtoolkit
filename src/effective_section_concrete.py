"""
pyengtoolkit - effective_section_composite
---
Calculation of effective width composite steel-concrete sections based on prescription of Eurocode 4.
"""

from dataclasses import dataclass
import pandas as pd
import numpy as np
from matplotlib.figure import Figure

def b_eff_i(l0: float, b_i: float):
    """
    A function that takes length_0 and b_i and returns the i-th b_eff.
    """
    if (0.2 * b_i + 0.1 * l0) <= 0.2 * l0:
        if (0.2 * b_i + 0.1 * l0) <= b_i:
            b_eff_i = (0.2 * b_i + 0.1 * l0)
        else:
            b_eff_i = b_i
        
    else:
        if 0.2 * l0 <= b_i:
            b_eff_i = 0.2 * l0
        else:
            b_eff_i = b_i
    return b_eff_i

def b_eff(b_w: float, b_eff_1: float, b_eff_2: float):
    """
    A function that takes b_w, b_ff_1 and b_eff_2 and returns the total b_eff.
    """
    b_eff = b_w + b_eff_1 + b_eff_2
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
    b_1: list
    b_2: list
    b_w: list
    
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
    def length_0(self):
        """
        A function that defines the zones in which the length_0 is constant.
        """
        if self.n_spans == 1:
            spans_coord_l0 = []
            spans_coord_l0.append(self.spans_length[0])
        if self.n_spans == 2:
            spans_coord_l0 = [0]
            if self.spans_type[0] == 1 and self.spans_type[1] == 1:
                first_coord  = 0.85 * self.spans_length[0]
                second_coord = 0.15 * self.spans_length[1]
            if self.spans_type[0] == 3 and self.spans_type[1] == 1:
                first_coord  = self.spans_length[0]
                second_coord = 0.15 * self.spans_length[1]
            if self.spans_type[0] == 1 and self.spans_type[1] == 3:
                first_coord  = 0.85 * self.spans_length[0]
                second_coord = 1.00 * self.spans_length[1]
            spans_coord_l0.append(first_coord)
            spans_coord_l0.append(self.spans_coord[1]+second_coord)
            spans_coord_l0.append(self.spans_coord[2])
        else:
            #Start span
            if self.spans_type[0] == 1:
                spans_coord_l0 = [0]
                start_span_coord  = 0.85 * self.spans_length[0]
                spans_coord_l0.append(start_span_coord)
            elif self.spans_type[0] == 2:
                raise ValueError(f"The span type for end span must be 1 or 3, not {self.spans_type[0]}") 
            elif self.spans_type[0] == 3:
                spans_coord_l0 = []
                start_span_coord  = 0
                spans_coord_l0.append(start_span_coord)
            #Central spans
            for count, span_type in enumerate(start=1, iterable=self.spans_type[1:-1:]):
                if span_type != 2:
                    raise ValueError(f"The span type for central span must be 2, not {span_type}") 
                else:
                    first_coord = 0.15 * self.spans_length[count]
                    second_coord = 0.70 * self.spans_length[count]
                    delta_length_1 = first_coord
                    span_coord_l0_1 = self.spans_coord[count] + delta_length_1
                    spans_coord_l0.append(span_coord_l0_1)
                    delta_length_2 = second_coord
                    span_coord_l0_2 = span_coord_l0_1 + delta_length_2
                    spans_coord_l0.append(span_coord_l0_2)
            #End span
            if self.spans_type[-1] == 1:
                end_span_coord  = 0.15 * self.spans_length[-1]
            elif self.spans_type[-1] == 2:
                raise ValueError(f"The span type for end span must be 1 or 3, not {self.spans_type[-1]}") 
            elif self.spans_type[-1] == 3:
                end_span_coord  = 0
            delta_length_end = end_span_coord 
            span_coord_l0_end = self.spans_coord[-2] + delta_length_end    
            spans_coord_l0.append(span_coord_l0_end)
            spans_coord_l0.append(self.spans_coord[-1])    
        return spans_coord_l0

    def b_eff_table(self):
        """
        A function that returns a table that contains all the informations related to effective width of section,
        along all the segments of the bridge.
        """
        #First the coordinates of spans, length_0 and total coordinates are defined
        x_span = self.spans_coord
        x_l0 = self.length_0
        x_tot = np.unique(x_span + x_l0).tolist()

        #Initial and final coordinates for every segment of spans coord
        x_span_in = []
        x_span_fin = []
        for count, coord in enumerate(x_span[0:(len(x_span)-1)]):
            x_span_in.append(x_span[count])
            x_span_fin.append(x_span[count+1])

        #Initial and final coordinates for every segment of length_0 segments
        x_l0_in = []
        x_l0_fin = []
        if self.n_spans == 1:
            x_l0_in = x_span_in
            x_l0_fin = x_span_fin
        else:
            for count, coord in enumerate(x_l0[0:(len(x_l0)-1)]):
                x_l0_in.append(x_l0[count])
                x_l0_fin.append(x_l0[count+1])

        #Initial and final coordinates for every segment of total segments
        x_tot_in = []
        x_tot_fin = []
        for count, coord in enumerate(x_tot[0:(len(x_tot)-1)]):
            x_tot_in.append(x_tot[count])
            x_tot_fin.append(x_tot[count+1])
        
        #Number of total segments
        l0_tot_segments = []
        for count in range(1,len(x_tot)):
            l0_tot_segments.append(count)

        #Values of length_0
        l0 = []
        for count, x in enumerate(x_l0_in):
            l0_i = x_l0_fin[count] - x_l0_in[count]
            l0.append(l0_i)

        #Values of length 0 for every segment of x_tot
        l0_tot = []
        #Values of l0 in x_tot
        for value1 in x_tot_fin:
            for count, value2 in enumerate(x_l0_fin):
                if value1 <= value2:
                    l0_tot.append(l0[count])
                    break
        
        #Values of b_w for every segment of x_tot
        b_w_tot = []
        #Values of length_0 in x_tot
        for value1 in x_tot_fin:
            for count, value2 in enumerate(x_span_fin):
                if value1 <= value2:
                    b_w_tot.append(self.b_w[count])
                    break

        b_1_tot = []
        b_2_tot = []
        #Values of b_1 in x_tot
        for value1 in x_tot_fin:
            for count, value2 in enumerate(x_span_fin):
                if value1 <= value2:
                    b_1_tot.append(self.b_1[count])
                    b_2_tot.append(self.b_2[count])
                    break

        b_eff_1 = []
        b_eff_2 = []
        for count, length in enumerate(l0_tot):
            b_eff_1_i = b_eff_i(length,b_1_tot[count])
            b_eff_2_i = b_eff_i(length,b_2_tot[count])
            b_eff_1.append(b_eff_1_i)
            b_eff_2.append(b_eff_2_i)
        
        b_eff_list = []
        for count, b in enumerate(b_eff_1):
            b_eff_list.append(b_eff(b_w_tot[count], b_eff_1[count], b_eff_2[count]))

        #Dictionary for the DataFrame is defined
        dict = {"Segment": l0_tot_segments,
                "x_in [m]": x_tot_in,
                "x_fin [m]": x_tot_fin,
                "l_0 [m]": l0_tot,
                "b_1 [m]": b_1_tot,
                "b_2 [m]": b_2_tot,
                "b_eff_1 [m]": b_eff_1,
                "b_eff_2 [m]": b_eff_2,
                "b_w [m]": b_w_tot,
                "b_eff [m]": b_eff_list
                }

        df = pd.DataFrame(data=dict)
        df_indexed = df.set_index(keys="Segment")
        return df_indexed #l0_tot_segments, x_l0_in, x_l0_fin, l0, l0_tot, b_1_tot, b_eff_1, b_eff_list

    def plot_b_eff(self):
        """
        A function that plot the values of b_eff throughout the axis of the bridge.
        """
        #First the data are defined
        df = self.b_eff_table()
        x_in = df["x_in [m]"].values.tolist()
        x_fin = df["x_fin [m]"].values.tolist()
        l0 = df["l_0 [m]"].values.tolist()
        b_eff = df["b_eff [m]"].values.tolist()


        x_axis = []
        y_axis = []
        for count, value in enumerate(x_in):
            x_axis.append(x_in[count])
            x_axis.append(x_fin[count])
            y_axis.append(b_eff[count])
            y_axis.append(b_eff[count])
        fig = Figure()
        ax = fig.gca()
        ax.step(x=x_axis, y=y_axis, color='firebrick')
        ax.set_xlabel("x [m]")
        ax.set_ylabel("Effective width b_eff [m]")
        ax.vlines(self.spans_coord, 0, max(y_axis), colors='lightgrey',linestyles='dashed')
        ax.grid(color='0.8', linestyle='--', linewidth=0.5)
        return fig