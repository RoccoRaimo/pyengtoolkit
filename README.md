# pyengtoolkit (Python Engineering Toolkit)

A kit of useful python tools for everyday calculations in **structural engineering** journey!


# How to install
~~You can install the library using pip~~ (WORK IN PROGRESS):

```python
pip install pyengtoolkit
```

You can choose the single structural module to use in your python files or jupyter notebooks importing with:

```python
from pyengtoolkit.src import 'name of the package that you want to use'
```


# List of tools available
- **effective_section_concrete**:
calculation of effective width for T, double-T and L shaped concrete sections, based on prescriptions of Eurocode 2 (EN 1992-1-1 §5.3.2.1);
- **effective_section_composite**:
calculation of effective width for steel-concrete sections, based on prescriptions of Eurocode 4 (EN 1994-1-1 §5.4.1.2);
- **concrete_shear**:
list of functions that calculate shear resistance (according to EC-NTC18-NBR) and render them in LaTeX with `handcalcs` library;
- **steel_beam**: 
list of functions that calculate resistance and stability for steel beams (according to EC) and render them in LaTeX with `handcalcs` library;
- **stay_cables**: calculation of stay cables elastic and geometric elongations;
- **PT_losses**: calculation of Post-Tension Losses for prestressing steel, according to EC.   
It is possible to consider:
    - Immediate elastic shortening of Concrete
    - Friction
    - Wedge draw-in of anchorage devices
    - Shrinkage of Concrete
    - Creep of Concrete
    - Relaxation of Prestressing Steel

# List of tools planned
- 

---

# Examples
You can find some examples of usage for every tool in examples folder!
