from sympy import *
from sympy.physics.quantum.spin import CG
import numpy as np
allok = True
def checky(a, b, log):
   good = np.allclose(a, float(b))
   global allok
   if not good:
       allok = False
   print(f'{a:.6f}', f'{float(b):.6f}', good, log)
