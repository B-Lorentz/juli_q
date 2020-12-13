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
checky(-0.6831300510639732, CG(  S(2)/S(1) ,  S(-1)/S(1) ,  S(4)/S(1) ,  S(4)/S(1) ,  S(3)/S(1) ,  S(3)/S(1) ).doit().evalf(),' ')
checky(0.5916079783099616, CG(  S(2)/S(1) ,  S(0)/S(1) ,  S(4)/S(1) ,  S(3)/S(1) ,  S(3)/S(1) ,  S(3)/S(1) ).doit().evalf(),' ')
checky(-0.3872983346207417, CG(  S(2)/S(1) ,  S(1)/S(1) ,  S(4)/S(1) ,  S(2)/S(1) ,  S(3)/S(1) ,  S(3)/S(1) ).doit().evalf(),' ')
checky(0.18257418583505536, CG(  S(2)/S(1) ,  S(2)/S(1) ,  S(4)/S(1) ,  S(1)/S(1) ,  S(3)/S(1) ,  S(3)/S(1) ).doit().evalf(),' ')
print('All tests passed: ' , allok)
