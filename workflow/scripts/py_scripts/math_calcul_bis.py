from math_calcul import Mathematica_1
import sys

init=Mathematica_1(sys.argv[1])
factor_a, factor_b=init.return_factors()
init.calcul_depth_min(sys.argv[2], factor_a, factor_b)