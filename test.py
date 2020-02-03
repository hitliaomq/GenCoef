#!python
# filename: test.py
# function: test for symmetry.py
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)
# 

import numpy as np
import symmetry

########################################################################
#Generate the coefficients for specific strain patterns and crystal type
CrystalType = 'c1'
Ord = 3
Strain = np.array([[1, 0, 0, 0, 0, 0]])

(Cijk, StrainModeCoef, StrainMode) = symmetry.CoefForSingleMode(CrystalType, Ord, Strain)
print(StrainMode)


########################################################################
#Generate the coefficient for specific crystal type
CrystalType = 'c1'
Ord = 3

(StrainModeCoef, StrainMode) = symmetry.gen_strain_mode(CrystalType, Ord)
print(StrainMode)
print(StrainModeCoef.coef3)
print(StrainModeCoef.coef2)
