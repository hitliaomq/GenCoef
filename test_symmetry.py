#!python
# filename: test_symmetry.py
# function: test for the coefficients(include A2 and A3) generation
#           in symmetry module.
#           
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)
# 

import numpy as np
import symmetry
import symmetry_stress

###########################################################################
#TEST FOR SYMMETRY MODULE
#Two main parts:
#    1. test for symmetry script
#    2. test for symmetry_stress script
#In each test of script, there are two main functions:
#    a. test for a given strain mode
#    b. test for a given symmetry

print("####################Start test for symmetry module#########################")

print("*************Start of the test for symmetry script*************")
print("1. This is the test for a specific strain mode.")
# Generate the coefficients for a given strain mode and specific crystal type
CrystalType = 'c1'
Ord = 3
Strain = np.array([[0, 1, 0, 0, 1, 0]])
(Cijk,StrainModeCoef,StrainMode)=symmetry.CoefForSingleMode(CrystalType, Ord, Strain)
print(StrainMode)
print("------End of test for a specific strain mode.------\n")

print("2. This is the test for a specific crystal symmetry.")
#Generate the coefficient for specific crystal type
CrystalType = 'c1'
Ord = 3
(StrainModeCoef, StrainMode) = symmetry.gen_strain_mode(CrystalType, Ord)
print(StrainMode)
print(StrainModeCoef.coef3)


print(StrainModeCoef.coef2)
print("------End of test for a given symmetry.------\n")
print("*************End of test for symmetry script*************")
print("\n\n")
#End of test for symmetry script


print("*************Start of test for symmetry_stress script*************")
print("1. This is the test for a specific strain mode.")
# Generate the coefficients for a given strain mode and specific crystal type
CrystalType = 'c1'
Ord = 3
#Strain = np.array([[1, 0.5, 0, 2, 2, 0]])
Strain = np.array([[1, 10, 0, 2, 2, 0]])
(Cijk,StrainModeCoef,StrainMode)=symmetry_stress.CoefForSingleMode(CrystalType, Ord, Strain)
print(StrainMode)
coef3_inv = np.linalg.inv(StrainModeCoef.coef3)
err_f = np.linalg.norm(StrainModeCoef.coef3) * np.linalg.norm(coef3_inv)
print(err_f)
print("------End of test for a specific strain mode.------\n")

print("2. This is the test for a specific crystal symmetry.")
#Generate the coefficient for specific crystal type
CrystalType = 'c1'
Ord = 3
(StrainModeCoef, StrainMode) = symmetry_stress.gen_strain_mode(CrystalType, Ord)
print(StrainMode)
print(StrainModeCoef.coef3)
print(StrainModeCoef.coef2)


print("------End of test for a given symmetry.------\n")
print("*************End of test for symmetry_stress script*************")
