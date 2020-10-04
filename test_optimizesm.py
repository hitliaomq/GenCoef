#!python
# filename: test_optimizesm.py
# function: test for optimize_sm.py
#           
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)
# 

from optimize_sm import OptimizeSM
CrystalType = 'n'
Ord = 3
flag_se = 'e'
pN = 20
max_iter = 50

dim = 336

OptimizeSM(dim=dim, pN=pN, max_iter=max_iter, CrystalType=CrystalType, Ord=Ord, flag_se=flag_se)