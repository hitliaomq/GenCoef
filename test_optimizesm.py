#!python
# filename: test_optimizesm.py
# function: test for optimize_sm.py
#           
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)
# 

from optimize_sm import OptimizeSM
CrystalType = 'c1'
Ord = 3
flag_se = 's'
pN = 10
max_iter = 20

dim = 6

OptimizeSM(dim=dim, pN=pN, max_iter=max_iter, CrystalType=CrystalType, Ord=Ord)