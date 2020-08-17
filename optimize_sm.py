#!python

import symmetry
import symmetry_stress
import matplotlib.pyplot as plt
import numpy as np
import random
import math
import os

def read_strainmode(STRAINMODE = "STRAINMODE"):
    strainmodetmp = np.loadtxt(STRAINMODE)
    m = strainmodetmp.shape
    if len(m) == 1:
        strainmode = np.zeros((1, 6))
        strainmode[0, :] = strainmodetmp
    else:
        strainmode = strainmodetmp
    return strainmode

def OptimizeSM(dim=6, pN=50, max_iter=200, CrystalType='c1', Ord=3, flag_se='s', flag_fig=1):
    my_pso = PSOFit(pN=pN, dim=dim, max_iter=max_iter, CrystalType=CrystalType, Ord=Ord, flag_se=flag_se)
    my_pso.init_Population()
    (fitness, gbest) = my_pso.iterator()

    print(gbest)
    print(fitness[-1])

    np.savetxt('RESULTSM', gbest.reshape((-1, 6)))

    if flag_fig:
        plt.title("Figure1")
        plt.xlabel("iterators", size=14)
        plt.ylabel("fitness", size=14)
        t = np.array([t for t in range(0,len(fitness))])
        fitness = np.array(fitness)
        plt.plot(t, fitness, color='b',linewidth=3)
        plt.show()

def fitfun(x, CrystalType='c1', Ord=3, flag_se='s'):
    #Ord = 0, default means Ord = len(x)
    x = np.array(x).reshape((-1, 6))
    if flag_se == 's':
        Cij_mode, coef_e, StrainMode = symmetry_stress.CoefForSingleMode(CrystalType, Ord, x)
    else:
        Cij_mode, coef_e, StrainMode = symmetry.CoefForSingleMode(CrystalType, Ord, x)
    if Ord == 2:
        coef = coef_e.coef2
    elif Ord == 3:
        coef = coef_e.coef3
    elif Ord == 4:
        coef = coef_e.coef4
    matrixrank = np.linalg.matrix_rank(coef)
    print('The rank of the coef matrix is {}'.format(matrixrank))
    if matrixrank == get_Nindepen(CrystalType=CrystalType, Ord=Ord):
        y = np.linalg.cond(coef)
    else:
        y = 1e10
    print('The condition number of the coef matrix is {}'.format(y))
    return y
    
#the part of PSO ref https://blog.csdn.net/ztf312/article/details/75669685
class PSOFit(object):
    """docstring for PSOFit"""
    def __init__(self, pN=100, dim=6, max_iter=1000, CrystalType='c1', Ord=3, flag_se='s'):
        super(PSOFit, self).__init__()
        self.w0 = 0.9
        self.w1 = 0.4
        self.c1 = 2
        self.c2 = 2
        self.r1 = 0.6
        self.r2 = 0.3
        self.pN = pN    #number of particle
        self.dim = dim  #dimension, according to the x
        self.max_iter = max_iter    #maximum iter
        self.X = np.random.rand(self.pN, self.dim)      #the positionn
        self.V = np.zeros((self.pN, self.dim))      #the v
        self.pbest = np.zeros((self.pN, self.dim))  #the position of best_particle
        self.gbest = np.zeros((1, self.dim))        #the position of best_group
        self.p_fit = np.zeros(self.pN)              #the best value of particle
        self.tor = 1e-5
        self.fit = 1e10
        self.CrystalType = CrystalType
        self.Ord = Ord
        self.flag_se = flag_se

    def init_Population(self, Nmin = -10, Nmax = 10, tor=0.5, InitMethod='tor', v_scale=10.0):
        if os.path.exists('INITSM'):
            xtmp = np.loadtxt('INITSM')
            if not xtmp.size == 0:
                x0 = xtmp.reshape((1, -1))[0]
        else:
            x0 = np.random.random((self.dim))
        for i in range(self.pN):
            for j in range(self.dim):
                if i == 0:
                    self.X[i][j] = x0[j]
                    self.V[i][j] = x0[j] * np.random.random() * v_scale
                else:
                    if InitMethod.lower().startswith('r'):
                        self.X[i][j] = random.uniform(Nmin, Nmax)
                        self.V[i][j] = random.uniform(Nmin, Nmax)
                    elif InitMethod.lower().startswith('t'):
                        self.X[i][j] = x0[j] + np.random.random() * tor
                        self.V[i][j] = x0[j] * np.random.random() * v_scale
            self.pbest[i] = self.X[i]
            tmp = fitfun(self.X[i], self.CrystalType, self.Ord, self.flag_se)
            self.p_fit[i] = tmp
            if(tmp < self.fit):
                self.fit = tmp
                self.gbest = self.X[i]

    def iterator(self):
        fitness = []
        for t in range(self.max_iter):
            w = (self.w0 - self.w1) * (self.max_iter - t)/self.max_iter + self.w1
            for i in range(self.pN):
               tmp = fitfun(self.X[i], self.CrystalType, self.Ord, self.flag_se)
               if(tmp<self.p_fit[i]):
                   self.p_fit[i] = tmp
                   self.pbest[i] = self.X[i]
                   if(self.p_fit[i] < self.fit):
                       self.gbest = self.X[i]
                       self.fit = self.p_fit[i]
            for i in range(self.pN):
                self.V[i] = w*self.V[i] + self.c1*self.r1*(self.pbest[i] - self.X[i]) + \
                            self.c2*self.r2*(self.gbest - self.X[i])
                self.X[i] = self.X[i] + self.V[i]
            fitness.append(self.fit)
            print(self.fit)
            #'''
            if t == 0:
                pass
            else:
                if np.max(self.p_fit) - np.min(self.p_fit) < self.tor:
                    break
            #'''
        return(fitness, self.gbest)

def simplify_SM(filename='RESULTSM', CrystalType='c1', Ord=3, flag_se='s'):
    x = read_strainmode(filename)
    y = fitfun(x, CrystalType=CrystalType, Ord=Ord, flag_se=flag_se)
    print(y)
    xabs = np.absolute(x)
    print(xabs)
    x1 = [[round(itemi/np.max(xabs), 2) for itemi in item] for item in x]
    y1 = fitfun(np.array(x1), CrystalType=CrystalType, Ord=Ord, flag_se=flag_se)
    print(y1)
    with open('FINALSM', 'w+') as fopen:
        for xi in x1:
            xi_str = [str(item) for item in xi]
            fopen.write('  '.join(xi_str) + '\n')
    print(x1)
    for i, x1i in enumerate(x1):
        for j, x1ij in enumerate(x1i):
            if abs(x1ij) < 0.1:
                x1[i][j] = 0.
    y2 = fitfun(np.array(x1), CrystalType=CrystalType, Ord=Ord, flag_se=flag_se)
    print(y2)
    print(x1)

def write_matrix(matrix, dec=3):
    '''
    Here matrix is list, not np.ndarray
    '''
    for xi in matrix:
        xi = [format(item, '.{}f'.format(dec)) for item in xi]
        print('  '.join(xi))

def trans_resultsm(filename='RESULTSM'):
    x = read_strainmode(filename)
    xmax = np.max(np.absolute(x))
    for xi in x:
        xi_m = symmetry.vec2matrix(xi.tolist())/xmax
        #xi_m = xi.tolist()/xmax
        write_matrix(xi_m.tolist())

def get_Nindepen(CrystalType='c1', Ord=3):
    Ctype = ['c1', 'c2', 'h1', 'h2', 't1', 't2', 'r1', 'r2', 'o', 'm', 'n']
    ind = Ctype.index(CrystalType)
    if Ord == 2:
        Nindepen = [3, 3, 5, 5, 6, 7, 6, 7, 9, 13, 21]
    elif Ord == 3:
        Nindepen = [6, 8, 10, 12, 12, 16, 14, 20, 20, 32, 56]
    elif Ord == 4:
        Nindepen = [11, 14, 19, 24, 25, 36, 28, 42, 42, 70, 126]
    return int(Nindepen[ind])

def get_dim(CrystalType='c1', Ord=3):
    N = get_Nindepen(CrystalType=CrystalType, Ord=Ord)
    if flag_se == 's':
        dim = math.ceil(N/6.) * 6
    else:
        dim = int(6)
    return dim
