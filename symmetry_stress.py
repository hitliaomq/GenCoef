#!python
#filename: symmetry_stress.py
#function: generate strain patterns for calculating third-order elastic constant 
#          using strain-stress method
#author: Mingqing Liao
#e-mail: liaomq1900127@163.com
#FGMS group @ Harbin Institute of Technology(HIT)
from itertools import combinations, product, permutations
import numpy as np
import symmetry
import symdata
import math

def strain2deformgrad(StrainMatrix):
    Y = 2 * StrainMatrix + 1 * np.eye(3)
    (D, V) = eigv(Y)
    VT = V.T
    DeformGrad = V.dot(np.sqrt(D)).dot(VT)
    return DeformGrad


def gen_strain(strain_sca):
    #function: generate the different strain patterns, containing 0, strain_sca and strain_sca/2
    #          if you want to change the strain patterns, modify the strainprod in line 19
    strain_sca = float(strain_sca)
    Strain = np.array([[0, 0, 0, 0, 0, 0]])
    for i in range(1, 7):
        straini = list(combinations([0, 1, 2, 3, 4, 5], i))
        for j in range(0, len(straini)):
            straintmp = np.array([[0, 0, 0, 0, 0, 0]], dtype = float)
            strainindex = np.array(straini[j], dtype = int)
            strainprod = list(product([strain_sca, strain_sca/2.], repeat=i))
            #strainprod = list(product([4.*strain_sca, 3.*strain_sca, 2.*strain_sca, strain_sca, strain_sca/2., strain_sca/4.], repeat=i))
            for k in range(0, len(strainprod)):
                strainprodk = np.array(strainprod[k])
                #print strainprodk
                straintmp[0, strainindex] = strainprodk         
                straintmp[0, 3:] = 2 * straintmp[0, 3:]
                Strain = np.append(Strain, straintmp, axis = 0)
    #Strain = np.append(Strain, np.array([[1, 0.5, 0, 0, 0, 0]]), axis = 0)
    Strain = np.delete(Strain, 0, 0)
    return Strain

def cauchy2pk(Cauchy, LagStrain):
    DeformGrd = strain2deformgrad(LagStrain)
    DeformGrdInv = np.linalg.inv(DeformGrd)
    PK = np.linalg.det(DeformGrd) * DeformGrdInv.dot(Cauchy).dot(DeformGrdInv)
    #print DeformGrd
    #print DeformGrdInv
    #print PK
    #PKvec = np.zeros((1, 6))
    PKvec = np.array([PK[0,0], PK[1,1], PK[2,2], PK[1,2], PK[0, 2], PK[0,1]])
    return (PK, PKvec)

def gen_strain_coef(StrainOrd, StrainOrdCoef, Ord = 3):
    #function: get the non-zero second- or third-(specified by the Ord parameter) order 
    #          elastic constant(CoefUniq) and corresponding coefficient(CoefUniqCoef)
    #          the input parameter can get by gen_strainord
    #          Ord the order, here 2 or 3
    n_strain = len(StrainOrd)
    N = int(n_strain ** (Ord - 1))
    #print N
    Coef = np.ones((6, N))
    Cijk = np.zeros((6, N))
    for i in range(0, N):
        IJK = symmetry.Num2IJK(i + 1, n_strain, int(Ord) - 1)
        for j in range(0, 6):
            coeftmp = np.zeros(int(Ord))
            coeftmp[int(Ord) - 1] = j + 1
            for i_ord in range(0, int(Ord) - 1):
                #print IJK[i_ord]
                coeftmp[i_ord] = StrainOrd[IJK[i_ord] - 1]
                Coef[j, i] = Coef[j, i] * StrainOrdCoef[IJK[i_ord] - 1]
            coeftmp = np.sort(coeftmp)
            for k in range(0, int(Ord)):
                Cijk[j, i] = Cijk[j, i] + coeftmp[k] * (10 ** (Ord - 1 - k))
    CoefUniqtmp, CoefIndextmp = np.unique(Cijk[0, :], return_index=True)
    n_uniq = len(CoefUniqtmp)
    CoefUniq = np.zeros((6, n_uniq))
    CoefUniqCoef = np.zeros((6, n_uniq))
    for i in range(0, 6):
        CoefUniqi, CoefIndex = np.unique(Cijk[i, :], return_index=True)
        CoefUniq[i, :] = CoefUniqi
        Coefi = Coef[i, CoefIndex]     
        for j in range(0, n_uniq):
            CoefUniqCoef[i, j] = sum(Cijk[i, :] == CoefUniqi[j])
        CoefUniqCoef[i, :] = CoefUniqCoef[i, :] * Coefi
    return (CoefUniq, CoefUniqCoef)

def gen_cijk_coef(CrystalType, Strain, Ord):
    #function: generate the coefficients of independent elastic constant according to 
    #          the crystal type(CrystalType) and the strain(Strain) and order(Ord)
    (StrainOrd, StrainOrdCoef) = symmetry.gen_strainord(Strain)
    CijkInd = symdata.coef_ind(CrystalType, Ord)
    CoefCoef = symdata.coef_crystal(CrystalType, Ord)              
    (Cijk, CijkCoef) = gen_strain_coef(StrainOrd, StrainOrdCoef, Ord)    
    CijkUniq = np.array(symdata.get_unique(CijkInd))
    CoefStrain = np.zeros((6, len(CijkUniq)))
    for k in range(0, 6):
        CijkNum = symmetry.cijk2num(Cijk[k, :], Ord)
        for i in range(0, len(CijkNum)):
            COrderi = CijkInd[CijkNum[i]]
            #print type(COrderi)
            #print CoefCoef[CijkNum[i]]
            if type(COrderi) is list:
                for j in range(0, len(COrderi)):
                    #print(i)
                    CoefOrderi = float(CoefCoef[CijkNum[i]][j]) * float(CijkCoef[k][i])
                    Index_CoefStrain = np.argwhere(CijkUniq == COrderi[j])
                    Index_CoefStrain = Index_CoefStrain[0][0]
                    CoefStrain[k, Index_CoefStrain] = CoefStrain[k, Index_CoefStrain] + CoefOrderi
            else:
                if COrderi != 0:
                    #print CoefCoef[CijkNum[i]]
                    CoefOrderi = float(CoefCoef[CijkNum[i]]) * float(CijkCoef[k, i])
                    #print CoefOrderi
                    Index_CoefStrain = np.argwhere(CijkUniq == COrderi)
                    Index_CoefStrain = Index_CoefStrain[0][0]
                    CoefStrain[k, Index_CoefStrain] = CoefStrain[k, Index_CoefStrain] + CoefOrderi
    return CoefStrain

def CoefForSingleMode(CrystalType, Ord, Strain):
    #function: generate the coefficient for a specified strain and crystal type
    #
    #Cijk = print_cijk(CrystalType, Ord)
    #print Cijk
    m = Strain.shape
    if len(m) == 1:
        m = 1
    else:
        m = m[0]
    StrainMode = np.zeros((m, 3, 3))
    StrainModeCoef = symmetry.CoefStr(Ord)
    for i in range(2, int(Ord) + 1):
        CijUniq = np.array(symdata.get_unique(symdata.coef_ind(CrystalType, i)))
        n_uniq = len(CijUniq)
        exec('StrainModeCoef.coef' + str(i) + ' = np.zeros((int(6*m), n_uniq), dtype = float)')
    #StrainModeCoef = CoefStr(Ord)
    Cijk = symmetry.CoefStr(Ord)
    for i in range(2, int(Ord) + 1):
        exec('Cijk.coef' + str(i) + ' = symmetry.print_cijk(CrystalType, i)')
        for j in range(0, m):
            StrainM = symmetry.vec2matrix(Strain[j, :])
            StrainMode[j, :, :] = StrainM
            #print StrainM
            exec('StrainModeCoef.coef' + str(i) + '[6*j:6*(j+1), :] = gen_cijk_coef(CrystalType, StrainM, i)')
        exec('print(StrainModeCoef.coef' + str(i) + ')')
    return (Cijk, StrainModeCoef, StrainMode)

def gen_fullrank_db(CrystalType, Ord):
    #function: generate full rank strain mode for 
    #          current crystal type(CrystalType) and order(Ord)
    #StrainAll = symmetry.gen_strain(1)
    StrainAll = gen_strain(1)
    n_strain = len(StrainAll)
    StrainMode_tmp = np.zeros((n_strain, 3, 3), dtype = float)
    count = 0
    for i in range(0, len(StrainAll)):
        Straini = symmetry.vec2matrix(StrainAll[i])
        StrainModei = gen_cijk_coef(CrystalType, Straini, Ord)
        if np.linalg.matrix_rank(StrainModei) == 6:
            # This value is to contral the strain need to be calculated
            StrainMode_tmp[count, :, :] = Straini
            #print Straini
            count = count + 1
    StrainMode = StrainMode_tmp[0:count, :, :]
    return StrainMode

def gen_alldiff_db():
    strainprod = list(permutations([2, 1, 3, 2, 2, 2], 6))
    n_strain = len(strainprod)
    StrainMode = np.zeros((n_strain, 3, 3), dtype = float)
    for i in range(0, n_strain):
        Straini = symmetry.vec2matrix(np.array(strainprod[i]))
        StrainMode[i, :, :] = Straini
    return StrainMode

def gen_strain_mode(CrystalType, Ord):
    #function: generate the strain mode for given crystal type(CrystalType) and order(Ord)
    CijkInd = symdata.coef_ind(CrystalType, Ord)
    CoefCoef = symdata.coef_crystal(CrystalType, Ord)
    CijkUniq = np.array(symdata.get_unique(CijkInd))
    n_uniq = len(CijkUniq)
    #print(n_uniq)
    if (Ord % 2.0 == 0) and (n_uniq > 5):
        #n_uniq > 5 keep doubt? to do
        n_row = int(math.ceil(n_uniq/6.)*6.) #+ 6
    else:
        n_row = int(math.ceil(n_uniq/6.)*6.)
    #print(n_row)
    Cijk = symmetry.num2cijk(CijkUniq, Ord)
    StrainMode = np.zeros((int(math.ceil(n_row/6.)), 3, 3), dtype = float)

    #StrainModeCoef = np.zeros((n_row, n_uniq), dtype = float)
    StrainModeCoef = symmetry.CoefStr(Ord)
    exec('StrainModeCoef.coef' + str(int(Ord)) + ' = np.zeros((n_row, n_uniq), dtype = float)')
    for i in range(2, int(Ord)):
        CijUniq = np.array(symdata.get_unique(symdata.coef_ind(CrystalType, i)))
        n_uniqi = len(CijUniq)
        exec('StrainModeCoef.coef' + str(i) + ' = np.zeros((n_row, n_uniqi), dtype = float)')

    if n_uniq < 6:
        #The ord must be 2
        StrainAll = gen_strain(1)
        for i in range(0, len(StrainAll)):
            StrainModei = symmetry.vec2matrix(StrainAll[i])
            StrainModeCoefi = gen_cijk_coef(CrystalType, StrainModei, Ord)
            SMRank = np.linalg.matrix_rank(StrainModeCoefi)
            if SMRank == n_uniq:
                StrainMode[0, :, :] = StrainModei
                StrainModeCoef = StrainModeCoefi
                break
    else:
        print("Begin of generate the full-rank database.\n")
        StrainAll = gen_fullrank_db(CrystalType, Ord)
        print("The full-rank database were generated.\n")
        #StrainAll = gen_alldiff_db()
        #print StrainAll
        n_strain = len(StrainAll)
        #print(n_strain)
        for j in range(0, n_strain):
            count = 0
            flag = 0
            for i in range(j, n_strain):
                Straini = StrainAll[i, :, :]
                StrainModei = gen_cijk_coef(CrystalType, Straini, Ord)

                if count == 0:
                    exec('StrainModeCoef.coef' + str(int(Ord)) + '[0:6, :] = StrainModei')
                    #StrainModeCoef[0:6, :] = StrainModei
                    StrainMode[count, :, :] = Straini
                    for k in range(2, int(Ord)):
                        exec('StrainModeCoef.coef' + str(int(k)) + '[0:6, :] = gen_cijk_coef(CrystalType, Straini, k)')
                    count = count + 1
                    #'''
                    if n_uniq == 6:
                        flag = 1
                        break
                    #'''
                else:
                    StrainModetmp = np.zeros((6*(count + 1), n_uniq), dtype = float)
                    StrainModetmp[0:6*count, :] = eval('StrainModeCoef.coef' + str(int(Ord)) + '[0:6*count, :]')
                    StrainModetmp[6*count:6*(count + 1), :] = StrainModei
                    SMRank = np.linalg.matrix_rank(StrainModetmp)
                    #print(SMRank)
                    #print(n_uniq)
                    #print(count)
                    ### NOTE: the default is SMRank == 6*(count + 1)
                    RnkLst = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120, 126]
                    #if SMRank == 6*(count + 1) or SMRank == n_uniq:
                    if SMRank > (RnkLst[count] - 1) or SMRank == n_uniq:
                    #if SMRank > 5*(count + 1) or SMRank == n_uniq:
                        #exec('print(StrainModeCoef.coef' + str(int(Ord)) + ')')
                        #print(StrainModei)
                        exec('StrainModeCoef.coef' + str(int(Ord)) + '[6*count:6*(count + 1), :] = StrainModei')
                        StrainMode[count, :, :] = Straini
                        for k in range(2, int(Ord)):
                            exec('StrainModeCoef.coef' + str(k) + '[6*count:6*(count + 1), :] = gen_cijk_coef(CrystalType, Straini, k)')
                        count = count + 1
                        if SMRank == n_uniq:
                            print(SMRank)
                            flag = 1
                            break
            if flag == 1:
                break
    return (StrainModeCoef, StrainMode)
