# README

## Introduction

This is a python script to generate the coefficients($A_2$ and $A_3$ in the following fomular) for calculating the Third-Order Elastic Constants(TOECs) using first principles calculations by strain-energy method or strain-stress method.

Strain-energy method

$\Delta E/V_0 = A_2\eta^2 + A_3\eta^3$ 

Here, $\Delta E $ is the energy difference between the deformed structure and undeformed structure

$V_0​$ is the volume of undeformed structure

$\eta​$ is the Lagrangeian strain

Strain-stress method

$t = A_2\eta + A_3\eta^2$

[Download from GitHub](https://github.com/hitliaomq/GenCoef)

## Usage

### Import the script

```python
import numpy as np
import symmetry
```

### Generate all the strain modes and corresponding coefficients for specific crystal type

```python
CrystalType = 'c1'
Ord = 3

(StrainModeCoef, StrainMode) = symmetry.gen_strain_mode(CrystalType, Ord)
#StrainMode is the strain modes, in 3 by 3 matrix
print(StrainMode)
#The StrainModeCoef is a new defined class in the script
#StrainModeCoef.coef3 is the coefficients for TOECs
#StrainModeCoef.coef2 is the coefficients for SOECs
#If the Ord=2, then there is only *.coef2
print(StrainModeCoef.coef3)
print(StrainModeCoef.coef2)
```

### Generate the coefficients for specific strain mode

```python
CrystalType = 'c1'
Ord = 3
Strain = np.array([[1, 0, 0, 2, 2, 0]])

(Cijk, StrainModeCoef, StrainMode) = symmetry.CoefForSingleMode(CrystalType, Ord, Strain)
#The Cijk is the independent TOECs
#Cijk.coef2 and Cijk.coef3
print(StrainMode)
```

For more usage, please ref. `test_symmetry.py` and `test_optimizesm.py`

## Reference

If you use this scripts, please ref the following reference

[1] Mingqing Liao, Yong Liu, Fei Zhou, Tianyi Han, et al. Selection of strain and fitting schemes for calculating higher-order elastic constants, submitted to *Physical Review Letters*

[2] Mingqing Liao, Yong Liu, Fei Zhou, Tianyi Han, et al. Comparison between methods for third-order elastic constants from first-principles, submitted to *Physical Review B*

In addition, when coding this script, the following reference are ref-ed.

The independent SOECs  and TOECs are taken from: 
[3] BRUGGER K. Pure Modes for Elastic Waves in Crystals. Journal of Applied Physics, 1965,36(3):759-768.

The method for calculating  the TOECs is taken from:  
[4] ZHAO J, WINEY J M, GUPTA Y M. First-principles calculations of second- and third-order elastic constants for single crystals of arbitrary symmetry. PHYSICAL REVIEW B, 2007,75(0941059):94105.

This code now has been a part of [Elastic3rd](https://github.com/hitliaomq/ELASTIC3rd), Please ref:
[5] Mingqing Liao, Yong Liu, Shun-Li Shang, Fei Zhou, et al., Elastic3rd: A tool for calculating third-order elastic constants from first-principles calculations, submitted to *Computer Physics Communication*

## Author information

Mingqing Liao

E-mail: liaomq1900127@163.com

FGMS group @ Harbin Institute of Technology(HIT)



