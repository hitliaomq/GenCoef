# README

## Introduction

This is a python script to generate the coefficients($A_2$ and $A_3$ in the following fomular) for calculating the Third-Order Elastic Constants(TOECs) using first principles calculations.

$\Delta E/V_0 = A_2\eta^2 + A_3\eta^3$ 

Here, $\Delta E $ is the energy difference between the deformed structure and undeformed structure

$V_0$ is the volume of undeformed structure

$\eta$ is the Lagrangeian strain

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

## Reference

If you use this scripts, please ref the following reference

Mingqing Liao, Yong Liu, Fei Zhou, Nan Qu et.al. Strain patterns for calculating third-order elastic constants of all crystal classes from first principles, submitted to Physical Review B

In addition, when coding this script, the following reference are ref-ed.

The independent SOECs  and TOECs are taken from: BRUGGER K. Pure Modes for Elastic Waves in Crystals. Journal of Applied Physics, 1965,36(3):759-768.

The method for calculating  the TOECs is taken from:  ZHAO J, WINEY J M, GUPTA Y M. First-principles calculations of second- and third-order elastic constants for single crystals of arbitrary symmetry. PHYSICAL REVIEW B, 2007,75(0941059):94105.

## Author information

Mingqing Liao

E-mail: liaomq1900127@163.com

FGMS group @ Harbin Institute of Technology(HIT)



