# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 17:26:07 2020

@author: EngenhariaCecal
"""
import numpy as np
import os
import pandas as pd

# =============================================================================
# # ===========================================================================
# # define data and result directory
# # ===========================================================================
# =============================================================================
maindirectory = os.getcwd()
resultdirectory=("results")
datadirectory=("data")
if not os.path.exists(resultdirectory):
    os.mkdir(resultdirectory)

# =============================================================================
# # ===========================================================================
# # Adding material to database based on elements
# # ===========================================================================
# =============================================================================
dfMoleMass=pd.read_csv(datadirectory+'/'+'MolecularMass.csv')
dfMoleMass['Type']=2 # this is added for dfdivide function to work. arbitrary 
dfMoleMass['CaO']=dfMoleMass['Ca']+dfMoleMass['O']
dfMoleMass['SiO2']=dfMoleMass['Si']+2*dfMoleMass['O']
dfMoleMass['MgO']=dfMoleMass['Mg']+dfMoleMass['O']
dfMoleMass['Al2O3']=2*dfMoleMass['Al']+3*dfMoleMass['O']
dfMoleMass['MnO']=dfMoleMass['Mn']+dfMoleMass['O']
dfMoleMass['P2O5']=2*dfMoleMass['P']+5*dfMoleMass['O']
dfMoleMass['Na2O']=2*dfMoleMass['Na']+dfMoleMass['O']
dfMoleMass['K2O']=2*dfMoleMass['K']+dfMoleMass['O']
dfMoleMass['TiO2']=dfMoleMass['Ti']+2*dfMoleMass['O']
dfMoleMass['FeO']=dfMoleMass['Fe']+dfMoleMass['O']
dfMoleMass['Fe3O4']=3*dfMoleMass['Fe']+4*dfMoleMass['O']
dfMoleMass['Fe2O3']=2*dfMoleMass['Fe']+3*dfMoleMass['O']
dfMoleMass['H2O']=2*dfMoleMass['H']+dfMoleMass['O']

# =============================================================================
# # ===========================================================================
# # defining functions
# # ===========================================================================
# =============================================================================

# =============================================================================
# Finction 1
# =============================================================================
def dfdivide(dfnumerator,dfdenominator):
    """Divides one data frame with several indices over another with 
    only one index; columns must be the same though"""
    dfnumeratorTemp=dfnumerator
    dfdenominatorTemp=dfdenominator
    dfnumeratorTemp['key'] = 1
    dfdenominatorTemp['key'] = 1
    newdf=pd.DataFrame(columns=dfnumeratorTemp.columns) 
    newdf=pd.merge(dfnumeratorTemp,dfdenominatorTemp,on="key").drop("key", 1) 
    for column in dfnumeratorTemp.columns.drop("key"):
        columnx = column + '_x'
        columny = column + '_y'
        if column not in ['Key']:
            newdf[column] = newdf[columnx]/newdf[columny]
            newdf = newdf.drop([columnx,columny], axis = 1)
        # elif column == 'Type':
        #     newdf[column] = newdf[columnx]
        #     newdf = newdf.drop([columnx,columny], axis = 1)
    newdf = newdf.set_index([pd.Index(dfnumeratorTemp.index)])
    return newdf
# =============================================================================
# Function 2
# =============================================================================
dfshamote=pd.read_csv(datadirectory+'/'+'Shamote.csv')
def enthalpy( material , temperature ):
    """This returns the enthalpy based on material and temperature, temperature
    must be in Kelvin and the enthalpy is given in SI unit"""
    dfvar = \
        dfshamote.query(" param == @material ").\
            query(" max_temp > @temperature >= min_temp ")
    temperature=temperature/1000
    coefs = dfvar[['A','B','C','D','E','F','G','H']].values.reshape([8,1])
    return coefs[0]*temperature+coefs[1]*temperature**2/2+\
                      coefs[2]*temperature**3/3+\
                    coefs[3]*temperature**4/4-coefs[4]/temperature+\
                        coefs[5]-coefs[7]
# =============================================================================
# Function 3
# =============================================================================
def dfdistibutioncalculator(dfDistribution,dfElement):
    y={}
    for indexDist in dfDistribution.index:
        list=[]
        for columnElement in dfElement:
            list.append(np.dot(dfDistribution.loc[indexDist].values,\
                                  dfElement[columnElement].values))
        y[indexDist]=list
    # display(y)
    dfNew=pd.DataFrame.from_dict\
        ( y,orient='index', columns=dfElement.columns )
    return dfNew


# =============================================================================
# # ===========================================================================
# # Getting and postprocessing the data 
# # ===========================================================================
# =============================================================================

# =============================================================================
# BAR - adding material to database based on elements
# =============================================================================
dfBARComposition=pd.read_csv(datadirectory+'/'+'BARComposition.csv', index_col=0, header=0)
# dfBARElementsMolarMass = dfMoleMass[dfBARComposition.columns]
# dfBARMole=dfdivide(dfBARComposition,dfBARElementsMolarMass)
# =============================================================================
# Fuel - adding material to database based on elements
# =============================================================================
dfFuelComposition=pd.read_csv(datadirectory+'/'+'FuelComposition.csv', index_col=0, header=0)
dfFuelComposition=dfFuelComposition.sort_index(axis=0)
dfFuelElementsMolarMass = dfMoleMass[dfFuelComposition.columns]
dfFuelElementsMole= dfdivide(dfFuelComposition,dfFuelElementsMolarMass)
# =============================================================================
# Boost fuel
# =============================================================================
dfBoostFuelDistribution=pd.read_csv(datadirectory+'/'+'BoostFuelDistribution.csv'\
                                    , index_col=0, header=0)
dfBoostFuelDistribution=dfBoostFuelDistribution.sort_index(axis=1)
dfBoostFuelMole = dfdistibutioncalculator(dfBoostFuelDistribution,dfFuelElementsMole)
dfBoostFuelComposition = dfdistibutioncalculator(dfBoostFuelDistribution,dfFuelComposition)
# =============================================================================
# Lateral Fuel    
# =============================================================================
dfLateralFuelDistribution=pd.read_csv(datadirectory+'/'+'LateralFuelDistribution.csv'\
                                , index_col=0, header=0)
dfLateralFuelDistribution=dfLateralFuelDistribution.sort_index(axis=1)
dfLateralFuelMole = dfdistibutioncalculator(dfLateralFuelDistribution,dfFuelElementsMole)
dfLateralFuelComposition = dfdistibutioncalculator(dfLateralFuelDistribution,dfFuelComposition)

# =============================================================================
# Ratio between boost and lateral fuel
# =============================================================================

dfFuelBoostLateralFlow=pd.read_csv\
    (datadirectory+'/'+'FuelBoostLateralFlow.csv', index_col=0, header=0)


# =============================================================================
# Scrap - adding material to database based on elements
# =============================================================================
dfScrapComposition=pd.read_csv(datadirectory+'/'+'ScrapComposition.csv', index_col=0, header=0)
dfScrapElementsMolarMass = dfMoleMass[dfScrapComposition.columns]
dfScrapElementsMole=dfdivide(dfScrapComposition,dfScrapElementsMolarMass)
dfScrapDistribution=pd.read_csv(datadirectory+'/'+'ScrapDistribution.csv', \
                                index_col=0, header=0 )
dfScrapMole = dfdistibutioncalculator(dfScrapDistribution,dfScrapElementsMole)
dfScrapTotalComposition = dfdistibutioncalculator(dfScrapDistribution,dfScrapComposition)

# display(dfScrapMole)
# =============================================================================
# Iron - adding material to database based on elements
# =============================================================================
dfIronComposition=pd.read_csv(datadirectory+'/'+'IronComposition.csv', index_col=0, header=0)
dfIronElementsMolarMass = dfMoleMass[dfIronComposition.columns]
dfIronElementsMole = dfIronComposition.\
    div(dfIronElementsMolarMass,axis="columns")

# =============================================================================
# Blast
# =============================================================================

# B2/22.4*1000*0.2
dfBlastFlow=pd.read_csv(datadirectory+'/'+'BlastFlow.csv')
dfBlastFlow['O2']=(dfBlastFlow['V1']+dfBlastFlow['V2'])/22.4*1000*0.21
dfBlastFlow['N2']=(dfBlastFlow['V1']+dfBlastFlow['V2'])/22.4*1000*0.79
# display(dfBlast)


# =============================================================================
# Fines
# =============================================================================

dfFineBlastRatio=pd.read_csv(datadirectory+'/'+'FineBlastRatio.csv', index_col=0, header=0)
dfFineComposition=pd.read_csv(datadirectory+'/'+'FineComposition.csv')

# =============================================================================
# Off Gas
# =============================================================================

dfOffGasComposition=pd.read_csv(datadirectory+'/'+'OffGasComposition.csv', index_col=0, header=0)

# =============================================================================
# =============================================================================
# # Matrices
# =============================================================================
# =============================================================================


# =============================================================================
# A matrix
# =============================================================================
Fe=dfMoleMass['Fe'].iloc[0]
O=dfMoleMass['O'].iloc[0]
C=dfMoleMass['C'].iloc[0]

FeO_BAR=dfBARComposition['FeO'].mean(axis=0)
Fe2O3_BAR=dfBARComposition['Fe2O3'].mean(axis=0)
Fe3O4_BAR=dfBARComposition['Fe3O4'].mean(axis=0)
Fe_Iron=dfIronComposition['Fe'].mean(axis=0)



A11 = -FeO_BAR*(Fe/(Fe+O)) -Fe2O3_BAR*(2*Fe/(2*Fe+3*O)) \
    -Fe3O4_BAR*(3*Fe/(3*Fe+4*O))
A12 = Fe_Iron
A13 = 0
A14 = 0

CO_OffGas=dfOffGasComposition['CO'].mean(axis=0)
CO2_OffGas=dfOffGasComposition['CO2'].mean(axis=0)

A21 = -FeO_BAR*(O/(Fe+O)) -Fe2O3_BAR*(3*O/(2*Fe+3*O)) \
    -Fe3O4_BAR*(4*O/(3*Fe+4*O))
A22 = 0
A23 = 0
A24 = CO_OffGas * (O/(C+O)) + CO2_OffGas * (2*O/(C+2*O))

C_BAR = dfBARComposition['C'].mean(axis=0)
C_Iron = dfIronComposition['C'].mean(axis=0)
C_BAR = dfFuelComposition['C'].mean(axis=0)
BoostRatio = dfFuelBoostLateralFlow['Boost'].mean(axis=0) / \
    (dfFuelBoostLateralFlow['Boost'].mean(axis=0) + \
     dfFuelBoostLateralFlow['Lateral'].mean(axis=0))
LateralRatio = dfFuelBoostLateralFlow['Lateral'].mean(axis=0) / \
    (dfFuelBoostLateralFlow['Boost'].mean(axis=0) + \
     dfFuelBoostLateralFlow['Lateral'].mean(axis=0))
C_BoostFuel = dfBoostFuelComposition['C'].mean(axis=0)
C_LateralFuel = dfLateralFuelComposition['C'].mean(axis=0)

A31 = - C_BAR
A32 = C_Iron
A33 = - (C_BoostFuel * BoostRatio + C_LateralFuel * LateralRatio)
A34 = CO_OffGas * (C/(C+O)) + CO2_OffGas * (C/(C+2*O))

A41 = 0
A42 = 0
A43 = - (C_BoostFuel * BoostRatio + C_LateralFuel * LateralRatio)
A44 = 0
 
A=np.matrix([\
            [A11, A12, A13, A14],\
                [A21, A22, A23, A24],\
                    [A31, A32, A33, A34],\
                        [A41, A42, A43, A44]\
                ])
    
Ainv = np.linalg.inv(A) 


Fe_Fine=dfFineComposition['Fe'].mean(axis=0)
Fine_Blast_Ratio= dfFineBlastRatio['Fine_Blast_Ratio'].mean(axis=0)
BlastFlow = (dfBlastFlow['V1']+dfBlastFlow['V2']).mean(axis=0) * 1.2
Y1 = - Fine_Blast_Ratio * Fe_Fine / 1000 * BlastFlow # gram * 1000 = kg

OxygenAirRatio = 0.2314 # from engineering toolbox
Y2= - BlastFlow * OxygenAirRatio

C_Fine = dfFineComposition['C'].mean(axis=0)

Y3 = - Fine_Blast_Ratio * C_Fine / 1000 * BlastFlow # gram * 1000 = kg

Y4 = - 0.3 * BlastFlow 
# =============================================================================
# Y vector
# =============================================================================
Y= np.matrix([\
            [Y1],\
                [Y2],\
                    [Y3],\
                        [Y4]\
                ])

# =============================================================================
# Calculations
# =============================================================================
X = Ainv * Y








