# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:38:53 2020

@author: EngenhariaCecal
"""
import matplotlib
import numpy
import matplotlib.pyplot as pyplot
import thermopy
from thermopy import nasa9polynomials as nasa9
from thermopy import burcat

db = nasa9.Database()

# print(compound.molecular_weight)
# print(compound.elements)
# print(compound.enthalpy(850+273))
# enthalpy kJ/kmol default

# equation 1 : O2+C -> CO2 + DeltaH

O2 = db.set_compound('o2')

C = db.set_compound('c')

CO2 = db.set_compound('co2')

DeltaHEquation1 = CO2.enthalpy(25+273) - C.enthalpy(25+273) - O2.enthalpy(25+273)

print("Delta H eq 1",float('%.0f' % (DeltaHEquation1)), "J/mol")

# print("Delta H",DeltaH*O2.molecular_weight, "(kJ/kmol) / (kg/kmol) = kJ/kg")

T=25+273
o2 = db.set_compound('o2')
c = db.set_compound('c')
co2 = db.set_compound('co2')
n2 = db.set_compound('n2')

reactants = (o2, c)
products = (co2, n2)
reactants_coefficients = (1, 1)
product_coefficients = (1, 0)
reaction1 = thermopy.nasa9polynomials.Reaction(T, reactants, products, reactants_coefficients, product_coefficients)
print(reaction1.enthalpy_reaction(T))

# equation 1 : O2+C -> CO2 + DeltaH

CO = db.set_compound('carbon monoxide')

DeltaHEquation2 = 2*CO.enthalpy(25+273) - C.enthalpy(25+273) - CO2.enthalpy(25+273)

print("Delta H eq 2",float('%.0f' % (DeltaHEquation2)), "J/mol")

# T=25+273
# c = db.set_compound('c')
# co2 = db.set_compound('co2')
# co = db.set_compound('Fe2O3(cr)')
# n2 = db.set_compound('n2')
# reactants = (co2, c)
# products = (co, n2) 
# reactants_coefficients = (1, 1)
# product_coefficients = (2, 0)
# reaction2 = thermopy.nasa9polynomials.Reaction(T, reactants, products, reactants_coefficients, product_coefficients)
# print(reaction2.enthalpy_reaction(T))
