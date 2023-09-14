#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: hamza"""
import pandas as pd 
import numpy as np



def extraction_of_element(composition):
    kayn =[]
    ch7al=[]
    for i in range(len(composition)):
      #  print(composition[i])
       
            if composition[i].isalpha():
                if composition[i+1].isalpha() and composition[i+1].islower() :
                #    print("dakhal ind+1",composition[i:i+2])
                    kayn.append(composition[i:i+2])           
                elif composition[i].isupper():
                    kayn.append(composition[i])
            elif composition[i].isnumeric():
                    if (i+1)<len(composition) and composition[i+1].isnumeric() :
                        ch7al.append(int(composition[i:i+2]))           
                    elif composition[i-1].isnumeric() :
                        pass
                    else :
                         ch7al.append(int(composition[i]))

    return kayn,ch7al


def  H_soreted_function(speceis,composition ):
    ind=speceis.index('H')
    if ind==0:
        pass 
    else :
        n=composition[0]
        x=speceis[0]
        composition[0]=composition[ind]
        speceis[0]=speceis[ind]
        composition[ind]=n
        speceis[ind]=x
    return  speceis,composition
#  We use the librarie periodictable to obtaine the mass_molar of eche atom in our data 

def Mass_mol(atom_symbol):
    from periodictable import elements
    Symbol_mass={}
    for element in elements:
        # Get the element symbol, name, and molar mass
        symbol = element.symbol
        molar_mass = element.mass
        Symbol_mass[symbol]= molar_mass
    mass_molaire=[]
    Molar_mass_for_Hydrogen= Symbol_mass["H"]
    for atom in range(len(atom_symbol)):
        if atom_symbol[atom] in Symbol_mass.keys():
            Molar_mass=Symbol_mass[atom_symbol[atom]]
            mass_molaire.append(Molar_mass)
        else :  
            mass_molaire.append(0)
    # molar mass and mass_molaire is the same, the deference help in memory storage 
    return mass_molaire,Molar_mass_for_Hydrogen


#The function's purpose is to extract the chemical elements and their corresponding counts from a compound formula
def Extract_Element_count(name):
	name=name
	element_count={}
	Hyd_count={}
	container=[]
	op=0
	if '(' in name or ')' in name:
		indix_left_paranthesis=name.index("(")
		indix_right_paranthesis=name.index(")")
		subname=name[indix_left_paranthesis+1:indix_right_paranthesis]
		elemenets_left_ofthe_left_paranthesis=name[:indix_left_paranthesis]
		coeficient_right_ofthe_right_paranthesis=int(name[indix_right_paranthesis+1:])
		k=2
		strange_name=[elemenets_left_ofthe_left_paranthesis,subname]
	else:
		k=1
	while(op<k):
		if k==2:
			name=strange_name[op]
			#element_count={}
			
		else :
			name=name
		for i in range(len(name)):
		    ch=name[i]
		    if ((ch>='A') & (ch<='Z')):    #  X
		        if i+1 < len(name):
		            chp1=name[i+1]
		            if ((chp1>='A') & (chp1<='Z')):   #XY
		                if ch=="H":
		                    Hyd_count["H"]=1
		                else:
		                    element_count[ch]=1
		            elif ((chp1>='a') & (chp1<='z')):   #  Xy
		                if (i+2)<len(name):
		                    chp2=name[i+2]
		                    if ((chp2>='A') & (chp2<='Z')):
		                        element_count[name[i:i+2]]=1
		                    else :
		                        if (i+3) < len(name):

		                            chp3=name[i+3]
		                            if ((chp3>='0')&(chp3 <= '9')):
		                                if (i+4)< len(name):
		                                    element_count[name[i:i+2]]=int(name[i+2:i+4])
		                                else: 
		                                    element_count[name[i:i+2]]=int(name[i+2:])
		                            else :
		                                element_count[name[i:i+2]]=int(name[i+2])
		                        else :
		                                element_count[name[i:i+2]]=int(name[i+2])
		                else: 
		                    element_count[name[i:]]=1
		            else :       #  Xa
		                if i+2 < len(name):
		                    chp2=name[i+2]
		                    if ((chp2>='0')&(chp2<='9')): #Xab                    
		                        if i+3 < len(name):
		                            if ch=="H":
		                                Hyd_count["H"]=int(name[i+1:i+3])
		                            else:
		                                element_count[ch]=int(name[i+1:i+3])
		                                
		                        else : 
		                            if ch=="H":
		                                Hyd_count["H"]=int(name[i+1:])
		                            else:
		                                element_count[ch]=int(name[i+1:])

		                    else :
		                        if ch=="H":
		                            Hyd_count["H"]=int(name[i+1])
		                        else:
		                            element_count[ch]=int(name[i+1])
		                else :
		                    if ch=="H":
		                        Hyd_count["H"]=int(name[i+1])
		                    else:
		                        element_count[ch]=int(name[i+1])

		        else :
		            if ch=="H":
		                Hyd_count["H"]=1
		            else:
		                element_count[ch]=1   
		    else : 
		        pass
		op+=1
		#container.append(element_count)


	if k==2:
		element_count[list(element_count.keys())[1]]=list(element_count.values())[1]*coeficient_right_ofthe_right_paranthesis
		Hyd_count["H"]=list(Hyd_count.values())[0]*coeficient_right_ofthe_right_paranthesis
	return element_count,Hyd_count



def Apply_Extract_Element_count(Series_name,Series_ntypes):
    list_of_Dict=[]
    for i in range (len(Series_name)):
        d1,d_H=Extract_Element_count(Series_name[i])
        list_keys=list(d1.keys())
        list_values=list(d1.values())
        list_values_H=list(d_H.values())
        ntype=Series_ntypes[i]
        Dict2={}
        if ntype>2:
            Dict2["A"]=list_keys[0]
            Dict2["nA"]=list_values[0]
            Dict2["B"]=list_keys[1]
            Dict2["nB"]=list_values[1]
            Dict2["H"]=list_values_H[0]

        else:       
            Dict2["A"]=list_keys[0]
            Dict2["nA"]=list_values[0]  
            Dict2["B"]="NaN"
            Dict2["nB"]=0
            Dict2["H"]=list_values_H[0]  


        list_of_Dict.append(Dict2)
    return list_of_Dict
