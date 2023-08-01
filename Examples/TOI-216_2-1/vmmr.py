
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pathlib
import subprocess
import sys
import os
from statistics import mean

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction
import vplot as vpl
from scipy.interpolate import interp2d
from tqdm import tqdm

import vplanet


"""
d = float (double in C)
c = string (character in C)
f = fraction
i = integer (or any single lower-case letter)
ax = an array of values of type(x)
l = list
s = dictionary
v = function
q = object

"""

# Path hacks
path = pathlib.Path(__file__).parents[0].absolute()
sys.path.insert(1, str(path.parents[0]))

iNumLibrations = 0
dInnerEccentricity = "nan"

####FUNCTIONS START HERE####

# A function that prints out a Default value if the user does not customize a requested value.
def vStateDefault(dTmp):
    return print("The default value is " + str(dTmp))

# A function that requests the user to customize a certain option. The result becomes a float.
def vInputRequest(dDefault):
    print("Press enter if you want to skip this customization step.")
    dValue = 0.0
    i = 0
    while i == 0:
        dValue = input()
        if dValue == "":
            dValue = dDefault #Default value
            print("VALUE SET TO THE ITS DEFAULT:", str(dValue))
            print("")
            i += 1
        else:
            try:
                float(dValue)
                dValue = float(dValue)
                i += 1
                print("VALUE SET TO:", str(dValue))
            except ValueError:
                print("")
                print("Undetected value. Specify value or skip by pressing enter.")
    return dValue

# Function below produces an integer output of a certain value's user option while modifying the domain of the input option
# iMinInt sets lowest int input. Default is 0
# iMaxInt sets highest int input. Default is infinity
# iValue is interpreted as the default value
def vInputRequestInt(iValue, iMinInt = "-inf", iMaxInt = "+inf"):  
    iOutput = int(vInputRequest(iValue))
    i = 0
    if isinstance(iMinInt, int) == True and isinstance(iMaxInt, int) == True and iMaxInt >= iMinInt:
        while i == 0:
            if iOutput < iMinInt or iOutput > iMaxInt:
                print("Must input an integer between " + 
                      str(iMinInt) + " and " + 
                      str(iMaxInt) + ". ")
                iOutput = int(vInputRequest(iValue))
            else:
                i += 1
    elif iMinInt == "-inf" and isinstance(iMaxInt, int) == True:
        while i == 0:
            if iOutput > iMaxInt:
                print("Must input an integer that is equal to " +  
                      "or lower than " + str(iMaxInt) + ". ")
                iOutput = int(vInputRequest(iValue))
            else:
                i += 1
    elif isinstance(iMinInt, int) == True and iMaxInt == "+inf":
        while i == 0:
            if iOutput < iMinInt:
                print("Must input an integer that is equal to " + 
                      "or higher than " + str(iMinInt) + ". ")
                iOutput = int(vInputRequest(iValue))
            else:
                i += 1
    else:
        print("Revise your integer inputs.")
    print("")
    return iOutput

#Function below produces a dictionary of each satellite's initial Orbital Period
def vFindInitOrbPeriods(aTmpBodies):
    #Array that assigns names
    acBodiesNames = [str(aTmpBodies[iBody]).split(" ")[-1][:-1] 
                        for iBody in range(len(aTmpBodies)) 
                        if hasattr(aTmpBodies[iBody], 'OrbPeriod')]
    #Array that finds the initial orbtial periods
    adTmpInitOrbPeriods = [aTmpBodies[iBody].OrbPeriod[0] 
                            for iBody in range(len(aTmpBodies))
                            if hasattr(aTmpBodies[iBody], 'OrbPeriod')]
    #Creating a dictionary for satellite name with their initial value
    sTmpInitOrbPeriods = {}
    for iBody in range(len(acBodiesNames)):
        sTmpInitOrbPeriods[acBodiesNames[iBody]] = adTmpInitOrbPeriods[iBody]
    
    #Sorting Initial Orbital Periods from highest to lowest
    adTmpInitOrbPeriods.sort(reverse = True)
    
    lTmpInitOrbPeriods = sorted(sTmpInitOrbPeriods.items(), 
                                key=lambda x:x[1], reverse = True)
    sTmpInitOrbPeriods = dict(lTmpInitOrbPeriods)
    return sTmpInitOrbPeriods

#Function below produces all potential integer:integer Orbital Period Ratios (OPRs) with their respected Order of Resonance (OR)
#The result is an array of dictionaries containing the satellite pair, OPR numerator and denominator, and the OR
def vFindOPRsAndORs(adInitialOrbitalPeriods, dMaxOrbPTolerance, dInnerEccentricity = "nan"):
    asResInfo = [] #States which satellites may be in resonance with additional info such as Orbital Period Ratio and the Order of Resonance
    #We iterate over the orbital period dictionary twice in two nested for loops to create Orbital Period Ratios (OPRs)
    for dKey1 in adInitialOrbitalPeriods:
        dTmpOrbP1 = adInitialOrbitalPeriods[dKey1]
        for dKey2 in adInitialOrbitalPeriods:
            dTmpOrbP2 = adInitialOrbitalPeriods[dKey2]                  
            # We want distinct ratios by making the first iterated value being more than the second 
            if dTmpOrbP1 > dTmpOrbP2:
                adOrbResStrengths = [] # Collecting Orbital Strengths
                for iDenom in range(1, iMaxDenominator + 1):
                    dTmpOrbPRatio = dTmpOrbP1 / dTmpOrbP2
                    fTmpOrbPFrac = Fraction(dTmpOrbPRatio).limit_denominator(iDenom)
                    iTmpNumerator = fTmpOrbPFrac.numerator
                    iTmpDenominator = fTmpOrbPFrac.denominator
                    iTmpk = iTmpNumerator - iTmpDenominator # Order of Resonance
                    if iDenom == 1:
                        iNumerator = iTmpNumerator
                        iDenominator = iTmpDenominator
                        k = iTmpk
                    if iTmpk <= kMax:
                        #dInnerEccentricity = "nan"# 0.027
                        if isinstance(dInnerEccentricity, float) == True:
                            dStrength = dInnerEccentricity ** abs(iTmpk)
                            adOrbResStrengths.append(dStrength)
                            # Need to make distinct k values
                            if iDenom == 1:
                                #kTmpPrevious = 0 #Treated as a dummy var
                                i = 0 # Keeps track if anything EVER changed
                                #dOrbPRatio = dTmpOrbPRatio
                                iNumerator = iTmpNumerator
                                iDenominator = iTmpDenominator
                                k = iTmpk
                                #fOrbPFrac = fTmpOrbPFrac
                            #Need to iterate at least twice before creating strength comparisons
                            elif iDenom > 1 and kTmpPrevious != iTmpk: 
                                dSCompare = adOrbResStrengths[-2] / adOrbResStrengths[-1]
                                i += 1
                                if dSCompare < dStrengthCompare:                                    
                                    print("dSCompare:", dSCompare)
                                    #dOrbPRatio = dTmpOrbPRatio
                                    iNumerator = iTmpNumerator
                                    iDenominator = iTmpDenominator
                                    k = iTmpk
                                    #fOrbPFrac = fTmpOrbPFrac
                            elif iDenom == iMaxDenominator and i == 0:
                                #dOrbPRatio = dTmpOrbPRatio
                                iNumerator = iTmpNumerator
                                iDenominator = iTmpDenominator
                                k = iTmpk
                                #fOrbPFrac = fTmpOrbPFrac
                                
                        elif dInnerEccentricity == "nan":
                            if iDenom == iMaxDenominator:
                                print("No eccentricity detected.")
                            #dOrbPRatio = dTmpOrbPRatio
                            iNumerator = iTmpNumerator
                            iDenominator = iTmpDenominator
                            k = iTmpk
                            #fOrbPFrac = fTmpOrbPFrac
                            
                    kTmpPrevious = iTmpk
                    
                    dExactOrbRes = iNumerator / iDenominator
                    print("iNumerator:", iNumerator)
                    #dTmpOrbPTolerance = abs(dTmpOrbPRatio - dExactOrbRes)
                    dTmpOrbPTolerance = abs(dTmpOrbPRatio - dExactOrbRes) / dExactOrbRes
                if dTmpOrbPTolerance < dMaxOrbPTolerance and k <= kMax:
                    print("Orbital Resonance candidate found!")
                    print("Satellite " + str(dKey1) + "'s Orbital Period:", dTmpOrbP1)
                    print("Satellite " + str(dKey2) + "'s Orbital Period:", dTmpOrbP2)
                    print("Orbital Period Ratio in decimal form:", dTmpOrbPRatio)
                    print("Expected Orbital Resonance:", 
                            str(iNumerator) + ":" + str(iDenominator))
                    print("")
                    print("Order of Resonance:", k)
                    print("Orbital Period Tolerance:", dTmpOrbPTolerance, 
                            "(" + str(round(dTmpOrbPTolerance*100, 2)) + "%)")
                    print("")
                    #Orbital Period of SatelliteOuter is always longer than SatelliteInner
                    sTmpResInfo = {"SatelliteOuter": dKey1,
                                    "SatelliteInner": dKey2,
                                    "OrbPerOuter": dTmpOrbP1,
                                    "OrbPerInner": dTmpOrbP2,
                                    "Numerator": iNumerator,
                                    "Denominator": iDenominator, 
                                    "ResOrder": k}
                    asResInfo.append(sTmpResInfo)
    return asResInfo

def vResArgCoefficients(iResOrder, j, cResType = "ecc"): # j is the numerator of the Orbital Period Ratio
    # cResType = <"ecc" | "inc" > # Chooses either eccentricity or inclination resonances.
    aaResArgCoeffs = []
    for i1 in range(iResOrder + 1):
        for i2 in range(iResOrder + 1):
            if i1 + i2 == iResOrder:
                if cResType == "ecc":
                    aiResArgCoeffs = [j, iResOrder - j]
                    aiResArgCoeffs.append(-i1)
                    aiResArgCoeffs.append(-i2)
                    aiResArgCoeffs.append(0)
                    aiResArgCoeffs.append(0)
                    aaResArgCoeffs.append(aiResArgCoeffs)
                if cResType == "inc":
                    aiResArgCoeffs = [j, iResOrder - j]
                    aiResArgCoeffs.append(0)
                    aiResArgCoeffs.append(0)                    
                    aiResArgCoeffs.append(-i1)
                    aiResArgCoeffs.append(-i2)
                    aaResArgCoeffs.append(aiResArgCoeffs)
                break
    return aaResArgCoeffs

# Finds the Orbital Elements of the Satellite which matches the Initial Orbital Period input
def vFindOrbElemsSatellite(aTmpBodies, dTmpOrbPerInit):
    saTmpAngles = {"MeanLongitude":'', "LongPericenter":'', "LongAscNode":''}
    for iTmpBody in range(len(aTmpBodies)):
        if hasattr(aTmpBodies[iTmpBody], 'OrbPeriod') and aTmpBodies[iTmpBody].OrbPeriod[0] == dTmpOrbPerInit:
            saTmpAngles["MeanLongitude"] = aTmpBodies[iTmpBody].MeanLongitude
            saTmpAngles["LongPericenter"] = aTmpBodies[iTmpBody].LongP
            saTmpAngles["LongAscNode"] = aTmpBodies[iTmpBody].LongA
            break
    return saTmpAngles

def vResonantArguments(aiCoefficients,
                        adTmpMeanLongitudeOuter, adTmpLongPericenterOuter, adTmpLongAscNodeOuter, 
                        adTmpMeanLongitudeInner, adTmpLongPericenterInner, adTmpLongAscNodeInner):
    #We must multiply the coefficients by their respected parameter
    aaAngles = [adTmpMeanLongitudeOuter, adTmpMeanLongitudeInner, 
                adTmpLongPericenterOuter, adTmpLongPericenterInner, 
                adTmpLongAscNodeOuter, adTmpLongAscNodeInner]
                
    # The next array stores all Resonant Arguments for a certain period of time within the simulation
    adResonantArguments = []
    """
    # The for loop follows the equation of the Resonant Argument
    dResonantArgument = (j1*MeanLongitudeOuter + j2*MeanLongitudeInner + 
                        j3*LongPericenterOuter + j4LongPericenterInner +
                        j5*LongAscNodeOuter + j6*LongAscNodeInner)
    """
    for iAngleElem in range(len(aaAngles[0])):
        dResonantArgument = sum(aiCoefficients[i]*aaAngles[i][iAngleElem] for i in range(len(aiCoefficients)))
        # Appends each index of a Resonant Argument into a final array
        # We also take the modulo to restrict angles between 0 and 360 degrees
        adResonantArguments.append(dResonantArgument % 360)

    # The final desired result of the function
    return adResonantArguments

####FUNCTIONS END HERE####

dDeftOrbP = 0.1 # Default Orbital Period Ratio Tolerance
dDeftStrengthCompare = 1000.0 # Default Strength comparison
# Highest possible order of resonance is second order k = 2
iDeftkMax = 3 # Default order of Resonance
iDeftMaxDenominator = 10 # Default maximum integer in Orb resonances denominator

print("Would you like to customize the tolerance of the resonance, " + 
      "the orbital resonance strength comparison, " + 
      "the highest order of resonance, " + 
      "or the highest denominator of an orbital resonance? " + 
      "More details in the README.txt.")
print("Options: 'y' or 'yes' + press enter. Only press enter to skip the entire customization.")
cCustomize = str(input())
if cCustomize == "yes" or cCustomize == "y":
    print("Specify the 'Tolerance of the Resonance'. \n" 
          "Example: a tolerance of 0.05 (5%) means that " +
          "the Orbital Period Ratio tested for a 2:1 " +
          "resonance is restricted between 1.90 and 2.10. " +
          "This example would convert an orbital period " +
          "ratio of 2.09 or 1.91 to a resonance of 2:1 " +
          "instead of a 21:10 or 19:10. ")
    vStateDefault(dDeftOrbP) 
    dOrbPTolerance = vInputRequest(dDeftOrbP)
    
    print("Specify the 'Orbital Resonance Strength Comparison'. \n" + 
          "Strengh is the inner satellite's eccentricity to the " +
          "power of the absolute value of the resonance order. " + 
          "Example: an inner satellite with ecc = 0.027 and an " +
          "resonance of 21:10 would have a strength of 10^-18 " +
          "while the 2:1 resonance has a strength of 10^-2 " + 
          "thus the 2:1 resonance would be chosen since their " +
          "strength comparison would be 10^15 > 10^3. ")
    vStateDefault(dDeftStrengthCompare)
    dStrengthCompare = vInputRequest(dDeftStrengthCompare)
    
    print("Specify the 'Highest Resonance Order'. \n" +  
          "Your current options are between 1 and 3" + 
          ". Decimal values will round to the nearest integer. ")
    vStateDefault(iDeftkMax)
    kMax = vInputRequestInt(iDeftkMax, 1)
    
    print("Specify the 'highest positive integer in the Orbital Resonances denominator'. \n" + 
          "Example: inputting a value 9 would convert an orbital period ratio of " + 
          "2.10 to a value of 2:1 rather than 21:10 since 9 < 10. \n")
    vStateDefault(iDeftMaxDenominator)
    iMaxDenominator = vInputRequestInt(iDeftMaxDenominator, 1)
else:    #We iterate over the orbital period dictionary twice to create Orbital Period Ratios (OPRs)
    dOrbPTolerance = dDeftOrbP
    dStrengthCompare = dDeftStrengthCompare
    kMax = iDeftkMax
    iMaxDenominator = iDeftMaxDenominator
    

    
    
    
"""Running Single Simulation"""

if not ((path / "vspace.in").exists()):
    # Defining Object obtained from VPLanet
    qOutput = vplanet.run(path / "vpl.in", quiet=True, units=False)
    aBodies = qOutput.bodies

    #Function that creates dictionary of satellites with their Initial Orbital Periods
    adInitOrbPeriods = vFindInitOrbPeriods(aBodies) 
    print("Satellites with their respected Orbital Periods listed below.") 
    print(adInitOrbPeriods)
    print("")
    
    print("Searching for Orbital Resonance candidates...")
    #Now we find the Orbital period ratios (OPR) of all satellites
    #We will also produce the potential Order of Resonances (ORs) per satellite pair in a potential Resonance
    #Note there may be more than one potential Resonance in a system, so we produce all possible Resonances.
    asPotentialResPairs = vFindOPRsAndORs(adInitOrbPeriods, dOrbPTolerance)
    print(asPotentialResPairs)
    print("")
    for sPair in asPotentialResPairs:
        dOrbPerInitOuter = sPair["OrbPerOuter"]
        dOrbPerInitInner = sPair["OrbPerInner"]
        cSatelliteOuter = sPair["SatelliteOuter"]
        cSatelliteInner = sPair["SatelliteInner"]
        print("Coefficients of the Resonant Arguments for the satellite pair", 
                cSatelliteOuter, "and", cSatelliteInner)
        print(vResArgCoefficients(sPair["ResOrder"], sPair["Numerator"]))
        #Finds the coefficients of every Resonant Argument
        aaAllResArgCoeff = vResArgCoefficients(sPair["ResOrder"], sPair["Numerator"])

        #The bottom two lines finds the related Orbital Elements to the outer and inner satellite, respectively
        saOrbElemsOuter = vFindOrbElemsSatellite(aBodies, dOrbPerInitOuter) 
        saOrbElemsInner = vFindOrbElemsSatellite(aBodies, dOrbPerInitInner)

        aaResonantArguments = []
        
        for aiResArgCoeff in aaAllResArgCoeff:
            aaResonantArguments.append(vResonantArguments(aiResArgCoeff, 
                saOrbElemsOuter["MeanLongitude"], saOrbElemsOuter["LongPericenter"], saOrbElemsOuter["LongAscNode"], 
                saOrbElemsInner["MeanLongitude"], saOrbElemsInner["LongPericenter"], saOrbElemsInner["LongAscNode"]))

        def vMakePlot(aTmpBodies, aaTmpAllResArgCoeff, aaTmpResonantArguments, cTmpSatelliteOuter, cTmpSatelliteInner):
            # Number of Resonant Arguments >= 2, unless if there are no potential ResArgs detected.
            # Counting the number of Resonant Arguments to use for our plot
            iNumResArgs = int(len(aaTmpResonantArguments))
            print("Number of Resonant Arguments for the satellite pair " + 
                cTmpSatelliteOuter + "/" + cTmpSatelliteInner + ":", iNumResArgs)
            cPlotName = "ResArgPair_" + cTmpSatelliteOuter + "_" + cTmpSatelliteInner
            bMakePlot = False
            for iTmpBody in range(len(aTmpBodies)):
                if hasattr(aTmpBodies[iTmpBody], 'Time') and iNumResArgs > 0:
                    bMakePlot = True
                    aTime = aTmpBodies[iTmpBody].Time
                    break
            # Determining the appropriate units to plot of the time axis.
            def vReadVplInFile(file_path):
                cUnitTime = "year"     # Default value if not found in the file
                dStopTime = 1.0        # Default value if not found in the file

                with open(file_path, 'r') as file:
                    for line in file:
                        if "sUnitTime" in line:
                            cUnitTime = str(line.split()[1].strip())
                        if "dStopTime" in line:
                            dStopTime = float(line.split()[1])

                return cUnitTime, dStopTime
            # Here are the desired units to calculate within SpiNBody
            sdTimeConversions = {
                'second': 1,
                'minute': 60,
                'hour': 3600,
                'day': 86400,
                'year': 31536000
            }
            # Here is a function that converts the time array to the appropriate units
            def vConvertTime(adTimeValues, cCurrentUnit, cTargetUnit): # 1st arg can be float too
                dConversionFactor = sdTimeConversions[cCurrentUnit] / sdTimeConversions[cTargetUnit]
                return adTimeValues * dConversionFactor
            # Read vpl.in and get the time unit and total length of time
            cVpl_in_File = "vpl.in"  # You may need to specify the correct path here
            cTimeUnit, dTotalTime = vReadVplInFile(cVpl_in_File)
            # Default target unit is set to the default time unit
            cTargetUnit = cTimeUnit
            # Determine the appropriate time unit based on the total length of time
            if cTimeUnit != "year":
                for cKey in sdTimeConversions:
                    dConversionUnit = sdTimeConversions[cKey]
                    if dConversionUnit >= sdTimeConversions[cTimeUnit]:
                        dMaxTime = vConvertTime(dTotalTime, cTimeUnit, cKey)
                        cTargetUnit = cKey
                        if dMaxTime <= 1000:
                            break
                # Converting the time array to the target unit
                if cTimeUnit != cTargetUnit:
                    aTime = vConvertTime(aTime, cTimeUnit, cTargetUnit)
            # continuing on towards potentially higher units.
            # Determining the prefix units of Time in order to label properly the Horizontal axis
            #k = kilo [1e3], M = Mega [1e6], G = Giga [1e9], T = Tera [1e12] (I doubt you will ever reach a Terayear)
            acUnitPrefixes = ["", "k", "M", "G", "T"]    
            iTimeUnitDetermine = 0
            while max(aTime) >= 1000.0:
                aTime = [dTime / 1000.0 for dTime in aTime]
                iTimeUnitDetermine += 1
                if iTimeUnitDetermine == len(acUnitPrefixes):
                    break
            # We have obtained the unit prefix
            cTimePrefix = acUnitPrefixes[iTimeUnitDetermine]
            def vSubscript(cName, iSubscriptLength = 3):
                cSubscript = cName
                if len(cName) > iSubscriptLength:
                    cSubscript = cName[:iSubscriptLength:]
                return cSubscript
            cSatOuterSubscript = vSubscript(cTmpSatelliteOuter)
            cSatInnerSubscript = vSubscript(cTmpSatelliteInner)
            # Determining the Resonant Argument Labels on the vertical axis
            # Preparing to set up the labels here
            def vCreateLabel(cSubOuter, cSubInner, aiTmpCoefficients):
                acSymbols = ["lambda_{{{}}}".format(cSubOuter), "lambda_{{{}}}".format(cSubInner), 
                                "varpi_{{{}}}".format(cSubOuter), "varpi_{{{}}}".format(cSubInner), 
                                "Omega_{{{}}}".format(cSubOuter), "Omega_{{{}}}".format(cSubInner)]
                
                # Defining an empty label to input the string
                cTmpLabel = ""
                for i in range(len(aiTmpCoefficients)):
                    # Only add symbol when the coefficient is nonzero
                    if aiTmpCoefficients[i] != 0:
                        if aiTmpCoefficients[i] == -1:
                            # Omitting the 1 and only including the '-' sign
                            cTmpLabel += "-" + "\\" + acSymbols[i]
                        else:
                            # Including the coefficient next to the associated symbol
                            cTmpLabel += str(aiTmpCoefficients[i]) + "\\" + acSymbols[i]
                            
                # Appending the Resonant Argument label while removing excessive '$' signs for clean Latex notation 
                cTmpLabel = "$" + cTmpLabel + "$"
                return cTmpLabel           


            # Creating size of the plot. For now we leave it with numbers, but I may generalize these numbers later for n ResArgs
            if bMakePlot == True:
                mpl.rcParams['figure.figsize'] = (6.5,6.5)
                mpl.rcParams['font.size'] = 11.0
                cColor = "black"
                iPointSize = 10.0
                cTimeLabel = "Time [" + cTimePrefix + cTargetUnit + "s]"
                if iNumResArgs > 0 and iNumResArgs <= 3:
                    fig, axes = plt.subplots(ncols = 1, nrows = iNumResArgs, sharey=False)
                    for iSim in range(iNumResArgs):
                        cLabel = vCreateLabel(cSatOuterSubscript, cSatInnerSubscript, aaTmpAllResArgCoeff[iSim])
                        print("Resonant Argument", str(iSim) + ":", cLabel)
                        axes[iSim].scatter(aTime, aaTmpResonantArguments[iSim], color = cColor, s = iPointSize)
                        axes[iSim].set_xlabel(cTimeLabel, fontsize=14)
                        axes[iSim].set_ylabel(cLabel + " [$^\circ$]", fontsize=16)
                        axes[iSim].set_xlim(min(aTime), max(aTime))
                        axes[iSim].set_ylim(0.0, 360.0)
                        
                if iNumResArgs >= 4:
                    iNumRows = int(iNumResArgs / 2)
                    fig, axes = plt.subplots(ncols = 2, nrows = iNumRows, sharey=False)
                    for iSim in range(iNumResArgs):
                        iRow = 0; iCol = iSim % iNumRows
                        if iSim >= iNumRows:
                            iRow += 1
                        cLabel = vCreateLabel(cSatOuterSubscript, cSatInnerSubscript, aaTmpAllResArgCoeff[iSim])
                        print("Resonant Argument", str(iSim) + ":", cLabel)
                        axes[iRow][iCol].scatter(aTime, aaTmpResonantArguments[iSim], color = cColor, s = iPointSize)
                        axes[iRow][iCol].set_xlabel(cTimeLabel, fontsize=12)
                        axes[iRow][iCol].set_ylabel(cLabel + " [$^\circ$]", fontsize=14)
                        axes[iRow][iCol].set_xlim(min(aTime), max(aTime))
                        axes[iRow][iCol].set_ylim(0.0, 360.0)
                fig.tight_layout()
                if (sys.argv[1] == 'pdf'):
                    plt.savefig(path / (cPlotName + ".pdf"), bbox_inches="tight", dpi=400)
                elif (sys.argv[1] == 'png'):
                    plt.savefig(path / (cPlotName + ".png"), bbox_inches="tight", dpi=400)
                else:
                    print("Your first argument needs to be either png or pdf")
                plt.show()
            return "Plotting Complete.\n"
        
        print(vMakePlot(aBodies, aaAllResArgCoeff, aaResonantArguments, 
                        cSatelliteOuter, cSatelliteInner))
    print("BAAAAAAM!!!\n")
                
"""Work on this after finishing single file case. 
    This section of code only for multiple simulations.
    Which means we have vspace.in in directory."""
print("checking to see if a vspace.in file exists...")
if ((path / "vspace.in").exists()):
    print("A vspace.in file exists.")
    #Locating necessary names of Destination Folder and Each file's Trial Name
    acVspaceLines = ["sDestFolder", "sTrialName"]
    acSourceNames = [] #Identifies names of sDestFolder and Trial Name
    with open(path / "vspace.in", "r") as file:
        iFoundNames = 0
        
        for line_number, line in enumerate(file, start=1):
            for cVspaceName in acVspaceLines:  
                if cVspaceName in line:
                    acSourceNames.append(line[len(cVspaceName) + 1:])
                    iFoundNames += 1
                if iFoundNames == len(acVspaceLines):
                    break
    cDestFolderName = acSourceNames[0].strip()
    cTrialName = acSourceNames[1].strip()
    cVspaceDir = path / cDestFolderName
    cTrailNamesDir = cTrialName + "*"
    print("Checking to see if a vspace directory exists...")
    if not cVspaceDir.exists():
        print("A directory does not exist. Using Vspace to create simulation directories and files...")
        subprocess.check_output(["vspace", "vspace.in"], cwd=str(path))

        lDirs = glob.glob(str(path / cDestFolderName / cTrailNamesDir))
        print("Created a directory called", cDestFolderName, "within the directory", path)
        print("Created", len(lDirs), "subdirectories inside", cVspaceDir)
    elif cVspaceDir.exists():
        lDirs = glob.glob(str(path / cDestFolderName / cTrailNamesDir))
        print("The vspace directory", cDestFolderName, "exists and contains", 
                len(lDirs), "subdirectories inside", cVspaceDir)
else:
    print("No vspace.in file found.")            


