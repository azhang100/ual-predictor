import visualizer

# references:
# DissociationCurve = https://derangedphysiology.com/main/cicm-primary-exam/
#                     respiratory-system/Chapter-112/oxyhaemoglobin-dissociation-curve
# AdvectionDiffusion = https://pmc.ncbi.nlm.nih.gov/articles/PMC2584207/pdf/nihms77385.pdf
# MorePromises = 

# ======== user input ========

MICROMOL_O2_PER_ML_BOUND = 7.5 # (μM of O2 / mL of blood)
MICROMOL_O2_PER_ML_FREE = 0.0139 # (μM of O2 / mL of blood / mmHg)

MICROMOL_PER_ML_TO_SO2 = 1.0 / MICROMOL_O2_PER_ML_BOUND * 100

DIFFUSIVITY_O2 = 0.0000162 # (cm^2 / sec)
# diffusionMicroMol = DIFFUSIVITY_O2 * MICROMOL_O2_PER_ML * ()

# PDMS
MEMBRANE_PERMEABILITY = 0.000000236 # (2.676E-7 μM of O2 * cm /cm2 /sec /mmHg)

# ==================== more ===================

# bloodFlow is mL/sec
# membraneThickness, capillaryHeight is cm
# gasExchangeSurfaceArea is cm^2
# inputSO2 is percent (0.00 to 100.0)
def calculate(capillaryHeight, membraneThickness, bloodFlow, inputSO2, gasExchangeSurfaceArea,
              gasPO2 = 760,
              timeTickSteps=64, bloodArrayLength=64, useTwoSidedDiffusion=True):
    
  # calculate a few things here so it is not repeatedly calculated
  dwellTime = gasExchangeSurfaceArea * capillaryHeight / bloodFlow
  timeTick = dwellTime/timeTickSteps
  
  bloodDiffusionDistance = capillaryHeight/bloodArrayLength
  bloodAdmittance = timeTick * DIFFUSIVITY_O2 * MICROMOL_O2_PER_ML_FREE * bloodDiffusionDistance**(-2) * MICROMOL_PER_ML_TO_SO2
  print("bloodAdmittance:", bloodAdmittance)
  membraneAdmittance = timeTick * MEMBRANE_PERMEABILITY / membraneThickness * MICROMOL_PER_ML_TO_SO2 / bloodDiffusionDistance
  print("timeTick:", timeTick)
  print("membraneThickness:", membraneThickness)
  print("membraneAdmittance:", membraneAdmittance)
  
  print("dwellTime:",dwellTime)
  
  bloodArray = [inputSO2]*bloodArrayLength
  bloodArrayNew = [inputSO2]*bloodArrayLength
  for t in range(0, timeTickSteps):
    print("time step:", t, "current SO2:", sum(bloodArray)/len(bloodArray))
    bloodArray[0                 ] = membraneDiffuse(bloodArray[0                 ],gasPO2,membraneAdmittance)
    bloodArray[bloodArrayLength-1] = membraneDiffuse(bloodArray[bloodArrayLength-1],gasPO2,membraneAdmittance*useTwoSidedDiffusion)
    # visualizer.draw(bloodArray)

    copyList(bloodArrayNew, bloodArray)
    for i in range(0,len(bloodArray)-1):
      diffusionSat = bloodDiffuse(bloodArray[i], bloodArray[i+1], bloodAdmittance)
      bloodArrayNew[i  ] += diffusionSat
      bloodArrayNew[i+1] -= diffusionSat
    copyList(bloodArray, bloodArrayNew)
      
  visualizer.draw(bloodArray)
      
  return sum(bloodArray)/len(bloodArray)

# diffuse from cellA to cellB
def bloodDiffuse(cellB, cellA, admittance):
  pressureB = saturationToPressure(cellB)
  pressureA = saturationToPressure(cellA)
  diffusionSat = admittance*(pressureA-pressureB)
  # diffusionMicroMol = timeTick * DIFFUSIVITY_O2 * MICROMOL_O2_PER_ML_FREE * (pressureA-pressureB) *
  #                     (bloodDiffusionDistance)^(-2)
  # diffusionSat = diffusionMicroMol / MICROMOL_O2_PER_ML_BOUND * 100
  return diffusionSat

def membraneDiffuse(saturationBlood, pressureGas, admittance):
  pressureB = saturationToPressure(saturationBlood)
  diffusionSat = admittance*(pressureGas-pressureB)
  # diffusionMicroMol = timeTick * (pressureGas-pressureB) / 
  #                     (MEMBRANE_RESISTANCE_PER_DISTANCE * membraneThickness)
  # diffusionSat = diffusionMicroMol / MICROMOL_O2_PER_ML_BOUND * 100
  #                      / bloodDiffusionDistance
  saturationBlood += diffusionSat
  saturationBlood = max(min(saturationBlood,100),0)
  return saturationBlood

# ================== MORE HELPER FUNCTIONS =========================

from bisect import bisect_left, bisect_right

# source: DissociationCurve
# TODO: Adjust dissocationCurveSaO2 to "effective Saturation" which includes freely dissolved O2
"""
dissociationCurvePaO2 = [    0,    1,    2,    4,    6,    8,   10,   12,
                            14,   16,   18,   20,   22,   24,   26,   28,   
                            30,   32,   34,   36,   38,   40,   42,   44,   
                            46,   48,   50,   52,   54,   56,   58,   60,   
                            65,   70,   75,   80,   85,   90,   95,  100,  
                           110,  120,  130,  140,  150,  175,  200,  225,  
                           250,  300,  400,  500,  760] # (mmHg)
dissociationCurveSaO2 = [ 0.00, 0.60, 1.19, 2.56, 4.37, 6.68, 9.58,12.96,
                         16.89,21.40,26.50,32.12,37.60,43.14,48.27,53.16,
                         57.54,61.69,65.16,68.63,71.94,74.69,77.29,79.55,
                         81.71,83.52,85.08,86.59,87.70,88.93,89.95,90.85,
                         92.73,94.06,95.10,95.84,96.42,96.88,97.25,97.49,
                         97.91,98.21,98.44,98.62,98.77,99.03,99.20,99.32,
                         99.41,99.53,99.65,99.72,100.0] # (%)
"""
dissociationCurvePaO2 = [    0,    1,    2,    4,    6,    8,   10,   12,
                            14,   16,   18,   20,   22,   24,   26,   28,   
                            30,   32,   34,   36,   38,   40,   42,   44,   
                            46,   48,   50,   52,   54,   56,   58,   60,   
                            65,   70,   75,   80,   85,   90,   95,  100,  
                           110,  120,  130,  140,  150,  175,  200] # (mmHg)
dissociationCurveSaO2 = [ 0.00, 0.60, 1.19, 2.56, 4.37, 6.68, 9.58,12.96,
                         16.89,21.40,26.50,32.12,37.60,43.14,48.27,53.16,
                         57.54,61.69,65.16,68.63,71.94,74.69,77.29,79.55,
                         81.71,83.52,85.08,86.59,87.70,88.93,89.95,90.85,
                         92.73,94.06,95.10,95.84,96.42,96.88,97.25,97.49,
                         97.91,98.21,98.44,98.62,98.77,100.0,200.0] # (%)
def saturationToPressure(saturation):
  # assert(saturation >= 0) # assert(saturation <= 100)
  indexL = bisect_right(dissociationCurveSaO2, saturation)-1
  indexR = bisect_left (dissociationCurveSaO2, saturation)
  return map_range(saturation,dissociationCurveSaO2[indexL],dissociationCurveSaO2[indexR],
                              dissociationCurvePaO2[indexL],dissociationCurvePaO2[indexR])

def pressureToSaturation(pressure):
  # assert(pressure >= 0)   # assert(pressure <= 760)
  indexL = bisect_right(dissociationCurvePaO2, saturation)-1
  indexR = bisect_left (dissociationCurvePaO2, saturation)
  return map_range(saturation,dissociationCurvePaO2[indexL],dissociationCurvePaO2[indexR],
                              dissociationCurveSaO2[indexL],dissociationCurveSaO2[indexR])

# https://stackoverflow.com/questions/70643627/python-equivalent-for-arduinos-map-function
def map_range(x, in_min, in_max, out_min, out_max):
  if (out_min == out_max):
    return out_min
  return (x - in_min) * (out_max - out_min) // (in_max - in_min) + out_min

def copyList(toList, fromList):
  assert(len(toList)==len(fromList))
  for i in range(0,len(fromList)):
    toList[i] = fromList[i]


