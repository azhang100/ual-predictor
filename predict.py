import visualizer

# references:
# DissociationCurve = https://derangedphysiology.com/main/cicm-primary-exam/
#                     respiratory-system/Chapter-112/oxyhaemoglobin-dissociation-curve
# AdvectionDiffusion = https://pmc.ncbi.nlm.nih.gov/articles/PMC2584207/pdf/nihms77385.pdf
# MorePromises = 

# ======== user input ========

MICROMOL_O2_PER_ML_BOUND = 7.5 # (μM of O2 / mL of blood)
MICROMOL_O2_PER_ML_FREE = 0.0139 # (μM of O2 / mL of blood / mmHg)

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
  bloodAdmittance = timeTick * DIFFUSIVITY_O2 * MICROMOL_O2_PER_ML_FREE * bloodDiffusionDistance**(-2)
  print("bloodAdmittance:", bloodAdmittance)
  membraneAdmittance = timeTick * MEMBRANE_PERMEABILITY / membraneThickness / bloodDiffusionDistance
  print("timeTick:", timeTick)
  print("membraneThickness:", membraneThickness)
  print("membraneAdmittance:", membraneAdmittance)
  
  print("dwellTime:",dwellTime)
  
  bloodArray    = [convertO2(inputSO2,dcSaO2,dcUMol)]*bloodArrayLength
  bloodArrayNew = [convertO2(inputSO2,dcSaO2,dcUMol)]*bloodArrayLength
  for t in range(0, timeTickSteps):
    print("time step:", t, "current SO2:", convertO2(sum(bloodArray)/len(bloodArray), dcUMol, dcSaO2))
    bloodArray[0                 ] = membraneDiffuse(bloodArray[0                 ],gasPO2,membraneAdmittance)
    bloodArray[bloodArrayLength-1] = membraneDiffuse(bloodArray[bloodArrayLength-1],gasPO2,membraneAdmittance*useTwoSidedDiffusion)
    # draw(bloodArray)

    copyList(bloodArrayNew, bloodArray)
    for i in range(0,len(bloodArray)-1):
      (bloodArrayNew[i],bloodArrayNew[i+1]) = bloodDiffuse(bloodArray[i], bloodArray[i+1], bloodAdmittance)
    copyList(bloodArray, bloodArrayNew)
      
    if (t%10==0):
      draw(bloodArray, flush=True)
      
  return convertO2(sum(bloodArray)/len(bloodArray), dcUMol, dcSaO2)

# diffuse from cellA to cellB
def bloodDiffuse(cellB, cellA, admittance):
  pressureB = convertO2(cellB, dcUMol, dcPaO2)
  pressureA = convertO2(cellA, dcUMol, dcPaO2)
  diffusionDiff = admittance*(pressureA-pressureB)
  # diffusionMicroDiff = timeTick * DIFFUSIVITY_O2 * MICROMOL_O2_PER_ML_FREE * (pressureA-pressureB) *
  #                     (bloodDiffusionDistance)^(-2)

  diffusionDiffMax = (cellA-cellB)/2
  if (abs(diffusionDiff) > abs(diffusionDiffMax)):
    print("WARN: blood diffusionDiff maxxed out, limiting")
    diffusionDiff /= 2
  
  newCellB = cellB + diffusionDiff
  newCellA = cellA - diffusionDiff
  
  return (newCellB, newCellA)

# returns bloodUMol
def membraneDiffuse(bloodUMol, pressureGas, admittance):
  pressureB = convertO2(bloodUMol, dcUMol, dcPaO2)
  diffusionDiff = admittance*(pressureGas-pressureB)
  # diffusionMicroMol = timeTick * (pressureGas-pressureB) / 
  #                     (MEMBRANE_RESISTANCE_PER_DISTANCE * membraneThickness) / bloodDiffusionDistance
  result = bloodUMol+diffusionDiff

  # limit maximum transfer so that blood pp remains less than gas pp
  bloodUMolMax = convertO2(pressureGas, dcPaO2, dcUMol)
  if (result>bloodUMolMax):
    print("WARN: membrane bloodUMol maxxed out, limiting")
    result = bloodUMolMax
  
  return result

# ================== MORE HELPER FUNCTIONS =========================

from bisect import bisect_left, bisect_right

# source: DissociationCurve

dcPaO2 = [  0.0,  1.0,  2.0,  4.0,  6.0,  8.0, 10.0, 12.0,
           14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 
           30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 
           46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 
           65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0,100.0,
          110.0,120.0,130.0,140.0,150.0,175.0,200.0,225.0,
          250.0,300.0,400.0,500.0,     1000.0,     2000.0] # (mmHg) # the last values would have been 900 if this was linear
dcSaO2 = [ 0.00, 0.60, 1.19, 2.56, 4.37, 6.68, 9.58,12.96,
          16.89,21.40,26.50,32.12,37.60,43.14,48.27,53.16,
          57.54,61.69,65.16,68.63,71.94,74.69,77.29,79.55,
          81.71,83.52,85.08,86.59,87.70,88.93,89.95,90.85,
          92.73,94.06,95.10,95.84,96.42,96.88,97.25,97.49,
          97.91,98.21,98.44,98.62,98.77,99.03,99.20,99.32,
          99.41,99.53,99.65,99.72,      100.0,      100.5] # (%)
"""
dcPaO2 = [    0,    1,    2,    4,    6,    8,   10,   12,
              14,   16,   18,   20,   22,   24,   26,   28,   
              30,   32,   34,   36,   38,   40,   42,   44,   
              46,   48,   50,   52,   54,   56,   58,   60,   
              65,   70,   75,   80,   85,   90,   95,  100,  
             110,  120,  130,  140,  150,  175,  200] # (mmHg)
dcSaO2 = [ 0.00, 0.60, 1.19, 2.56, 4.37, 6.68, 9.58,12.96,
           16.89,21.40,26.50,32.12,37.60,43.14,48.27,53.16,
           57.54,61.69,65.16,68.63,71.94,74.69,77.29,79.55,
           81.71,83.52,85.08,86.59,87.70,88.93,89.95,90.85,
           92.73,94.06,95.10,95.84,96.42,96.88,97.25,97.49,
           97.91,98.21,98.44,98.62,98.77,100.0,200.0] # (%)
"""
dcUMol = []
for i in range(len(dcPaO2)):
  microMolOxygen = dcPaO2[i]*MICROMOL_O2_PER_ML_FREE + dcSaO2[i]/100.0*MICROMOL_O2_PER_ML_BOUND
  dcUMol.append(microMolOxygen)
assert len(dcUMol) == len(dcPaO2)
assert len(dcUMol) == len(dcSaO2)

def convertO2(value,dcFrom, dcTo):
  indexL = bisect_right(dcFrom, value)-1
  indexR = bisect_left (dcFrom, value)
  
  assert (indexR < len(dcFrom)), "dc table out of bounds"
  assert (dcFrom[indexL] <= value)
  assert (value <= dcFrom[indexR])
  
  return map_range(value,dcFrom[indexL],dcFrom[indexR],
                           dcTo[indexL],  dcTo[indexR])

# https://stackoverflow.com/questions/70643627/python-equivalent-for-arduinos-map-function
def map_range(x, in_min, in_max, out_min, out_max):
  if (out_min == out_max):
    return out_min
  ratio = (x - in_min)/(in_max - in_min)
  return out_min + (out_max - out_min)*ratio
  # return (x - in_min) * (out_max - out_min) // (in_max - in_min) + out_min

#======== utility =========

def copyList(toList, fromList):
  assert(len(toList)==len(fromList))
  for i in range(0,len(fromList)):
    toList[i] = fromList[i]
    
def draw(bloodArray, flush=False):
  for uMol in bloodArray:
    bloodSat = convertO2(uMol, dcUMol, dcSaO2)
    text = "{:5.2f}".format(bloodSat).zfill(5)
    if (len(text) > 5):
      text = text[:-1]
    draw.staticText += "|"+text
  draw.staticText += "|\n"

  if (flush):
    print(draw.staticText)
    draw.staticText = ""
  
draw.staticText = ""

def warnIf(condition, text):
  if (condition):
    print(text)

#===========test=====
"""
print(dcUMol)

SaO2 = 65.0
print(" 65%:", SaO2)
uMol = convertO2(SaO2,dcSaO2, dcUMol)
print("UMol:", uMol)
SaO2 = convertO2(uMol,dcUMol, dcSaO2)
print("SaO2:", SaO2)
"""
