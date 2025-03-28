from predict import*

# ======== computational constants ========

BLOOD_ARRAY_LENGTH = 32
TIME_TICK_STEPS = 512 # 128 # must be an order of magnitude larger than BLOOD_ARRAY_LENGTH for small blood capillary heights
USE_TWO_SIDED_DIFFUSION = True
GAS_PO2 = 760.0 # (mmHg)

inputSO2 = 65.0 # (%)
membraneThickness = 0.0050 # (50 um = 0.0050 cm)
capillaryHeight = 0.0200 # (200 um = 0.0200 cm)
gasExchangeSurfaceArea = 0.0297*1.8 # (cm^2) # 300um x 10mm
bloodFlow = 164.0/10875/60 # (164 mL/min across 10875 capillaries)

# ======== should be 95 ========

outputSO2 = calculate(capillaryHeight, membraneThickness, bloodFlow, inputSO2, gasExchangeSurfaceArea,
                      gasPO2 = GAS_PO2,
                      timeTickSteps=TIME_TICK_STEPS, bloodArrayLength=BLOOD_ARRAY_LENGTH,
                      useTwoSidedDiffusion=USE_TWO_SIDED_DIFFUSION)
print("outputSO2:", outputSO2)

#========= ISOLATE =========

bloodFlow = 500.0/10875/60 # (164 mL/min across 10875 capillaries)

# ======== should be 49 mmHg =========

outputSO2 = calculate(capillaryHeight, membraneThickness, bloodFlow, inputSO2, gasExchangeSurfaceArea,
                      gasPO2 = GAS_PO2,
                      timeTickSteps=TIME_TICK_STEPS, bloodArrayLength=BLOOD_ARRAY_LENGTH,
                      useTwoSidedDiffusion=USE_TWO_SIDED_DIFFUSION)
print("outputSO2:", outputSO2)

# ======== should be 52 mmHg (no membrane) ========

membraneThickness = 0.0000001

outputSO2 = calculate(capillaryHeight, membraneThickness, bloodFlow, inputSO2, gasExchangeSurfaceArea,
                      gasPO2 = GAS_PO2,
                      timeTickSteps=TIME_TICK_STEPS, bloodArrayLength=BLOOD_ARRAY_LENGTH,
                      useTwoSidedDiffusion=USE_TWO_SIDED_DIFFUSION)
print("outputSO2:", outputSO2)

# ======== should be 115 mmHg (no blood) ========

# TODO: Manually change DIFFUSIVITY_O2 to 128 or some other high number
