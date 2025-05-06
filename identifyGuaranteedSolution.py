"""
    *** isGuaranteedSolution specifies whether a minimum flux distribution of a node X can be the maximized
    *** In this function all of the desired chemicals and their threshold are defined and set.
"""
from cobra.flux_analysis import flux_variability_analysis

def isGuaranteedSolution(X, model):

    #
    # As mentioned in the 4th step of README.md file, the desired chemical and biomass (and their threshold) reactions should be defined in the following lines.
    chemical = "r_1761" # ethanol
    biomass = "r_2111"
    Th_chemical = 8.0
    Th_biomass = 0.1
    #
    #
   
   
    fva_result = flux_variability_analysis(model, None, loopless = False, fraction_of_optimum=0.99, pfba_factor = 1.1, processes = 8)

    if (fva_result['minimum'][chemical] >Th_chemical and fva_result['minimum'][biomass] > Th_biomass):
            X.chemical = fva_result['minimum'][chemical]
            X.biomass = fva_result['minimum'][biomass]

    return 
