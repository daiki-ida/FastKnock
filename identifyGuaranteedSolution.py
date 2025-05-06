"""
    *** isGuaranteedSolution specifies whether a minimum flux distribution of a node X can be the maximized
    *** In this function all of the desired chemicals and their threshold are defined and set.
"""
from cobra.flux_analysis.fastcc import fastcc
from cobra.flux_analysis import flux_variability_analysis
from cobra.exceptions import OptimizationError

def isGuaranteedSolution(X, model):

    # 1. FASTCC で一貫性モデルを取得
    try:
        consistent_model = fastcc(model, zero_cutoff=1e-6)
    except Exception as e:
        raise RuntimeError(f"FASTCC failed: {e}")

    #
    # As mentioned in the 4th step of README.md file, the desired chemical and biomass (and their threshold) reactions should be defined in the following lines.
    chemical = "r_1761" # ethanol
    biomass = "r_2111"
    Th_chemical = 20.0
    Th_biomass = 0.2
    #
    #
   
    # Loopless FVA
    try:
        fva_result = flux_variability_analysis(model = consistent_model, 
                                               reaction_list = None, 
                                               loopless = True, 
                                               fraction_of_optimum=0.99, 
                                               pfba_factor = 1.1, 
                                               processes = 8)
    except OptimizationError as oe:
        # unbounded や infeasible が発生した場合は例外をログに出して終了
        print(f"[Warning] FVA failed on consistent_model: {oe}")
        return

    # 閾値チェックとXへの格納
    if (fva_result['minimum'][chemical] > Th_chemical and fva_result['minimum'][biomass] > Th_biomass):
            X.chemical = fva_result['minimum'][chemical]
            X.biomass = fva_result['minimum'][biomass]

    return 
