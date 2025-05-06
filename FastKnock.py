"""
    *** FastKnock.py is the main method of the algorithm
    *** Processes are defined in this method for traversing the tree branches.
    *** Each Process writes its solutions into the separate files.
"""

# 実行は FastKnock の仮想環境に入った後 (. .venv/bin/activate), uv run FastKnock.py

from Node import Node
import PreProcessing as PreProc
from identifyTargetSpace import identifyTargetSpace
from constructSubTree import constructSubTree
from traverseTree import procTraverseTree
import sys
import multiprocessing
from findingCoKnockOut import findingCoKnockoutRxns
import mergeFile 
import cobra



def construct_data_structure (number):
    data_struct = [[] for i in range(int(number)+1)]
    return data_struct

"""
    * medium culture is determined and 
    * by preprocessing the model Removable and
    * coKnockoutRxns sets are obtained in this method
"""
def model_preparation (model_name):
    ori_model = cobra.io.load_matlab_model(model_name)
    model = ori_model.copy()
    model.solver = "cplex"
    
    #
    # As mentioned in the third step of README.md file, the medium culture should be specified in the following lines.
    # すべての可逆型取り込み反応 (uptake) をいったん閉じる
    for ex_rxn in model.exchanges:
        if "_REV" in ex_rxn.id:
            ex_rxn.bounds = 0.0, 0.0
    # 各種排出反応を調整
    model.reactions.r_1992.bound = 0.0, 0.0 # oxygen
    model.reactions.r_1714.bound = 0.0, 0.0 # glucose
    model.reactions.r_1718.bound = 0.0, 0.0 # xylose
    model.reactions.r_1761.upper_bound = 1000.0 # ethanol
    # 各種反応の制約を定義
    model.reactions.r_2111.bounds = 0.0, 1000.0 # biomass reacion
    model.reactions.r_1714_REV.bounds = 0.0, 10.0 # glucose exchange
    model.reactions.r_1718_REV.bounds = 0.0, 10.0 # xylose exchange
    model.reactions.r_1093No1.bounds = 0.0, 1000.0 # xylose reductase
    model.reactions.r_2104_REV.bounds = 0.0, 0.0 # xylitol exchange
    model.reactions.r_1092No1.bounds = 0.0, 1000.0 # xylitol dehydrogenase
    model.reactions.r_4490.bounds = 0.0, 1000.0 # xylulose⇒xylitol
    model.reactions.r_1654_REV.bounds = 0.0, 1000.0 # ammonium exchange
    model.reactions.r_1992_REV.bounds = 0.0, 10.0 # oxygen exchange
    model.reactions.r_1832_REV.bounds = 0.0, 1000.0 # proton exchange
    model.reactions.r_1861_REV.bounds = 0.0, 1000.0 # iron(2+) exchange
    model.reactions.r_2005_REV.bounds = 0.0, 1000.0 # phosphate exchange
    model.reactions.r_2020_REV.bounds = 0.0, 1000.0 # potassium exchange
    model.reactions.r_2049_REV.bounds = 0.0, 1000.0 # sodium exchange
    model.reactions.r_2060_REV.bounds = 0.0, 1000.0 # sulphate exchange
    model.reactions.r_2100_REV.bounds = 0.0, 1000.0 # water exchange
    model.reactions.r_4593_REV.bounds = 0.0, 1000.0 # chloride change
    model.reactions.r_4594_REV.bounds = 0.0, 1000.0 # Cu2(+) exchange
    model.reactions.r_4595_REV.bounds = 0.0, 1000.0 # Mn(2+) exchange
    model.reactions.r_4596_REV.bounds = 0.0, 1000.0 # Zn(2+) exchange
    model.reactions.r_4597_REV.bounds = 0.0, 1000.0 # Mg(2+) exchange
    model.reactions.r_4600_REV.bounds = 0.0, 1000.0 # Ca(2+) exchange
      

    coKnockoutRxns = findingCoKnockoutRxns(model)
    PreProc.identifying_eliminate_list(model)

    eliminate_list = []
    with open("eliminate_list.txt" , 'r') as f:
        for line in f:
            if line not in ['\n', '\r\n', '\t']:
                line = line.replace('\n', '')
                eliminate_list.append(line)
    model = PreProc.remove_blocked_rxn(model)

    Removable = []
    for i in model.reactions:
        if i.id not in eliminate_list:
            Removable.append(i.id)

    return model, Removable, coKnockoutRxns


def print_node(X):
    print (X.level)
    print (X.deleted_rxns)
    print (X.target_space)
    print (X.flux_dist)


""" The solutions are written in the file by this method"""
def writeInFile (level, solution):

    file_name = str(level)+".txt"
    with open(file_name , 'w') as f:
        for X in solution:
            sol = [X.deleted_rxns, X.biomass, X.chemical]
            f.write(str(sol))
            f.write('\n')
    f.close()
    
    

def main():
    
    target_level = input("Enter your desired level (maximum number of simultaneous reaction knockout ):  ")

    guaranteed_flag = 0
    if ( int(input("If you want to use FastKnock in guaranteed mood, press 1 othrwise press 0: ")) == 1):
        guaranteed_flag = 1

    num_of_processes = input("Enter number of proccessor cores: ")

    queue = construct_data_structure (target_level)
    checked = construct_data_structure (target_level)
    solution = construct_data_structure (target_level)
    guaranteed_solution = construct_data_structure (target_level)

    #  
    #As mentioned in the second step of README.md file, the name of model should be replaced in the next line.
    model, Removable, coKnockoutRxns = model_preparation("ecYeastGEM_batch (CellFactory-ecYeastGEM).mat")

    root = Node (0,[] , [] , [],0,0)
    root = identifyTargetSpace(root, model, Removable,coKnockoutRxns)

    all_fba_call = [[] for i in range(int(target_level)+1)]
    for i in range (int(target_level)+1):
        all_fba_call[i] = 0
    
    level_one, queue, checked, solution, all_fba_call = constructSubTree(root, target_level, checked, queue, solution, model, Removable, all_fba_call, coKnockoutRxns, guaranteed_flag, guaranteed_solution )



    queue_index = 0
    elements = (len(queue[int(level_one)]) - 3)/(int(num_of_processes) - 2)

    for j in range (int(num_of_processes)):
        globals()['queue_p%s' % (j+1)] = construct_data_structure (target_level)
        globals()['checked_p%s' % (j+1)] = construct_data_structure (target_level)
        globals()['solution_p%s' % (j+1)] = construct_data_structure (target_level)
        globals()['guaranteed_solution_p%s' % (j+1)] = construct_data_structure (target_level)

        globals()['all_fba_call_p%s' % (j+1)] = construct_data_structure (target_level)
        for i in range (int(target_level)+1):
            globals()['all_fba_call_p%s' % (j+1)][i] = 0

        if (j == 0):
            globals()['queue_p%s' % (j+1)][int(level_one)] = queue[int(level_one)][0:1]
        
        elif (j == 1):
            globals()['queue_p%s' % (j+1)][int(level_one)] = queue[int(level_one)][1:3]
            globals()['checked_p%s' % (j+1)][int(level_one)] = queue[int(level_one)][0:1]

             
        else:
            globals()['queue_p%s' % (j+1)][int(level_one)] = queue[int(level_one)][int(elements)*(j-1):int(elements)*j]

        if(j > 1):
            globals()['checked_p%s' % (j+1)][int(level_one)] = queue[int(level_one)][0:int(elements)*(j-1)]

        
        globals()[f"p{j+1}"] = multiprocessing.Process(target = procTraverseTree, args=(f"p{j+1}", level_one, globals()['queue_p%s' % (j+1)], globals()['checked_p%s' % (j+1)], globals()['solution_p%s' % (j+1)], target_level, model, Removable, globals()['all_fba_call_p%s' % (j+1)], coKnockoutRxns , guaranteed_flag, globals()['guaranteed_solution_p%s' % (j+1)] ))



    for j in range (int(num_of_processes)):
       globals()['p%s' % (j+1)].start()

    for j in range (int(num_of_processes)):
       globals()['p%s' % (j+1)].join()

    mergeFile.merge(int(num_of_processes), int(target_level))

if __name__ =='__main__':
    main()
