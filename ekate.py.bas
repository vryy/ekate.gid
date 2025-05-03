##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025-2026 by Hoang-Giang Bui (UoB)            #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on __DaTe__
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import rEpLaCeMeNtStRiNg_include as simulation_include
from rEpLaCeMeNtStRiNg_include import **
model_name_ = 'rEpLaCeMeNtStRiNg'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging)
    model.InitializeModel()

*if(strcmp(GenData(Simulation_Script),"standard")==0)
    time = 0.0
    delta_time = *GenData(time_step_length)

    for step in range( 0, *GenData(time_steps) ):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        model.WriteOutput( time )
        print("###################")
        print("step "+str(time)+" done.")
        print("###################")
    print("Calculation done")
*else
    # ============================================ #
    # |       USER CALCULATION SCRIPT            | #
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv #
*endif


    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    main(logging=False, output=False)

def tag():
    return "unknown"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
