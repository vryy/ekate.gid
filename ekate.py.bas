##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on __DaTe__ __TiMe__
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./rEpLaCeMeNtStRiNg.gid')
import rEpLaCeMeNtStRiNg_include
from rEpLaCeMeNtStRiNg_include import **
# calculate insitu-stress for geology_virgin.gid
model = rEpLaCeMeNtStRiNg_include.Model('rEpLaCeMeNtStRiNg',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

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
print "Calculation done"
*else
*# user-defined script is used (will be appended automatically)
# =====================
# | USER SCRIPT FOR CALCULATION OF EKATE.GID |
# vvvvvvvvvvvvvvvvvvvvv
*endif


##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
