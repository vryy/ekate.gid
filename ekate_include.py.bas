##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
import time as time_module
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import **
from KratosMultiphysics.StructuralApplication import **
from KratosMultiphysics.ExternalSolversApplication import **
from KratosMultiphysics.MKLSolversApplication import **
from KratosMultiphysics.BRepApplication import **
*if(strcmp(GenData(Ekate_Auxiliary_Application),"1")==0)
from KratosMultiphysics.EkateAuxiliaryApplication import **
from KratosMultiphysics.ExternalConstitutiveLawsApplication import **
*endif
*if(strcmp(GenData(Enable_Mortar_Contact),"1")==0||strcmp(GenData(Mortar_Application),"1")==0)
from KratosMultiphysics.MortarApplication import **
*endif
*if(strcmp(GenData(Layer_Application),"1")==0||strcmp(GenData(VTK_Output),"1")==0)
from KratosMultiphysics.LayerApplication import **
*endif
*if(strcmp(GenData(Finite_Cell_Application),"1")==0)
from KratosMultiphysics.FiniteCellApplication import **
*endif
*if(strcmp(GenData(Finite_Cell_Structural_Application),"1")==0)
from KratosMultiphysics.FiniteCellStructuralApplication import **
*endif
*if(strcmp(GenData(Soil_Mechanics_Application),"1")==0)
from KratosMultiphysics.SoilMechanicsApplication import **
*if(strcmp(GenData(Finite_Cell_Application),"1")==0)
from KratosMultiphysics.FiniteCellSoilMechanicsApplication import **
*endif
*endif
*if(strcmp(GenData(Plate_And_Shell_Application),"1")==0)
from KratosMultiphysics.PlateAndShellApplication import **
*endif
*if(strcmp(GenData(P4est_Application),"1")==0)
from KratosMultiphysics.mpi import **
from KratosMultiphysics.P4estApplication import **
*endif
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, results_path ):
        #setting the domain size for the problem to be solved
*set var domainsize(int)=GenData(Dimension,int)
        self.domain_size = *domainsize
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart(problem_name)
        self.path = path
        self.results_path = results_path+os.sep
        self.problem_name = problem_name
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################
        # reading simulation parameters
*set var ntimesteps(int)=GenData(time_steps,int)
        number_of_time_steps = *ntimesteps
        self.analysis_parameters = {}
        # content of analysis_parameters:
        # perform_contact_analysis_flag
        # penalty value for normal contact
        # maximum number of uzawa iterations
        # friction coefficient
        # penalty value for frictional contact
        # contact_double_check_flag
        # contact_ramp_penalties_flag
        # maximum penalty value for normal contact
        # ramp criterion for normal contact
        # ramp factor for normal contact
        # maximum penalty value for frictional contact
        # ramp criterion for frictional contact
        # ramp factor for frictional contact
        # analysis type: static (0), quasi-static (1) or dynamic (2)
*if(strcmp(GenData(analysis_type),"static")==0)
*set var analysistype(int)=0
*endif
*if(strcmp(GenData(analysis_type),"quasi-static")==0)
*set var analysistype(int)=1
*endif
*if(strcmp(GenData(analysis_type),"dynamic")==0)
*set var analysistype(int)=2
*endif
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
        perform_contact_analysis_flag = True
        # performing contact analysis: reading contact parameters
*set var penalty(real)=GenData(Penalty_Value,real)
        penalty = *penalty
*set var maxuzawa(int)=GenData(Max_Uzawa_Iterations,int)
        maxuzawa = *maxuzawa
*set var friction(real)=GenData(Friction_Coefficient,real)
        friction = *friction
*set var frictionpenalty(real)=GenData(Friction_Penalty_Value,real)
        frictionpenalty = *frictionpenalty
*if(strcmp(GenData(Bilateral_Contact),"1")==0)
        contact_double_check_flag = True
*else
        contact_double_check_flag = False
*endif
*if(strcmp(GenData(Ramp_Penalties),"1")==0)
        contact_ramp_penalties_flag = True
*set var maxpenalty(real)=GenData(Maximum_Penalty,real)
        maxpenalty = *maxpenalty
*set var rampcriterion(real)=GenData(Ramp_Criterion,real)
        rampcriterion = *rampcriterion
*set var rampfactor(real)=GenData(Ramp_Factor,real)
        rampfactor = *rampfactor
*set var fricmaxpenalty(real)=GenData(Friction_Maximum_Penalty,real)
        fricmaxpenalty = *fricmaxpenalty
*set var fricrampcriterion(real)=GenData(Friction_Ramp_Criterion,real)
        fricrampcriterion = *fricrampcriterion
*set var fricrampfactor(real)=GenData(Friction_Ramp_Factor,real)
        fricrampfactor = *fricrampfactor
*else
        contact_ramp_penalties_flag = False
        maxpenalty = penalty
        rampcriterion = 0.0
        rampfactor = 0.0
        fricmaxpenalty = penalty
        fricrampcriterion = 0.0
        fricrampfactor = 0.0
*endif
*else
        perform_contact_analysis_flag = False
        penalty = 0.0
        maxuzawa = 0.0
        friction = 0.0
        frictionpenalty = 0.0
        contact_double_check_flag = False
        contact_ramp_penalties_flag = False
        maxpenalty = 0.0
        rampcriterion = 0.0
        rampfactor = 0.0
        fricmaxpenalty = 0.0
        fricrampcriterion = 0.0
        fricrampfactor = 0.0
*endif
        self.analysis_parameters['dimension'] = self.domain_size # *domainsize
        self.analysis_parameters['perform_contact_analysis_flag'] = perform_contact_analysis_flag
        self.analysis_parameters['penalty'] = penalty
        self.analysis_parameters['maxuzawa'] = maxuzawa
        self.analysis_parameters['friction'] = friction
        self.analysis_parameters['frictionpenalty'] = frictionpenalty
        self.analysis_parameters['contact_double_check_flag'] = contact_double_check_flag
        self.analysis_parameters['contact_ramp_penalties_flag'] = contact_ramp_penalties_flag
        self.analysis_parameters['maxpenalty'] = maxpenalty
        self.analysis_parameters['rampcriterion'] = rampcriterion
        self.analysis_parameters['rampfactor'] = rampfactor
        self.analysis_parameters['fricmaxpenalty'] = fricmaxpenalty
        self.analysis_parameters['fricrampcriterion'] = fricrampcriterion
        self.analysis_parameters['fricrampfactor'] = fricrampfactor
*if(strcmp(GenData(Plot_Matrix_Structure),"0")==0)
        self.analysis_parameters['print_sparsity_info_flag'] = False
*else
        self.analysis_parameters['print_sparsity_info_flag'] = True
*endif
        self.analysis_parameters['analysis_type'] = *analysistype
        self.analysis_parameters['dissipation_radius'] = 0.1
        self.analysis_parameters['decouple_build_and_solve'] = True
        self.analysis_parameters['solving_scheme'] = 'monolithic'
*if(strcmp(GenData(Stop_Newton_Raphson_If_Not_Converged),"1")==0)
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
*else
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = False
*endif
        self.analysis_parameters['list_dof'] = True

*if(strcmp(GenData(Absolute_Tolerance),"custom")==0)
*set var abstol(real)=GenData(Custom_Absolute_Tolerance,real)
*else
*set var abstol(real)=GenData(Absolute_Tolerance,real)
*endif
        self.abs_tol = *abstol
*if(strcmp(GenData(Relative_Tolerance),"custom")==0)
*set var reltol(real)=GenData(Custom_Relative_Tolerance,real)
*else
*set var reltol(real)=GenData(Relative_Tolerance,real)
*endif
        self.rel_tol = *reltol

        ## generating solver
*if(strcmp(GenData(Enable_Mortar_Contact),"1")==0)
        import mortar_gpts_contact_strategy
*if(strcmp(GenData(analysis_type),"static")==0)
        self.solver = mortar_gpts_contact_strategy.SampleSolver(self.model_part, self.abs_tol, self.rel_tol, self.analysis_parameters)
*else
        self.solver = mortar_gpts_contact_strategy.SampleSolverEkateQuasiStatic(self.model_part, self.abs_tol, self.rel_tol, self.analysis_parameters)
*endif
        mortar_gpts_contact_strategy.AddVariables( self.model_part )
*else
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        #import ekate_solver_parallel
        #self.solver = ekate_solver_parallel.EkateSolver( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
        #ekate_solver_parallel.AddVariables( self.model_part )
*endif
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
*if(strcmp(GenData(Calculate_Reactions),"1")==0)
        (self.solver).CalculateReactionFlag = True
*else
        (self.solver).CalculateReactionFlag = False
*endif
*if(strcmp(GenData(Read_Deactivation_File),"1")==0)
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        self.element_assignments = {}
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
                self.element_assignments[int(val_set[1])] = int(val_set[2])
        print("input data read OK")
        #print "+++++++++++++++++++++++++++++++++++++++"
        #for node in self.model_part.Nodes:
        #    print(node)
        #print("+++++++++++++++++++++++++++++++++++++++")
*endif
        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## POST_PROCESSING ###############################################
        ##################################################################

        self.write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        self.write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
*if(strcmp(GenData(Output_Format),"ASCII")==0)
        self.post_mode = GiDPostMode.GiD_PostAscii
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.multi_file_flag = MultiFileFlag.MultipleFiles
*else
        self.multi_file_flag = MultiFileFlag.SingleFile
*endif
*else
        self.post_mode = GiDPostMode.GiD_PostBinary
        self.multi_file_flag = MultiFileFlag.MultipleFiles
*endif
*if(strcmp(GenData(Layer_Application),"1")==0)
        self.gid_io = SDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )
*else
        self.gid_io = StructuralGidIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )
*endif
*if(strcmp(GenData(VTK_Output),"1")==0)
*if(strcmp(GenData(VTK_Output_Format),"ASCII")==0)
        vtk_post_mode = VTKPostMode.VTK_PostAscii
*elseif(strcmp(GenData(VTK_Output_Format),"Binary")==0)
        vtk_post_mode = VTKPostMode.VTK_PostBinary
*endif
*if(strcmp(GenData(VTK_Output_Type),"VTU")==0)
        self.vtk_io = VtkVTUIO( self.results_path+"vtk_results/"+self.problem_name, vtk_post_mode )
*elseif(strcmp(GenData(VTK_Output_Type),"VTM")==0)
        self.vtk_io = VtkVTMIO( self.results_path+"vtk_results/"+self.problem_name, vtk_post_mode )
*endif
        if not os.path.isdir(self.results_path+"vtk_results"):
            os.mkdir(self.results_path+"vtk_results")
*endif

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
*if(strcmp(GenData(Solver),"BiCGStabLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        #preconditioner = ILU0Preconditioner()
        preconditioner = Preconditioner()
        plinear_solver = BICGSTABSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"CGLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        #preconditioner = ILU0Preconditioner()
        preconditioner = Preconditioner()
        #plinear_solver = DeflatedCGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner,1)
        plinear_solver = CGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"SuperLU")==0)
        plinear_solver = SuperLUSolver()
*endif
*if(strcmp(GenData(Solver),"SkylineLUFactorizationSolver")==0)
        plinear_solver = SkylineLUFactorizationSolver()
*endif
*if(strcmp(GenData(Solver),"GMRESSolver")==0)
        #preconditioner = Preconditioner()
        #plinear_solver = GMRESSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
        plinear_solver = MKLGMRESSolver()
*endif
*if(strcmp(GenData(Solver),"Pardiso")==0)
        plinear_solver = MKLPardisoSolver()
*endif
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)
        (self.solver.solver).max_iter = 10 #control the maximum iterations of Newton Raphson loop
*if(strcmp(GenData(Move_Mesh_Flag),"1")==0)
        (self.solver.solver).MoveMeshFlag = True
*else
        (self.solver.solver).MoveMeshFlag = False
*endif

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )

*if(strcmp(GenData(Topology_Optimization),"1")==0)
        ##################################################################
        ## TOPOLOGY OPTIMIZATION SETUP ####################################
        ##################################################################
*set var mvolfrac(real)=GenData(Volume_Fraction,real)
*set var mrmin(real)=GenData(Search_Radius,real)
        volfrac = *mvolfrac
        rmin = *mrmin
*if(strcmp(GenData(Filter_Type),"Sensitivity_Filter")==0)
        ft = 1
*elseif(strcmp(GenData(Filter_Type),"Density_Filter")==0)
        ft = 2
*else
        ft = 0
*endif
        self.topopt_proc = TopologyUpdateProcess(self.model_part, volfrac, rmin, ft)
        self.topopt_proc.SetBinSize(100)
        self.solver.solver.attached_processes.append(self.topopt_proc)
*endif

    def AddDofsForNodes(self, nodes):
*if(strcmp(GenData(Perform_MultiFlow_Analysis),"1")==0)
        for node in nodes:
            node.AddDof( WATER_PRESSURE )
*if(strcmp(GenData(Perform_ThreePhase_Analysis),"1")==0)
            node.AddDof( AIR_PRESSURE )
*endif
*endif
*if(strcmp(GenData(Enable_Mortar_Contact),"1")==0)
        import mortar_gpts_contact_strategy
        mortar_gpts_contact_strategy.AddDofsForNodes( nodes )
*if(strcmp(GenData(Perform_MultiFlow_Analysis),"1")==0)
        mortar_gpts_contact_strategy.AddFluidDofsForNodes( nodes )
*endif
*else
        import structural_solver_advanced
        structural_solver_advanced.AddDofsForNodes( nodes )
*endif

    def AddDofs(self, model_part):
*if(strcmp(GenData(Perform_MultiFlow_Analysis),"1")==0)
        for node in model_part.Nodes:
            node.AddDof( WATER_PRESSURE )
*if(strcmp(GenData(Perform_ThreePhase_Analysis),"1")==0)
            node.AddDof( AIR_PRESSURE )
*endif
*endif
*if(strcmp(GenData(Enable_Mortar_Contact),"1")==0)
        import mortar_gpts_contact_strategy
        mortar_gpts_contact_strategy.AddDofs( model_part )
*if(strcmp(GenData(Perform_MultiFlow_Analysis),"1")==0)
        mortar_gpts_contact_strategy.AddFluidDofs( model_part )
*endif
*else
        import structural_solver_advanced
        structural_solver_advanced.AddDofs( model_part )
*endif

    def SetModelPart(self, model_part):
        self.model_part = model_part
        number_of_time_steps = *ntimesteps

        ## generating solver
*if(strcmp(GenData(Enable_Mortar_Contact),"1")==0)
        import mortar_gpts_contact_strategy
*if(strcmp(GenData(analysis_type),"static")==0)
        self.solver = mortar_gpts_contact_strategy.SampleSolver(self.model_part, self.abs_tol, self.rel_tol, self.analysis_parameters)
*else
        self.solver = mortar_gpts_contact_strategy.SampleSolverEkateQuasiStatic(self.model_part, self.abs_tol, self.rel_tol, self.analysis_parameters)
*endif
*else
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
*endif
*if(strcmp(GenData(Calculate_Reactions),"1")==0)
        (self.solver).CalculateReactionFlag = True
*else
        (self.solver).CalculateReactionFlag = False
*endif
        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
*if(strcmp(GenData(Solver),"BiCGStabLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        #preconditioner = ILU0Preconditioner()
        preconditioner = Preconditioner()
        plinear_solver = BICGSTABSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"CGLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        #preconditioner = ILU0Preconditioner()
        preconditioner = Preconditioner()
        #plinear_solver = DeflatedCGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner,1)
        plinear_solver = CGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"SuperLU")==0)
        plinear_solver = SuperLUSolver()
*endif
*if(strcmp(GenData(Solver),"SkylineLUFactorizationSolver")==0)
        plinear_solver = SkylineLUFactorizationSolver()
*endif
*if(strcmp(GenData(Solver),"GMRESSolver")==0)
        #preconditioner = Preconditioner()
        #plinear_solver = GMRESSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
        plinear_solver = MKLGMRESSolver()
*endif
*if(strcmp(GenData(Solver),"Pardiso")==0)
        plinear_solver = MKLPardisoSolver()
*endif
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)

    def SetOutputPath(self, path):
        self.results_path = path
*if(strcmp(GenData(Layer_Application),"1")==0)
        self.gid_io = SDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )
*else
        self.gid_io = StructuralGidIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )
*endif

    def FixPressureNodes( self, free_node_list_water, free_node_list_air):
        for node in self.model_part.Nodes:
            if (node.IsFixed(WATER_PRESSURE)==0):
                node.Fix(WATER_PRESSURE)
                free_node_list_water.append(node)
            if (node.IsFixed(AIR_PRESSURE)==0):
                node.Fix(AIR_PRESSURE)
                free_node_list_air.append(node)

    def ApplyInsituWaterPressure( self, free_node_list_water, free_node_list_air, z_zero, gravity_z):
        water_density=1000.0;
        for node in self.model_part.Nodes:
            water_pressure= water_density**gravity_z**(z_zero-(node.Z-node.GetSolutionStepValue(DISPLACEMENT_Z,0)))
            if( water_pressure < 1.0 ):
                water_pressure = 1.0
            node.SetSolutionStepValue(WATER_PRESSURE, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_EINS, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_NULL, water_pressure)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

    def ApplyInsituWaterPressure2D( self, free_node_list_water, free_node_list_air, y_zero, gravity_y):
        water_density=1000.0;
        for node in self.model_part.Nodes:
            water_pressure= water_density**gravity_y**(y_zero-(node.Y-node.GetSolutionStepValue(DISPLACEMENT_Y,0)))
            if( water_pressure < 1.0 ):
                water_pressure = 1.0
            node.SetSolutionStepValue(WATER_PRESSURE, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_EINS, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_NULL, water_pressure)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

    def SetReferenceWaterPressure( self ):
        self.SetReferenceWaterPressureForElements(self.model_part.Elements)

    def SetReferenceWaterPressureForElements( self, elements ):
        for element in elements:
            self.SetReferenceWaterPressureForElement(element)

    def SetReferenceWaterPressureForElement( self, element ):
        water_pressures = element.GetValuesOnIntegrationPoints( WATER_PRESSURE, self.model_part.ProcessInfo )
        pressure_list = []
        for item in water_pressures:
            pressure_list.append( item[0] )
        element.SetValuesOnIntegrationPoints( REFERENCE_WATER_PRESSURE, pressure_list, self.model_part.ProcessInfo )

    def FreePressureNodes(self,free_node_list_water, free_node_list_air):
        for item in free_node_list_water:
            #self.model_part.Nodes[item].Free(WATER_PRESSURE)
            item.Free(WATER_PRESSURE)
        for item in free_node_list_air:
            #self.model_part.Nodes[item].Free(AIR_PRESSURE)
            item.Free(AIR_PRESSURE)

    def WriteMaterialParameters( self, time, indices ):
        self.gid_io.OpenResultFile( self.results_path+self.problem_name, GiDPostMode.GiD_PostBinary)
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        #self.gid_io.ChangeOutputName( self.results_path+self.problem_name +str(time), GiDPostMode.GiD_PostBinary )
*endif
        for index in indices:
            self.gid_io.SuperPrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, index)
        self.gid_io.CloseResultFile()

    def WriteMonitoringSectionResults( self, time ):
        outfile = open("step_"+str(time)+".dat",'w')
        outfile.write("ekate result file for step "+str(time)+"\n")
        outfile.close()

    def WriteOutput( self, time ):
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.gid_io.InitializeMesh( time )
        mesh = self.model_part.GetMesh()
        print("mesh at time " + str(time) + " is ready for printing")
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
*else
        if( self.meshWritten == False ):
            self.gid_io.InitializeMesh( 0.0 )
            mesh = self.model_part.GetMesh()
            self.gid_io.WriteMesh( mesh )
            self.meshWritten = True
            self.gid_io.FinalizeMesh()
*endif
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
*else
        self.gid_io.InitializeResults( 0.0, self.model_part.GetMesh() )
*endif
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 0)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 1)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 2)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 3)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 4)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 5)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 6)
*if(strcmp(GenData(Displacements),"1")==0)
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, time, 0)
        print("nodal DISPLACEMENT written")
*endif
*if(strcmp(GenData(Reactions),"1")==0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part.Nodes, time, 0)
        print("nodal REACTION written")
*endif
*if(strcmp(GenData(PK2_Stresses),"1")==0)
        self.gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR, self.model_part, time)
*endif
*if(strcmp(GenData(Green_Lagrange_Strains),"1")==0)
        self.gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, self.model_part, time)
*endif
*if(strcmp(GenData(Stresses),"1")==0)
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        print("gauss point STRESSES written")
*endif
*if(strcmp(GenData(Jack_Forces),"1")==0)
        self.gid_io.PrintOnGaussPoints(JACK_FORCE, self.model_part, time)
*endif
*if(strcmp(GenData(Insitu_Stress),"1")==0)
        self.gid_io.PrintOnGaussPoints(PRESTRESS, self.model_part, time)
        print("gauss point PRESTRESS written")
*endif
*if(strcmp(GenData(Plastic_Strains),"1")==0)
        self.gid_io.PrintOnGaussPoints(PLASTIC_STRAIN_VECTOR, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(PLASTICITY_INDICATOR, self.model_part, time)
        print("gauss point PLASTICITY_INDICATOR written")
*endif
*if(strcmp(GenData(Internal_Variables),"1")==0)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 0)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 1)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 2)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 3)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 4)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 5)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 6)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 7)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 8)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 9)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 10)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 11)
        self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, 12)

*endif
*if(strcmp(GenData(Air_Pressure),"1")==0)
        self.gid_io.PrintOnGaussPoints(AIR_PRESSURE, self.model_part, time)
*endif
*if(strcmp(GenData(Water_Pressure),"1")==0)
        #self.gid_io.PrintOnGaussPoints(WATER_PRESSURE, self.model_part, time)
        #print("gauss point WATER_PRESSURE written")
        self.gid_io.WriteNodalResults(WATER_PRESSURE, self.model_part.Nodes, time, 0)
        print("nodal WATER_PRESSURE written")
        self.gid_io.PrintOnGaussPoints(EXCESS_PORE_WATER_PRESSURE, self.model_part, time)
        print("gauss point EXCESS_PORE_WATER_PRESSURE written")
*endif
*if(strcmp(GenData(Saturation),"1")==0)
        self.gid_io.PrintOnGaussPoints(SATURATION, self.model_part, time)
        print("gauss point SATURATION written")
*endif
*if(strcmp(GenData(Face_Load),"1")==0)
        self.gid_io.WriteNodalResults(FACE_LOAD, self.model_part.Nodes, time, 0)
        print("nodal FACE_LOAD written")
*endif
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
        self.gid_io.PrintOnGaussPoints(CONTACT_PENETRATION,self.model_part, time, 0)
        self.gid_io.PrintOnGaussPoints(NORMAL_CONTACT_STRESS,self.model_part, time, 0)
*endif
*if(strcmp(GenData(Topology_Optimization),"1")==0)
        self.gid_io.PrintOnGaussPoints(MATERIAL_DENSITY, self.model_part, time)
        print("gauss point MATERIAL_DENSITY written")
*endif
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.gid_io.FinalizeResults()
        self.gid_io.Reset()
*endif
*if(strcmp(GenData(VTK_Output),"1")==0)
*if(strcmp(GenData(VTK_Output_Type),"VTU")==0)
        print("write results to vtu")
*elseif(strcmp(GenData(VTK_Output_Type),"VTM")==0)
        print("write results to vtm")
*endif
        self.vtk_io.Initialize(time, mesh)
*if(strcmp(GenData(Displacements),"1")==0)
        print("write nodal displacements to vtu")
        self.vtk_io.RegisterNodalResults(DISPLACEMENT, 0)
*endif
*if(strcmp(GenData(Reactions),"1")==0)
        print("write nodal reactions to vtu")
        self.vtk_io.RegisterNodalResults(REACTION, 0)
*endif
*if(strcmp(GenData(Water_Pressure),"1")==0)
        print("write nodal water pressure to vtu")
        self.vtk_io.RegisterNodalResults(WATER_PRESSURE, 0)
*endif
*if(strcmp(GenData(Air_Pressure),"1")==0)
        print("write nodal air pressure to vtu")
        self.vtk_io.RegisterNodalResults(AIR_PRESSURE, 0)
*endif
        self.vtk_io.Finalize()
*endif

    def InitializeModel( self ):
        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        model_layers = __import__(self.problem_name+"_layers")
        ## ELEMENTS on layers ############################################
        self.layer_sets = model_layers.ReadLayerSets()
        ## NODES on layers ###############################################
        self.layer_nodes_sets = model_layers.ReadLayerNodesSets()
        ## CONTACT MASTER NODES ##########################################
        #self.contact_master_nodes = model_layers.ReadContactMasterNodes()
        ## CONTACT SLAVE NODES ###########################################
        #self.contact_slave_nodes = model_layers.ReadContactSlaveNodes()
        ##################################################################
        print("layer sets stored")
        ##################################################################
        ## STORE NODES ON GROUND SURFACE #################################
        ##################################################################
        self.top_surface_nodes = model_layers.ReadTopSurfaceNodes()
        print("nodes on ground surface stored")
        ##################################################################
        ## STORE NODES ON SIDE ###########################################
        ##################################################################
        self.boundary_nodes = model_layers.ReadBoundaryNodes()
        print("nodes on side surface stored")
        ##################################################################
        ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        ##################################################################
        self.node_groups = model_layers.ReadNodeGroups()
        print("node groups stored")
        ##################################################################
        ## EXTRACT CONDITIONS FROM NODE GROUPS ###########################
        ##################################################################
        start_time = time_module.time()
        self.layer_cond_sets = {}
        for layer, node_group in self.node_groups.iteritems():
            self.layer_cond_sets[layer] = []
        for layer, node_group in self.node_groups.iteritems():
            for cond in self.model_part.Conditions:
                in_group = True
                for node in cond.GetNodes():
                    if node.Id not in node_group:
                        in_group = False
                        break
                if in_group:
                    self.layer_cond_sets[layer].append(cond.Id)
        end_time = time_module.time()
        print("conditions in node groups stored, time = " + str(end_time - start_time) + "s")
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
*loop materials
*if(strcmp(MatProp(ConstitutiveLaw),"Isotropic3D")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(THICKNESS, 1.0 )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
        print("Linear elastic model selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"PlaneStrain")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(THICKNESS, *MatProp(Thickness,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
        print("Linear elastic model selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"PlaneStress")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(THICKNESS, *MatProp(Thickness,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, PlaneStress() )
        print("Linear elastic model selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"TutorialDamageModel")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(DAMAGE_E0, *MatProp(E0,real) )
        self.model_part.Properties[*MatNum].SetValue(DAMAGE_EF, *MatProp(Ef,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, TutorialDamageModel() )
        print("Tutorial damage model selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"GroutingMortar")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(PRIMARY_HYDRATION_TIME, *MatProp(prim_hyd_time,real) )
        self.model_part.Properties[*MatNum].SetValue(PRIMARY_HYDRATION_TIME_GRADIENT, *MatProp(gradient_prim_hyd_time,real) )
        self.model_part.Properties[*MatNum].SetValue(STIFFNESS_RATIO, *MatProp(E_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, GroutingMortar() )
        print("Grouting Mortar material selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"DruckerPrager")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(COHESION, *MatProp(Cohesion,real) )
        self.model_part.Properties[*MatNum].SetValue(INTERNAL_FRICTION_ANGLE, *MatProp(Friction_angle,real) )
        self.model_part.Properties[*MatNum].SetValue(ISOTROPIC_HARDENING_MODULUS, *MatProp(Isotropic_hardening_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, DruckerPrager() )
        print("Drucker Prager material selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"IsotropicDamage3D")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Compressive_Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONCRETE_YOUNG_MODULUS_C, *MatProp(Compressive_Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONCRETE_YOUNG_MODULUS_T, *MatProp(Tensile_Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(FC, *MatProp(Compressive_strength,real) )
        self.model_part.Properties[*MatNum].SetValue(FT, *MatProp(Tensile_strength,real) )
        self.model_part.Properties[*MatNum].SetValue(YIELD_STRESS, *MatProp(Yield_stress,real) )
        self.model_part.Properties[*MatNum].SetValue(FRACTURE_ENERGY, *MatProp(Fracture_energy,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        print("Isotropic Damage material selected for *MatProp(0), description: *MatProp(Description)")
*if(strcmp(MatProp(HardeningLaw),"Linear")==0)
        HardeningLaw_*MatNum = LinearSoftening()
*elseif(strcmp(MatProp(HardeningLaw),"Exponential")==0)
        HardeningLaw_*MatNum = ExponentialSoftening()
*endif
*if(strcmp(MatProp(FlowRule),"VonMises")==0)
        self.model_part.Properties[*MatNum].SetValue(CRUSHING_ENERGY, *MatProp(Crushing_energy,real) )
        FlowRule_*MatNum = VonMisesYieldFunction( State.Tri_D, HardeningLaw_*MatNum )
*elseif(strcmp(MatProp(FlowRule),"Energy")==0)
        FlowRule_*MatNum = EnergyYieldFunction( State.Tri_D )
*elseif(strcmp(MatProp(FlowRule),"MohrCoulomb")==0)
        self.model_part.Properties[*MatNum].SetValue(CRUSHING_ENERGY, *MatProp(Crushing_energy,real) )
        self.model_part.Properties[*MatNum].SetValue(INTERNAL_FRICTION_ANGLE, *MatProp(Friction_angle,real) )
        self.model_part.Properties[*MatNum].SetValue(COHESION, *MatProp(Cohesion,real) )
        self.model_part.Properties[*MatNum].SetValue(DILATANCY_ANGLE, *MatProp(Dilatancy_angle,real) )
        FlowRule_*MatNum = MorhCoulombYieldFunction( HardeningLaw_*MatNum, HardeningLaw_*MatNum, HardeningLaw_*MatNum, State.Tri_D, myPotencialPlastic.Associated )
*elseif(strcmp(MatProp(FlowRule),"MohrCoulombNonAssociated")==0)
        self.model_part.Properties[*MatNum].SetValue(CRUSHING_ENERGY, *MatProp(Crushing_energy,real) )
        self.model_part.Properties[*MatNum].SetValue(INTERNAL_FRICTION_ANGLE, *MatProp(Friction_angle,real) )
        self.model_part.Properties[*MatNum].SetValue(COHESION, *MatProp(Cohesion,real) )
        self.model_part.Properties[*MatNum].SetValue(DILATANCY_ANGLE, *MatProp(Dilatancy_angle,real) )
        MohrCoulomb_*MatNum = MorhCoulombYieldFunction( HardeningLaw_*MatNum, HardeningLaw_*MatNum, HardeningLaw_*MatNum, State.Tri_D, myPotencialPlastic.Not_Associated )
        Rankine_*MatNum = IsotropicRankineYieldFunction( HardeningLaw_*MatNum, State.Tri_D )
        FlowRule_*MatNum = ModifiedMorhCoulombYieldFunction( State.Tri_D, MohrCoulomb_*MatNum, Rankine_*MatNum )
*elseif(strcmp(MatProp(FlowRule),"Rankine")==0)
        FlowRule_*MatNum = IsotropicRankineYieldFunction( HardeningLaw_*MatNum, State.Tri_D )
*endif
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, IsotropicDamage3D( FlowRule_*MatNum , HardeningLaw_*MatNum , self.model_part.Properties[*MatNum] ) )
*elseif(strcmp(MatProp(ConstitutiveLaw),"UserDefined")==0)
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, DummyConstitutiveLaw() )
        print("User-defined material selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"TrussMaterial")==0)
        print("Truss material selected for *MatProp(0), description: *MatProp(Description)")
*elseif(strcmp(MatProp(ConstitutiveLaw),"BeamMaterial")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        print("Beam material selected for *MatProp(0), description: *MatProp(Description)")
*else
        print("Material *MatProp(0) *MatProp(ConstitutiveLaw) *MatProp(Description)")
*endif
*if(strcmp(GenData(Topology_Optimization),"1")==0)
*set var penalizationfactor(real)=GenData(Penalization_Factor,real)
*set var minimummodulus(real)=GenData(Minimum_Modulus,real)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS_MIN, *minimummodulus )
        self.model_part.Properties[*MatNum].SetValue(PENALIZATION_FACTOR, *penalizationfactor )
*endif
*end materials
*if(strcmp(GenData(Enable_Mortar_Contact),"1")==0||strcmp(GenData(Mortar_Application),"1")==0)
        ##################################################################
        ## MORTAR TYING/CONTACT ##########################################
        ##################################################################
        ## set the correct master and slave index set
        ## in addition, add additional condition set
        for cond in self.model_part.Conditions:
            if cond.Has(MASTER_INDEX):
                master_index_set = IntegerVector(1)
                master_index_set[0] = cond.GetValue(MASTER_INDEX)
                cond.SetValue(MASTER_INDEX_SET, master_index_set)
                layer = "mortar_master_" + str(master_index_set[0])
                if layer not in self.layer_cond_sets:
                    self.layer_cond_sets[layer] = []
                self.layer_cond_sets[layer].append(cond.Id)
            if cond.Has(SLAVE_INDEX):
                slave_index_set = IntegerVector(1)
                slave_index_set[0] = cond.GetValue(SLAVE_INDEX)
                cond.SetValue(SLAVE_INDEX_SET, slave_index_set)
                layer = "mortar_slave_" + str(slave_index_set[0])
                if layer not in self.layer_cond_sets:
                    self.layer_cond_sets[layer] = []
                self.layer_cond_sets[layer].append(cond.Id)
        print("Mortar slave/master index set are assigned")
        print("Mortar slave/master condition set are created")
*endif
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
        self.solver.solver.Initialize()
        self.model_part.Check( self.model_part.ProcessInfo )
        print("activation utility initialized")
        print("model successfully initialized")

    def WriteRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time)
        serializer = Serializer(fn)
        serializer.Save("ModelPart", self.model_part)
        serializer = 0
        print("Write restart data to " + fn + ".rest completed")

    def LoadRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time)
        serializer = Serializer(fn)
        serializer.Load("ModelPart", self.model_part)
        serializer = 0
        print("Load restart data from " + fn + ".rest completed")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()

    # solve with deactivation/reactivation
    # element/condition with nonzero ACTIVATION_LEVEL in [from_deac, to_deac] will be deactivated
    # element/condition with negative ACTIVATION_LEVEL will also be deactivated
    # element/condition with ACTIVATION_LEVEL in [from_reac, to_reac] will be re-activated
    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()

    # solve nothing (good for debugging) with deactivation/reactivation
    def DrySolve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)

    # solve without deactivation
    def SolveModel(self, time):
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()

    # solve nothing without deactivation
    def DrySolveModel(self, time):
        self.model_part.CloneTimeStep(time)

##################################################################
