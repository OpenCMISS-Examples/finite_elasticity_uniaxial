!Main program
PROGRAM UniaxialExtension

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif
    
  !Test program parameters
    
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP
  
  REAL(CMISSRP), PARAMETER :: C1=0.5_CMISSRP
  REAL(CMISSRP), PARAMETER :: C2=0.1_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MATERIAL_FIELD_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: DERIVED_FIELD_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program types
 
  !Program variables
  INTEGER(CMISSIntg) :: argumentLength,componentIdx,computationalNodeNumber,decompositionIndex,dimensionIdx1,dimensionIdx2, &
    & equationsSetIndex,err,interpolationType,nodeDomain,nodeNumber,numberOfArguments,numberOfComputationalNodes, &
    & numberOfDimensions,numberOfGaussXi,numberOfGLobalXNodes,numberOfGlobalYNodes,numberOfGlobalZNodes,numberOfNodes, &
    & numberOfNodesPerElement,numberOfTensorComponents,numberOfXi,status,xNodeIdx,yNodeIdx,zNodeIdx
  REAL(CMISSRP) :: alpha,analAlpha,errorAlpha,beta,analBeta,errorBeta,gamma,analGamma,errorGamma,J,analJ,errorJ,p,analP,errorP, &
    & F(3,3),analF(3,3),errorF(3,3),C(3,3),analC(3,3),errorC(3,3),E(3,3),analE(3,3),errorE(3,3),S(3,3),analS(3,3),errorS(3,3), &
    & sigma(3,3),analSigma(3,3),errorSigma(3,3),initialBeta,initialP,xi(3),yValue,zValue
  LOGICAL :: directoryExists=.FALSE.
  CHARACTER(LEN=255) :: commandArgument
  
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: derivedField,dependentField,equationsSetField,geometricField,materialsField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver,linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup
 
  !Setup problem 
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 1) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberOfDimensions
    IF(numberOfDimensions<2.OR.numberOfDimensions>3) CALL HandleError("Invalid number of dimensions.")
    IF(numberOfArguments>=2) THEN
      CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
      IF(status>0) CALL HandleError("Error for command argument 2.")
      READ(commandArgument(1:argumentLength),*) alpha
      IF(alpha<0.000001_CMISSRP.OR.alpha>0.999999_CMISSRP) CALL HandleError("Invalid alpha specification.")
      IF(numberOfArguments>=3) THEN
        CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
        IF(status>0) CALL HandleError("Error for command argument 3.")
        READ(commandArgument(1:argumentLength),*) interpolationType
        IF(interpolationType<1.OR.interpolationType>9) CALL HandleError("Invalid interpolation specification.")
      ELSE
        interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
      ENDIF
    ELSE
      alpha=0.1_CMISSRP
      interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    ENDIF
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfDimensions=2
    alpha=0.1_CMISSRP
    interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  ENDIF

  SELECT CASE(interpolationType)
  CASE(1,4,7)
    numberOfNodesPerElement=2
  CASE(2,8)
    numberOfNodesPerElement=3
  CASE(3,9)
    numberOfNodesPerElement=4
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")    
  END SELECT
  
  SELECT CASE(numberOfDimensions)
  CASE(1)
    numberOfGlobalXNodes=numberOfNodesPerElement
    numberOfGlobalYNodes=1
    numberOfGlobalZNodes=1
    IF(interpolationType==3.OR.interpolationType==4) THEN
      numberOfGaussXi=3
    ELSE
      numberOfGaussXi=2
    ENDIF
    numberOfXi=1
    xi(1)=0.5_CMISSRP
    numberOfTensorComponents=1
    CALL HandleError("Not implemented.")    
  CASE(2)
    numberOfGlobalXNodes=numberOfNodesPerElement
    numberOfGlobalYNodes=numberOfNodesPerElement
    numberOfGlobalZNodes=1
    IF(interpolationType==3.OR.interpolationType==4) THEN
      numberOfGaussXi=3
    ELSE
      numberOfGaussXi=2
    ENDIF
    numberOfXi=2
    xi(1)=0.5_CMISSRP
    xi(2)=0.5_CMISSRP
    numberOfTensorComponents=3
  CASE(3)
    numberOfGlobalXNodes=numberOfNodesPerElement
    numberOfGlobalYNodes=numberOfNodesPerElement
    numberOfGlobalZNodes=numberOfNodesPerElement
    numberOfGaussXi=3
    numberOfXi=3
    xi(1)=0.5_CMISSRP
    xi(2)=0.5_CMISSRP
    xi(3)=0.5_CMISSRP
    numberOfTensorComponents=6
  CASE DEFAULT
    CALL HandleError("Invalid number of dimensions.")    
  END SELECT
  numberOfNodes=numberOfGlobalXNodes*numberOfGlobalYNodes*numberOfGlobalZNodes
  
  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"],err)
  !Set output
  CALL cmfe_OutputSetOn("Uniaxial",err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)

  
  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)  

  !Create a rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"Region",err)
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL cmfe_Region_CreateFinish(region,err)

  !Define geometric basis
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  IF(numberOfDimensions==2) THEN
    CALL cmfe_Basis_NumberOfXiSet(basis,2,err)
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    CALL cmfe_Basis_NumberOfXiSet(basis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF
  CALL cmfe_Basis_CreateFinish(basis,err)
  
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  IF(numberOfDimensions==2) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[1,1],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[1,1,1],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL cmfe_Field_ScalingTypeSet(geometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
  CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(dependentField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(dependentField,geometricField,err)
  CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,err)
  CALL cmfe_Field_NumberOfVariablesSet(dependentField,2,err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfDimensions+1,err)
  CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,numberOfDimensions+1,err)
  !Set the pressure to be a constant element
  CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,numberOfDimensions+1, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  CALL cmfe_Field_CreateFinish(dependentField,err)
  
  !Create the material field
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_Field_CreateStart(MATERIAL_FIELD_USER_NUMBER,region,materialsField,err)
  CALL cmfe_Field_TypeSet(materialsField,CMFE_FIELD_MATERIAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(materialsField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(materialsField,geometricField,err)
  CALL cmfe_Field_NumberOfVariablesSet(materialsField,1,err)
  CALL cmfe_Field_VariableLabelSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL cmfe_Field_NumberOfComponentsSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,2,err)
  CALL cmfe_Field_CreateFinish(materialsField,err)

  !Create the derived fields
  CALL cmfe_Field_Initialise(derivedField,err)
  CALL cmfe_Field_CreateStart(DERIVED_FIELD_USER_NUMBER,region,derivedField,err)
  CALL cmfe_Field_TypeSet(derivedField,CMFE_FIELD_GENERAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(derivedField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(derivedField,geometricField,err)
  CALL cmfe_Field_DependentTypeSet(derivedField,CMFE_FIELD_DEPENDENT_TYPE,err)
  CALL cmfe_Field_NumberOfVariablesSet(derivedField,5,err)
  CALL cmfe_Field_VariableTypesSet(derivedField,[CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_U2_VARIABLE_TYPE, &
    & CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_U5_VARIABLE_TYPE],err)
  CALL cmfe_Field_VariableLabelSet(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,"DeformationGradient",err)
  CALL cmfe_Field_VariableLabelSet(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,"RightCauchyGreen",err)
  CALL cmfe_Field_VariableLabelSet(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,"GreenLagrange",err)
  CALL cmfe_Field_VariableLabelSet(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,"SecondPiolaKirchoff",err)
  CALL cmfe_Field_VariableLabelSet(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,"Cauchy",err)
  CALL cmfe_Field_NumberOfComponentsSet(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL cmfe_Field_NumberOfComponentsSet(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL cmfe_Field_NumberOfComponentsSet(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL cmfe_Field_NumberOfComponentsSet(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL cmfe_Field_NumberOfComponentsSet(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,numberOfTensorComponents,err)
  DO componentIdx=1,numberOfTensorComponents
    CALL cmfe_Field_ComponentInterpolationSet(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,componentIdx, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,componentIdx, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,componentIdx, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,componentIdx, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,componentIdx, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  ENDDO !componentIdx
  CALL cmfe_Field_CreateFinish(derivedField,err)
  
  !Create the equations_set
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the equations set material field 
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIAL_FIELD_USER_NUMBER,materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01.
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C1,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C2,err)

  !Create the derived fields
  CALL cmfe_EquationsSet_DerivedCreateStart(equationsSet,DERIVED_FIELD_USER_NUMBER,derivedField,err)
  CALL cmfe_EquationsSet_DerivedVariableSet(equationsSet,CMFE_EQUATIONS_SET_DERIVED_DEFORMATION_GRADIENT, &
    & CMFE_FIELD_U1_VARIABLE_TYPE,err)
  CALL cmfe_EquationsSet_DerivedVariableSet(equationsSet,CMFE_EQUATIONS_SET_DERIVED_R_CAUCHY_GREEN_DEFORMATION, &
    & CMFE_FIELD_U2_VARIABLE_TYPE,err)
  CALL cmfe_EquationsSet_DerivedVariableSet(equationsSet,CMFE_EQUATIONS_SET_DERIVED_GREEN_LAGRANGE_STRAIN, &
    & CMFE_FIELD_U3_VARIABLE_TYPE,err)
  CALL cmfe_EquationsSet_DerivedVariableSet(equationsSet,CMFE_EQUATIONS_SET_DERIVED_SECOND_PK_STRESS, &
    & CMFE_FIELD_U4_VARIABLE_TYPE,err)
  CALL cmfe_EquationsSet_DerivedVariableSet(equationsSet,CMFE_EQUATIONS_SET_DERIVED_CAUCHY_STRESS, &
    & CMFE_FIELD_U5_VARIABLE_TYPE,err)
  CALL cmfe_EquationsSet_DerivedCreateFinish(equationsSet,err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and displacement bcs
  DO dimensionIdx1=1,numberOfDimensions
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & dimensionIdx1,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,dimensionIdx1,err)
  ENDDO !dimensionIdx1

  !Define the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  CALL cmfe_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  
  !Set x=0 nodes to no displacement in the x-direction
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO yNodeIdx=1,numberOfGlobalYNodes
      nodeNumber=1+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
      ENDIF
    ENDDO !yNodeIdx
  ENDDO !zNodeIdx
  !Set y=0 nodes to no displacement in the y-direction
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO xNodeIdx=1,numberOfGlobalXNodes
      nodeNumber=xNodeIdx+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !y-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
      ENDIF
    ENDDO !yNodeIdx
  ENDDO !zNodeIdx
  IF(numberOfDimensions==3) THEN
    !Set z=0 nodes to no displacement in the z-direction
    DO yNodeIdx=1,numberOfGlobalYNodes
      DO xNodeIdx=1,numberOfGlobalXNodes
        nodeNumber=xNodeIdx+(yNodeIdx-1)*numberOfGlobalXNodes
        CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !z-direction
          CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        ENDIF
      ENDDO !yNodeIdx
    ENDDO !zNodeIdx
  ENDIF
  !Fix the x=LENGTH nodes to alpha% x-displacement
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO yNodeIdx=1,numberOfGlobalYNodes
      nodeNumber=numberOfGlobalXNodes+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMFE_BOUNDARY_CONDITION_FIXED,alpha,err)
      ENDIF
    ENDDO !yNodeIdx
  ENDDO !zNodeIdx
   
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Set initial values for the dependent field
  SELECT CASE(numberOfDimensions)
  CASE(2)
    initialBeta=0.09090957777_CMISSRP
    initialP  =-0.19177686982_CMISSRP
    !Set initial values for the top
    DO zNodeIdx=1,numberOfGlobalZNodes
      DO xNodeIdx=1,numberOfGlobalXNodes
        nodeNumber=xNodeIdx+(numberOfGlobalYNodes-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
        CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !y-direction
          CALL cmfe_Field_ParameterSetUpdateNode(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
            & nodeNumber,2,1.0_CMISSRP-initialBeta,err)
        ENDIF
      ENDDO !yNodeIdx
    ENDDO !zNodeIdx
    CALL cmfe_Field_ParameterSetUpdateElement(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,numberOfDimensions+1,initialP,err)
  CASE(3)
  CASE DEFAULT
    CALL HandleError("Invalid number of dimensions.")
  END SELECT

  !Solve problem
  CALL cmfe_Problem_Solve(problem,err)

  !Calculate the derived fields
  CALL cmfe_EquationsSet_DerivedVariableCalculate(equationsSet,CMFE_EQUATIONS_SET_DERIVED_DEFORMATION_GRADIENT,err)
  CALL cmfe_EquationsSet_DerivedVariableCalculate(equationsSet,CMFE_EQUATIONS_SET_DERIVED_R_CAUCHY_GREEN_DEFORMATION,err)
  CALL cmfe_EquationsSet_DerivedVariableCalculate(equationsSet,CMFE_EQUATIONS_SET_DERIVED_GREEN_LAGRANGE_STRAIN,err)
  CALL cmfe_EquationsSet_DerivedVariableCalculate(equationsSet,CMFE_EQUATIONS_SET_DERIVED_SECOND_PK_STRESS,err)
  CALL cmfe_EquationsSet_DerivedVariableCalculate(equationsSet,CMFE_EQUATIONS_SET_DERIVED_CAUCHY_STRESS,err)
  !Construct tensors
  analAlpha=alpha
  CALL ComputeError(alpha,analAlpha,errorAlpha)
  CALL cmfe_Field_ParameterSetGetNode(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,numberOfNodes,2, &
    & yValue,err)
  beta=1-yValue
  analBeta=beta
  CALL ComputeError(beta,analBeta,errorBeta)
  IF(numberOfDimensions==3) THEN
    CALL cmfe_Field_ParameterSetGetNode(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,numberOfNodes,3, &
      & zValue,err)
    gamma=1-zValue
    analGamma=gamma
    CALL ComputeError(gamma,analGamma,errorGamma)
  ELSE
    zValue=0.0_CMISSRP
    gamma=0.0_CMISSRP
    errorGamma=0.0_CMISSRP
  ENDIF
  CALL cmfe_Field_ParameterSetGetElement(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,numberOfDimensions+1,p,err)
  analP=p
  CALL ComputeError(p,analP,errorP)
  F=0.0_CMISSRP
  analF=0.0_CMISSRP
  errorF=0.0_CMISSRP
  C=0.0_CMISSRP
  analC=0.0_CMISSRP
  errorC=0.0_CMISSRP
  E=0.0_CMISSRP
  analE=0.0_CMISSRP
  errorE=0.0_CMISSRP
  S=0.0_CMISSRP
  analS=0.0_CMISSRP
  errorS=0.0_CMISSRP
  sigma=0.0_CMISSRP
  analSigma=0.0_CMISSRP
  errorSigma=0.0_CMISSRP
  SELECT CASE(numberOfDimensions)
  CASE(2)
    J=(1.0_CMISSRP+alpha)*(1.0_CMISSRP-beta)
    analJ=(1.0_CMISSRP+analAlpha)*(1.0_CMISSRP-analBeta)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,F(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,F(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,F(1,2),err)
    F(2,1)=F(1,2)
    analF(1,1)=1.0_CMISSRP+analAlpha
    analF(2,2)=1.0_CMISSRP-analBeta
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,C(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,C(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,C(1,2),err)
    C(2,1)=C(1,2)
    analC(1,1)=(1.0_CMISSRP+analAlpha)/(1.0_CMISSRP-analBeta)
    analC(2,2)=(1.0_CMISSRP-analBeta)/(1.0_CMISSRP+analAlpha)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,E(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,E(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,E(1,2),err)
    E(2,1)=E(1,2)
    analE(1,1)=(analAlpha+analBeta)/(2.0_CMISSRP*(1.0_CMISSRP-analBeta))
    analE(2,2)=-1.0_CMISSRP*(analAlpha+analBeta)/(2.0_CMISSRP*(1.0_CMISSRP+analAlpha))
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,S(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,S(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,S(1,2),err)
    S(2,1)=S(1,2)
    analS(1,1)=2.0_CMISSRP*(C1+C2*(1.0_CMISSRP-analBeta)/(1.0_CMISSRP+analAlpha))
    analS(2,2)=2.0_CMISSRP*(C1+C2*(1.0_CMISSRP+analAlpha)/(1.0_CMISSRP-analBeta))
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,sigma(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,sigma(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,sigma(1,2),err)
    sigma(2,1)=sigma(1,2)
    analSigma(1,1)=analJ**(-2.5_CMISSRP)*C1*(analAlpha*(2.0_CMISSRP+analAlpha)-analBeta*(2.0_CMISSRP-analBeta))+analP
    analSigma(2,2)=-1.0*CMISSRP*analJ**(-2.5_CMISSRP)*C1*(analAlpha*(2.0_CMISSRP+analAlpha)-analBeta*(2.0_CMISSRP-analBeta))+analP
  CASE(3)
    J=(1.0_CMISSRP+alpha)*(1.0_CMISSRP-beta)*(1.0_CMISSRP-gamma)
    analJ=(1.0_CMISSRP+analAlpha)*(1.0_CMISSRP-analBeta)*(1.0_CMISSRP-analGamma)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,F(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,F(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,F(3,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,4,F(2,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,F(1,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,6,F(1,2),err)
    F(2,1)=F(1,2)
    F(3,1)=F(1,3)
    F(3,2)=F(2,3)
    analF(1,1)=1.0_CMISSRP+analAlpha
    analF(2,2)=1.0_CMISSRP-analBeta
    analF(3,3)=1.0_CMISSRP-analGamma
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,C(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,C(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,C(3,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,4,C(2,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,C(1,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,6,C(1,2),err)
    C(2,1)=C(1,2)
    C(3,1)=C(1,3)
    C(3,2)=C(2,3)
    analC(1,1)=(((1.0_CMISSRP+analAlpha)**4.0_CMISSRP)**(1.0_CMISSRP/3.0_CMISSRP))/ &
      & ((((1.0_CMISSRP-analBeta)**2.0_CMISSRP)*((1.0_CMISSRP-analGamma)**2.0_CMISSRP))**(1.0_CMISSRP/3.0_CMISSRP))
    analC(2,2)=(((1.0_CMISSRP-analBeta)**4.0_CMISSRP)**(1.0_CMISSRP/3.0_CMISSRP))/ &
      & ((((1.0_CMISSRP+analAlpha)**2.0_CMISSRP)*((1.0_CMISSRP-analGamma)**2.0_CMISSRP))**(1.0_CMISSRP/3.0_CMISSRP))
    analC(3,3)=(((1.0_CMISSRP-analGamma)**4.0_CMISSRP)**(1.0_CMISSRP/3.0_CMISSRP))/ &
      & ((((1.0_CMISSRP+analAlpha)**2.0_CMISSRP)*((1.0_CMISSRP-analBeta)**2.0_CMISSRP))**(1.0_CMISSRP/3.0_CMISSRP))
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,E(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,E(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,E(3,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,4,E(2,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,E(1,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,6,E(1,2),err)
    E(2,1)=E(1,2)
    E(3,1)=E(1,3)
    E(3,2)=E(2,3)
    analE(1,1)=0.5_CMISSRP*(analC(1,1)-1.0_CMISSRP)
    analE(2,2)=0.5_CMISSRP*(analC(2,2)-1.0_CMISSRP)
    analE(3,3)=0.5_CMISSRP*(analC(3,3)-1.0_CMISSRP)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,S(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,S(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,S(3,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,4,S(2,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,S(1,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U4_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,6,S(1,2),err)
    S(2,1)=S(1,2)
    S(3,1)=S(1,3)
    S(3,2)=S(2,3)
    analS(1,1)=2.0_CMISSRP*(C1+(C2*((1.0_CMISSRP-analBeta)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)**2.0_CMISSRP))/ &
      & (((1.0_CMISSRP+analAlpha)**2.0_CMISSRP*(1.0_CMISSRP-analBeta)**2.0_CMISSRP*(1.0_CMISSRP-analGamma)**2.0_CMISSRP)** &
      & (1.0_CMISSRP/3.0_CMISSRP)))
    analS(2,2)=2.0_CMISSRP*(C1+(C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)**2.0_CMISSRP))/ &
      & (((1.0_CMISSRP+analAlpha)**2.0_CMISSRP*(1.0_CMISSRP-analBeta)**2.0_CMISSRP*(1.0_CMISSRP-analGamma)**2.0_CMISSRP)** &
      & (1.0_CMISSRP/3.0_CMISSRP)))
    analS(3,3)=2.0_CMISSRP*(C1+(C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analBeta)**2.0_CMISSRP))/ &
      & (((1.0_CMISSRP+analAlpha)**2.0_CMISSRP*(1.0_CMISSRP-analBeta)**2.0_CMISSRP*(1.0_CMISSRP-analGamma)**2.0_CMISSRP)** &
      & (1.0_CMISSRP/3.0_CMISSRP)))
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,sigma(1,1),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,sigma(2,2),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,sigma(3,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,4,sigma(2,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,sigma(1,3),err)
    CALL cmfe_Field_ParameterSetGetElement(derivedField,CMFE_FIELD_U5_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,6,sigma(1,2),err)
    sigma(2,1)=sigma(1,2)
    sigma(3,1)=sigma(1,3)
    sigma(3,2)=sigma(2,3)
    analSigma(1,1)=(2.0_CMISSRP/(3.0_CMISSRP*analJ**3.0_CMISSRP))*((2.0_CMISSRP*(1.0_CMISSRP+analAlpha)**2.0_CMISSRP* &
      & (C1+C2*((1.0_CMISSRP-analBeta)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)**2.0_CMISSRP))) &
      & -((1.0_CMISSRP-analBeta)**2.0_CMISSRP*(C1+C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)**2.0_CMISSRP))) &
      & -((1.0_CMISSRP-analGamma)**2.0_CMISSRP*(C1+C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analBeta)**2.0_CMISSRP))) &
      & )+analP
    analSigma(2,2)=(2.0_CMISSRP/(3.0_CMISSRP*analJ**3.0_CMISSRP))*((-1.0_CMISSRP*(1.0_CMISSRP+analAlpha)**2.0_CMISSRP* &
      & (C1+C2*((1.0_CMISSRP+analBeta)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)**2.0_CMISSRP))) &
      & +(2.0_CMISSRP*(1.0_CMISSRP-analBeta)**2.0_CMISSRP*(C1+C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)** &
      & 2.0_CMISSRP))) &
      & -((1.0_CMISSRP-analGamma)**2.0_CMISSRP*(C1+C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analBeta)**2.0_CMISSRP))) &
      & )+analP
    analSigma(3,3)=(2.0_CMISSRP/(3.0_CMISSRP*analJ**3.0_CMISSRP))*((-1.0_CMISSRP*(1.0_CMISSRP+analAlpha)**2.0_CMISSRP* &
      & (C1+C2*((1.0_CMISSRP+analBeta)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)**2.0_CMISSRP))) &
      & -((1.0_CMISSRP-analBeta)**2.0_CMISSRP*(C1+C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analGamma)** &
      & 2.0_CMISSRP))) &
      & +(2.0_CMISSRP*(1.0_CMISSRP-analGamma)**2.0_CMISSRP*(C1+C2*((1.0_CMISSRP+analAlpha)**2.0_CMISSRP+(1.0_CMISSRP-analBeta)** &
      & 2.0_CMISSRP))))+analP
  CASE DEFAULT
    CALL HandleError("Invalid number of dimensions.")    
  END SELECT
  CALL ComputeError(J,analJ,errorJ)
  DO dimensionIdx1=1,numberOfDimensions
    DO dimensionIdx2=1,numberOfDimensions
      CALL ComputeError(F(dimensionIdx1,dimensionIdx2),analF(dimensionIdx1,dimensionIdx2),errorF(dimensionIdx1,dimensionIdx2))
      CALL ComputeError(C(dimensionIdx1,dimensionIdx2),analC(dimensionIdx1,dimensionIdx2),errorC(dimensionIdx1,dimensionIdx2))
      CALL ComputeError(E(dimensionIdx1,dimensionIdx2),analE(dimensionIdx1,dimensionIdx2),errorE(dimensionIdx1,dimensionIdx2))
      CALL ComputeError(S(dimensionIdx1,dimensionIdx2),analS(dimensionIdx1,dimensionIdx2),errorS(dimensionIdx1,dimensionIdx2))
      CALL ComputeError(sigma(dimensionIdx1,dimensionIdx2),analSigma(dimensionIdx1,dimensionIdx2), &
        & errorSigma(dimensionIdx1,dimensionIdx2))
    ENDDO !dimensionIdx2
  ENDDO !dimensionIdx

  CALL WriteScalarError(1,numberOfDimensions,"alpha",alpha,analAlpha,errorAlpha)
  CALL WriteScalarError(2,numberOfDimensions,"beta",beta,analBeta,errorBeta)
  CALL WriteScalarError(2,numberOfDimensions,"gamma",gamma,analGamma,errorGamma)
  CALL WriteScalarError(2,numberOfDimensions,"p",p,analP,errorP)
  CALL WriteScalarError(2,numberOfDimensions,"J",J,analJ,errorJ)
  CALL WriteTensorError(1,numberOfDimensions,"F",F,analF,errorF)
  CALL WriteTensorError(2,numberOfDimensions,"C",C,analC,errorC)
  CALL WriteTensorError(2,numberOfDimensions,"E",E,analE,errorE)
  CALL WriteTensorError(2,numberOfDimensions,"S",S,analS,errorS)
  CALL WriteTensorError(2,numberOfDimensions,"sigma",sigma,analSigma,errorSigma)

  INQUIRE(FILE="./results",EXIST=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"./results/Uniaxial","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"./results/Uniaxial","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    
    STOP
  END SUBROUTINE HandleError

  SUBROUTINE ComputeError(actual,analytic,error)
    
    REAL(CMISSRP), INTENT(IN) :: actual
    REAL(CMISSRP), INTENT(IN) :: analytic
    REAL(CMISSRP), INTENT(OUT) :: error
    
    IF(ABS(analytic)<=5.0*CMISSRP*EPSILON(0.0_CMISSRP)) THEN
      error=(actual-analytic)
    ELSE
      error=(actual-analytic)*100.0_CMISSRP/analytic
    ENDIF
    
    RETURN
  END SUBROUTINE ComputeError

  SUBROUTINE WriteScalarError(outputType,numberOfDimensions,symbol,actual,analytic,error)
    
    INTEGER(CMISSIntg), INTENT(IN) :: outputType
    INTEGER(CMISSIntg), INTENT(IN) :: numberOfDimensions
    CHARACTER(LEN=*), INTENT(IN) :: symbol
    REAL(CMISSRP), INTENT(IN) :: actual
    REAL(CMISSRP), INTENT(IN) :: analytic
    REAL(CMISSRP), INTENT(IN) :: error

    CHARACTER(LEN=10) :: outputSymbol
    
    IF(outputType==1) THEN
      WRITE(*,*)
      SELECT CASE(numberOfDimensions)
      CASE(2)
        WRITE(*,'(33X,"Actual",20X,"Analytic",22X,"%Error")')
      CASE(3)
        WRITE(*,'(46X,"Actual",33X,"Analytic",35X,"%Error")')
      CASE DEFAULT
        CALL HandleError("Invalid number of dimensions.")
      END SELECT
    ENDIF
    WRITE(*,*)
    WRITE(outputSymbol,'(A)') ADJUSTL(symbol(1:LEN_TRIM(symbol)))
    SELECT CASE(numberOfDimensions)
    CASE(2)
      WRITE(*,'(A," = ",14X,E12.5,2X,14X,E12.5,2X,14X,E12.5)') outputSymbol,actual,analytic,error
    CASE(3)
      WRITE(*,'(A," = ",27X,E12.5,2X,27X,E12.5,2X,27X,E12.5)') outputSymbol,actual,analytic,error
    CASE DEFAULT
      CALL HandleError("Invalid number of dimensions.")
    END SELECT
    
    RETURN
  END SUBROUTINE WriteScalarError
  
  SUBROUTINE WriteTensorError(outputType,numberOfDimensions,symbol,actual,analytic,error)
    
    INTEGER(CMISSIntg), INTENT(IN) :: outputType
    INTEGER(CMISSIntg), INTENT(IN) :: numberOfDimensions
    CHARACTER(LEN=*), INTENT(IN) :: symbol
    REAL(CMISSRP), INTENT(IN) :: actual(:,:)
    REAL(CMISSRP), INTENT(IN) :: analytic(:,:)
    REAL(CMISSRP), INTENT(IN) :: error(:,:)
    
    INTEGER(CMISSIntg) :: dimensionIdx1,dimensionIdx2
    CHARACTER(LEN=5) :: outputSymbol

    IF(outputType==1) THEN
      WRITE(*,*)
      SELECT CASE(numberOfDimensions)
      CASE(2)
        WRITE(*,'(33X,"Actual",20X,"Analytic",22X,"%Error")')
      CASE(3)
        WRITE(*,'(46X,"Actual",33X,"Analytic",35X,"%Error")')
      CASE DEFAULT
        CALL HandleError("Invalid number of dimensions.")
      END SELECT
    ENDIF
    WRITE(*,*)
    WRITE(outputSymbol,'(A)') ADJUSTL(symbol(1:LEN_TRIM(symbol)))
    SELECT CASE(numberOfDimensions)
    CASE(2)
      DO dimensionIdx1=1,numberOfDimensions
        WRITE(*,'(A,"(",I1,",:) = ",2(X,E12.5),2X,2(X,E12.5),2X,2(X,E12.5))') outputSymbol,dimensionIdx1, &
          & (actual(dimensionIdx1,dimensionIdx2),dimensionIdx2=1,numberOfDimensions), &
          & (analytic(dimensionIdx1,dimensionIdx2),dimensionIdx2=1,numberOfDimensions), &
          & (error(dimensionIdx1,dimensionIdx2),dimensionIdx2=1,numberOfDimensions)
      ENDDO !dimensionIdx1
    CASE(3)
      DO dimensionIdx1=1,numberOfDimensions
        WRITE(*,'(A,"(",I1,",:) = "3(X,E12.5),2X,3(X,E12.5),2X,3(X,E12.5))') outputSymbol,dimensionIdx1, &
          & (actual(dimensionIdx1,dimensionIdx2),dimensionIdx2=1,numberOfDimensions), &
          & (analytic(dimensionIdx1,dimensionIdx2),dimensionIdx2=1,numberOfDimensions), &
          & (error(dimensionIdx1,dimensionIdx2),dimensionIdx2=1,numberOfDimensions)
      ENDDO !dimensionIdx1
    CASE DEFAULT
      CALL HandleError("Invalid number of dimensions.")
    END SELECT
    
    RETURN
  END SUBROUTINE WriteTensorError

END PROGRAM UniaxialExtension
