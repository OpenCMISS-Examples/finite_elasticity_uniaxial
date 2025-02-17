!Main program
PROGRAM UniaxialExtension

  USE OpenCMISS

  IMPLICIT NONE
    
  !Test program parameters
    
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP !x-extent
  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP !y-extent
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP  !z-extent
  
  REAL(OC_RP), PARAMETER :: ALPHA=0.1_OC_RP
  
  REAL(OC_RP), PARAMETER :: C1=0.5_OC_RP
  REAL(OC_RP), PARAMETER :: C2=0.1_OC_RP
  
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MATERIAL_FIELD_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: DERIVED_FIELD_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program types
 
  !Program variables
  INTEGER(OC_Intg) :: argumentLength,componentIdx,computationalNodeNumber,decompositionIndex,dimensionIdx1,dimensionIdx2, &
    & equationsSetIndex,err,interpolationType,nodeDomain,nodeNumber,numberOfArguments,numberOfComputationalNodes, &
    & numberOfDimensions,numberOfGaussXi,numberOfGLobalXNodes,numberOfGlobalYNodes,numberOfGlobalZNodes,numberOfNodes, &
    & numberOfNodesPerElement,numberOfTensorComponents,numberOfXi,status,xNodeIdx,yNodeIdx,zNodeIdx
  REAL(OC_RP) :: analAlpha,errorAlpha,analBeta,beta,errorBeta,analGamma,errorGamma,gamma,J,analJ,errorJ,p,analP,errorP, &
    & F(3,3),analF(3,3),errorF(3,3),C(3,3),analC(3,3),errorC(3,3),E(3,3),analE(3,3),errorE(3,3),S(3,3),analS(3,3),errorS(3,3), &
    & sigma(3,3),analSigma(3,3),errorSigma(3,3),initialBeta,initialP,xi(3),yValue,zValue
  LOGICAL :: directoryExists=.FALSE.
  CHARACTER(LEN=255) :: commandArgument
  
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: derivedField,dependentField,equationsSetField,geometricField,materialsField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup
 
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
      READ(commandArgument(1:argumentLength),*) interpolationType
      IF(interpolationType<1.OR.interpolationType>9) CALL HandleError("Invalid interpolation specification.")
    ELSE
      interpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    ENDIF
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfDimensions=3
    interpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
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
    xi(1)=0.5_OC_RP
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
    xi(1)=0.5_OC_RP
    xi(2)=0.5_OC_RP
    numberOfTensorComponents=3
  CASE(3)
    numberOfGlobalXNodes=numberOfNodesPerElement
    numberOfGlobalYNodes=numberOfNodesPerElement
    numberOfGlobalZNodes=numberOfNodesPerElement
    numberOfGaussXi=3
    numberOfXi=3
    xi(1)=0.5_OC_RP
    xi(2)=0.5_OC_RP
    xi(3)=0.5_OC_RP
    numberOfTensorComponents=6
  CASE DEFAULT
    CALL HandleError("Invalid number of dimensions.")    
  END SELECT
  numberOfNodes=numberOfGlobalXNodes*numberOfGlobalYNodes*numberOfGlobalZNodes
  
  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  !CALL OC_DiagnosticsSetOn(OC_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"],err)
  !Set output
  CALL OC_OutputSetOn("Uniaxial",err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  
  !Get the number of computational nodes and this computational node number
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)  

  !Create a rectangular cartesian coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"Region",err)
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL OC_Region_CreateFinish(region,err)

  !Define geometric basis
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  IF(numberOfDimensions==2) THEN
    CALL OC_Basis_NumberOfXiSet(basis,2,err)
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    CALL OC_Basis_NumberOfXiSet(basis,3,err)
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF
  CALL OC_Basis_CreateFinish(basis,err)
  
  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  IF(numberOfDimensions==2) THEN
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[1,1],err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH,HEIGHT,WIDTH],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[1,1,1],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)
  
  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL OC_Field_ScalingTypeSet(geometricField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
  CALL OC_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create the dependent field
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
  CALL OC_Field_TypeSet(dependentField,OC_FIELD_GEOMETRIC_GENERAL_TYPE,err)
  CALL OC_Field_DecompositionSet(dependentField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(dependentField,geometricField,err)
  CALL OC_Field_DependentTypeSet(dependentField,OC_FIELD_DEPENDENT_TYPE,err)
  CALL OC_Field_NumberOfVariablesSet(dependentField,2,err)
  CALL OC_Field_VariableTypesSet(dependentField,[OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_T_VARIABLE_TYPE],err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,numberOfDimensions+1,err)
  CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,numberOfDimensions+1,err)
  !Set the pressure to be a constant element
  CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
    & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,numberOfDimensions+1, &
    & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  CALL OC_Field_CreateFinish(dependentField,err)
  
  !Create the material field
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_Field_CreateStart(MATERIAL_FIELD_USER_NUMBER,region,materialsField,err)
  CALL OC_Field_TypeSet(materialsField,OC_FIELD_MATERIAL_TYPE,err)
  CALL OC_Field_DecompositionSet(materialsField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(materialsField,geometricField,err)
  CALL OC_Field_NumberOfVariablesSet(materialsField,1,err)
  CALL OC_Field_VariableLabelSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL OC_Field_NumberOfComponentsSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,2,err)
  CALL OC_Field_CreateFinish(materialsField,err)

  !Create the derived fields
  CALL OC_Field_Initialise(derivedField,err)
  CALL OC_Field_CreateStart(DERIVED_FIELD_USER_NUMBER,region,derivedField,err)
  CALL OC_Field_TypeSet(derivedField,OC_FIELD_GENERAL_TYPE,err)
  CALL OC_Field_DecompositionSet(derivedField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(derivedField,geometricField,err)
  CALL OC_Field_DependentTypeSet(derivedField,OC_FIELD_DEPENDENT_TYPE,err)
  CALL OC_Field_NumberOfVariablesSet(derivedField,5,err)
  CALL OC_Field_VariableTypesSet(derivedField,[OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_U2_VARIABLE_TYPE, &
    & OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_U5_VARIABLE_TYPE],err)
  CALL OC_Field_VariableLabelSet(derivedField,OC_FIELD_U1_VARIABLE_TYPE,"DeformationGradient",err)
  CALL OC_Field_VariableLabelSet(derivedField,OC_FIELD_U2_VARIABLE_TYPE,"RightCauchyGreen",err)
  CALL OC_Field_VariableLabelSet(derivedField,OC_FIELD_U3_VARIABLE_TYPE,"GreenLagrange",err)
  CALL OC_Field_VariableLabelSet(derivedField,OC_FIELD_U4_VARIABLE_TYPE,"SecondPiolaKirchoff",err)
  CALL OC_Field_VariableLabelSet(derivedField,OC_FIELD_U5_VARIABLE_TYPE,"Cauchy",err)
  CALL OC_Field_NumberOfComponentsSet(derivedField,OC_FIELD_U1_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL OC_Field_NumberOfComponentsSet(derivedField,OC_FIELD_U2_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL OC_Field_NumberOfComponentsSet(derivedField,OC_FIELD_U3_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL OC_Field_NumberOfComponentsSet(derivedField,OC_FIELD_U4_VARIABLE_TYPE,numberOfTensorComponents,err)
  CALL OC_Field_NumberOfComponentsSet(derivedField,OC_FIELD_U5_VARIABLE_TYPE,numberOfTensorComponents,err)
  DO componentIdx=1,numberOfTensorComponents
    CALL OC_Field_ComponentInterpolationSet(derivedField,OC_FIELD_U1_VARIABLE_TYPE,componentIdx, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(derivedField,OC_FIELD_U2_VARIABLE_TYPE,componentIdx, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(derivedField,OC_FIELD_U3_VARIABLE_TYPE,componentIdx, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(derivedField,OC_FIELD_U4_VARIABLE_TYPE,componentIdx, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(derivedField,OC_FIELD_U5_VARIABLE_TYPE,componentIdx, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  ENDDO !componentIdx
  CALL OC_Field_CreateFinish(derivedField,err)
  
  !Create the equations_set
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the equations set material field 
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIAL_FIELD_USER_NUMBER,materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01.
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,C1,err)
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,C2,err)

  !Create the derived fields
  CALL OC_EquationsSet_DerivedCreateStart(equationsSet,DERIVED_FIELD_USER_NUMBER,derivedField,err)
  CALL OC_EquationsSet_DerivedVariableSet(equationsSet,OC_EQUATIONS_SET_DERIVED_DEFORMATION_GRADIENT, &
    & OC_FIELD_U1_VARIABLE_TYPE,err)
  CALL OC_EquationsSet_DerivedVariableSet(equationsSet,OC_EQUATIONS_SET_DERIVED_R_CAUCHY_GREEN_DEFORMATION, &
    & OC_FIELD_U2_VARIABLE_TYPE,err)
  CALL OC_EquationsSet_DerivedVariableSet(equationsSet,OC_EQUATIONS_SET_DERIVED_GREEN_LAGRANGE_STRAIN, &
    & OC_FIELD_U3_VARIABLE_TYPE,err)
  CALL OC_EquationsSet_DerivedVariableSet(equationsSet,OC_EQUATIONS_SET_DERIVED_SECOND_PK_STRESS, &
    & OC_FIELD_U4_VARIABLE_TYPE,err)
  CALL OC_EquationsSet_DerivedVariableSet(equationsSet,OC_EQUATIONS_SET_DERIVED_CAUCHY_STRESS, &
    & OC_FIELD_U5_VARIABLE_TYPE,err)
  CALL OC_EquationsSet_DerivedCreateFinish(equationsSet,err)

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and displacement bcs
  DO dimensionIdx1=1,numberOfDimensions
    CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & dimensionIdx1,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,dimensionIdx1,err)
  ENDDO !dimensionIdx1

  !Define the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL OC_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)
  !CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  CALL OC_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  
  !Set x=0 nodes to no displacement in the x-direction
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO yNodeIdx=1,numberOfGlobalYNodes
      nodeNumber=1+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
      ENDIF
    ENDDO !yNodeIdx
  ENDDO !zNodeIdx
  !Set y=0 nodes to no displacement in the y-direction
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO xNodeIdx=1,numberOfGlobalXNodes
      nodeNumber=xNodeIdx+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !y-direction
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
      ENDIF
    ENDDO !yNodeIdx
  ENDDO !zNodeIdx
  IF(numberOfDimensions==3) THEN
    !Set z=0 nodes to no displacement in the z-direction
    DO yNodeIdx=1,numberOfGlobalYNodes
      DO xNodeIdx=1,numberOfGlobalXNodes
        nodeNumber=xNodeIdx+(yNodeIdx-1)*numberOfGlobalXNodes
        CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !z-direction
          CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
            & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        ENDIF
      ENDDO !yNodeIdx
    ENDDO !zNodeIdx
  ENDIF
  !Fix the x=LENGTH nodes to 1+ALPHA*LENGTH
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO yNodeIdx=1,numberOfGlobalYNodes
      nodeNumber=numberOfGlobalXNodes+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        WRITE(*,'("Setting x boundary condition for node ",I0)') nodeNumber
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & OC_BOUNDARY_CONDITION_FIXED,ALPHA*LENGTH,err)
      ENDIF
    ENDDO !yNodeIdx
  ENDDO !zNodeIdx
  ! !Fix the y=HEIGHT nodes to 1-BETA*HEIGHT
  ! DO zNodeIdx=1,numberOfGlobalZNodes
  !   DO xNodeIdx=1,numberOfGlobalXNodes
  !     nodeNumber=numberOfGlobalXNodes*(numberOfGlobalYNodes-1)+xNodeIdx+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
  !     CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
  !     IF(nodeDomain==computationalNodeNumber) THEN
  !       !y-direction
  !       WRITE(*,'("Setting y boundary condition for node ",I0)') nodeNumber
  !       CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
  !         & OC_BOUNDARY_CONDITION_FIXED,-BETA*HEIGHT,err)
  !     ENDIF
  !   ENDDO !xNodeIdx
  ! ENDDO !zNodeIdx
  ! IF(numberOfDimensions==3) THEN
  !   !Fix the z=WIDTH nodes to 1-GAMMA*WIDTH
  !   DO yNodeIdx=1,numberOfGlobalYNodes
  !     DO xNodeIdx=1,numberOfGlobalXNodes
  !       nodeNumber=xNodeIdx+(yNodeIdx-1)*numberOfGlobalXNodes+(numberOfGlobalZNodes-1)*numberOfGlobalXNodes*numberOfGlobalYNodes 
  !       CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
  !       IF(nodeDomain==computationalNodeNumber) THEN
  !         !z-direction
  !         WRITE(*,'("Setting z boundary condition for node ",I0)') nodeNumber
  !         CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
  !           & OC_BOUNDARY_CONDITION_FIXED,-GAMMA*WIDTH,err)
  !       ENDIF
  !     ENDDO !yNodeIdx
  !   ENDDO !yNodeIdx
  ! ENDIF
  
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Set initial values for the dependent field
  initialP=0.1_OC_RP !Compressive
  CALL OC_Field_ParameterSetUpdateElement(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,numberOfDimensions+1,initialP,err)

  !Solve problem
  CALL OC_Problem_Solve(problem,err)

  !Calculate the derived fields
  CALL OC_EquationsSet_DerivedVariableCalculate(equationsSet,OC_EQUATIONS_SET_DERIVED_DEFORMATION_GRADIENT,err)
  CALL OC_EquationsSet_DerivedVariableCalculate(equationsSet,OC_EQUATIONS_SET_DERIVED_R_CAUCHY_GREEN_DEFORMATION,err)
  CALL OC_EquationsSet_DerivedVariableCalculate(equationsSet,OC_EQUATIONS_SET_DERIVED_GREEN_LAGRANGE_STRAIN,err)
  CALL OC_EquationsSet_DerivedVariableCalculate(equationsSet,OC_EQUATIONS_SET_DERIVED_SECOND_PK_STRESS,err)
  CALL OC_EquationsSet_DerivedVariableCalculate(equationsSet,OC_EQUATIONS_SET_DERIVED_CAUCHY_STRESS,err)
  !Construct tensors
  analAlpha=alpha
  analBeta=SQRT(1.0_OC_RP/(1.0_OC_RP+ALPHA))
  analGamma=analBeta
  CALL ComputeError(alpha,analAlpha,errorAlpha)
  CALL OC_Field_ParameterSetGetNode(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,numberOfNodes,2, &
    & yValue,err)
  beta=yValue
  CALL ComputeError(beta,analBeta,errorBeta)
  IF(numberOfDimensions==3) THEN
    CALL OC_Field_ParameterSetGetNode(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,numberOfNodes,3, &
      & zValue,err)
    gamma=zValue
   CALL ComputeError(gamma,analGamma,errorGamma)
  ELSE
    zValue=0.0_OC_RP
    errorGamma=0.0_OC_RP
  ENDIF
  CALL OC_Field_ParameterSetGetElement(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,numberOfDimensions+1,p,err)
  analP=p
  CALL ComputeError(p,analP,errorP)
  F=0.0_OC_RP
  analF=0.0_OC_RP
  errorF=0.0_OC_RP
  C=0.0_OC_RP
  analC=0.0_OC_RP
  errorC=0.0_OC_RP
  E=0.0_OC_RP
  analE=0.0_OC_RP
  errorE=0.0_OC_RP
  S=0.0_OC_RP
  analS=0.0_OC_RP
  errorS=0.0_OC_RP
  sigma=0.0_OC_RP
  analSigma=0.0_OC_RP
  errorSigma=0.0_OC_RP
  SELECT CASE(numberOfDimensions)
  CASE(2)
    J=(1.0_OC_RP+alpha)*(1.0_OC_RP-beta)
    analJ=(1.0_OC_RP+analAlpha)*(1.0_OC_RP-analBeta)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,F(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,F(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,F(1,2),err)
    F(2,1)=F(1,2)
    analF(1,1)=1.0_OC_RP+analAlpha
    analF(2,2)=1.0_OC_RP-analBeta
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,C(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,C(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,C(1,2),err)
    C(2,1)=C(1,2)
    analC(1,1)=(1.0_OC_RP+analAlpha)/(1.0_OC_RP-analBeta)
    analC(2,2)=(1.0_OC_RP-analBeta)/(1.0_OC_RP+analAlpha)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,E(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,E(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,E(1,2),err)
    E(2,1)=E(1,2)
    analE(1,1)=(analAlpha+analBeta)/(2.0_OC_RP*(1.0_OC_RP-analBeta))
    analE(2,2)=-1.0_OC_RP*(analAlpha+analBeta)/(2.0_OC_RP*(1.0_OC_RP+analAlpha))
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,S(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,S(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,S(1,2),err)
    S(2,1)=S(1,2)
    analS(1,1)=2.0_OC_RP*(C1+C2*(1.0_OC_RP-analBeta)/(1.0_OC_RP+analAlpha))
    analS(2,2)=2.0_OC_RP*(C1+C2*(1.0_OC_RP+analAlpha)/(1.0_OC_RP-analBeta))
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,sigma(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,sigma(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,sigma(1,2),err)
    sigma(2,1)=sigma(1,2)
    analSigma(1,1)=analJ**(-2.5_OC_RP)*C1*(analAlpha*(2.0_OC_RP+analAlpha)-analBeta*(2.0_OC_RP-analBeta))+analP
    analSigma(2,2)=-1.0*OC_RP*analJ**(-2.5_OC_RP)*C1*(analAlpha*(2.0_OC_RP+analAlpha)-analBeta*(2.0_OC_RP-analBeta))+analP
  CASE(3)
    J=(1.0_OC_RP+alpha)*(1.0_OC_RP-beta)*(1.0_OC_RP-gamma)
    analJ=(1.0_OC_RP+analAlpha)*(1.0_OC_RP-analBeta)*(1.0_OC_RP-analGamma)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,F(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,F(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,F(3,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,4,F(2,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,5,F(1,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,6,F(1,2),err)
    F(2,1)=F(1,2)
    F(3,1)=F(1,3)
    F(3,2)=F(2,3)
    analF(1,1)=1.0_OC_RP+analAlpha
    analF(2,2)=1.0_OC_RP-analBeta
    analF(3,3)=1.0_OC_RP-analGamma
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,C(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,C(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,C(3,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,4,C(2,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,5,C(1,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,6,C(1,2),err)
    C(2,1)=C(1,2)
    C(3,1)=C(1,3)
    C(3,2)=C(2,3)
    analC(1,1)=(((1.0_OC_RP+analAlpha)**4.0_OC_RP)**(1.0_OC_RP/3.0_OC_RP))/ &
      & ((((1.0_OC_RP-analBeta)**2.0_OC_RP)*((1.0_OC_RP-analGamma)**2.0_OC_RP))**(1.0_OC_RP/3.0_OC_RP))
    analC(2,2)=(((1.0_OC_RP-analBeta)**4.0_OC_RP)**(1.0_OC_RP/3.0_OC_RP))/ &
      & ((((1.0_OC_RP+analAlpha)**2.0_OC_RP)*((1.0_OC_RP-analGamma)**2.0_OC_RP))**(1.0_OC_RP/3.0_OC_RP))
    analC(3,3)=(((1.0_OC_RP-analGamma)**4.0_OC_RP)**(1.0_OC_RP/3.0_OC_RP))/ &
      & ((((1.0_OC_RP+analAlpha)**2.0_OC_RP)*((1.0_OC_RP-analBeta)**2.0_OC_RP))**(1.0_OC_RP/3.0_OC_RP))
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,E(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,E(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,E(3,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,4,E(2,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,5,E(1,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U3_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,6,E(1,2),err)
    E(2,1)=E(1,2)
    E(3,1)=E(1,3)
    E(3,2)=E(2,3)
    analE(1,1)=0.5_OC_RP*(analC(1,1)-1.0_OC_RP)
    analE(2,2)=0.5_OC_RP*(analC(2,2)-1.0_OC_RP)
    analE(3,3)=0.5_OC_RP*(analC(3,3)-1.0_OC_RP)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,S(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,S(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,S(3,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,4,S(2,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,5,S(1,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U4_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,6,S(1,2),err)
    S(2,1)=S(1,2)
    S(3,1)=S(1,3)
    S(3,2)=S(2,3)
    analS(1,1)=2.0_OC_RP*(C1+(C2*((1.0_OC_RP-analBeta)**2.0_OC_RP+(1.0_OC_RP-analGamma)**2.0_OC_RP))/ &
      & (((1.0_OC_RP+analAlpha)**2.0_OC_RP*(1.0_OC_RP-analBeta)**2.0_OC_RP*(1.0_OC_RP-analGamma)**2.0_OC_RP)** &
      & (1.0_OC_RP/3.0_OC_RP)))
    analS(2,2)=2.0_OC_RP*(C1+(C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analGamma)**2.0_OC_RP))/ &
      & (((1.0_OC_RP+analAlpha)**2.0_OC_RP*(1.0_OC_RP-analBeta)**2.0_OC_RP*(1.0_OC_RP-analGamma)**2.0_OC_RP)** &
      & (1.0_OC_RP/3.0_OC_RP)))
    analS(3,3)=2.0_OC_RP*(C1+(C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analBeta)**2.0_OC_RP))/ &
      & (((1.0_OC_RP+analAlpha)**2.0_OC_RP*(1.0_OC_RP-analBeta)**2.0_OC_RP*(1.0_OC_RP-analGamma)**2.0_OC_RP)** &
      & (1.0_OC_RP/3.0_OC_RP)))
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,sigma(1,1),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2,sigma(2,2),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3,sigma(3,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,4,sigma(2,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,5,sigma(1,3),err)
    CALL OC_Field_ParameterSetGetElement(derivedField,OC_FIELD_U5_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,6,sigma(1,2),err)
    sigma(2,1)=sigma(1,2)
    sigma(3,1)=sigma(1,3)
    sigma(3,2)=sigma(2,3)
    analSigma(1,1)=(2.0_OC_RP/(3.0_OC_RP*analJ**3.0_OC_RP))*((2.0_OC_RP*(1.0_OC_RP+analAlpha)**2.0_OC_RP* &
      & (C1+C2*((1.0_OC_RP-analBeta)**2.0_OC_RP+(1.0_OC_RP-analGamma)**2.0_OC_RP))) &
      & -((1.0_OC_RP-analBeta)**2.0_OC_RP*(C1+C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analGamma)**2.0_OC_RP))) &
      & -((1.0_OC_RP-analGamma)**2.0_OC_RP*(C1+C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analBeta)**2.0_OC_RP))) &
      & )+analP
    analSigma(2,2)=(2.0_OC_RP/(3.0_OC_RP*analJ**3.0_OC_RP))*((-1.0_OC_RP*(1.0_OC_RP+analAlpha)**2.0_OC_RP* &
      & (C1+C2*((1.0_OC_RP+analBeta)**2.0_OC_RP+(1.0_OC_RP-analGamma)**2.0_OC_RP))) &
      & +(2.0_OC_RP*(1.0_OC_RP-analBeta)**2.0_OC_RP*(C1+C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analGamma)** &
      & 2.0_OC_RP))) &
      & -((1.0_OC_RP-analGamma)**2.0_OC_RP*(C1+C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analBeta)**2.0_OC_RP))) &
      & )+analP
    analSigma(3,3)=(2.0_OC_RP/(3.0_OC_RP*analJ**3.0_OC_RP))*((-1.0_OC_RP*(1.0_OC_RP+analAlpha)**2.0_OC_RP* &
      & (C1+C2*((1.0_OC_RP+analBeta)**2.0_OC_RP+(1.0_OC_RP-analGamma)**2.0_OC_RP))) &
      & -((1.0_OC_RP-analBeta)**2.0_OC_RP*(C1+C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analGamma)** &
      & 2.0_OC_RP))) &
      & +(2.0_OC_RP*(1.0_OC_RP-analGamma)**2.0_OC_RP*(C1+C2*((1.0_OC_RP+analAlpha)**2.0_OC_RP+(1.0_OC_RP-analBeta)** &
      & 2.0_OC_RP))))+analP
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
  CALL OC_Fields_Initialise(fields,err)
  CALL OC_Fields_Create(region,fields,err)
  CALL OC_Fields_NodesExport(fields,"./results/Uniaxial","FORTRAN",err)
  CALL OC_Fields_ElementsExport(fields,"./results/Uniaxial","FORTRAN",err)
  CALL OC_Fields_Finalise(fields,err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    
    STOP
  END SUBROUTINE HandleError

  SUBROUTINE ComputeError(actual,analytic,error)
    
    REAL(OC_RP), INTENT(IN) :: actual
    REAL(OC_RP), INTENT(IN) :: analytic
    REAL(OC_RP), INTENT(OUT) :: error
    
    IF(ABS(analytic)<=5.0*OC_RP*EPSILON(0.0_OC_RP)) THEN
      error=(actual-analytic)
    ELSE
      error=(actual-analytic)*100.0_OC_RP/analytic
    ENDIF
    
    RETURN
  END SUBROUTINE ComputeError

  SUBROUTINE WriteScalarError(outputType,numberOfDimensions,symbol,actual,analytic,error)
    
    INTEGER(OC_Intg), INTENT(IN) :: outputType
    INTEGER(OC_Intg), INTENT(IN) :: numberOfDimensions
    CHARACTER(LEN=*), INTENT(IN) :: symbol
    REAL(OC_RP), INTENT(IN) :: actual
    REAL(OC_RP), INTENT(IN) :: analytic
    REAL(OC_RP), INTENT(IN) :: error

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
    
    INTEGER(OC_Intg), INTENT(IN) :: outputType
    INTEGER(OC_Intg), INTENT(IN) :: numberOfDimensions
    CHARACTER(LEN=*), INTENT(IN) :: symbol
    REAL(OC_RP), INTENT(IN) :: actual(:,:)
    REAL(OC_RP), INTENT(IN) :: analytic(:,:)
    REAL(OC_RP), INTENT(IN) :: error(:,:)
    
    INTEGER(OC_Intg) :: dimensionIdx1,dimensionIdx2
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
