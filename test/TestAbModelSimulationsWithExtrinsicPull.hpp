#ifndef TESTABMODELSIMULATIONS_HPP_
#define TESTABMODELSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"



#include "BernoulliTrialCellCycleModel.hpp"
#include "CellSrnModel.hpp"
#include "FarhadifarForce.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "AbEdgeSrnModel.hpp"
#include "AbEdgeTrackingModifier.hpp"

// from Sidekick
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "FakePetscSetup.hpp"
#include "ExtrinsicPullModifier.hpp"
#include "SidekickBoundaryCondition.hpp"
#include "ForceForScenario4.hpp"


/*
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "CellSrnModel.hpp"
#include "FakePetscSetup.hpp"
#include "FarhadifarForce.hpp"
#include "AbEdgeSrnModel.hpp"
#include "AbEdgeTrackingModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"


// from Sidekick
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
//#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
//#include "NoCellCycleModel.hpp"
//#include "ConstantTargetAreaModifier.hpp"
//#include "OffLatticeSimulation.hpp"
//#include "FakePetscSetup.hpp"
#include "ExtrinsicPullModifier.hpp"
#include "SidekickBoundaryCondition.hpp"
#include "ForceForScenario4.hpp"

*/

static const double M_DT = 0.1;
static const double M_RELAXATION_TIME = 10; 
static const double M_EXTENSION_TIME = 10; 
static const double M_VIS_TIME_STEP = 1;
static const unsigned M_NUM_CELLS_WIDE = 14;
static const unsigned M_NUM_CELLS_HIGH = 20;
// Specify mechanical parameter values
    // -0.259,0.172
double k = 1.0; //1.0
double lambda_bar = 0.05; //0.05
double gamma_bar = 0.04; //0.04
double heterotypic_line_tension_multiplier = 3.0; //2.0
double supercontractile_line_tension_multiplier = 2.0; //2.0


class TestAbModelSimulations : public AbstractCellBasedTestSuite
{
public:


    void TestAbModelWithoutCellDivision()
    {
        // Create regular hexagonal mesh object
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();

        // Create cell objects, one for each hexagonal element of the mesh
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);
        for (unsigned index = 0; index < p_mesh->GetNumElements(); ++index)
        {
            // Omit cell division, so use 'no cell cycle' model
            NoCellCycleModel* p_cc_model = new NoCellCycleModel();
            p_cc_model->SetDimension(2);

            // Initial conditions on each edge
            std::vector<double> initial_conditions;
            initial_conditions.push_back(1.0); // A
            initial_conditions.push_back(1.0); // B
            initial_conditions.push_back(0.0); // C
            initial_conditions.push_back(0.0); // neigh_C

            // Create cell SRN with associated edge SRNs
            auto p_cell_srn = new CellSrnModel();
            auto p_element = p_mesh->GetElement(index);
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                MAKE_PTR(AbEdgeSrnModel, p_edge_srn);
                p_edge_srn->SetInitialConditions(initial_conditions);
                p_cell_srn->AddEdgeSrnModel(p_edge_srn);
            }

            // Create cell object
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(-1.0);
            cells.push_back(p_cell);
        }

        // Create vertex-based cell population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create cell-based simulation object
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestAbModelWithoutCellDivisionWithExtrinsicPull");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10.0);

        // Create and pass tracking modifier to simulation
        MAKE_PTR(AbEdgeTrackingModifier<2>, p_edge_modifier);
        simulator.AddSimulationModifier(p_edge_modifier);


         MAKE_PTR(ForceForScenario4<2>, p_force);
        //p_force->SetNumStripes(0);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.0));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.0));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.0));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*lambda_bar*pow(k,1.0));
        p_force->SetUseCombinedInterfacesForLineTension(false);

          
        simulator.AddForce(p_force);

        

        MAKE_PTR(ExtrinsicPullModifier, p_pull_modifier);
        p_pull_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_pull_modifier->SetSpeed(0.06);
        simulator.AddSimulationModifier(p_pull_modifier);

       // MAKE_PTR(NagaiHondaForce<2>, p_Nforce);
       // simulator.AddForce(p_Nforce);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0.01);
        simulator.AddSimulationModifier(p_growth_modifier);
        //TS_ASSERT_THROWS_NOTHING(simulator.Solve());


        simulator.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);

        // Run the simulation by calling Solve()
        simulator.Solve();
    }


    void TestAbModelWithCellDivision()
    {
        // Create regular hexagonal mesh object
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();

        // Create cell objects, one for each hexagonal element of the mesh
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        for (unsigned index = 0; index < p_mesh->GetNumElements(); ++index)
        {
            // Create 'Bernoulli trial' based cell cycle model
            BernoulliTrialCellCycleModel* p_cc_model = new BernoulliTrialCellCycleModel();
            p_cc_model->SetDimension(2);
            p_cc_model->SetBirthTime(-1.0);
            p_cc_model->SetDivisionProbability(0.001);

            // Initial conditions on each edge
            std::vector<double> initial_conditions;
            initial_conditions.push_back(1.0); // A
            initial_conditions.push_back(1.0); // B
            initial_conditions.push_back(0.0); // C
            initial_conditions.push_back(0.0); // neigh_C

            // Create cell SRN with associated edge SRNs
            auto p_cell_srn = new CellSrnModel();
            auto p_element = p_mesh->GetElement(index);
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                MAKE_PTR(AbEdgeSrnModel, p_edge_srn);
                p_edge_srn->SetInitialConditions(initial_conditions);
                p_cell_srn->AddEdgeSrnModel(p_edge_srn);
            }

            // Create cell object
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn));
            p_cell->SetCellProliferativeType(p_transit_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create vertex-based cell population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create cell-based simulation object
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestAbModelWithCellDivisionWithExtrinsicPull");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(50.0);

        // Create and pass tracking modifiers to simulation
        MAKE_PTR(AbEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create and pass force law to simulation
        //MAKE_PTR(FarhadifarForce<2>, p_force);
        //simulator.AddForce(p_force);



        // Create and pass cell division rule to simulation
        boost::shared_ptr<ShortAxisVertexBasedDivisionRule<2> > p_division_rule(new ShortAxisVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule);


         MAKE_PTR(ForceForScenario4<2>, p_force);
        //p_force->SetNumStripes(0);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.0));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.0));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.0));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*lambda_bar*pow(k,1.0));
        p_force->SetUseCombinedInterfacesForLineTension(false);

          
        simulator.AddForce(p_force);

        

        MAKE_PTR(ExtrinsicPullModifier, p_pull_modifier);
        p_pull_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_pull_modifier->SetSpeed(0.06);
        simulator.AddSimulationModifier(p_pull_modifier);

       // MAKE_PTR(NagaiHondaForce<2>, p_Nforce);
       // simulator.AddForce(p_Nforce);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0.01);
        simulator.AddSimulationModifier(p_growth_modifier);
        //TS_ASSERT_THROWS_NOTHING(simulator.Solve());


        simulator.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);

        // Run the simulation by calling Solve()
        simulator.Solve();
    }
};

#endif /*TESTABMODELSIMULATIONS_HPP_*/