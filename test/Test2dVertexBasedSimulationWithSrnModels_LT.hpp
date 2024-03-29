#ifndef TESTDELTANOTCHEDGEINTERIORODESIMULATION_HPP_
#define TESTDELTANOTCHEDGEINTERIORODESIMULATION_HPP_

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

#include "CellSrnModel.hpp"

#include "FtDsEdgeSrnModel.hpp"
#include "FtDsEdgeTrackingModifier.hpp"

#include "FtDsInteriorSrnModel.hpp"
#include "FtDsEdgeInteriorTrackingModifier.hpp"


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

// for Voronoi mesh generator
#include "VoronoiVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "Warnings.hpp"

static const double M_DT = 0.1;
static const double M_RELAXATION_TIME = 10; 
static const double M_EXTENSION_TIME = 10; 
static const double M_VIS_TIME_STEP = 1;
static const unsigned M_NUM_CELLS_WIDE = 14;
static const unsigned M_NUM_CELLS_HIGH = 20;
// Specify mechanical parameter values
    // -0.259,0.172
double k = 1.0;
double lambda_bar = 0.05;
double gamma_bar = 0.04;
double heterotypic_line_tension_multiplier = 2.0;
double supercontractile_line_tension_multiplier = 2.0;

/**
 * These tests check and demonstrate simulation of vertex based models with edge and interior Srn models
 */
class TestFtDsEdgeOnlyODESimulation_LT : public AbstractCellBasedTestSuite
{
public:
    /*
     * Test vertex based simulations when both edge AND interior SRN models are specified
     */
    void TestRunningMultiODECellWithEdgesAndInterior()
    {
        EXIT_IF_PARALLEL;
        // /* First we create a regular vertex mesh. */
        // HoneycombVertexMeshGenerator generator(6, 6);
        // MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

       #if BOOST_VERSION >= 105200

        // Generate a mesh that is 20 cells wide, 12 high, with 4 Lloyd's relaxation steps and target average element area 1.23
        VoronoiVertexMeshGenerator generator(20, 12, 4, 1.23);
        MutableVertexMesh<2,2>* p_mesh_a = generator.GetMesh();
        MutableVertexMesh<2,2>* p_mesh_b = generator.GetMeshAfterReMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            auto p_element = p_mesh->GetElement(elem_index);
            /* Initialise edge based SRN */
            auto p_cell_edge_srn_model = new CellSrnModel();
            /* We choose to initialise the total concentrations to random levels */
            auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
            auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

            double total_edge_length = 0.0;
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                auto p_elem_edge = p_element->GetEdge(i);
                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch vary depending on the edge length */
                initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
                initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

                MAKE_PTR(FtDsEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }
            //Add interior SRN models to cells
            MAKE_PTR(FtDsInteriorSrnModel, p_cell_srn_model);
            std::vector<double> zero_conditions(2);
            p_cell_srn_model->SetInitialConditions(zero_conditions);
            p_cell_edge_srn_model->SetInteriorSrnModel(p_cell_srn_model);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFtDsEdgeInteriorODESimulation_LT");
        //simulator.SetSamplingTimestepMultiple(10);
        //simulator.SetEndTime(10.0);

        /* Update CellData and CellEdgeData so that SRN simulations can run properly */
        MAKE_PTR(FtDsEdgeInteriorTrackingModifier<2>, p_cell_modifier);
        simulator.AddSimulationModifier(p_cell_modifier);
        MAKE_PTR(FtDsEdgeTrackingModifier<2>, p_edge_modifier);
        simulator.AddSimulationModifier(p_edge_modifier);

        MAKE_PTR(ForceForScenario4<2>, p_force);
        //p_force->SetNumStripes(0);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*lambda_bar*pow(k,1.5));
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
        simulator.AddSimulationModifier(p_growth_modifier);
        //TS_ASSERT_THROWS_NOTHING(simulator.Solve());


        simulator.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulator.Solve();
    }

    /*
     * Test whether vertex based models can run when ONLY edge SRN models are specified
     */
    void TestRunningMultiODECellWithEdges()
    {
        // /* First we create a regular vertex mesh. */
        // HoneycombVertexMeshGenerator generator(6, 6);
        // MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

       #if BOOST_VERSION >= 105200

        // Generate a mesh that is 20 cells wide, 12 high, with 4 Lloyd's relaxation steps and target average element area 1.23
        VoronoiVertexMeshGenerator generator(20, 12, 4, 1.23);
        MutableVertexMesh<2,2>* p_mesh_a = generator.GetMesh();
        MutableVertexMesh<2,2>* p_mesh_b = generator.GetMeshAfterReMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            auto p_element = p_mesh->GetElement(elem_index);

            /* Initialise edge based SRN */
            auto p_cell_edge_srn_model = new CellSrnModel();

            /* We choose to initialise the total concentrations to random levels */
            auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
            auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

            double total_edge_length = 0.0;
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                auto p_elem_edge = p_element->GetEdge(i);
                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch vary depending on the edge length */
                initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
                initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

                MAKE_PTR(FtDsEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFtDsEdgeOnlyODESimulation_LT");
        //simulator.SetSamplingTimestepMultiple(10);
        //simulator.SetEndTime(10.0);

        /* Update CellEdgeData so that SRN simulations can run properly */
        MAKE_PTR(FtDsEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);


        MAKE_PTR(ForceForScenario4<2>, p_force);
        //p_force->SetNumStripes(0);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);


        simulator.AddForce(p_force);

       // MAKE_PTR(NagaiHondaForce<2>, p_Nforce);
       // simulator.AddForce(p_Nforce);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        //TS_ASSERT_THROWS_NOTHING(simulator.Solve());


        

        MAKE_PTR(ExtrinsicPullModifier, p_pull_modifier);
        p_pull_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_pull_modifier->SetSpeed(0.06);
        simulator.AddSimulationModifier(p_pull_modifier);

       

        simulator.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulator.Solve();

    }
};


#endif /*TESTDELTANOTCHEDGEINTERIORODESIMULATION_HPP_*/
