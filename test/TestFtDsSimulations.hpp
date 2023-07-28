#ifndef TESTFTDSSIMULATIONS_HPP_
#define TESTFTDSSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "CellSrnModel.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DeltaNotchInteriorSrnModel.hpp" ///\todo This should be replaced with a new class FtDsInteriorSrnModel
#include "DeltaNotchEdgeInteriorTrackingModifier.hpp" ///\todo This should be replaced with a new class FtDsEdgeInteriorTrackingModifier
#include "FtDsEdgeSrnModel.hpp"
#include "FtDsEdgeTrackingModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

class TestFtDsSimulations : public AbstractCellBasedTestSuite
{
public:

    void TestFtDsPolarisationWithoutCellDivision()
    {
        EXIT_IF_PARALLEL;

        // Create regular heaxgonal mesh object
        HoneycombVertexMeshGenerator generator(5, 3);
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

            // Initialise edge-based SRN
            auto p_cell_edge_srn_model = new CellSrnModel();
<<<<<<< HEAD:test/TestFtDsSimulations.hpp
            auto centroid_x = p_mesh->GetCentroidOfElement(index)[0];
            auto centroid_y = p_mesh->GetCentroidOfElement(index)[1];

            // Initialise protein concentrations
            ///\todo why are notch and delta specified?
=======
            auto centroid_x = p_mesh->GetCentroidOfElement(elem_index)[0];
            auto centroid_y = p_mesh->GetCentroidOfElement(elem_index)[1];
            /* We choose to initialise the total concentrations of DsP and FtP to reflect a radial Fj gradient, and Ds and Ft to oppose it */
>>>>>>> d011ac11d1325cef62adf074e44bee35c91f5f16:test/Test2dVertexBasedSimulationWithSrnModels_LT_ODE.hpp
            auto notch_concentration = -((5+centroid_x * centroid_x) * (5+(centroid_y -1.5) * (centroid_y -1.5))) +1.5+216; // 50*(1-((80 + 10 * centroid)/ (80 + 10 * 6)));
            auto delta_concentration = -((5+centroid_x * centroid_x) * (5+(centroid_y -1.5) * (centroid_y -1.5))) +1.5+216; // 50*(1-((80 + 10 * elem_index)) / (80 + 10 * 3));
            auto DsP_concentration = (5+centroid_x * centroid_x) * (5+(centroid_y -1.5) * (centroid_y -1.5)) - 1.5; //50*((80 + 10 * centroid) / (80 + 10 * 6));
            auto FtP_concentration = (5+centroid_x * centroid_x) * (5+(centroid_y -1.5) * (centroid_y -1.5)) - 1.5; //50*((80 + 10 * centroid) / (80 + 10 * 6));
            auto A_concentration = 0;
            auto B_concentration = 0;
            auto C_concentration = 0;
            auto D_concentration = 0;
            auto A__concentration = 0; ///\todo how is this different to A_concentration?
            auto B__concentration = 0; ///\todo how is this different to B_concentration?
            auto C__concentration = 0; ///\todo how is this different to C_concentration?
            auto D__concentration = 0; ///\todo how is this different to D_concentration?

            double total_edge_length = 0.0;
            auto p_element = p_mesh->GetElement(index);
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            // Create SRN model for each edge
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                auto p_elem_edge = p_element->GetEdge(i);
                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                // Initial concentration of delta and notch vary depending on the edge length
                initial_conditions.push_back(p_edge_length / total_edge_length * notch_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * delta_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * DsP_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * FtP_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * A_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * B_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * C_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * D_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * A__concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * B__concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * C__concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * D__concentration);

                MAKE_PTR(FtDsEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

            // Create interior SRN model
            MAKE_PTR(DeltaNotchInteriorSrnModel, p_cell_srn_model);
            std::vector<double> zero_conditions(2);
            p_cell_srn_model->SetInitialConditions(zero_conditions);
            p_cell_edge_srn_model->SetInteriorSrnModel(p_cell_srn_model);

            // Create cell object
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(-1.0);
            cells.push_back(p_cell);
        }

        // Create vertex-based cell population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create cell-based simulation object
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFtDsPolarisationWithoutCellDivision");        
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetEndTime(1.0);

        // Create and pass tracking modifiers to simulation
        MAKE_PTR(DeltaNotchEdgeInteriorTrackingModifier<2>, p_cell_modifier);
        simulator.AddSimulationModifier(p_cell_modifier);
        MAKE_PTR(FtDsEdgeTrackingModifier<2>, p_edge_modifier);
        simulator.AddSimulationModifier(p_edge_modifier);

<<<<<<< HEAD:test/TestFtDsSimulations.hpp
        // Run the simulation by calling Solve()
=======
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.SetOutputDivisionLocations(true);
        // Add division rule
        boost::shared_ptr<ShortAxisVertexBasedDivisionRule<2> > p_division_rule(new ShortAxisVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        
        simulator.SetSamplingTimestepMultiple(1.0);
        //simulator.SetDt(0.1);
        simulator.SetEndTime(1.0); // error when this is larger than 1.0

>>>>>>> d011ac11d1325cef62adf074e44bee35c91f5f16:test/Test2dVertexBasedSimulationWithSrnModels_LT_ODE.hpp
        simulator.Solve();
    }

    void noTestFtDsPolarisationWithCellDivision()
    {
        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);

        for (unsigned elem_index = 0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            BernoulliTrialCellCycleModel* p_cc_model = new BernoulliTrialCellCycleModel();
            p_cc_model->SetDimension(2);
            p_cc_model->SetBirthTime(-1.0);
            p_cc_model->SetDivisionProbability(0.1);

            auto p_element = p_mesh->GetElement(elem_index);

            /* Initialise edge based SRN */
            auto p_cell_edge_srn_model = new CellSrnModel();

            /* We choose to initialise the total concentrations as in Eman's initial conditions */

            auto notch_concentration = 50*(1-((80 + 10 * elem_index) / (80 + 10 * 6)));
            auto delta_concentration = 50*(1-((80 + 10 * elem_index) / (80 + 10 * 6)));
            auto DsP_concentration =  50*((80 + 10 * elem_index) / (80 + 10 * 6));// [(6 + elem_index * elem_index) / (6+(elem_index-3)*(elem_index-3))]-3; 
            auto FtP_concentration =  50*((80 + 10 * elem_index) / (80 + 10 * 6));// [(6 + elem_index * elem_index) / (6+(elem_index-3)*(elem_index-3))]-3;
            auto A_concentration = 0;
            auto B_concentration = 0;
            auto C_concentration = 0;
            auto D_concentration = 0;
            auto A__concentration = 0;
            auto B__concentration = 0;
            auto C__concentration = 0;
            auto D__concentration = 0;

            // auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
            // auto notch_concentration = RandomNumberGenerator::Instance()->ranf();
            // auto C_concentration = RandomNumberGenerator::Instance()->ranf();

            double total_edge_length = 0.0;
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                auto p_elem_edge = p_element->GetEdge(i);
                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch vary depending on the edge length */

                initial_conditions.push_back(p_edge_length / total_edge_length * notch_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * delta_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * DsP_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * FtP_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * A_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * B_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * C_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * D_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * A__concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * B__concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * C__concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * D__concentration);

                MAKE_PTR(FtDsEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create vertex-based cell population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create cell-based simulation object
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFtDsPolarisationWithCellDivision");
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create and pass tracking modifiers to simulation
        MAKE_PTR(FtDsEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create and pass target area growth rule to simulation
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create and pass cell division rule to simulation
        boost::shared_ptr<ShortAxisVertexBasedDivisionRule<2> > p_division_rule(new ShortAxisVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule);

        // Run the simulation by calling Solve()
        simulator.Solve();
    }
};

#endif /*TESTFTDSSIMULATIONS_HPP_*/
