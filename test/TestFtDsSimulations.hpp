#ifndef TESTFTDSSIMULATIONS_HPP_
#define TESTFTDSSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "CellSrnModel.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "FakePetscSetup.hpp"
#include "FarhadifarForce.hpp"
#include "FtDsEdgeSrnModel.hpp"
#include "FtDsEdgeTrackingModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
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
        // Create regular hexagonal mesh object
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
            double centroid_x = p_mesh->GetCentroidOfElement(index)[0];
            double centroid_y = p_mesh->GetCentroidOfElement(index)[1];

            /*
             * Initialise the total concentrations of DsP and FtP to reflect a 
             * radial Fj gradient, and Ds and Ft to oppose it.
             */
            double Ds_concentration = -((5 + centroid_x * centroid_x) * (5 + (centroid_y -1.5) * (centroid_y -1.5))) + 1.5 + 216;
            double Ft_concentration = -((5 + centroid_x * centroid_x) * (5 + (centroid_y -1.5) * (centroid_y -1.5))) + 1.5 + 216;
            double DsP_concentration = (5 + centroid_x * centroid_x) * (5 + (centroid_y -1.5) * (centroid_y -1.5)) - 1.5;
            double FtP_concentration = (5 + centroid_x * centroid_x) * (5 + (centroid_y -1.5) * (centroid_y -1.5)) - 1.5;
            double A_concentration = 0;
            double B_concentration = 0;
            double C_concentration = 0;
            double D_concentration = 0;
            double neigh_A_concentration = 0;
            double neigh_B_concentration = 0;
            double neigh_C_concentration = 0;
            double neigh_D_concentration = 0;

            double total_edge_length = 0.0;
            auto p_element = p_mesh->GetElement(index);
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            // Create SRN model for each edge
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                // Initial concentration of Ds and Ft depend on edge length
                double factor = p_element->GetEdge(i)->rGetLength() / total_edge_length;
                
                std::vector<double> initial_conditions;
                initial_conditions.push_back(factor * Ds_concentration);
                initial_conditions.push_back(factor * Ft_concentration);
                initial_conditions.push_back(factor * DsP_concentration);
                initial_conditions.push_back(factor * FtP_concentration);
                initial_conditions.push_back(factor * A_concentration);
                initial_conditions.push_back(factor * B_concentration);
                initial_conditions.push_back(factor * C_concentration);
                initial_conditions.push_back(factor * D_concentration);
                initial_conditions.push_back(factor * neigh_A_concentration);
                initial_conditions.push_back(factor * neigh_B_concentration);
                initial_conditions.push_back(factor * neigh_C_concentration);
                initial_conditions.push_back(factor * neigh_D_concentration);

                MAKE_PTR(FtDsEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

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
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(10.0);

        // Create and pass tracking modifier to simulation
        MAKE_PTR(FtDsEdgeTrackingModifier<2>, p_edge_modifier);
        simulator.AddSimulationModifier(p_edge_modifier);

        // Run the simulation by calling Solve()
        simulator.Solve();
    }

    void TestFtDsPolarisationWithCellDivision()
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
            p_cc_model->SetDivisionProbability(0.1);

            // Initialise edge-based SRN
            auto p_cell_edge_srn_model = new CellSrnModel();

            // Initialise the total concentrations as in Eman's work
            double Ds_concentration = 50*(1 - ((80 + 10 * index) / (80 + 10 * 6)));
            double Ft_concentration = 50*(1 - ((80 + 10 * index) / (80 + 10 * 6)));
            double DsP_concentration = 50*((80 + 10 * index) / (80 + 10 * 6));
            double FtP_concentration = 50*((80 + 10 * index) / (80 + 10 * 6));
            double A_concentration = 0;
            double B_concentration = 0;
            double C_concentration = 0;
            double D_concentration = 0;
            double neigh_A_concentration = 0;
            double neigh_B_concentration = 0;
            double neigh_C_concentration = 0;
            double neigh_D_concentration = 0;

            double total_edge_length = 0.0;
            auto p_element = p_mesh->GetElement(index);
            for (unsigned i = 0; i < p_element->GetNumEdges(); ++i)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            // Create SRN model for each edge
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                // Initial concentration of Ds and Ft depend on edge length
                double factor = p_element->GetEdge(i)->rGetLength() / total_edge_length;

                std::vector<double> initial_conditions;
                initial_conditions.push_back(factor * Ds_concentration);
                initial_conditions.push_back(factor * Ft_concentration);
                initial_conditions.push_back(factor * DsP_concentration);
                initial_conditions.push_back(factor * FtP_concentration);
                initial_conditions.push_back(factor * A_concentration);
                initial_conditions.push_back(factor * B_concentration);
                initial_conditions.push_back(factor * C_concentration);
                initial_conditions.push_back(factor * D_concentration);
                initial_conditions.push_back(factor * neigh_A_concentration);
                initial_conditions.push_back(factor * neigh_B_concentration);
                initial_conditions.push_back(factor * neigh_C_concentration);
                initial_conditions.push_back(factor * neigh_D_concentration);

                MAKE_PTR(FtDsEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create vertex-based cell population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create cell-based simulation object
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFtDsPolarisationWithCellDivision");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(10.0);

        // Create and pass tracking modifiers to simulation
        MAKE_PTR(FtDsEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create and pass force law to simulation
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);

        // Create and pass target area growth rule to simulation
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(10.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create and pass cell division rule to simulation
        boost::shared_ptr<ShortAxisVertexBasedDivisionRule<2> > p_division_rule(new ShortAxisVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule);

        // Run the simulation by calling Solve()
        simulator.Solve();
    }
};

#endif /*TESTFTDSSIMULATIONS_HPP_*/
