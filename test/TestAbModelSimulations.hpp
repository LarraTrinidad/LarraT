#ifndef TESTABMODELSIMULATIONS_HPP_
#define TESTABMODELSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "CellSrnModel.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "FakePetscSetup.hpp"
#include "FarhadifarForce.hpp"
#include "AbEdgeSrnModel.hpp"
#include "AbEdgeTrackingModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"



class TestAbModelSimulations : public AbstractCellBasedTestSuite
{
public:

    void TestAbModelWithoutCellDivision()
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
        simulator.SetOutputDirectory("TestAbModelWithoutCellDivision");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10.0);

        // Create and pass tracking modifier to simulation
        MAKE_PTR(AbEdgeTrackingModifier<2>, p_edge_modifier);
        simulator.AddSimulationModifier(p_edge_modifier);

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
        simulator.SetOutputDirectory("TestAbModelWithCellDivision");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(10.0); //10

        // Create and pass tracking modifiers to simulation
        MAKE_PTR(AbEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create and pass force law to simulation
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);

        // Create and pass target area growth rule to simulation
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0.01);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create and pass cell division rule to simulation
        boost::shared_ptr<ShortAxisVertexBasedDivisionRule<2> > p_division_rule(new ShortAxisVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule);

        // Run the simulation by calling Solve()
        simulator.Solve();
    }
};

#endif /*TESTABMODELSIMULATIONS_HPP_*/
