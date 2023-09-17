#include "AbEdgeTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellSrnModel.hpp"
#include "AbEdgeSrnModel.hpp"

template<unsigned DIM>
AbEdgeTrackingModifier<DIM>::AbEdgeTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mUnboundProteinDiffusionCoefficient(0.03)
{
}

template<unsigned DIM>
AbEdgeTrackingModifier<DIM>::~AbEdgeTrackingModifier()
{
}

template<unsigned DIM>
void AbEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AbEdgeTrackingModifier<DIM>::SetupSolve(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation,
    std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AbEdgeTrackingModifier<DIM>::UpdateCellData(
    AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    double D = mUnboundProteinDiffusionCoefficient;

    // For ease, store a static cast of the vertex-based cell population
    assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
    auto p_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over cells
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get this cell's cell-level SRN model and number of edges
        assert(dynamic_cast<CellSrnModel*>(cell_iter->GetSrnModel()));
        auto p_cell_srn = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        unsigned num_edges = p_cell_srn->GetNumEdgeSrn();

        // Cell's edge data
        std::vector<double> A_old;
        std::vector<double> B_old;
        std::vector<double> C_old;
        std::vector<double> edge_lengths;

        
        //\ looks like the new cell edges (from junctional rearrangements after cell division) don't accumulate protein - new junctions aren't being taken into account in this for loop?
        for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {  
            // Store the current unbound protein concentrations on this edge
            auto p_edge_srn = boost::static_pointer_cast<AbEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
            A_old.push_back(p_edge_srn->GetA());
            B_old.push_back(p_edge_srn->GetB());
            C_old.push_back(p_edge_srn->GetC());

            // Store this edge's length
            auto p_element = p_population->GetElementCorrespondingToCell(*cell_iter);
            double edge_length = p_element->GetEdge(edge_index)->rGetLength();
            edge_lengths.push_back(edge_length);
        }

        /*
         * Update unbound protein concentrations based on a linear diffusive 
         * flux between neighbouring edges.
         */
        std::vector<double> A_new(num_edges);
        std::vector<double> B_new(num_edges);
        std::vector<double> C_new(num_edges);

           for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {
           auto p_edge_srn = boost::static_pointer_cast<AbEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
           unsigned prev_index = (edge_index == 0) ? num_edges - 1 : edge_index - 1;
           unsigned next_index = (edge_index == num_edges - 1) ? 0 : edge_index + 1;

           ///\todo consider validity of diffusive flux expression
           double dx = 0.5 * (edge_lengths[prev_index] + edge_lengths[edge_index]);
           double dt = SimulationTime::Instance()->GetTimeStep();
           double time_elapsed = SimulationTime::Instance()->GetTimeStepsElapsed();
           /* This if statement must be used as UpdateCellData is called for setup solve
             and at the end of each time step. Thus if no time has elapsed diffusion should
             not be occuring */
             if (time_elapsed == 0){
              A_new[edge_index] = A_old[edge_index];
              B_new[edge_index] = B_old[edge_index];
              C_new[edge_index] = C_old[edge_index];
             } else {
              A_new[edge_index] = A_old[edge_index] + D*(A_old[prev_index] - 2*A_old[edge_index] + A_old[next_index])*dt/(dx*dx);
              B_new[edge_index] = B_old[edge_index] + D*(B_old[prev_index] - 2*B_old[edge_index] + B_old[next_index])*dt/(dx*dx);
              C_new[edge_index] = C_old[edge_index];
             }
             
             p_edge_srn->SetA(A_new[edge_index]);
             p_edge_srn->SetB(B_new[edge_index]);
             p_edge_srn->SetC(C_new[edge_index]);
        }
        

        // Note: state variables must be in the same order as in AbOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge A", A_new);
        cell_iter->GetCellEdgeData()->SetItem("edge B", B_new);
        cell_iter->GetCellEdgeData()->SetItem("edge C", C_new);
    }

    // After the edge data is filled, fill the edge neighbour data

    // Iterate over cells again
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get this cell's cell-level SRN model and number of edges
        auto p_cell_srn = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        unsigned num_edges = p_cell_srn->GetNumEdgeSrn();

        std::vector<double> neigh_mean_A(num_edges);
        std::vector<double> neigh_mean_B(num_edges);

        for (unsigned edge_index = 0; edge_index < num_edges; ++edge_index)
        {
            // Get neighbouring cell's values
            ///\todo is this correct?
            auto elem_neighbours = p_population->GetNeighbouringEdgeIndices(*cell_iter, edge_index);
            for (auto neighbour : elem_neighbours)
            {
                auto p_cell = p_population->GetCellUsingLocationIndex(neighbour.first);
                auto p_data = p_cell->GetCellEdgeData();
                std::vector<double> neighbour_A_vec = p_data->GetItem("edge A");
                std::vector<double> neighbour_B_vec = p_data->GetItem("edge B");
                neigh_mean_A[edge_index] += neighbour_A_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_B[edge_index] += neighbour_B_vec[neighbour.second] / elem_neighbours.size();
            }
        }

        cell_iter->GetCellEdgeData()->SetItem("neighbour A", neigh_mean_A);
        cell_iter->GetCellEdgeData()->SetItem("neighbour B", neigh_mean_B);
    }
}

template<unsigned DIM>
double AbEdgeTrackingModifier<DIM>::GetUnboundProteinDiffusionCoefficient()
{
    return mUnboundProteinDiffusionCoefficient;
}

template<unsigned DIM>
void AbEdgeTrackingModifier<DIM>::SetUnboundProteinDiffusionCoefficient(
    double unboundProteinDiffusionCoefficient)
{
    mUnboundProteinDiffusionCoefficient = unboundProteinDiffusionCoefficient;
}

template<unsigned DIM>
void AbEdgeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AbEdgeTrackingModifier<1>;
template class AbEdgeTrackingModifier<2>;
template class AbEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(AbEdgeTrackingModifier)