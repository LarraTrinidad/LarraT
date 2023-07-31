#include "FtDsEdgeTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellSrnModel.hpp"
#include "FtDsEdgeSrnModel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FtDsEdgeTrackingModifier<DIM>::FtDsEdgeTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mUnboundProteinDiffusionCoefficient(0.03)
{
}

template<unsigned DIM>
FtDsEdgeTrackingModifier<DIM>::~FtDsEdgeTrackingModifier()
{
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::SetupSolve(
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
void FtDsEdgeTrackingModifier<DIM>::UpdateCellData(
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
        std::vector<double> Ds_old;
        std::vector<double> Ft_old;
        std::vector<double> DsP_old;
        std::vector<double> FtP_old;
        std::vector<double> edge_lengths;

        for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {  
            // Store the current unbound protein concentrations on this edge
            auto p_edge_srn = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
            Ds_old.push_back(p_edge_srn->GetDs());
            Ft_old.push_back(p_edge_srn->GetFt());
            DsP_old.push_back(p_edge_srn->GetDsP());
            FtP_old.push_back(p_edge_srn->GetFtP());

            // Store this edge's length
            auto p_element = p_population->GetElementCorrespondingToCell(*cell_iter);
            double edge_length = p_element->GetEdge(edge_index)->rGetLength();
            edge_lengths.push_back(edge_length);
        }

        /*
         * Update unbound protein concentrations based on a linear diffusive 
         * flux between neighbouring edges.
         */
        std::vector<double> Ds_new(num_edges);
        std::vector<double> Ft_new(num_edges);
        std::vector<double> DsP_new(num_edges);
        std::vector<double> FtP_new(num_edges);
        for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {
           unsigned prev_index = (edge_index == 0) ? num_edges - 1 : edge_index - 1;
           unsigned next_index = (edge_index == num_edges - 1) ? 0 : edge_index + 1;

           ///\todo consider validity of diffusive flux expression
           double dx = 0.5 * (edge_lengths[prev_index] + edge_lengths[edge_index]);
           double dt = SimulationTime::Instance()->GetTimeStep();
           Ds_new[edge_index] = Ds_old[edge_index] + D*(Ds_old[prev_index] - 2*Ds_old[edge_index] + Ds_old[next_index])*dt/(dx*dx);
           Ft_new[edge_index] = Ft_old[edge_index] + D*(Ft_old[prev_index] - 2*Ft_old[edge_index] + Ft_old[next_index])*dt/(dx*dx);
           DsP_new[edge_index] = DsP_old[edge_index] + D*(DsP_old[prev_index] - 2*DsP_old[edge_index] + DsP_old[next_index])*dt/(dx*dx);
           FtP_new[edge_index] = FtP_old[edge_index] + D*(FtP_old[prev_index] - 2*FtP_old[edge_index] + FtP_old[next_index])*dt/(dx*dx);

        }

        // Note: state variables must be in the same order as in FtDsOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge Ds", Ds_new);
        cell_iter->GetCellEdgeData()->SetItem("edge Ft", Ft_new);
        cell_iter->GetCellEdgeData()->SetItem("edge DsP", DsP_new);
        cell_iter->GetCellEdgeData()->SetItem("edge FtP", FtP_new);
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

        std::vector<double> neigh_mean_Ds(num_edges);
        std::vector<double> neigh_mean_Ft(num_edges);
        std::vector<double> neigh_mean_DsP(num_edges);
        std::vector<double> neigh_mean_FtP(num_edges);

        for (unsigned edge_index = 0; edge_index < num_edges; ++edge_index)
        {
            // Get neighbouring cell's values
            ///\todo is this correct?
            auto elem_neighbours = p_population->GetNeighbouringEdgeIndices(*cell_iter, edge_index);
            for (auto neighbour : elem_neighbours)
            {
                auto p_cell = p_population->GetCellUsingLocationIndex(neighbour.first);
                auto p_data = p_cell->GetCellEdgeData();
                std::vector<double> neighbour_Ds_vec = p_data->GetItem("edge Ds");
                std::vector<double> neighbour_Ft_vec = p_data->GetItem("edge Ft");
                std::vector<double> neighbour_DsP_vec = p_data->GetItem("edge DsP");
                std::vector<double> neighbour_FtP_vec = p_data->GetItem("edge FtP");
                
                neigh_mean_Ds[edge_index] += neighbour_Ds_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_Ft[edge_index] += neighbour_Ft_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_DsP[edge_index] += neighbour_DsP_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_FtP[edge_index] += neighbour_FtP_vec[neighbour.second] / elem_neighbours.size();
            }
        }

        cell_iter->GetCellEdgeData()->SetItem("neighbour Ds", neigh_mean_Ds);
        cell_iter->GetCellEdgeData()->SetItem("neighbour Ft", neigh_mean_Ft);
        cell_iter->GetCellEdgeData()->SetItem("neighbour DsP", neigh_mean_DsP);
        cell_iter->GetCellEdgeData()->SetItem("neighbour FtP", neigh_mean_FtP);
    }
}

template<unsigned DIM>
double FtDsEdgeTrackingModifier<DIM>::GetUnboundProteinDiffusionCoefficient()
{
    return mUnboundProteinDiffusionCoefficient;
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::SetUnboundProteinDiffusionCoefficient(
    double unboundProteinDiffusionCoefficient)
{
    mUnboundProteinDiffusionCoefficient = unboundProteinDiffusionCoefficient;
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class FtDsEdgeTrackingModifier<1>;
template class FtDsEdgeTrackingModifier<2>;
template class FtDsEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(FtDsEdgeTrackingModifier)