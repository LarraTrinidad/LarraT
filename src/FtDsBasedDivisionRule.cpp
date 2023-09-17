/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "FtDsBasedDivisionRule.hpp"
#include "FtDsEdgeTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellSrnModel.hpp"
#include "FtDsEdgeSrnModel.hpp"
#include "Debug.hpp"

/*template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> FtDsBasedDivisionRule<SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation)
{
    c_vector<double, SPACE_DIM> random_vector;
    double random_angle = 2.0*M_PI*0;
    random_vector(0) = cos(random_angle);
    random_vector(1) = sin(random_angle);
    return random_vector;
}*/
/*
template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::UpdateCellData(
    AbstractCellPopulation<DIM, DIM>& rCellPopulation)
    */


template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> FtDsBasedDivisionRule<SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation)
{
  c_vector<double, SPACE_DIM> division_vector;
  // For ease, store a static cast of the vertex-based cell population
    assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation));
    auto p_population = static_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);

    std::vector<double> Ft_levels;
    std::vector<double> edge_midpoint;

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
        ///std::vector<double> Ft_levels;
        ///std::vector<double> edge_midpoint;

        for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {  
            // Store the current Ft (bound and unbound) concentrations on this edge
            auto p_edge_srn = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
            Ft_levels.push_back(p_edge_srn->GetFt()+p_edge_srn->GetFtP()+(p_edge_srn->GetA()) + (p_edge_srn->GetB())+(p_edge_srn->GetC()) + (p_edge_srn->GetD()));

          // Store this edge's centroid
          auto p_element = p_population->GetElementCorrespondingToCell(*cell_iter);
          auto edge_midpoint = p_element->GetEdge(edge_index)->rGetCentreLocation();
          //edge_midpoint.push_back(edge_midpoint);
        }
        /*
        // Get the maximum level of Ft and which edge
        auto max_Ft_level = 0;
        for (unsigned edge = 0; edge < num_edges; ++edge)
        {
          if (Ft_levels(edge) > max_Ft_level)
          {
            max_Ft_level = Ft_levels(edge);
            double max_Ft_edge = edge;
          } 
        // edge's midpoint
        std::vector<double> max_Ft_edge_midpoint;
        max_Ft_edge_midpoint = edge_midpoint[edge];

        // cell's centroid
        std::vector<double< cell_centroid;
        cell_centroid = 
            */

        
          }

        double division_angle = 2.0*M_PI*0;
        division_vector(0) = cos(division_angle);
        division_vector(1) = sin(division_angle);
        return division_vector;


    }
    
  
// Explicit instantiation
template class FtDsBasedDivisionRule<1>;
template class FtDsBasedDivisionRule<2>;
template class FtDsBasedDivisionRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FtDsBasedDivisionRule)
