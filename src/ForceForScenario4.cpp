#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "ForceForScenario4.hpp"
#include "StressTensor.hpp"

template<unsigned DIM>
ForceForScenario4<DIM>::ForceForScenario4()
  : FarhadifarForce<DIM>(),
    mHomotypicLineTensionParameter(1.0),
    mHeterotypicLineTensionParameter(1.0),
    mSupercontractileLineTensionParameter(1.0),
    mNumStripes(4),
    mUseCombinedInterfacesForLineTension(false),
    mUseDistinctStripeMismatchesForCombinedInterfaces(false),
    mGradualContractilityIncreaseTime(1.0)
    {
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
double ForceForScenario4<DIM>::GetCombinedInterfaceLength(Node<DIM>* pNode,
                                                          unsigned elemIndex,
                                                          unsigned cell1StripeIdentity,
                                                          unsigned cell2StripeIdentity,
                                                          VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
  
  VertexElement<DIM, DIM>* p_elem = rVertexCellPopulation.GetElement(elemIndex);
  unsigned num_nodes = p_elem->GetNumNodes();
  
  Node<DIM>* p_prev_node = pNode;
  std::set<unsigned> prev_node_elems = p_prev_node->rGetContainingElementIndices();
  unsigned prev_local_idx = p_elem->GetNodeLocalIndex(pNode->GetIndex());
  
  double combined_interface_length = 0.0;
  bool part_of_combined_interface = true;
  
  // We include this condition in the while loop to avoid an infinite loop in
  // the case of a cell surrounded by cells of different stripe identities
  unsigned num_prev_nodes_considered = 0;
  
  while (part_of_combined_interface && (num_prev_nodes_considered < num_nodes))
  {
    unsigned this_node_index = p_prev_node->GetIndex();
    
    std::set<unsigned> this_node_elems = prev_node_elems;
    prev_local_idx = (prev_local_idx - 1 + num_nodes)%num_nodes;
    
    p_prev_node = p_elem->GetNode(prev_local_idx);
    prev_node_elems = p_prev_node->rGetContainingElementIndices();
    
    part_of_combined_interface = false;
    
    std::set<unsigned> shared_elems;
    std::set_intersection(prev_node_elems.begin(), prev_node_elems.end(), this_node_elems.begin(), this_node_elems.end(), std::inserter(shared_elems, shared_elems.begin()));
    
    if (shared_elems.size() == 2)
    {
      CellPtr p_cell_1 = rVertexCellPopulation.GetCellUsingLocationIndex(*(shared_elems.begin()));
      CellPtr p_cell_2 = rVertexCellPopulation.GetCellUsingLocationIndex(*(++(shared_elems.begin())));
      unsigned cell_1_stripe_identity = p_cell_1->GetCellData()->GetItem("stripe");
      unsigned cell_2_stripe_identity = p_cell_2->GetCellData()->GetItem("stripe");
      
      if (mUseDistinctStripeMismatchesForCombinedInterfaces)
      {
        // Only continue if the neighbouring cells have the same stripes (e.g. 1 and 2) as the interface for which we are calling this method
        if (   ((cell_1_stripe_identity == cell1StripeIdentity) && (cell_2_stripe_identity == cell2StripeIdentity))
                 || ((cell_1_stripe_identity == cell2StripeIdentity) && (cell_2_stripe_identity == cell1StripeIdentity)) )
        {
          part_of_combined_interface = true;
          combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(p_prev_node->GetIndex(), this_node_index);
        }
      }
      else if (cell_1_stripe_identity != cell_2_stripe_identity)
      {
        part_of_combined_interface = true;
        combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(p_prev_node->GetIndex(), this_node_index);
      }
    }
    
    num_prev_nodes_considered++;
  }
  
  Node<DIM>* p_next_node = pNode;
  std::set<unsigned> next_node_elems = p_next_node->rGetContainingElementIndices();
  unsigned next_local_idx = p_elem->GetNodeLocalIndex(pNode->GetIndex());
  
  part_of_combined_interface = true;
  
  // We include this condition in the while loop to avoid an infinite loop in
  // the case of a cell surrounded by cells of different stripe identities
  unsigned num_next_nodes_considered = 0;
  
  while (part_of_combined_interface && (num_next_nodes_considered < num_nodes))
  {
    unsigned this_node_index = p_next_node->GetIndex();
    std::set<unsigned> this_node_elems = next_node_elems;
    
    next_local_idx = (next_local_idx + 1)%num_nodes;
    
    p_next_node = p_elem->GetNode(next_local_idx);
    next_node_elems = p_next_node->rGetContainingElementIndices();
    
    part_of_combined_interface = false;
    
    std::set<unsigned> shared_elems;
    std::set_intersection(this_node_elems.begin(), this_node_elems.end(), next_node_elems.begin(), next_node_elems.end(), std::inserter(shared_elems, shared_elems.begin()));
    
    if (shared_elems.size() == 2)
    {
      CellPtr p_cell_1 = rVertexCellPopulation.GetCellUsingLocationIndex(*(shared_elems.begin()));
      CellPtr p_cell_2 = rVertexCellPopulation.GetCellUsingLocationIndex(*(++(shared_elems.begin())));
      unsigned cell_1_stripe_identity = p_cell_1->GetCellData()->GetItem("stripe");
      unsigned cell_2_stripe_identity = p_cell_2->GetCellData()->GetItem("stripe");
      
      if (mUseDistinctStripeMismatchesForCombinedInterfaces)
      {
        // Only continue if the neighbouring cells have the same stripes (e.g. 1 and 2) as the interface for which we are calling this method
        if (   ((cell_1_stripe_identity == cell1StripeIdentity) && (cell_2_stripe_identity == cell2StripeIdentity))
                 || ((cell_1_stripe_identity == cell2StripeIdentity) && (cell_2_stripe_identity == cell1StripeIdentity)) )
        {
          part_of_combined_interface = true;
          combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(this_node_index, p_next_node->GetIndex());
        }
      }
      else if (cell_1_stripe_identity != cell_2_stripe_identity)
      {
        part_of_combined_interface = true;
        combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(this_node_index, p_next_node->GetIndex());
      }
    }
    num_next_nodes_considered++;
  }
  
  return combined_interface_length;
}

template<unsigned DIM>
double ForceForScenario4<DIM>::GetCombinedInterfaceScaleFactor(Node<DIM>* pNodeA,
                                                               Node<DIM>* pNodeB,
                                                               unsigned element1Index,
                                                               unsigned element2Index,
                                                               unsigned cell1StripeIdentity,
                                                               unsigned cell2StripeIdentity,
                                                               VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
  // Compute the total length of the 'boundary' interface that contains this edge, from each cell's perspective
  
  // Cell 1
  double interface_length_in_cell_1 = GetCombinedInterfaceLength(pNodeA, element1Index, cell1StripeIdentity, cell2StripeIdentity, rVertexCellPopulation);
  
  // Cell 2
  double interface_length_in_cell_2 = GetCombinedInterfaceLength(pNodeA, element2Index, cell1StripeIdentity, cell2StripeIdentity, rVertexCellPopulation);
  
  assert(interface_length_in_cell_1 > 0);
  assert(interface_length_in_cell_2 > 0);
  
  // Now compute the ratio of the longer to the short 'boundary' interfaces
  double shorter_interface = std::min(interface_length_in_cell_1, interface_length_in_cell_2);
  double scale_factor = 1.0/shorter_interface;
  
  return scale_factor;
}

template<unsigned DIM>
double ForceForScenario4<DIM>::GetLineTensionParameter(Node<DIM>* pNodeA,
                                                       Node<DIM>* pNodeB,
                                                       VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
  double line_tension = 0.0;
  
  // If either node is a boundary node, then use mBoundaryLineTensionParameter...
  if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
  {
    line_tension = this->mBoundaryLineTensionParameter;
  }
  else
  {
    // ...otherwise find the elements shared by these nodes
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(), elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(), elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    assert(!shared_elements.empty());
    
    // If the nodes share a single edge, then it must be on the boundary, so use mBoundaryLineTensionParameter...
    if (shared_elements.size() == 1)
    {
      line_tension = this->mBoundaryLineTensionParameter;
    }
    else
    {
      // ...otherwise proceed
      assert(shared_elements.size() == 2);
      unsigned elem_1_index = *(shared_elements.begin());
      unsigned elem_2_index = *(++(shared_elements.begin()));
      
      // Find the stripe identities on either side of the edge
      CellPtr p_cell_1 = rVertexCellPopulation.GetCellUsingLocationIndex(elem_1_index);
      CellPtr p_cell_2 = rVertexCellPopulation.GetCellUsingLocationIndex(elem_2_index);
      
      unsigned cell_1_stripe_identity = p_cell_1->GetCellData()->GetItem("stripe");
      unsigned cell_2_stripe_identity = p_cell_2->GetCellData()->GetItem("stripe");
      
      // If the edge is within a stripe (an 'internal' interface), then use mHomotypicLineTensionParameter
      if (cell_1_stripe_identity == cell_2_stripe_identity)
      {
        // The factor of 0.5 arises because each internal interface is visited twice in our iteration
        line_tension = 0.5*mHomotypicLineTensionParameter;
      }
      else
      {
        // Label numbers wrap around, so check to find smallest difference in stripe identities
        unsigned mismatch = abs(int(cell_1_stripe_identity - cell_2_stripe_identity));
        if (mismatch > mNumStripes/2)
        {
          mismatch = mNumStripes - mismatch;
        }
        
        // If the edge is between stripes whose identities differ by 1, then use mHeterotypicLineTensionParameter
        if (mismatch == 1)
        {
          // The factor of 0.5 arises because each internal interface is visited twice in our iteration
          line_tension = 0.5*mHeterotypicLineTensionParameter;
        }
        else
        {
          // The factor of 0.5 arises because each internal interface is visited twice in our iteration
          line_tension = 0.5*mSupercontractileLineTensionParameter;
        }
        
        if (mUseCombinedInterfacesForLineTension)
        {
          double scale_factor = GetCombinedInterfaceScaleFactor(pNodeA, pNodeB, elem_1_index, elem_2_index, cell_1_stripe_identity, cell_2_stripe_identity, rVertexCellPopulation);
          
          line_tension *= scale_factor;
        }
        
        // Scale scale_factor relative to simulation time (to simulate contractility gradually increasing)
        double sim_time = SimulationTime::Instance()->GetTime();
        if (sim_time < mGradualContractilityIncreaseTime)
        {
          double gradual_scaling = sim_time / mGradualContractilityIncreaseTime;
          line_tension *= gradual_scaling;
        }
      }
    }
  }
  
  return line_tension;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetGradualIncreaseContractilityTime(double gradualContractilityIncreaseTime)
{
  mGradualContractilityIncreaseTime = gradualContractilityIncreaseTime;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetHomotypicLineTensionParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
  mHomotypicLineTensionParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetHeterotypicLineTensionParameter(double labelledCellCellAdhesionEnergyParameter)
{
  mHeterotypicLineTensionParameter = labelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetSupercontractileLineTensionParameter(double supercontractileLineTensionParameter)
{
  mSupercontractileLineTensionParameter = supercontractileLineTensionParameter;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetNumStripes(unsigned numStripes)
{
  mNumStripes = numStripes;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetUseCombinedInterfacesForLineTension(bool useCombinedInterfaceLineTension)
{
  mUseCombinedInterfacesForLineTension = useCombinedInterfaceLineTension;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::SetUseDistinctStripeMismatchesForCombinedInterfaces(bool useDistinctStripeMismatchesForCombinedInterfaces)
{
  mUseDistinctStripeMismatchesForCombinedInterfaces = useDistinctStripeMismatchesForCombinedInterfaces;
}

template<unsigned DIM>
void ForceForScenario4<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
  // Output member variables then call method on direct parent class
  *rParamsFile << "\t\t\t<HomotypicLineTensionParameter>" << mHomotypicLineTensionParameter << "</HomotypicLineTensionParameter> \n";
  *rParamsFile << "\t\t\t<HeterotypicLineTensionParameter>" << mHeterotypicLineTensionParameter << "</HeterotypicLineTensionParameter> \n";
  *rParamsFile << "\t\t\t<NumStripes>" << mNumStripes << "</NumStripes> \n";
  *rParamsFile << "\t\t\t<UseCombinedInterfacesForLineTension>" << mUseCombinedInterfacesForLineTension << "</UseCombinedInterfacesForLineTension> \n";
  
  FarhadifarForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ForceForScenario4<1>;
template class ForceForScenario4<2>;
template class ForceForScenario4<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ForceForScenario4)