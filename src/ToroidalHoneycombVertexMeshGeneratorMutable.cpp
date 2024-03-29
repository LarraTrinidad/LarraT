#include "ToroidalHoneycombVertexMeshGeneratorMutable.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"

ToroidalHoneycombVertexMeshGeneratorMutable::ToroidalHoneycombVertexMeshGeneratorMutable(
    unsigned numElementsAcross,
    unsigned numElementsUp,
    double cellRearrangementThreshold,
    double t2Threshold,
    double elementArea)
{
    // numElementsAcross and numElementsUp must be even for toroidal meshes
    assert(numElementsAcross > 1);
    assert(numElementsUp > 1);
    assert(numElementsAcross%2 == 0); ///\todo This should be an exception
    assert(numElementsUp%2 == 0); ///\todo This should be an exception
    
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);
    
    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;
    
    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index;
    
    // Create the nodes
    for (unsigned j=0; j<2*numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
            double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3.0);
            
            Node<2>* p_node = new Node<2>(node_index, false , x_coord, y_coord);
            nodes.push_back(p_node);
            node_index++;
        }
    }
    
    /*
    * Create the elements. The array node_indices contains the global node 
    indices from bottom, going anticlockwise.
    */
    for (unsigned j = 0; j < numElementsUp; ++j)
    {
        for (unsigned i = 0; i < numElementsAcross; ++i)
        {
            element_index = j*numElementsAcross + i;
            
            node_indices[0] = 2*j*numElementsAcross + i + 1*(j%2==1);
            node_indices[1] = node_indices[0] + numElementsAcross + 1*(j%2==0);
            node_indices[2] = node_indices[0] + 2*numElementsAcross + 1*(j%2==0);
            node_indices[3] = node_indices[0] + 3*numElementsAcross;
            node_indices[4] = node_indices[0] + 2*numElementsAcross - 1*(j%2==1);
            node_indices[5] = node_indices[0] + numElementsAcross - 1*(j%2==1);
            
            if (i == numElementsAcross-1) // on far right
            {
                node_indices[0] -= numElementsAcross*(j%2==1);
                node_indices[1] -= numElementsAcross;
                node_indices[2] -= numElementsAcross;
                node_indices[3] -= numElementsAcross*(j%2==1);
            }
            if (j == numElementsUp-1) // on far top
            {
                node_indices[2] -= 2*numElementsAcross*numElementsUp;
                node_indices[3] -= 2*numElementsAcross*numElementsUp;
                node_indices[4] -= 2*numElementsAcross*numElementsUp;
            }
            
            std::vector<Node<2>*> element_nodes;
            for (unsigned k = 0; k < 6; ++k)
            {
               element_nodes.push_back(nodes[node_indices[k]]);
            }
            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            elements.push_back(p_element);
        }
    }
    
    double mesh_width = numElementsAcross;
    double mesh_height = 1.5*numElementsUp/sqrt(3.0);
    
    mpMesh = new Toroidal2dVertexMeshWithMutableSize(mesh_width, mesh_height, nodes, elements, cellRearrangementThreshold, t2Threshold);
    
    // Scale the mesh so that each element's area takes the value elementArea
    if (elementArea != 0)
    {
        assert(elementArea > 0);
        
        mpMesh->Scale(sqrt(elementArea*2.0/sqrt(3.0)), sqrt(elementArea*2.0/sqrt(3.0)));
        
        Toroidal2dVertexMeshWithMutableSize* p_static_cast_mesh_toroidal = static_cast<Toroidal2dVertexMeshWithMutableSize*>(mpMesh);
        double height = p_static_cast_mesh_toroidal->GetWidth(1);
        double width = p_static_cast_mesh_toroidal->GetWidth(0);
        p_static_cast_mesh_toroidal->SetWidth(0, width * sqrt(elementArea*2.0/sqrt(3.0)));
        p_static_cast_mesh_toroidal->SetWidth(1, height * sqrt(elementArea*2.0/sqrt(3.0)));
        p_static_cast_mesh_toroidal->SetBoxCoords(1, width * sqrt(elementArea*2.0/sqrt(3.0)));
        p_static_cast_mesh_toroidal->SetBoxCoords(3, height * sqrt(elementArea*2.0/sqrt(3.0)));
    }
}

MutableVertexMesh<2,2>* ToroidalHoneycombVertexMeshGeneratorMutable::GetMesh()
{
  EXCEPTION("A toroidal mesh was created but a normal mesh is being requested.");
  return mpMesh; // Not really
}

Toroidal2dVertexMeshWithMutableSize* ToroidalHoneycombVertexMeshGeneratorMutable::GetMutableToroidalMesh()
{
  return (Toroidal2dVertexMeshWithMutableSize*) mpMesh;
}