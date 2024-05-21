#ifndef SIDEKICKBOUNDARYCONDITION_HPP_
#define SIDEKICKBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template<unsigned DIM>
class SidekickBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:
  
  /** Needed for serialization. */
  friend class boost::serialization::access;
  /**
   * Serialize the object.
   *
   * @param archive the archive
   * @param version the current version of this class
   */
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version)
  {
    archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
  }
  
public:
  
  /**
   * Constructor.
   *
   * @param pCellPopulation pointer to the cell population
   */
  explicit SidekickBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation);
  
  /**
   * Overridden ImposeBoundaryCondition() method.
   *
   * @param rOldLocations the node locations before any boundary conditions are applied
   */
  void ImposeBoundaryCondition(
      const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations) override;
  
  /**
   * Overridden VerifyBoundaryCondition() method.
   *
   * @return whether the boundary conditions are satisfied.
   */
  bool VerifyBoundaryCondition() override;
  
  /**
   * Overridden OutputCellPopulationBoundaryConditionParameters() method.
   *
   * @param rParamsFile the file stream to which the parameters are output
   */
  void OutputCellPopulationBoundaryConditionParameters(
      out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SidekickBoundaryCondition)
  
  namespace boost
  {
  namespace serialization
  {
  /**
   * Serialize information required to construct a SidekickBoundaryCondition.
   */
  template<class Archive, unsigned DIM>
  inline void save_construct_data(
      Archive & ar, const SidekickBoundaryCondition<DIM>* t, const unsigned int file_version)
  {
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
  }
  
  /**
   * De-serialize constructor parameters and initialize a SidekickBoundaryCondition.
   */
  template<class Archive, unsigned DIM>
  inline void load_construct_data(
      Archive & ar, SidekickBoundaryCondition<DIM>* t, const unsigned int file_version)
  {
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    
    // Invoke inplace constructor to initialise instance
    ::new(t)SidekickBoundaryCondition<DIM>(p_cell_population);
  }
  }
  } // namespace ...

#endif /*SIDEKICKBOUNDARYCONDITION_HPP_*/