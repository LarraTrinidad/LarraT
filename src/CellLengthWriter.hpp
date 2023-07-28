#ifndef CELLLENGTHWRITER_HPP_
#define CELLLENGTHWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing cell lengths to file.
 *
 * The output file is called celllength.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored 
 * inthe VTK cell data "Cell Length" by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellLengthWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
  /** Needed for serialization. */
  friend class boost::serialization::access;
  /**
   * Serialize the object and its member variables.
   *
   * @param archive the archive
   * @param version the current version of this class
   */
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version)
  {
    archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
  }
  
public:
  
  /**
   * Default constructor.
   */
  CellLengthWriter();
  
  /**
   * Overridden GetCellDataForVtkOutput() method.
   *
   * @param pCell a cell
   * @param pCellPopulation a pointer to the cell population owning the cell
   *
   * @return data associated with the cell
   */
  c_vector<double, SPACE_DIM> GetCellDataForVtkOutputAsVector(
      CellPtr pCell,
      AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
  
  /**
   * Overridden VisitCell() method.
   * 
   * @param pCell a cell
   * @param pCellPopulation a pointer to the cell population owning the cell
   */
  virtual void VisitCell(
      CellPtr pCell,
      AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellLengthWriter)
  
#endif /* CELLLENGTHWRITER_HPP_ */