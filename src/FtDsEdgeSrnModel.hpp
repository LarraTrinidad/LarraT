#ifndef FTDSEDGESRNMODEL_HPP_
#define FTDSEDGESRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "FtDsEdgeOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

class FtDsEdgeSrnModel : public AbstractOdeSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:

    /**
     * Protected copy-constructor for use by CreateSrnModel(). The only way for 
     * external code to create a copy of a SRN model is by calling that method, 
     * to ensure that a model of the correct subclass is created. This 
     * copy-constructor helps subclasses to ensure that all member variables are 
     * correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a 
     * daughter cell upon cell division. Note that the parent SRN model will 
     * have had ResetForDivision() called just before CreateSrnModel() is 
     * called, so performing an exact copy of the parent is suitable behaviour. 
     * Any daughter-cell-specific initialisation can be done in 
     * InitialiseDaughterCell().
     *
     * @param rModel  the SRN model to copy.
     */
    FtDsEdgeSrnModel(const FtDsEdgeSrnModel& rModel);

public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver 
     *     object (allows the use of different ODE solvers)
     */
    FtDsEdgeSrnModel(
        boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    virtual AbstractSrnModel* CreateSrnModel() override;

    /**
     * Overridden Initialise() method to initialise the SRN model at the start 
     * of a simulation.
     */
    virtual void Initialise() override;

    /**
     * This method is called when a new edge is created (e.g. after cell 
     * division or T1 swap)
     */
    virtual void InitialiseDaughterCell() override;

    /**
     * Overridden SimulateToTime() method for custom behaviour.
     */
    virtual void SimulateToCurrentTime() override;

    /**
     * Update the protein concentrations at neighbouring edge sensed by this 
     * edge; that is, fetch neighbour values from CellEdgeData object, storing 
     * the sensed information, into this model.
     */
    void UpdateFtDs();

    /** @return the current Ds concentration at this edge. */
    double GetDs();

    /**
     * Set the Ds concentration at this edge.
     * 
     * @param value
     */
    void SetDs(double value);

    /** @return the current Ft concentration at this edge. */
    double GetFt();

    /**
     * Set the Ft concentration at this edge.
     * 
     * @param value
     */
    void SetFt(double value);

    /** @return the current DsP concentration at this edge. */
    double GetDsP();

    /**
     * Set the DsP concentration at this edge.
     * 
     * @param value
     */
    void SetDsP(double value);

    /** @return the current FtP concentration at this edge. */
    double GetFtP();

    /**
     * Set the FtP concentration at this edge.
     * 
     * @param value
     */
    void SetFtP(double value);

    /** @return the current complex A concentration at this edge. */
    double GetA();

    /**
     * Set the complex A concentration at this edge.
     * 
     * @param value
     */
    void SetA(double value);

    /** @return the current complex B concentration at this edge. */
    double GetB();

    /**
     * Set the complex B concentration at this edge.
     * 
     * @param value
     */
    void SetB(double value);

    /** @return the current complex C concentration at this edge. */
    double GetC();

    /**
     * Set the complex C concentration at this edge.
     * 
     * @param value
     */
    void SetC(double value);

    /** @return the current complex D concentration at this edge. */
   double GetD();

    /**
     * Set the complex D concentration at this edge.
     * 
     * @param value
     */
    void SetD(double value);

    /**
     * @return the current complex A concentration at the neighbouring cell's 
     *     edge.
     */
    double GetNeighA();

    /**
     * Set the complex A concentration at the neighbouring cell's edge.
     * 
     * @param value
     */
    void SetNeighA(double value);

    /**
     * @return the current complex B concentration at the neighbouring cell's 
     *     edge.
     */
    double GetNeighB();

    /**
     * Set the complex B concentration at the neighbouring cell's edge.
     * 
     * @param value
     */
    void SetNeighB(double value);

    /**
     * @return the current complex C concentration at the neighbouring cell's 
     *     edge.
     */
    double GetNeighC();

    /**
     * Set the complex C concentration at the neighbouring cell's edge.
     * 
     * @param value
     */
    void SetNeighC(double value);

    /**
     * @return the current complex D concentration at the neighbouring cell's 
     *     edge.
     */
   double GetNeighD();

    /**
     * Set the complex D concentration at the neighbouring cell's edge.
     * 
     * @param value
     */
    void SetNeighD(double value);

    /**
     * Overridden OutputSrnModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile) override;

    /**
     * Overridden AddSrnQuantities() method.
     *
     * @param pOtherSrn
     * @param scale
     */
    virtual void AddSrnQuantities(AbstractSrnModel* pOtherSrn,
                                  const double scale = 1.0) override;

    /**
     * Overridden AddShrunkEdgeSrn() method.
     * 
     * @param pShrunkEdgeSrn
     */
    virtual void AddShrunkEdgeSrn(AbstractSrnModel* pShrunkEdgeSrn) override;

    /**
     * Overridden AddMergedEdgeSrn() method.
     * 
     * @param pMergedEdgeSrn
     */
    virtual void AddMergedEdgeSrn(AbstractSrnModel* pMergedEdgeSrn) override;

    /**
     * Overridden SplitEdgeSrn() method.
     * 
     * @param relativePosition
     */
    virtual void SplitEdgeSrn(const double relativePosition) override;
};

typedef boost::shared_ptr<FtDsEdgeSrnModel> FtDsEdgeSrnModelPtr;

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(FtDsEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(FtDsEdgeSrnModel)

#endif  /* FTDSEDGESRNMODEL_HPP_ */
