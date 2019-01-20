//
// Created by gq on 1/15/19.
//

#ifndef CELLCORSONWRITER_HPP_
#define CELLCORSONWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing to file, for each cell,
 * the level of delta, notch and mean level of delta among neighbouring cells.
 *
 * The output file is called celldeltanotch.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Cell delta" by default.
 *
 * Note: if you use a DeltaNotchSrnModel then the delta and notch levels are
 * stored in CellData, and thus (if VTK is switched on) will be output as VTK cell
 * data already.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellCorsonWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
    CellCorsonWriter();

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get a double associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden VisitCell() method.
     *
     * Visit a cell and write its delta and notch data.
     *
     * Outputs a line of space-separated values of the form:
     * ...[location index] [cell id] [x-pos] [y-pos] [z-pos] [delta] [notch] [mean neighbouring delta]...
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively.
     *
     * This is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCorsonWriter)
#endif //CELLCORSONWRITER_HPP_
