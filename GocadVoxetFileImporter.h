//---------------------------------------------------------------------------
#ifndef GocadVoxetFileImporterH
#define GocadVoxetFileImporterH
//---------------------------------------------------------------------------

#include "../IO/Gocad/GocadFileImporter.h"

#include "3DVector.h"

#include <vcl.h>

using namespace std;

namespace InSight
{

// Insight blockmesh index iteration order is (z,x,y)
enum GridAxis { W_AXIS = 0, U_AXIS, V_AXIS };

struct IGocadGrid
{
    IGocadGrid();
    I3DVector origin; // this is the origin of the voxet coord system

    // vectors defining the voxet coord system
    I3DVector axis_u;
    I3DVector axis_v;
    I3DVector axis_w;

    // beginning and end of the voxet extents as fraction of above coord system
    // ({0,0,0} to {1,1,1})
    I3DVector axis_min;
    I3DVector axis_max;

    int nu, nv, nw;

    // verifies if the axes allow for import, and stores information necessary for
    // correct indexing
    bool fixAxes();

    // NOTE: the information gathered from a fixAxes call is used in the following functions
    // converts the ijk of the GOCAD grid and returns the correct blockmesh index
    int actualIJK(int u, int v, int w) const;
    I3DVector actualOrigin() const;
    I3DVector actualAxis(GridAxis axis) const;
    int actualN(GridAxis axis) const;

private:
    GridAxis _u_rank;
    GridAxis _v_rank;
    GridAxis _w_rank;
};

class IGocadVoxetFileImporterParams : public IGocadFileImporterParams
{
public:
    IGocadVoxetFileImporterParams();
    virtual ~IGocadVoxetFileImporterParams() {}

    virtual IParameterList* clone() const;

protected:
	 virtual void populateParams(const IContentsInfo* contents);
//	 virtual void setColorTable(const IGocadProperty& cur_prop);
};

struct IGocadVoxetFileContentsInfo : public IGocadFileContentsInfo
{
    virtual ~IGocadVoxetFileContentsInfo() {}

    IGocadGrid mGrid;

    /**
     * Returns true if the file contents are able to be imported (at least in part).
     */ // TODO implement instead of putting all the code in the data reading function
    //virtual bool areValid() const;
};

class IGocadVoxetFileImporter : public IGocadFileImporter
{
public:
    virtual ~IGocadVoxetFileImporter(){}

    virtual string getProgressCaption() const;
    const IGocadVoxetFileContentsInfo* getContentsInfo();
    const IGocadVoxetFileImporterParams* getParameters() const;

	void setGridGeometry(IMesh* mesh);

protected:
	IGocadVoxetFileImporterParams* getParameters();

    virtual void readContentsInfo();
    virtual void readData();
	virtual vector<string> dataStartStrings();
	virtual void setColorTable(const IGocadProperty& cur_prop);

};

class TActionImportGocadVoxet : public TBasicAction
{
public:
   DYNAMIC bool __fastcall Execute(void);
   __fastcall TActionImportGocadVoxet(TComponent* AOwner) : TBasicAction(AOwner) {}
};


} // namespace InSight

#endif

