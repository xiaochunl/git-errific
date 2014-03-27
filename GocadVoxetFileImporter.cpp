//---------------------------------------------------------------------------


#pragma hdrstop

#include <vcl.h>
#include "pch.h"
#include "GocadVoxetFileImporter.h"

#include "Universe.h"
#include "World.h"
#include "DataTypeList.h"
#include "DataType.h"
#include "Group.h"
#include "Mesh.h"
#include "Cell.h"
#include "BaseData.h"
#include "Vert.h"
#include "BlockMeshMesh.h"

#include <Types.h>
#include <ParsingAPI.h>
#include <ExceptionInterface.h>
#include "../Colors/ColorTable.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/format.hpp>
#include <algorithm>

//---------------------------------------------------------------------------

#pragma package(smart_init)

extern IUniverse *pUniverse;
using namespace std;

namespace InSight
{

// (9' .')9 ---------------[ IGocadGrid ]--------------- o-('. '9)

IGocadGrid::IGocadGrid()
: nu(0), nv(0), nw(0), _u_rank(U_AXIS), _v_rank(V_AXIS), _w_rank(W_AXIS) 
{}

bool IGocadGrid::fixAxes()
{
	I3DVector u_vec = axis_u;
    I3DVector v_vec = axis_v;
    I3DVector w_vec = axis_w;
    u_vec.Normalize();
    v_vec.Normalize();
    w_vec.Normalize();

	const float fZero = 0.01;
    const float PI = 3.14159265359;

    // find W axis    
	int nb_z_varying_axes = 0;
	if( abs(u_vec.z) > fZero)
	{
		nb_z_varying_axes++;
		_u_rank = W_AXIS;

        // look at angle between the remaining two to find U_AXIS and V_AXIS
        I3DVector v_vec_rotated = v_vec;
        v_vec_rotated.Rotate(I3DVector(0,0,-1), PI/2.0);
		if((v_vec_rotated-w_vec).Length() < fZero)
        {
            _v_rank = V_AXIS;
            _w_rank = U_AXIS;
		}
        else
		{
            I3DVector w_vec_rotated = w_vec;
			w_vec_rotated.Rotate(I3DVector(0,0,-1), PI/2.0);
            if((w_vec_rotated-v_vec).Length() < fZero)
            {
                _w_rank = V_AXIS;
                _v_rank = U_AXIS;
			}
			else
				return false;
		}
    }
	if( abs(v_vec.z) > fZero)
	{
		nb_z_varying_axes++;
		_v_rank = W_AXIS;

		I3DVector w_vec_rotated = w_vec;
		w_vec_rotated.Rotate(I3DVector(0,0,-1), PI/2.0);
		if((w_vec_rotated-u_vec).Length() < fZero)
		{
			_w_rank = V_AXIS;
			_u_rank = U_AXIS;
		}
		else
		{
			I3DVector u_vec_rotated = u_vec;
			u_vec_rotated.Rotate(I3DVector(0,0,-1), PI/2.0);
			if((u_vec_rotated-w_vec).Length() < fZero)
			{
                _u_rank = V_AXIS;
				_w_rank = U_AXIS;
			}
			else
				return false;
		}
    }
	if( abs(w_vec.z) > fZero)
    {
		nb_z_varying_axes++;
        _w_rank = W_AXIS;

        I3DVector v_vec_rotated = v_vec;
		v_vec_rotated.Rotate(I3DVector(0,0,-1), PI/2.0);
		if((v_vec_rotated-u_vec).Length() < fZero)
        {
			_v_rank = V_AXIS;
			_u_rank = U_AXIS;
		}
		else
		{
			I3DVector u_vec_rotated = u_vec;
			u_vec_rotated.Rotate(I3DVector(0,0,-1), PI/2.0);
            if((u_vec_rotated-v_vec).Length() < fZero)
			{
				_u_rank = V_AXIS;
				_v_rank = U_AXIS;
			}
			else
				return false;
        }
    }

	if(nb_z_varying_axes != 1)
        return false; // only one axis can vary vertically

    return true;
}

int IGocadGrid::actualIJK(int u, int v, int w) const
{
	int ijk = 0;
	if(_u_rank==W_AXIS)
		ijk += u;
	else if(_u_rank==U_AXIS)
        ijk += u*actualN(W_AXIS);
    else // V_AXIS
        ijk += u*actualN(W_AXIS)*actualN(U_AXIS);

	if(_v_rank==W_AXIS)
		ijk += v;
    else if(_v_rank==U_AXIS)
        ijk += v*actualN(W_AXIS);
    else // V_AXIS
        ijk += v*actualN(W_AXIS)*actualN(U_AXIS);

	if(_w_rank==W_AXIS)
		ijk += w;
    else if(_w_rank==U_AXIS)
        ijk += w*actualN(W_AXIS);
    else // V_AXIS
        ijk += w*actualN(W_AXIS)*actualN(U_AXIS);

	return ijk;
}

I3DVector IGocadGrid::actualOrigin() const
{
    I3DVector actual_o = origin;
	actual_o += axis_u * axis_min.x;
    actual_o += axis_v * axis_min.y;
    actual_o += axis_w * axis_min.z;
    return actual_o;
}

I3DVector IGocadGrid::actualAxis(GridAxis axis) const
{
    if(_w_rank == axis)
        return axis_w * (axis_max.z-axis_min.z);
    else if(_v_rank == axis)
        return axis_v * (axis_max.y-axis_min.y);
    else
        return axis_u * (axis_max.x-axis_min.x);
}

int IGocadGrid::actualN(GridAxis axis) const
{
    if(_w_rank == axis)
        return nw;
    else if(_v_rank == axis)
        return nv;
    else
        return nu;
}


// ----------------- Parameters functions -----------------

IGocadVoxetFileImporterParams::IGocadVoxetFileImporterParams()
: IGocadFileImporterParams()
{
    IFileParameter* file_name_param = getParameter<IFileParameter>("File path");
    if(file_name_param != NULL)
    {
        file_name_param->allowExtension("*.vo", "GOCAD Voxet files");
        file_name_param->allowAllFiles();
    }
}

IParameterList* IGocadVoxetFileImporterParams::clone() const
{
    // no special treatment for now (no member pointers)
    IGocadVoxetFileImporterParams* params = new IGocadVoxetFileImporterParams(*this);

    return params;
}



//void IGocadVoxetFileImporterParams::setColorTable(const IGocadProperty& cur_prop)
void IGocadVoxetFileImporter::setColorTable(const IGocadProperty& cur_prop)
{

	IDataTypeList *dataTypes = pUniverse->ThisWorld()->DataTypeList();
	IDataType *dt = dataTypes->DataType(cur_prop.class_name);
	if (!dt) return;
	if (dt->PrimitiveType() != T_REFERENCED) return;
	if( cur_prop.cmap.size()==0) return;
	IColorTable& pdataColorTable = dt->ColorTable();
	IColorTable::EntryTable m_table;
	m_table.clear();
	for (int i = 0;	i < cur_prop.cmap.size(); i++)
	{
		IColorRecord c;
		unsigned int refIndex, cR, cG, cB;
		float fVal;
		if(i%4 == 0){
			 refIndex = atoi(cur_prop.cmap[i].c_str());
		}else if (i%4 == 1) {
			 fVal = atof(cur_prop.cmap[i].c_str());
			 cR = fVal*255;
		}else if (i%4 == 2) {
			 fVal = atof(cur_prop.cmap[i].c_str());
			 cG = fVal*255;
		}else if (i%4 == 3) {
			 fVal = atof(cur_prop.cmap[i].c_str());
			 cB = fVal*255;
			 int x=1;
			IColorRecord c = IColorRecord(RGB(cR, cG, cB));
			m_table.push_back(IColorTable::Entry(refIndex, c));
		}
	}
	pdataColorTable.SetTable(m_table);

}


void IGocadVoxetFileImporterParams::populateParams(const IContentsInfo* contents)
{
    const IGocadVoxetFileContentsInfo* info = dynamic_cast<const IGocadVoxetFileContentsInfo*>(contents);
    if(!info)
        throw IProcException("Invalid GOCAD Voxet file contents.");

	INewMeshParameter* mesh_param = new INewMeshParameter("BlockMesh");
	mesh_param->setMeshType(IBlockMeshMeshTypeID);

	addParameter(mesh_param);
	mesh_param->setGroup(getParentGroup());

	// use object name if specified in file attributes
    map<string,string>::const_iterator it = info->attributes.find("name");
	string mesh_name;
    if(it != info->attributes.end())
        mesh_name = it->second;

    // otherwise just use the file name
    if(mesh_name.empty())
		mesh_name = IParsingAPI::getBaseFileName(getFilePath());

	mesh_param->setMeshName(mesh_name);

	IIndex pi=1;
	for(IIndex i=0; i<info->properties.size(); i++)
	{
		const IGocadProperty& cur_prop = info->properties[i];

		boost::format fmt = boost::format("Property %s") % (pi++);
		INewDataSetParameter* cur_dsp = new INewDataSetParameter( fmt.str(),
		  A_CELL,
		  T_FLOAT,
		  M_BLOCK,
		  mesh_param );

        if(cur_prop.discrete)
			cur_dsp->setPrimitiveType(T_REFERENCED);

		string cur_prop_name = cur_prop.name;
		cur_dsp->setDatasetName(cur_prop_name);
		cur_dsp->setDataType(cur_prop.class_name);
//		setColorTable(cur_prop);
		addParameter(cur_dsp);
	}
}

// ----------------- Procedure functions ------------------

string IGocadVoxetFileImporter::getProgressCaption() const
{
	string caption = string("Importing GOCAD Voxet File \"") +
        IParsingAPI::getBaseFileName(getParameters()->getFilePath()) + "\"...";
    return caption;
}

const IGocadVoxetFileContentsInfo* IGocadVoxetFileImporter::getContentsInfo()
{
    return dynamic_cast<const IGocadVoxetFileContentsInfo*>(IImporter::getContentsInfo());
}

const IGocadVoxetFileImporterParams* IGocadVoxetFileImporter::getParameters() const
{
    return IProcedure::getParameters<IGocadVoxetFileImporterParams>();
}

IGocadVoxetFileImporterParams* IGocadVoxetFileImporter::getParameters()
{
    return IProcedure::getParameters<IGocadVoxetFileImporterParams>();
}

void SwapDW(void *dwords, int n)
{
	   unsigned *tmp;
	   unsigned tmp1;
	   tmp = (unsigned*)dwords;

	   for(int i = 0; i < n; i++)
	   {
			  tmp1 = tmp[i];
			  tmp[i] = (tmp1&0xff)<<24 |
						   (tmp1&0xff00)<<8 |
						   (tmp1&0xff0000)>>8 |
						   (tmp1&0xff000000)>>24;
	   }
}

void SwapIW(void *dwords, int n)
{
	   short *tmp;
	   short tmp1;
	   tmp = (short*)dwords;

	   for(int i = 0; i < n; i++)
	   {
			  tmp1 = tmp[i];
			  tmp[i] = (tmp1&0xff)<<8 |(tmp1&0xff00)>>8;
	   }
}

void readVector(const vector<string>& words, I3DVector& vec)
{
    if(words.size() < 4) return;
    try
    {
        vec.x = boost::lexical_cast<double>(words[1]);
        vec.y = boost::lexical_cast<double>(words[2]);
        vec.z = boost::lexical_cast<double>(words[3]);
    }
    catch(boost::bad_lexical_cast&) {}
}


void IGocadVoxetFileImporter::readContentsInfo()
{
    // CODECOPY (some code copied from IGocadFileImporter)

    // create correct ContentsInfo object (every override should do this)
    IGocadVoxetFileContentsInfo* contents_info = dynamic_cast<IGocadVoxetFileContentsInfo*>(_contents_info);
    if(contents_info==NULL)
    {
        delete _contents_info; // make sure we start with a valid contents info object
        contents_info = new IGocadVoxetFileContentsInfo();
        _contents_info = contents_info;
    }

    // read general file information
    IFileImporter::readContentsInfo();

    string cur_line;
    streampos pos = _file_stream.tellg();
    vector<string> data_start_strings = dataStartStrings();

    const int props_offset = 3;
    // load values from msh file
	std::string sProperty,  sPropertyFile, sPropertyType, sValue;
    int ncx = 0;
    int ncy = 0;
	int ncz = 0;
	vector<float> NDVvalues;
	std::string cmap_name="";
	int cmap_size=0;

    // read file header
    while(getline(_file_stream, cur_line))
    {
        vector<string> words = IParsingAPI::cleanSplit(cur_line);
        if( words.size() == 0 )
            continue;

        if( boost::starts_with(cur_line, "GOCAD") )
        {
            if(words.size() > 1)
                contents_info->object_type = words[1];
        }
        else if( boost::starts_with(cur_line, "HEADER") )
        {
            // read attributes until closing bracket
            while(getline(_file_stream, cur_line) && !boost::starts_with(cur_line, "}"))
            {
                words = IParsingAPI::cleanSplit(cur_line, ":");
                if(words.size() == 2)
                    contents_info->attributes[ words[0] ] = words[1];
            }
        }
        else if( boost::starts_with(cur_line, "HDR") )
        {
            // read single attribute line
            if(cur_line.size() > 4)
            {
                string rest_of_line = cur_line.substr(4);
                words = IParsingAPI::cleanSplit(rest_of_line, ":");
                if(words.size() == 2)
                    contents_info->attributes[ words[0] ] = words[1];
            }
        }
        else if( boost::starts_with(cur_line, "AXIS_O"))
        {
            readVector(words, contents_info->mGrid.origin);
        }
        else if( boost::starts_with(cur_line, "AXIS_U "))
        {
            readVector(words, contents_info->mGrid.axis_u);
        }
        else if( boost::starts_with(cur_line, "AXIS_V "))
        {
            readVector(words, contents_info->mGrid.axis_v);
        }
        else if( boost::starts_with(cur_line, "AXIS_W "))
        {
            readVector(words, contents_info->mGrid.axis_w);
        }
        else if( boost::starts_with(cur_line, "AXIS_MIN "))
        {
            readVector(words, contents_info->mGrid.axis_min);
        }
        else if( boost::starts_with(cur_line, "AXIS_MAX "))
        {
            readVector(words, contents_info->mGrid.axis_max);
        }
        else if( boost::starts_with(cur_line, "AXIS_N "))
        {
            if(words.size() < 4)
                continue;

            try
            {
                contents_info->mGrid.nu = boost::lexical_cast<int>(words[1]);
                contents_info->mGrid.nv = boost::lexical_cast<int>(words[2]);
                contents_info->mGrid.nw = boost::lexical_cast<int>(words[3]);
            }
            catch(boost::bad_lexical_cast&) {}
        }
		else if( boost::starts_with(cur_line, "PROPERTY "))
        {
            if(words.size() < 3) continue;
			IGocadProperty cur_prop;
            sValue = words[2];
			//remove quotation marks:
            sValue.erase(remove( sValue.begin(),
                sValue.end(), '\"' ),sValue.end());

			cur_prop.name = sValue;
			cmap_name = "";
			cmap_size = 0;
			std::string cmap_title="*colormap*";

            // read property_file
			while(getline(_file_stream, cur_line))
            {
                words = IParsingAPI::cleanSplit(cur_line);
//                if(words.size() < 3)
//                    continue;

				if( boost::starts_with(cur_line, "PROPERTY_CLASS "))
				{
					sValue = words[2];
					//remove quotation marks:
					sValue.erase(remove( sValue.begin(), sValue.end(), '\"' ),sValue.end());
					cur_prop.class_name = sValue;
				}
				else if( boost::starts_with(cur_line, "colormap:"))
				{
					sValue = words[0];
					sValue =  sValue.substr(sValue.find_first_of(":")+1,sValue.length()-1);
					cmap_name = sValue;
					cmap_title = "*colormap*" + cmap_name + "*colors:";
				}
				else if( boost::starts_with(cur_line, "*colormap*nbcolors:"))
				{
					sValue = words[0];
					sValue =  sValue.substr(sValue.find_first_of(":")+1,sValue.length()-1);
					cmap_size = atoi(sValue.c_str());
				}
				else if( boost::starts_with(cur_line, cmap_title))
				{
					sValue = words[0];
					sValue =  sValue.substr(sValue.find_first_of(":")+1,sValue.length()-1);
					words[0] = sValue;
					cur_prop.cmap = words;
//					cur_prop.cmap_name = sValue;
				}
				else if( boost::starts_with(cur_line, "PROPERTY_SUBCLASS "))
                {
                    cur_prop.discrete = (words[3] != "Float");
                }
                else if( boost::starts_with(cur_line, "PROP_NO_DATA_VALUE "))
                {
                    try
                    {
                        cur_prop.no_data_value = boost::lexical_cast<double>(words[2]);
                    } catch(boost::bad_lexical_cast&)   {}
                }
                else if( boost::starts_with(cur_line, "PROP_ESIZE "))
                {
                    try
                    {
                        cur_prop.bytes =  boost::lexical_cast<int>(words[2]);
                    } catch(boost::bad_lexical_cast&)   {}
                }
                else if( boost::starts_with(cur_line, "PROP_FILE "))
                {
                    if(words.size() < props_offset) continue;
					cur_prop.file_name = words[2];
                    break;
                }
            }

            contents_info->properties.push_back(cur_prop);
        }
    }

    contents_info->mGrid.fixAxes();
    contents_info->data_start_pos = _file_stream.tellg();
}

void IGocadVoxetFileImporter::readData()
{
	// gather parameters and contents information
	const IGocadVoxetFileContentsInfo* contents_info = getContentsInfo();
	IGocadVoxetFileImporterParams* params = getParameters();

     // create groups/objects
	INewMeshParameter* mesh_param = params->getParameter<INewMeshParameter>("BlockMesh");
	IMesh* mesh = mesh_param->createMesh();
    if(!mesh)
		throw IProcException("Unable to create mesh!");

    // create properties
	vector<IBaseData*> datasets = IProcedureAPI::createAllDataSets(params, mesh_param);

	// read data: populate cells and vertices

	_file_stream.seekg(contents_info->data_start_pos);

	int ncx = 0;
	int ncy = 0;
	int ncz = 0;

	setGridGeometry(mesh);
	vector<IGocadProperty> properties = contents_info->properties;

	std::vector<LPBaseData> DataChannel;
	DataChannel.reserve( properties.size() );

    // use the original nu/nv/nw for reading
	ncx = contents_info->mGrid.nu;
	ncy = contents_info->mGrid.nv;
	ncz = contents_info->mGrid.nw;


	for( unsigned int i=0; i < datasets.size(); i++ )
	{
        // [aw] Create a data object and store it for later channel processing ...
		LPBaseData pData = datasets[i];
//		setColorTable(&pData, properties[i]);
		DataChannel.push_back( (LPBaseDataFloat)pData );

		// [aw] Resize space for all blocks (since the order of Insight Blockmesh input is different from
		//      that of GOCAD output):
		if( pData->isBlockManaged() )
		{
			IBaseManageBlock* ibmb = dynamic_cast<IBaseManageBlock*>(pData);
			ibmb->Resize(ncx*ncy*ncz);
		}
	}

	IString sval;
	string strval;

    setProgressCount(properties.size()*ncz*ncy*ncx);
	// Populate Records
	for(IIndex ch=0; ch<properties.size(); ch++)
	{
		string file_name = params->getFilePath();
		int iEnd = file_name.find_last_of("/\\")+1;
		if(iEnd != string::npos){
			file_name = file_name.substr(0, iEnd);
		}

		setColorTable(properties[ch]);

		file_name = file_name + properties[ch].file_name;
		FILE *in = fopen( file_name.c_str() , "rb" );

		if(in == NULL) continue;

		ICellID icell = 0;

		for( IIndex iz= 0; iz < ncz; iz++ )
		{
			for( IIndex iy= 0; iy < ncy; iy++ )
			{
				int y = iy;
				for( IIndex ix = 0; ix < ncx; ix++)
				{
					icell = contents_info->mGrid.actualIJK(ix,iy,iz);
					setProgressPos(ix + iy*ncx + iz*ncx*ncy + ch*ncx*ncy*ncz);
					if(properties[ch].bytes==1)
					{
						unsigned char val1;
						size_t result = fread( &val1, properties[ch].bytes, 1, in);
						unsigned short int val = (unsigned short int) val1;

						if( val!=properties[ch].no_data_value){
							sval = IntToStr(val);
							setDataValue(icell, DataChannel[ch], sval.c_str());
						} else{
							sval = 0;
						}
					}
					else if(properties[ch].bytes==2)
                    {
						short val;
						size_t result = fread( &val, properties[ch].bytes, 1, in);
						SwapIW(&val,1);
						if( val!=properties[ch].no_data_value){
							sval = IntToStr(val);
                            setDataValue(icell, DataChannel[ch], sval.c_str());
						}
					}
					else
                    {
						float val;
						size_t result = fread( &val, properties[ch].bytes, 1, in);
						SwapDW(&val,1);
						if( val!=properties[ch].no_data_value){
							sval = FloatToStr(val);
                            setDataValue(icell, DataChannel[ch], sval.c_str());
						}
					}
				}
			}
		}

		fclose(in);
	}
}

void IGocadVoxetFileImporter::setGridGeometry(IMesh* mesh)
{
    const IGocadVoxetFileContentsInfo* contents_info = getContentsInfo();
	const IGocadGrid& mGrid = contents_info->mGrid;
	int ncx = mGrid.actualN(U_AXIS);
	int ncy = mGrid.actualN(V_AXIS);
	int ncz = mGrid.actualN(W_AXIS);

	// initialize special "Grid" property that contains the whole geometry
	LPBaseDataFloat pGridData = (LPBaseDataFloat)mesh->AddData(
		"Grid", "Grid", T_FLOAT, M_BLOCK, A_MESH);
    I3DVector orig = mGrid.actualOrigin();
	(*pGridData)[0] = orig.x;
	(*pGridData)[1] = orig.y;
	(*pGridData)[2] = orig.z;
	(*pGridData)[3] = ncx+1;
	(*pGridData)[4] = ncy+1;
	(*pGridData)[5] = ncz+1;
	(*pGridData)[6] = mGrid.actualAxis(U_AXIS).dipdir();

	float cell_size_x = mGrid.actualAxis(U_AXIS).Length()/ncx;
	float cell_size_y = mGrid.actualAxis(V_AXIS).Length()/ncy;
    // reading z directly, since must be upright and can be negative
	float cell_size_z = mGrid.actualAxis(W_AXIS).z/ncz;
	vector<float> xoffsets(ncx+1);
	vector<float> yoffsets(ncy+1);
	vector<float> zoffsets(ncz+1);

	xoffsets[0] = 0;
	for (int i = 1; i <= ncx; i++)
	{
		xoffsets[i] = cell_size_x + xoffsets[i-1];
	}

	yoffsets[0] = 0;
	for (int i = 1; i <= ncy; i++)
	{
		yoffsets[i] = cell_size_y + yoffsets[i-1];
	}

	zoffsets[0] = 0;
	for (int i = 1; i <= ncz; i++)
	{
		zoffsets[i] = cell_size_z + zoffsets[i-1];
	}

	std::vector<float>::iterator it_float;
	IIndex icell = BLOCKMESH_GRID_CELLS_INDEX_OFFSET;


	for( it_float=xoffsets.begin(); it_float!=xoffsets.end(); it_float++ )
		(*pGridData)[icell++] = *it_float;

	for( it_float=yoffsets.begin(); it_float!=yoffsets.end(); it_float++ )
		(*pGridData)[icell++] = *it_float;

	for( it_float=zoffsets.begin(); it_float!=zoffsets.end(); it_float++ )
		(*pGridData)[icell++] = *it_float;

	// create cells
	int ncell = ncx*ncy*ncz;
	mesh->ReserveCells(ncell);
}


vector<string> IGocadVoxetFileImporter::dataStartStrings()
{
	vector<string> data_start_strings;

    data_start_strings.push_back("ILINE");
    data_start_strings.push_back("VRTX");
    data_start_strings.push_back("PVRTX");
	data_start_strings.push_back("SEG");

    return data_start_strings;
}

//-----------------------------------------------------------------------------
bool __fastcall TActionImportGocadVoxet::Execute(void)
{
    InSightExceptions::Clear();

    // TODO : build dialog around IFileParameter
    boost::shared_ptr<TOpenDialog> dialog(new TOpenDialog(Application));

	dialog->Filter = "GOCAD Voxet files (*.vo)|*.vo|"
        "All files (*.*)|*.*|";
    dialog->FilterIndex = 0;
    dialog->Title = "Import GOCAD Voxet Data";

    // configure the dialog options
    TOpenOptions options = dialog->Options;
    options << ofAllowMultiSelect << ofFileMustExist << ofPathMustExist
        >> ofHideReadOnly;
    dialog->Options = options;

    if (!dialog->Execute())
    {
        return false;
    }

    for (int i = 0; i < dialog->Files->Count; i++)
    {
        IString filename(dialog->Files->Strings[i]);

        IGocadVoxetFileImporterParams* params = new IGocadVoxetFileImporterParams();
		params->setFilePath(filename);

		IGocadVoxetFileImporter importer;
		importer.call(params);
    }

    return true;
}

}
