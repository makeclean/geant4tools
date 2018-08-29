#include <iostream>
#include "moab/Core.hpp"
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"
#include "uwuw.hpp"
#include "xercesc/util/XMLString.hpp"
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include "xercesc/dom/DOM.hpp"

// ---------------------------------------------------------------------------
//  This is a simple class that lets us do easy (though not terribly efficient)
//  trancoding of char* data to
// ---------------------------------------------------------------------------
class XStr {
 public :
  // -----------------------------------------------------------------------
  //  Constructors and Destructor
  // -----------------------------------------------------------------------
  XStr(const char* const toTranscode) {
    // Call the private transcoding method
    fUnicodeForm = xercesc::XMLString::transcode(toTranscode);
  }
 
  ~XStr() {
    xercesc::XMLString::release(&fUnicodeForm);
  }
  
  // -----------------------------------------------------------------------
  //  Getter methods
  // -----------------------------------------------------------------------
  const XMLCh* unicodeForm() const {
    return fUnicodeForm;
  }
  private :
  // -----------------------------------------------------------------------
  //  Private data members
  //
  //  fUnicodeForm
  //      This is the Unicode XMLCh format of the string.
  // -----------------------------------------------------------------------
  XMLCh*   fUnicodeForm;
};

#define X(str) XStr(str).unicodeForm()

// globals
xercesc::DOMDocument* gdml;
moab::DagMC *dag;

// convert double to string including scientific
// formatting
std::string to_string(double value) {
  std::ostringstream sstream;
  sstream.precision(16);
  sstream << std::scientific;
  sstream << value;
  std::string string_double = sstream.str();
  return string_double;
}

void write_xml(std::string output_filename) {
  xercesc::DOMImplementation *impl = xercesc::DOMImplementationRegistry::getDOMImplementation(X("LS"));
  if(impl == NULL) {
    std::cout << "Couldnt access DOM LS writer" << std::endl;
    return;
  }
  xercesc::DOMLSSerializer *serialiser = ((xercesc::DOMImplementationLS*)impl)->createLSSerializer();
  //
  xercesc::DOMConfiguration* config = serialiser->getDomConfig();
  config->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);

  xercesc::XMLFormatTarget *file = new xercesc::LocalFileFormatTarget(X(output_filename.c_str()));
  // Write the serialized output to the target.
  xercesc::DOMLSOutput* lsoutput = ((xercesc::DOMImplementationLS*)impl)->createLSOutput();
  lsoutput->setByteStream(file);

  serialiser->write(gdml,lsoutput);
}

// make the materials
void make_materials(std::string filename) {
  xercesc::DOMElement* rootElem = gdml->getDocumentElement();
  xercesc::DOMElement* materialsElem = gdml->createElement(X("materials"));
  rootElem->appendChild(materialsElem);
 
  //
  UWUW *uwuw = new UWUW(filename);
  // get the materials used in the problem
  std::map<std::string,pyne::Material> mats = uwuw->material_library;
  std::map<std::string,pyne::Material>::iterator it;
  pyne::Material mat;
  // first make the nuclides
  for ( it = mats.begin() ; it != mats.end() ; it++ ) {
    mat = mat + it->second;
  }
  
  // this summed material has every nuclide present in the
  // problem
  std::map<int,double> nuclides = mat.comp;
  std::map<int,double>::iterator nuc;
  // loop over the consituents and write each to the xml file
  for ( nuc = nuclides.begin() ; nuc != nuclides.end() ; nuc++ ) {
    xercesc::DOMElement* nucElem = gdml->createElement(X("isotope"));
    int nucid = nuc->first;
    int N = pyne::nucname::anum(nucid);
    int Z = pyne::nucname::znum(nucid);
    std::string name = pyne::nucname::name(nucid)+"_i";
    nucElem->setAttribute(X("N"),X(std::to_string(N).c_str()));
    nucElem->setAttribute(X("Z"),X(std::to_string(Z).c_str()));
    nucElem->setAttribute(X("name"),X(name.c_str()));
    double comp = nuc->second;
    materialsElem->appendChild(nucElem);

    xercesc::DOMElement* atomElem = gdml->createElement(X("atom"));
    double mass = pyne::atomic_mass(nucid);
    atomElem->setAttribute(X("unit"),X("g/mole"));
    atomElem->setAttribute(X("value"),X(std::to_string(mass).c_str()));
    nucElem->appendChild(atomElem);
  }
  
  // nuclides all done now define an element for each nuclide 
  // loop over the consituents and write each to the xml file
  for ( nuc = nuclides.begin() ; nuc != nuclides.end() ; nuc++ ) {
    // make an element entry for each 
    xercesc::DOMElement* nucElement = gdml->createElement(X("element"));
    int nucid = nuc->first;
    std::string name = pyne::nucname::name(nucid);
    materialsElem->appendChild(nucElement);
    nucElement->setAttribute(X("name"),X(name.c_str()));
    xercesc::DOMElement* nuclide = gdml->createElement(X("fraction"));
    nuclide->setAttribute(X("n"),X("1.0"));
    std::string ref = name + "_i";
    nuclide->setAttribute(X("ref"),X(ref.c_str()));
    nucElement->appendChild(nuclide);
  }  
   
  // loop over the materials
  for ( it = mats.begin() ; it != mats.end() ; it++ ) {
    pyne::Material mat = it->second;
    double dDensity = mat.density;
    // convert to atom fraction - this is what gdml needs
    mat = mat.to_atom_frac();
    // make the material
    xercesc::DOMElement* material = gdml->createElement(X("material"));
    material->setAttribute(X("name"), X((it->first).c_str()));
    material->setAttribute(X("state"),X("solid"));
    materialsElem->appendChild(material);

    xercesc::DOMElement* density = gdml->createElement(X("D"));
    density->setAttribute(X("unit"),X("g/cm3"));
    density->setAttribute(X("value"), X(std::to_string(dDensity).c_str()));
    material->appendChild(density);
    
    std::map<int,double> nuclides = mat.comp;
    std::map<int,double>::iterator nuc;
    // loop over the consituents and write each to the xml file
    for ( nuc = nuclides.begin() ; nuc != nuclides.end() ; nuc++ ) {
      int nucid = nuc->first;
      double comp = nuc->second;
      std::string name = pyne::nucname::name(nucid);
      if ( comp > 0.0 ) {
     	  xercesc::DOMElement* fraction = gdml->createElement(X("fraction"));
	   	  fraction->setAttribute(X("n"),X(to_string(comp).c_str()));
	   	  fraction->setAttribute(X("ref"), X(name.c_str()));
	   	  material->appendChild(fraction);
      }
    }
  }
  
  // make a vacuum material
  xercesc::DOMElement* material = gdml->createElement(X("material"));
  material->setAttribute(X("name"), X("mat:Vacuum"));
  material->setAttribute(X("state"),X("gas"));
  materialsElem->appendChild(material);
  // set the density to be 1e-30 g/cc
  xercesc::DOMElement* density = gdml->createElement(X("D"));
  density->setAttribute(X("unit"),X("g/cm3"));
  density->setAttribute(X("value"), X("1.e-25"));
  material->appendChild(density);
  // set the material - need to do something more clever here
  xercesc::DOMElement* fraction = gdml->createElement(X("fraction"));
  fraction->setAttribute(X("n"),X("1.0"));
  fraction->setAttribute(X("ref"), X("H1"));
  material->appendChild(fraction);
  
  
  // materials all done
  delete uwuw;
}

// write the geometry into the gdml file
void make_vertices(moab::Core *mbi, std::string filename) {
  // get the vertices
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> vertices;
  rval = mbi->get_entities_by_type(0,moab::MBVERTEX,vertices);
  double coord[3];
  
  xercesc::DOMElement* rootElem = gdml->getDocumentElement();
  xercesc::DOMElement* define = gdml->createElement(X("define"));
  rootElem->appendChild(define);
  
  for ( moab::EntityHandle vertex_handle : vertices ) {
    xercesc::DOMElement* position = gdml->createElement(X("position"));
    rval = mbi->get_coords(&vertex_handle,1,coord);
    position->setAttribute(X("name"),X(("v_"+std::to_string(vertex_handle)).c_str()));
    position->setAttribute(X("unit"),X("cm"));
    position->setAttribute(X("x"),X(to_string(coord[0]).c_str()));
    position->setAttribute(X("y"),X(to_string(coord[1]).c_str()));
    position->setAttribute(X("z"),X(to_string(coord[2]).c_str()));
    define->appendChild(position);
  }
  
}

// get the bounds of a given volume
void get_bounds(moab::Core *mbi, moab::EntityHandle vol, double box[3]) {
  moab::Range surfaces; // child surfaces
  moab::ErrorCode rval = mbi->get_child_meshsets(vol,surfaces);
  //
  moab::Range vertices;
  moab::Range triangles;
  //
  for ( moab::EntityHandle surface : surfaces ) {
    std::cout << triangles.size() << std::endl;
    rval = mbi->get_entities_by_type(surface, moab::MBTRI, triangles);
  }
  //
  rval = mbi->get_connectivity(triangles,vertices);
  //
  double coord[3];
  for ( moab::EntityHandle vertex : vertices ) {
    std::cout << "wo" << std::endl;
    rval = mbi->get_coords(&vertex,1,coord);
    box[0] = std::max(std::abs(coord[0]),std::abs(box[0]));
    box[1] = std::max(std::abs(coord[1]),std::abs(box[1]));
    box[2] = std::max(std::abs(coord[2]),std::abs(box[2]));
  }
  std::cout << box[0] << " " << box[1] << " " << box[2] << std::endl;
  return;
}

// write the geometry into the gdml file
void make_solids(moab::Core *mbi, std::string filename) {
  moab::ErrorCode rval;
  dagmcMetaData *dmd = new dagmcMetaData(dag);
  dmd->load_property_data();

  xercesc::DOMElement* rootElem = gdml->getDocumentElement();
  xercesc::DOMElement* solids = gdml->createElement(X("solids"));
  rootElem->appendChild(solids);

  // for each volume
  int dim = 3;
  const void *val = &dim;
  moab::Tag geom_dim_tag, id_tag;
  rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geom_dim_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  rval = mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  // get the volumes
  moab::Range volumes;
  rval = mbi->get_entities_by_type_and_tag(0,moab::MBENTITYSET ,&geom_dim_tag, &val, 1, volumes);

  // need to make the world volume
  // by looping over the set of volumes find the graveyard
  // and use that to make a box of the right size
  for ( moab::EntityHandle vol : volumes ) {
    std::string mat = dmd->volume_material_data_eh[vol];
    if (mat == "Graveyard") {
      double max[3];
      get_bounds(mbi,vol,max);
      std::cout << max[0] << " " << max[1] << " " << max[2] << std::endl;
      double x = max[0];
      double y = max[1];
      double z = max[2];
      xercesc::DOMElement* box = gdml->createElement(X("box"));
      solids->appendChild(box);
      box->setAttribute(X("name"),X("world_box"));
      box->setAttribute(X("x"),X(to_string(x).c_str()));
      box->setAttribute(X("y"),X(to_string(y).c_str()));
      box->setAttribute(X("z"),X(to_string(z).c_str()));
    }
  }
  
  // this loop makes the tessellated solids
  // for each vol
  for ( moab::EntityHandle vol : volumes ) {
    std::string mat = dmd->volume_material_data_eh[vol];
    if ( mat == "Vacuum" ) {
      std::cout << "Volume with material Vacuum being ignored" << std::endl;
      continue;
    }
    if ( mat == "Graveyard" ) {
      std::cout << "Volume with material Graveyard being ignored" << std::endl;
      continue;      
    }
    // get the child surfaces
    moab::Range surfaces;    
    xercesc::DOMElement* tessellated = gdml->createElement(X("tessellated"));
    solids->appendChild(tessellated);
    tessellated->setAttribute(X("aunit"),X("degrees"));
    tessellated->setAttribute(X("lunit"),X("cm"));
    // id of the volume
    int id;
    std::vector<moab::EntityHandle> volume;
    volume.push_back(vol);
    rval = mbi->tag_get_data(id_tag,&volume[0],1,&(id));
    std::string vol_name = "volume_"+std::to_string(id);
    tessellated->setAttribute(X("name"),X(vol_name.c_str()));
    rval = mbi->get_child_meshsets(vol,surfaces);
    // for each surface
    for ( moab::EntityHandle surface : surfaces ) {
      moab::Range triangles;
      int sense;
      rval = dag->surface_sense(vol, surface, sense);
      
      rval = mbi->get_entities_by_type(surface, moab::MBTRI, triangles);
      // for each triangle
      for ( moab::EntityHandle triangle : triangles ) {
	moab::Range tri;      
	moab::Range vertices;
	tri.insert(triangle);
	rval = mbi->get_vertices(tri,vertices);
	//std::cout << vertices.size() << std::endl;
	assert(vertices.size() == 3);
	xercesc::DOMElement* triangular = gdml->createElement(X("triangular"));
	if ( sense != 1 ) { // forward sense
	  triangular->setAttribute(X("vertex1"),X(("v_"+std::to_string(vertices[0])).c_str()));
	  triangular->setAttribute(X("vertex2"),X(("v_"+std::to_string(vertices[1])).c_str()));
	  triangular->setAttribute(X("vertex3"),X(("v_"+std::to_string(vertices[2])).c_str()));
	} else { // reverse sense - flip normals
	  triangular->setAttribute(X("vertex1"),X(("v_"+std::to_string(vertices[2])).c_str()));
	  triangular->setAttribute(X("vertex2"),X(("v_"+std::to_string(vertices[1])).c_str()));
	  triangular->setAttribute(X("vertex3"),X(("v_"+std::to_string(vertices[0])).c_str()));
	}
	tessellated->appendChild(triangular);
	tri.clear();
	vertices.clear();
      }
    }
  }

  xercesc::DOMElement* structure = gdml->createElement(X("structure"));
  rootElem->appendChild(structure);
  
  // this loop makes the structure section
  for ( moab::EntityHandle vol : volumes ) {
    std::string mat = dmd->volume_material_data_eh[vol];
    if ( mat == "Vacuum" || mat == "Graveyard" ) {
      continue;      
    }
    xercesc::DOMElement* volume = gdml->createElement(X("volume"));
    structure->appendChild(volume);

    // get the id to use as the name
    int id;
    std::vector<moab::EntityHandle> v;
    v.push_back(vol);
    // 
    rval = mbi->tag_get_data(id_tag,&v[0],1,&(id));
    std::string vol_name = "logical_"+std::to_string(id);
    volume->setAttribute(X("name"),X(vol_name.c_str()));
    // volume 
    xercesc::DOMElement* materialref = gdml->createElement(X("materialref"));
    std::string mat_name = "mat:"+mat;
    materialref->setAttribute(X("ref"),X(mat_name.c_str()));
    volume->appendChild(materialref);
    xercesc::DOMElement* solidref = gdml->createElement(X("solidref"));
    std::string name = "volume_" + std::to_string(id);
    solidref->setAttribute(X("ref"),X(name.c_str()));
    volume->appendChild(solidref);
  }

  // make the world volume
  xercesc::DOMElement* volume = gdml->createElement(X("volume"));
  structure->appendChild(volume);
  volume->setAttribute(X("name"),X("WorldVolume"));
  xercesc::DOMElement* material_ref = gdml->createElement(X("materialref"));
  volume->appendChild(material_ref);
  material_ref->setAttribute(X("ref"),X("mat:Vacuum"));
  xercesc::DOMElement* solid_ref = gdml->createElement(X("solidref"));
  solid_ref->setAttribute(X("ref"),X("world_box"));
  volume->appendChild(solid_ref);
  // now make all tesellated solids children of this

  // this loop makes the structure section
  for ( moab::EntityHandle vol : volumes ) {
    std::string mat = dmd->volume_material_data_eh[vol];
    if ( mat == "Vacuum" || mat == "Graveyard" ) {
      continue;      
    }
    xercesc::DOMElement* physvol = gdml->createElement(X("physvol"));
    volume->appendChild(physvol);

    xercesc::DOMElement* volumeref = gdml->createElement(X("volumeref"));
    physvol->appendChild(volumeref);
    
    // get the id to use as the name
    int id;
    std::vector<moab::EntityHandle> v;
    v.push_back(vol);
    // 
    rval = mbi->tag_get_data(id_tag,&v[0],1,&(id));
    std::string vol_name = "logical_"+std::to_string(id);
    volumeref->setAttribute(X("ref"),X(vol_name.c_str()));
  }

  
  // make the setup
  xercesc::DOMElement* setup = gdml->createElement(X("setup"));
  rootElem->appendChild(setup);
  setup->setAttribute(X("name"),X("Default"));
  setup->setAttribute(X("version"),X("1.0"));

  xercesc::DOMElement* world = gdml->createElement(X("world"));
  setup->appendChild(world);
  world->setAttribute(X("ref"),X("WorldVolume"));
}

int main(int argc, char* argv[]) {
  // start up xml
  xercesc::XMLPlatformUtils::Initialize();
  xercesc::DOMImplementation* impl =  xercesc::DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if ( impl == NULL) {
    std::cout << "Unable to get the requested xml implementation" << std::endl;
    return -1;
  }
  if ( argc != 2 ) return 0;
  std::string filename(argv[1]);
  
  // create the xml document
  gdml = impl->createDocument(0,X("gdml"),0);

  // load file
  moab::ErrorCode rval = moab::MB_FAILURE;
  moab::Core *mbi = new moab::Core();
  dag = new moab::DagMC(mbi);
  rval = dag->load_file(filename.c_str());
  rval = dag->setup_indices();
  /*
  return 0;

  moab::EntityHandle input_set;
  rval = mbi->create_meshset(moab::MESHSET_SET, input_set); 
  //MB_CHK_SET_ERR(rval, "failed to create meshset");
  std::cout << "Loading input file..." << std::endl;
  //rval = mbi->load_file(filename.c_str(), &input_set);  
  */
  
  make_vertices(mbi,filename);
  // load uwuw material
  make_materials(filename);
  // load tesellated solids
  make_solids(mbi,filename);


  // write the xml file

  // make the outptu file
  std::string output_file("test.gdml");
  write_xml(output_file);
  std::cout << "All done :) " << std::endl;
  
  // all done xml
  gdml->release();
  xercesc::XMLPlatformUtils::Terminate();
  delete mbi;
  return 0;
}
