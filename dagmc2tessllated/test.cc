#include <cmath>
#include <cassert>

#include <iostream>
#include "DagMC.hpp"
#include "DagSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4RandomDirection.hh"

int main(int argc, char* argv[]) {

  DagMC* dagmc = new moab::DagMC(); // create dag instance
  
  // dag_volumes
  const char* h5mfilename = argv[1];//"spheres.h5m";
  dagmc->load_file(h5mfilename);
  dagmc->init_OBBTree();

  G4VSolid *vol;
  bool dagsolid = std::atoi(argv[2]); //false;
  // new volume
  int id = 1;
  if(dagsolid) {
    vol = new DagSolid("vol_1", dagmc, id);
  } else {
    G4TessellatedSolid *gts = new G4TessellatedSolid();
    moab::Interface *mbi = dagmc->moab_instance();
    moab::EntityHandle volume = dagmc->entity_by_id(3,id);
    moab::Range child_surfs;
    moab::ErrorCode rval;
    rval = mbi->get_child_meshsets(volume,child_surfs);
    moab::Range::iterator it;
    // for each child surface get triangles
    moab::Range triangle_set, triangles;
    // we should probably check for the merge tag here, and invert normals
    // if sense needs to be reversed
    for ( it = child_surfs.begin() ; it != child_surfs.end() ; ++it ) {
      rval = mbi->get_entities_by_type(*it, moab::MBTRI, triangles);
      triangle_set.merge(triangles);
    }  
    for (moab::EntityHandle triangle : triangle_set) {
      std::vector<moab::EntityHandle> verts;
      rval = mbi->get_connectivity(&(triangle),1, verts);
      double pos1[3],pos2[3],pos3[3];
      rval = mbi->get_coords(&verts[0],1,&pos1[0]);
      rval = mbi->get_coords(&verts[1],1,&pos2[0]);
      rval = mbi->get_coords(&verts[2],1,&pos3[0]);
      G4TriangularFacet *facet = new G4TriangularFacet(
        G4ThreeVector(pos1[0],pos1[1],pos1[2]),
        G4ThreeVector(pos2[0],pos2[1],pos2[2]),
        G4ThreeVector(pos3[0],pos3[1],pos3[2]),
        ABSOLUTE);
      gts->AddFacet((G4VFacet*) facet);
    }
    gts->SetSolidClosed(true);
    vol = gts;
  }
  
  // start the clock
  auto t1 = std::chrono::high_resolution_clock::now();
  int count = 100000;
  for ( int i = 0 ; i < count ; i++ ) {
    G4ThreeVector direction = G4RandomDirection();
    int inside = vol->Inside(direction);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Time ";
  double time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << " total " << time << " us";
  std::cout << " per call " << time/(double)count;
  std::cout << " us" << std::endl;
  
  t1 = std::chrono::high_resolution_clock::now();
  G4ThreeVector pos = G4ThreeVector(0.,0.,0.);
  double distance = 0.;
  double distances = 0.;
  for ( int i = 0 ; i < count ; i++ ) {
    G4ThreeVector direction = G4RandomDirection();
    distance = vol->DistanceToOut(pos,direction);
    distances += distance;
  }

  std::cout << distances/(double)count << std::endl;;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "DisanceToIn" << std::endl;
  std::cout << "Time ";
  time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << " total " << time << " us";
  std::cout << " per call " << time/(double)count;
  std::cout << " us" << std::endl; 

  t1 = std::chrono::high_resolution_clock::now();
  G4ThreeVector wpos = G4ThreeVector(0.,0.,0.);
  distance = 0.;
  distances = 0.;
  for ( int i = 0 ; i < count ; i++ ) {
    G4ThreeVector direction = G4RandomDirection();
    distance = vol->DistanceToIn(pos,direction);
    distances += distance;
  }

  std::cout << distances/(double)count << std::endl;;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "DistanceToOut" << std::endl;
  std::cout << "Time ";
  time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << " total " << time << " us";
  std::cout << " per call " << time/(double)count;
  std::cout << " us" << std::endl; 

  t1 = std::chrono::high_resolution_clock::now();
  wpos = G4ThreeVector(0.,0.,0.);
  distance = 0.;
  distances = 0.;
  for ( int i = 0 ; i < count ; i++ ) {
    G4ThreeVector direction = G4RandomDirection();
    distance = vol->DistanceToIn(pos);
    distances += distance;
  }

  std::cout << distances/(double)count << std::endl;;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "DistanceToOut" << std::endl;
  std::cout << "Time ";
  time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << " total " << time << " us";
  std::cout << " per call " << time/(double)count;
  std::cout << " us" << std::endl; 
  t1 = std::chrono::high_resolution_clock::now();

  wpos = G4ThreeVector(0.,0.,0.);
  distance = 0.;
  distances = 0.;
  for ( int i = 0 ; i < count ; i++ ) {
    G4ThreeVector direction = G4RandomDirection();
    distance = vol->DistanceToOut(pos);
    distances += distance;
  }

  std::cout << distances/(double)count << std::endl;;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "DistanceToOut" << std::endl;
  std::cout << "Time ";
  time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << " total " << time << " us";
  std::cout << " per call " << time/(double)count;
  std::cout << " us" << std::endl; 
  
  return 0;
}
