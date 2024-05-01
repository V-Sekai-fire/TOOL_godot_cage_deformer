#include <iostream>
#include "include/somig.h"
#include "include/BasicIO.h"

using namespace std;

int main(int argc, char *argv[])
{
  std::cout << " Load cage " << std::endl;
  // Load cage:
  std::vector< point3d > cage_vertices;
  std::vector< std::vector< unsigned int > > cage_triangles;
  OBJIO::open("../models/cage.obj", cage_vertices, cage_triangles, true); // TRIANGULATE THE FACES OF 

  std::cout << " Load mesh " << std::endl;
  // Load mesh:
  std::vector< point3d > mesh_vertices;
  std::vector< std::vector< unsigned int > > mesh_triangles;
  OBJIO::open("../models/mesh.obj", mesh_vertices, mesh_triangles, true);

  somig_deformer_3 dfm;
  dfm.set_mesh(mesh_triangles, mesh_vertices);
  dfm.set_cage(cage_triangles, cage_vertices);

  bool sanity_check = true;
  dfm.precompute_somig_coords(sanity_check);
  
  std::cout << " Apply some deformation to the cage " << std::endl;
  // Apply some deformation to the cage:
  std::vector< point3d > cage_modified_vertices;
  OBJIO::open("../models/cage_deformed.obj", cage_modified_vertices);

  // Diffuse deformation
  std::vector< point3d > mesh_modified_vertices(mesh_vertices.size());
  dfm.deform(cage_modified_vertices, mesh_modified_vertices);
  
  std::cout << " Sava deformed cage" << std::endl;
  // Save deformed cage:
  OBJIO::save("cage_deformed.obj", cage_modified_vertices, cage_triangles);
  
  std::cout << " Save deformed mesh " << std::endl;
  // Save deformed mesh:
  OBJIO::save("mesh_deformed.obj", mesh_modified_vertices, mesh_triangles);

  cout << "[INFO] done" << endl;
  return 0;
}
