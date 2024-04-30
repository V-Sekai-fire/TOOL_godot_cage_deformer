
#include "include/BHC.h"
#include "include/BasicIO.h"






int main(int argc, char const* argv[])
{
    std::cout << " Load cage " << std::endl;
    // Load cage:
    std::vector< point3d > cage_vertices;
    std::vector< std::vector< unsigned int > > cage_triangles;
    OBJIO::open("../models/cage.obj", cage_vertices, cage_triangles, true); // TRIANGULATE THE FACES OF THE CAGE HERE! (obviously)



    std::cout << " Load mesh " << std::endl;
    // Load mesh:
    std::vector< point3d > mesh_vertices;
    std::vector< std::vector< unsigned int > > mesh_triangles;
    OBJIO::open("../models/mesh.obj", mesh_vertices, mesh_triangles, true);
     


    std::cout << " Coordinates " << std::endl;
    // Coordinates:
    std::vector< std::vector< double > >  BHC_h_phi, BHC_h_psi, BHC_bh_phi, BHC_bh_psi; // unconstrained
    std::vector< std::vector< double > >  BHConstrainedC_13_phi, BHConstrainedC_13_psi; // (1,3) version



    std::cout << " Compute (1,3)-regularized matrices " << std::endl;
    // Compute (1,3)-regularized matrices (depends only on the cage, not on the mesh):
    Eigen::MatrixXd ConstrainedBiH_13_C11;
    Eigen::MatrixXd ConstrainedBiH_13_C12;
    Eigen::MatrixXd ConstrainedBiH_13_C21;
    Eigen::MatrixXd ConstrainedBiH_13_C22;
    double gamma_D_13BC = 1.0; // you can play with this

    BiharmonicCoordinates3D::computeConstrainedBiharmonicMatrices_13
    (
        cage_triangles, cage_vertices,
        ConstrainedBiH_13_C11, ConstrainedBiH_13_C12, ConstrainedBiH_13_C21, ConstrainedBiH_13_C22,
        gamma_D_13BC
    );



    std::cout << " Compute unconstrained coordinates " << std::endl;
    // Compute unconstrained coordinates for the mesh vertices:
    {
        BHC_h_phi.resize(mesh_vertices.size());
        BHC_bh_phi.resize(mesh_vertices.size());
        BHC_h_psi.resize(mesh_vertices.size());
        BHC_bh_psi.resize(mesh_vertices.size());

#pragma omp parallel for 
        for (int p_idx = 0; p_idx < mesh_vertices.size(); ++p_idx)
        {
            BiharmonicCoordinates3D::computeCoordinates(mesh_vertices[p_idx],
                cage_triangles,
                cage_vertices,
                BHC_h_phi[p_idx], BHC_h_psi[p_idx], BHC_bh_phi[p_idx], BHC_bh_psi[p_idx]);

        }
    }



    std::cout << " Compute (1,3)-regularized BHC " << std::endl;
    // Compute (1,3)-regularized BHC:
    {
        BHConstrainedC_13_phi.resize(mesh_vertices.size());
        BHConstrainedC_13_psi.resize(mesh_vertices.size());
        for (int p_idx = 0; p_idx < mesh_vertices.size(); ++p_idx)
        {
            BiharmonicCoordinates3D::compute_13_blending_from_unconstrained_biharmonics(
                BHC_h_phi[p_idx], BHC_h_psi[p_idx], BHC_bh_phi[p_idx], BHC_bh_psi[p_idx],
                ConstrainedBiH_13_C11, ConstrainedBiH_13_C12, ConstrainedBiH_13_C21, ConstrainedBiH_13_C22,
                BHConstrainedC_13_phi[p_idx], BHConstrainedC_13_psi[p_idx]
            );
        }
    }


    std::cout << " Apply some deformation to the cage " << std::endl;
    // Apply some deformation to the cage:
    std::vector< point3d > cage_modified_vertices;
    OBJIO::open("../models/cage_deformed.obj", cage_modified_vertices);



    std::cout << " Compute cage triangle normals " << std::endl;
    // Compute cage triangle normals:
    std::vector< point3d > cage_triangle_normals(cage_triangles.size(), point3d(0, 0, 0));
    for (unsigned int tIt = 0; tIt < cage_triangles.size(); ++tIt) {
        auto& t = cage_triangles[tIt];
        cage_triangle_normals[tIt] = point3d::cross(cage_modified_vertices[t[1]] - cage_modified_vertices[t[0]], cage_modified_vertices[t[2]] - cage_modified_vertices[t[0]]).direction();
    }


    std::cout << " Update the mesh deformation, from the cage deformation " << std::endl;
    // Update the mesh deformation, from the cage deformation:
    std::vector< point3d > mesh_modified_vertices(mesh_vertices.size());
    {
#pragma omp parallel for
        for (int v = 0; v < mesh_vertices.size(); ++v)
        {
            point3d pos(0, 0, 0);
            for (unsigned int vc = 0; vc < cage_modified_vertices.size(); ++vc)
                pos += BHConstrainedC_13_phi[v][vc] * cage_modified_vertices[vc];
            for (unsigned int tc = 0; tc < cage_triangles.size(); ++tc) {
                pos += BHConstrainedC_13_psi[v][tc] * cage_triangle_normals[tc];
            }

            mesh_modified_vertices[v] = pos;
        }
    }


    std::cout << " Save deformed mesh " << std::endl;
    // Save deformed mesh:
    OBJIO::save("../models/mesh_deformed.obj", mesh_modified_vertices, mesh_triangles);

    return 0;
}