// Copyright (c) 2014-present K. S. Ernest (iFire) Lee.
// SPDX-License-Identifier: MIT
//
//------------------------------------------- DISCLAIMER ----------------------------------------------------//
// Implementation published as additional material                                                           //
// It comes without any guarantees. Use it at you own risks.                                                 //
// In case you use this software, please acknowledge it and cite:                                            //
//    Biharmonic Coordinates and their Derivatives for Triangular 3D Cages                                   //
//    Jean-Marc THIERY, Elie MICHEL, Jiong CHEN                                                              //
//    Siggraph 2024                                                                                          //
//-----------------------------------------------------------------------------------------------------------//

#ifndef MESH_MORPH_3D_H
#define MESH_MORPH_3D_H

#include <godot_cpp/godot.hpp>
#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/surface_tool.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/resource_loader.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include "point3.h"

namespace godot {

class MeshMorph3D : public MeshInstance3D {
    GDCLASS(MeshMorph3D, MeshInstance3D)
private:
    float gamma_D_13BC = 1.0;
    String cage_mesh_path = "triangle_3d_cage/art/cage.obj";
    String cage_deformed_path = "triangle_3d_cage/art/cage_deformed.obj";
    String mesh_path = "triangle_3d_cage/art/mesh.obj";
    bool deformation_switch = false;

	std::vector<point3d> convert_godot_array_to_vector(const Array &godot_array);
	std::vector<point3d> extract_vertices(Ref<ArrayMesh> mesh);
	const std::vector<std::vector<unsigned int>> extract_triangles(Ref<ArrayMesh> mesh);
	std::vector<point3d> calculate_normals(Ref<ArrayMesh> mesh);
	std::vector<point3d> apply_deformation(const std::vector<point3d> &original_vertices, const std::vector<std::vector<unsigned int>> &cage_triangles, const std::vector<point3d> &cage_vertices, const std::vector<point3d> &cage_modified_vertices, std::vector<point3d> &cage_triangle_normals);

protected:
	static void _bind_methods();
    
public:
    MeshMorph3D();
    ~MeshMorph3D();

    void _init();
    void apply_deformation_to_children();
    void set_gamma_D_13BC(float value) { gamma_D_13BC = value; }
    float get_gamma_D_13BC() const { return gamma_D_13BC; }
	void set_cage_mesh_path(String path);
	String get_cage_mesh_path() const;
	void set_cage_deformed_path(String path);
	String get_cage_deformed_path() const;
	void set_mesh_path(String path);
	String get_mesh_path() const;
	void set_deformation_switch(bool value);
	bool get_deformation_switch() const;
};
}

#endif // MESH_MORPH_3D_H
