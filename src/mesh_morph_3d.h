/**************************************************************************/
/*  mesh_morph_3d.h                                                       */
/**************************************************************************/
/*                         This file is part of:                          */
/*                             GODOT ENGINE                               */
/*                        https://godotengine.org                         */
/**************************************************************************/
/* Copyright (c) 2014-present Godot Engine contributors (see AUTHORS.md). */
/* Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur.                  */
/*                                                                        */
/* Permission is hereby granted, free of charge, to any person obtaining  */
/* a copy of this software and associated documentation files (the        */
/* "Software"), to deal in the Software without restriction, including    */
/* without limitation the rights to use, copy, modify, merge, publish,    */
/* distribute, sublicense, and/or sell copies of the Software, and to     */
/* permit persons to whom the Software is furnished to do so, subject to  */
/* the following conditions:                                              */
/*                                                                        */
/* The above copyright notice and this permission notice shall be         */
/* included in all copies or substantial portions of the Software.        */
/*                                                                        */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. */
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 */
/**************************************************************************/

#ifndef MESH_MORPH_3D_H
#define MESH_MORPH_3D_H
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

#include "point3.h"
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/resource_loader.hpp>
#include <godot_cpp/classes/surface_tool.hpp>
#include <godot_cpp/godot.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

class MeshMorph3D : public MeshInstance3D {
	GDCLASS(MeshMorph3D, MeshInstance3D)
private:
	float gamma_D_13BC = 1.0;
    NodePath cage_mesh_path;
    NodePath cage_deformed_path;
    NodePath source_path;
	bool deformation_switch = false;
	std::vector<point3d> convert_godot_array_to_vector(const Array &godot_array);
	std::vector<point3d> extract_vertices(Ref<ArrayMesh> mesh);
	const std::vector<std::vector<unsigned int>> extract_triangles(Ref<ArrayMesh> mesh);
	std::vector<point3d> calculate_normals(Ref<ArrayMesh> mesh);
	std::vector<point3d> apply_deformation(const std::vector<point3d> &original_vertices, const std::vector<std::vector<unsigned int>> &cage_triangles, const std::vector<point3d> &cage_vertices, const std::vector<point3d> &cage_modified_vertices, std::vector<point3d> &cage_triangle_normals);
	void extract_mesh_data(const Ref<ArrayMesh> mesh,
			std::vector<point3d> &vertices,
			std::vector<std::vector<unsigned int>> &triangles,
			bool include_triangles = true);
	void apply_deformation_to_children();
protected:
	static void _bind_methods();

public:
	MeshMorph3D();
	~MeshMorph3D();
	void set_gamma_D_13BC(float value) { gamma_D_13BC = value; }
	float get_gamma_D_13BC() const { return gamma_D_13BC; }
	void set_cage_mesh_from_path(NodePath path);
	NodePath get_cage_mesh_from_path() const;
	void set_cage_deformed_from_path(NodePath path);
	NodePath get_cage_deformed_from_path() const;
	void set_source_mesh_from_path(NodePath p_mesh);
	NodePath get_source_mesh_from_path() const;
	void set_deformation_switch(bool value);
	bool get_deformation_switch() const;
};
} //namespace godot

#endif // MESH_MORPH_3D_H
