/**************************************************************************/
/*  mesh_morph_3d.cpp                                                     */
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

#include "mesh_morph_3d.h"
#include "BHC.h"
#include "point3.h"
#include <Eigen/Dense>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/mesh_data_tool.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

using namespace godot;

void MeshMorph3D::_bind_methods() {
	ClassDB::bind_method(D_METHOD("apply_deformation_to_children"), &MeshMorph3D::apply_deformation_to_children);
	ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "deform_mesh"), "set_deform_mesh", "get_deform_mesh");
	ClassDB::bind_method(D_METHOD("set_gamma_D_13BC", "gamma"), &MeshMorph3D::set_gamma_D_13BC);
	ClassDB::bind_method(D_METHOD("get_gamma_D_13BC"), &MeshMorph3D::get_gamma_D_13BC);
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "gamma"), "set_gamma_D_13BC", "get_gamma_D_13BC");

	ClassDB::bind_method(D_METHOD("set_cage_mesh", "mesh"), &MeshMorph3D::set_cage_mesh);
    ClassDB::bind_method(D_METHOD("get_cage_mesh"), &MeshMorph3D::get_cage_mesh);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "cage_mesh", PROPERTY_HINT_RESOURCE_TYPE, "ArrayMesh"), "set_cage_mesh", "get_cage_mesh");

    ClassDB::bind_method(D_METHOD("set_cage_deformed", "path"), &MeshMorph3D::set_cage_deformed);
    ClassDB::bind_method(D_METHOD("get_cage_deformed"), &MeshMorph3D::get_cage_deformed);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "cage_deformed", PROPERTY_HINT_RESOURCE_TYPE, "ArrayMesh"), "set_cage_deformed", "get_cage_deformed");

    ClassDB::bind_method(D_METHOD("set_source_mesh", "mesh"), &MeshMorph3D::set_source_mesh);
    ClassDB::bind_method(D_METHOD("get_source_mesh"), &MeshMorph3D::get_source_mesh);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "source_mesh", PROPERTY_HINT_RESOURCE_TYPE, "ArrayMesh"), "set_source_mesh", "get_source_mesh");

	ClassDB::bind_method(D_METHOD("set_deformation_switch", "value"), &MeshMorph3D::set_deformation_switch);
	ClassDB::bind_method(D_METHOD("get_deformation_switch"), &MeshMorph3D::get_deformation_switch);
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "deformation_switch"), "set_deformation_switch", "get_deformation_switch");
}

MeshMorph3D::MeshMorph3D() {
}

MeshMorph3D::~MeshMorph3D() {
}

void MeshMorph3D::apply_deformation_to_children() {
	if (cage_mesh.is_null()) {
		UtilityFunctions::printerr("Cage mesh is not available.");
		return;
	}
	std::vector<point3d> cage_vertices;
	std::vector<std::vector<unsigned int>> cage_triangles;
    extract_mesh_data(cage_mesh,
		cage_vertices,
		cage_triangles,
		true);
	UtilityFunctions::print("Cage model loaded successfully.");
	std::vector<point3d> mesh_vertices;
	std::vector<std::vector<unsigned int>> mesh_triangles;
	extract_mesh_data(source_mesh,
		mesh_vertices,
		mesh_triangles,
		true);
	UtilityFunctions::print(String("Vertex count in original mesh: ") + String::num_int64(mesh_vertices.size()));

	{
		Ref<ArrayMesh> array_mesh = memnew(ArrayMesh);
		Array arrays;
		arrays.resize(Mesh::ARRAY_MAX);

		PackedVector3Array vertices;
		for (const auto &vertex : mesh_vertices) {
			vertices.push_back(Vector3(vertex[0], vertex[1], vertex[2]));
		}
		arrays[Mesh::ARRAY_VERTEX] = vertices;

		PackedInt32Array indices;
		for (const auto &triangle : mesh_triangles) {
			for (unsigned int index : triangle) {
				indices.push_back(index);
			}
		}
		arrays[Mesh::ARRAY_INDEX] = indices;

		array_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);

		std::vector<point3d> new_mesh_vertices;
		std::vector<std::vector<unsigned int>> new_mesh_triangles;

		if (array_mesh->get_surface_count() > 0) {
			Array new_arrays = array_mesh->surface_get_arrays(0);
			Array new_vertices = new_arrays[Mesh::ARRAY_VERTEX];
			Array new_indices = new_arrays[Mesh::ARRAY_INDEX];

			for (int i = 0; i < new_vertices.size(); ++i) {
				Vector3 v = new_vertices[i];
				new_mesh_vertices.push_back({ v.x, v.y, v.z });
			}

			for (int i = 0; i < new_indices.size(); i += 3) {
				std::vector<unsigned int> triangle = {
					static_cast<unsigned int>(new_indices[i]),
					static_cast<unsigned int>(new_indices[i + 1]),
					static_cast<unsigned int>(new_indices[i + 2])
				};
				new_mesh_triangles.push_back(triangle);
			}
		}
		mesh_vertices = new_mesh_vertices;
		mesh_triangles = new_mesh_triangles;
	}
	// Compute (1,3)-regularized matrices
	Eigen::MatrixXd ConstrainedBiH_13_C11;
	Eigen::MatrixXd ConstrainedBiH_13_C12;
	Eigen::MatrixXd ConstrainedBiH_13_C21;
	Eigen::MatrixXd ConstrainedBiH_13_C22;
	std::vector<std::vector<double>> BHConstrainedC_13_phi, BHConstrainedC_13_psi; // (1,3) version
	std::vector<std::vector<double>> BHC_h_phi, BHC_h_psi, BHC_bh_phi, BHC_bh_psi; // unconstrained
	BiharmonicCoordinates3D::computeConstrainedBiharmonicMatrices_13(
			cage_triangles, cage_vertices,
			ConstrainedBiH_13_C11, ConstrainedBiH_13_C12, ConstrainedBiH_13_C21, ConstrainedBiH_13_C22,
			gamma_D_13BC);
	UtilityFunctions::print("(1,3)-regularized matrices computed successfully.");
	BHC_h_phi.resize(mesh_vertices.size());
	BHC_bh_phi.resize(mesh_vertices.size());
	BHC_h_psi.resize(mesh_vertices.size());
	BHC_bh_psi.resize(mesh_vertices.size());

	// Compute unconstrained coordinates
	UtilityFunctions::print("Compute (1,3)-regularized BHC ");
	BHConstrainedC_13_phi.resize(mesh_vertices.size());
	BHConstrainedC_13_psi.resize(mesh_vertices.size());

#pragma omp parallel for
	for (int p_idx = 0; p_idx < mesh_vertices.size(); ++p_idx) {
		BiharmonicCoordinates3D::computeCoordinates(mesh_vertices[p_idx],
				cage_triangles, cage_vertices,
				BHC_h_phi[p_idx], BHC_h_psi[p_idx], BHC_bh_phi[p_idx], BHC_bh_psi[p_idx]);
	}
	UtilityFunctions::print("Unconstrained coordinates computed.");

	// Compute (1,3)-regularized blending from unconstrained biharmonics
	for (int p_idx = 0; p_idx < mesh_vertices.size(); ++p_idx) {
		BiharmonicCoordinates3D::compute_13_blending_from_unconstrained_biharmonics(
				BHC_h_phi[p_idx], BHC_h_psi[p_idx], BHC_bh_phi[p_idx], BHC_bh_psi[p_idx],
				ConstrainedBiH_13_C11, ConstrainedBiH_13_C12, ConstrainedBiH_13_C21, ConstrainedBiH_13_C22,
				BHConstrainedC_13_phi[p_idx], BHConstrainedC_13_psi[p_idx]);
	}
	std::vector<point3d> cage_modified_vertices;
    extract_mesh_data(cage_deformed,
		cage_modified_vertices,
		cage_triangles,
		false);
	// Compute cage triangle normals
	std::vector<point3d> cage_triangle_normals(cage_triangles.size(), point3d(0, 0, 0));
	for (unsigned int tIt = 0; tIt < cage_triangles.size(); ++tIt) {
		auto &t = cage_triangles[tIt];
		cage_triangle_normals[tIt] = point3d::cross(cage_modified_vertices[t[1]] - cage_modified_vertices[t[0]],
				cage_modified_vertices[t[2]] - cage_modified_vertices[t[0]])
											 .direction();
	}
	UtilityFunctions::print("Cage triangle normals computed.");

	// Apply deformation to mesh vertices based on computed coordinates and normals
#pragma omp parallel for
	for (int v = 0; v < mesh_vertices.size(); ++v) {
		point3d pos(0, 0, 0);
		for (unsigned int vc = 0; vc < cage_vertices.size(); ++vc) {
			pos += BHConstrainedC_13_phi[v][vc] * cage_modified_vertices[vc];
		}
		for (unsigned int tc = 0; tc < cage_triangles.size(); ++tc) {
			pos += BHConstrainedC_13_psi[v][tc] * cage_triangle_normals[tc];
		}
		mesh_vertices[v] = pos;
	}
    UtilityFunctions::print("Mesh deformation updated from cage deformation.");
	Ref<SurfaceTool> st = memnew(SurfaceTool);
	st->begin(Mesh::PRIMITIVE_TRIANGLES);
	for (const auto& vertex : mesh_vertices) {
		Vector3 godot_vertex(vertex.x(), vertex.z(), -vertex.y());
		st->add_vertex(godot_vertex);
	}
	for (const auto& triangle : mesh_triangles) {
		for (int i = triangle.size() - 1; i >= 0; --i) {
			st->add_index(triangle[i]);
		}
	}
	st->generate_normals();
	st->generate_tangents();
	set_mesh(st->commit());
}

std::vector<point3d> godot::MeshMorph3D::convert_godot_array_to_vector(const Array &godot_array) {
	std::vector<point3d> vec;
	vec.reserve(godot_array.size());
	for (int i = 0; i < godot_array.size(); ++i) {
		Vector3 v = godot_array[i];
		vec.emplace_back(point3d(v.x, v.y, v.z));
	}
	return vec;
}

std::vector<point3d> godot::MeshMorph3D::extract_vertices(Ref<ArrayMesh> mesh) {
	std::vector<point3d> vertices;
	Array arrays = mesh->surface_get_arrays(0);
	PackedVector3Array vertex_array = arrays[Mesh::ARRAY_VERTEX];
	PackedInt32Array index_array = arrays[Mesh::ARRAY_INDEX];
	if (index_array.is_empty()) {
		for (int i = 0; i < vertex_array.size(); ++i) {
			Vector3 godot_vertex = vertex_array[i];
			// Swapping Y and Z coordinates, and negating the new Y
			point3d converted_vertex(godot_vertex.x, -godot_vertex.z, godot_vertex.y);
			vertices.push_back(converted_vertex);
		}
	} else {
		for (int i = 0; i < index_array.size(); ++i) {
			int index = index_array[i];
			if (index >= 0 && index < vertex_array.size()) {
				Vector3 godot_vertex = vertex_array[index];
				// Swapping Y and Z coordinates, and negating the new Y
				point3d converted_vertex(godot_vertex.x, -godot_vertex.z, godot_vertex.y);
				vertices.push_back(converted_vertex);
			} else {
				UtilityFunctions::printerr("ERROR: Index out of bounds in the vertex array.");
				return std::vector<point3d>();
			}
		}
	}
	return vertices;
}

const std::vector<std::vector<unsigned int>> godot::MeshMorph3D::extract_triangles(Ref<ArrayMesh> mesh) {
	std::vector<std::vector<unsigned int>> triangles;
	// Assuming mesh has only one surface
	Array arrays = mesh->surface_get_arrays(0);
	PackedInt32Array index_array = arrays[Mesh::ARRAY_INDEX];
	for (int i = 0; i < index_array.size(); i += 3) {
		if (i + 2 < index_array.size()) {
			std::vector<unsigned int> triangle;
			triangle.push_back(index_array[i]);
			triangle.push_back(index_array[i + 1]);
			triangle.push_back(index_array[i + 2]);
			triangles.push_back(triangle);
		}
	}
	return triangles;
}

bool godot::MeshMorph3D::get_deformation_switch() const {
	return deformation_switch;
}

void godot::MeshMorph3D::set_deformation_switch(bool value) {
	deformation_switch = value;
	if (deformation_switch) {
		set_mesh(Ref<ArrayMesh>());
		apply_deformation_to_children();
	}
}

Ref<ArrayMesh> godot::MeshMorph3D::get_cage_deformed() const {
	return cage_deformed;
}

void godot::MeshMorph3D::set_cage_deformed(Ref<ArrayMesh> p_mesh) {
	cage_deformed = p_mesh;
}

Ref<ArrayMesh> godot::MeshMorph3D::get_cage_mesh() const {
	return cage_mesh;
}

void godot::MeshMorph3D::set_cage_mesh(Ref<ArrayMesh> p_mesh) {
	cage_mesh = p_mesh;
}

void godot::MeshMorph3D::extract_mesh_data(const Ref<ArrayMesh> mesh,
		std::vector<point3d> &vertices,
		std::vector<std::vector<unsigned int>> &triangles,
		bool include_triangles) {
	if (mesh.is_null()) {
		UtilityFunctions::printerr("Mesh is not available.", __FUNCTION__, __FILE__, __LINE__);
		return;
	}
	godot::Ref<godot::MeshDataTool> mdt = memnew(MeshDataTool);
	mdt->create_from_surface(mesh, 0);

	int vertex_count = mdt->get_vertex_count();
	for (int i = 0; i < vertex_count; ++i) {
		godot::Vector3 vertex = mdt->get_vertex(i);
		std::swap(vertex.y, vertex.z);
		vertex.y = -vertex.y;
		vertices.push_back(vertex);
	}

	if (include_triangles) {
		int face_count = mdt->get_face_count();
		for (int i = 0; i < face_count; ++i) {
			std::vector<unsigned int> triangle;
			triangle.push_back(mdt->get_face_vertex(i, 0));
			triangle.push_back(mdt->get_face_vertex(i, 2));
			triangle.push_back(mdt->get_face_vertex(i, 1));
			triangles.push_back(triangle);
		}
	}
	UtilityFunctions::print("Mesh data extracted successfully.");
}
