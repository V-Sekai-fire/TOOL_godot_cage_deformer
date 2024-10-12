/**************************************************************************/
/*  BHC.h                                                                 */
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

#ifndef BHC_H
#define BHC_H

//-----------------------------------   ASSOCIATED REFERENCES   ---------------------------------------------//
//                                                                                                           //
// "Green Coordinates". Lipman et al. 2008, ACM Siggraph:                                                    //
// this paper presents the use of Green coordinates for deforming 3D meshes                                  //
//                                                                                                           //
// "Variational harmonic maps for space deformation". Ben-Chen et al. 2009, ACM Siggraph                     //
// this paper contains an alternative computation of the coordinates, along with gradients and Hessians,     //
// based on a work of Urago (see annex of the paper)                                                         //
//                                                                                                           //
// You can find these expressions in URAGO 2000 :                                                            //
// "Analytical integrals of fundamental solution of three-dimensional Laplace equation and their gradients"  //
//-----------------------------------------------------------------------------------------------------------//

//------------------------------------------- DISCLAIMER ----------------------------------------------------//
// Implementation published as additional material                                                           //
// It comes without any guarantees. Use it at you own risks.                                                 //
// In case you use this software, please acknowledge it and cite                                             //
//    Biharmonic Coordinates and their Derivatives for Triangular 3D Cages                                   //
//    Jean-Marc THIERY, Elie MICHEL, Jiong CHEN                                                              //
//    Siggraph 2024                                                                                          //
//-----------------------------------------------------------------------------------------------------------//

#include "point3.h"
#include <cmath>
#include <vector>

#include <Eigen/Dense>

namespace BiharmonicCoordinates3D {

// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// FINAL FORMULATION:

double get_signed_solid_angle(point3d const &a, point3d const &b, point3d const &c) {
	typedef double T;
	T det = point3d::dot(a, point3d::cross(b, c));

	// TO IMPROVE :
	if (fabs(det) < 0.0000000001) // then you're on the limit case where you cover half the sphere
		return 0; // chances are that you are outside the triangle

	//	if (fabs(det) < 0.0000000001) // then you're on the limit case where you cover half the sphere
	//		return 2.0 * M_PI; // that is particularly shitty, because the sign is difficult to estimate...

	T al = a.norm(), bl = b.norm(), cl = c.norm();

	T div = al * bl * cl + point3d::dot(a, b) * cl + point3d::dot(a, c) * bl + point3d::dot(b, c) * al;
	T at = atan2(fabs(det), div);
	if (at < 0)
		at += M_PI; // If det>0 && div<0 atan2 returns < 0, so add pi.
	T omega = 2.0 * at;

	if (det > 0.0)
		return omega;
	return -omega;
}

point3d get_signed_solid_angle_gradient_subrouting(point3d const &eta, point3d const &v0, point3d const &v1) {
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm(), e = (v1 - v0).norm();
	return -(2.0 * (l0 + l1)) * point3d::cross(v1 - eta, v0 - eta) / ((l0 + l1 + e) * (l0 + l1 - e) * l0 * l1);
}
point3d get_signed_solid_angle_gradient(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	return (
			get_signed_solid_angle_gradient_subrouting(eta, v0, v1) +
			get_signed_solid_angle_gradient_subrouting(eta, v1, v2) +
			get_signed_solid_angle_gradient_subrouting(eta, v2, v0));
}

double square_t(double x) { return x * x; }
mat33d MT(point3d const &p, point3d const &q) {
	return mat33d::tensor(p, q) + mat33d::tensor(q, p); // \mathcal{T}( p ; q ) in the paper.
}
mat33d MT(mat33d const &p, mat33d const &q) {
	return p * q.getTranspose() + q * p.getTranspose(); // \mathcal{T}( p ; q ) in the paper.
}

// HARMONIC PART :
double h_psi(
		point3d const &eta,
		point3d const *tri_vertices // an array of 3 points
) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;

	double psi = 0.0;

	Nt /= NtNorm;

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.f * M_PI);
	T signed_volume = point3d::dot(point3d::cross(e[0], e[1]), e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = -3.0 * signed_omega * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);

	return psi;
}

void h_psi_with_derivatives(
		point3d const &eta,
		point3d const *tri_vertices, // an array of 3 points
		double &psi,
		point3d &psi_grad) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;

	psi = 0.0;

	Nt /= NtNorm;

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.f * M_PI);
	T signed_volume = point3d::dot(point3d::cross(e[0], e[1]), e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = -3.0 * signed_omega * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	psi_grad = -Pt;
}

void h_psi_with_derivatives(
		point3d const &eta,
		point3d const *tri_vertices, // an array of 3 points
		double &psi,
		point3d &psi_grad,
		mat33d &psi_H) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;

	psi = 0.0;

	Nt /= NtNorm;

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.f * M_PI);
	T signed_volume = point3d::dot(point3d::cross(e[0], e[1]), e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = -3.0 * signed_omega * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	psi_grad = -Pt;

	for (unsigned int v = 0; v < 3; ++v)
		C_tilde[v] = 1.0 / (2.0 * M_PI * (R[v] + d_norm[v]) * (R[v] - d_norm[v]));
	for (unsigned int v = 0; v < 3; ++v)
		signed_solid_angle_gradient += (C_tilde[v] * (1.0 / e_norm[(v + 1) % 3] + 1.0 / e_norm[(v + 2) % 3])) * J[v];

	mat33d const &Jacobian1 = -mat33d::getFromCols(
									  C_tilde[0] * (e_normalized[1] + e_normalized[2]),
									  C_tilde[1] * (e_normalized[0] + e_normalized[2]),
									  C_tilde[2] * (e_normalized[1] + e_normalized[0])) *
					mat33d::getFromRows(d[0], d[1], d[2]) *
					mat33d::vectorial(Nt).getTranspose() -
			mat33d::tensor(signed_solid_angle_gradient, Nt);

	psi_H = Jacobian1;
}

void h_psi_with_derivatives_inside_triangle(
		point3d const &eta,
		point3d const *tri_vertices, // an array of 3 points
		double &psi,
		point3d &psi_grad,
		mat33d &psi_H) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;

	psi = 0.0;

	Nt /= NtNorm;

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = 0.5;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = 0;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	psi_grad = -Pt;

	for (unsigned int v = 0; v < 3; ++v)
		C_tilde[v] = 1.0 / (2.0 * M_PI * (R[v] + d_norm[v]) * (R[v] - d_norm[v]));
	for (unsigned int v = 0; v < 3; ++v)
		signed_solid_angle_gradient += (C_tilde[v] * (1.0 / e_norm[(v + 1) % 3] + 1.0 / e_norm[(v + 2) % 3])) * J[v];

	mat33d const &Jacobian1 = -mat33d::getFromCols(
									  C_tilde[0] * (e_normalized[1] + e_normalized[2]),
									  C_tilde[1] * (e_normalized[0] + e_normalized[2]),
									  C_tilde[2] * (e_normalized[1] + e_normalized[0])) *
					mat33d::getFromRows(d[0], d[1], d[2]) *
					mat33d::vectorial(Nt).getTranspose() -
			mat33d::tensor(signed_solid_angle_gradient, Nt);

	psi_H = Jacobian1;
}

void h_coordinates(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi,
		double *phi) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.0 * M_PI);
	T signed_volume = point3d::dot(point3d::cross(e[0], e[1]), e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = -3.0 * signed_omega * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);

	for (unsigned int v = 0; v < 3; ++v) {
		phi[v] = point3d::dot(Pt, J[v]) / (2.0 * At);
	}
}

void h_coordinates_inside_triangle(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi,
		double *phi) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = 0.5;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = 0;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);

	for (unsigned int v = 0; v < 3; ++v) {
		phi[v] = point3d::dot(Pt, J[v]) / (2.0 * At);
	}
}

void h_coordinates_with_derivatives(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi,
		double *phi,
		point3d &psi_grad,
		point3d *phi_grad) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.0 * M_PI);
	T signed_volume = point3d::dot(point3d::cross(e[0], e[1]), e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = -3.0 * signed_omega * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);
	psi_grad = -Pt;

	for (unsigned int v = 0; v < 3; ++v) {
		phi[v] = point3d::dot(Pt, J[v]) / (2.0 * At);
		phi_grad[v] = point3d::cross(Pt, d[v]) / (2.0 * At);
	}
}

// HANDLE THIS AS A SPECIAL CASE TO MAKE SURE YOU HAVE ROBUST COMPUTATIONS THERE:
void h_coordinates_normal_derivatives_inside_triangle(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi_normal_derivative,
		double *phi_normal_derivative) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3]; // needed
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3]; // needed
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm(); // needed
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]); // needed

	point3d Pt(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);

	psi_normal_derivative = 0.5; // 2 * M_PI; // divided by 4 pi

	for (unsigned int v = 0; v < 3; ++v) {
		phi_normal_derivative[v] = point3d::dot(point3d::cross(Pt, d[v]) / (2.0 * At), Nt);
	}
}

// HANDLE THIS AS A SPECIAL CASE TO MAKE SURE YOU HAVE ROBUST COMPUTATIONS THERE:
void h_coordinates_on_triangle_corner(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi,
		double *phi,
		unsigned int cornerIt,
		point3d const &corner_safety_offset) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v) {
		if (v != cornerIt)
			e_normalized[v] = e[v] / e_norm[v];
	}
	e_normalized[cornerIt] = corner_safety_offset.direction();

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.0 * M_PI);

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3]; // Needed FOR cornerIt ONLY
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3]; // Needed FOR cornerIt ONLY
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm(); // Needed FOR cornerIt ONLY
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]); // Needed FOR cornerIt ONLY

	psi = 2 * At * C[cornerIt];

	phi[0] = phi[1] = phi[2] = 0;
	phi[cornerIt] = fabs(signed_omega); // CAREFUL ! it needs a special handling.
}

void h_coordinates_with_derivatives(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi,
		double *phi,
		point3d &psi_grad,
		point3d *phi_grad,
		mat33d &psi_H,
		mat33d *phi_H) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.0 * M_PI);
	T signed_volume = point3d::dot(point3d::cross(e[0], e[1]), e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = -3.0 * signed_omega * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);
	psi_grad = -Pt;

	for (unsigned int v = 0; v < 3; ++v) {
		phi[v] = point3d::dot(Pt, J[v]) / (2.0 * At);
		phi_grad[v] = point3d::cross(Pt, d[v]) / (2.0 * At);
	}

	for (unsigned int v = 0; v < 3; ++v)
		C_tilde[v] = 1.0 / (2.0 * M_PI * (R[v] + d_norm[v]) * (R[v] - d_norm[v]));
	for (unsigned int v = 0; v < 3; ++v)
		signed_solid_angle_gradient += (C_tilde[v] * (1.0 / e_norm[(v + 1) % 3] + 1.0 / e_norm[(v + 2) % 3])) * J[v];

	mat33d const &Jacobian1 = mat33d::getFromCols(
									  C_tilde[0] * (e_normalized[1] + e_normalized[2]),
									  C_tilde[1] * (e_normalized[0] + e_normalized[2]),
									  C_tilde[2] * (e_normalized[1] + e_normalized[0])) *
					mat33d::getFromRows(d[0], d[1], d[2]) *
					mat33d::vectorial(Nt).getTranspose() +
			mat33d::tensor(signed_solid_angle_gradient, Nt);

	psi_H = -1.0 * Jacobian1;

	for (unsigned int v = 0; v < 3; ++v)
		phi_H[v] = (-0.5 / At) * mat33d::vectorial(d[v]) * Jacobian1;
}

void h_coordinates_with_derivatives_inside_triangle(
		point3d const &eta,
		point3d const *tri_vertices,
		double &psi,
		double *phi,
		point3d &psi_grad,
		point3d *phi_grad,
		mat33d &psi_H,
		mat33d *phi_H) {
	typedef double T;
	point3d Nt = point3d::cross(tri_vertices[1] - tri_vertices[0], tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt.normalize();

	point3d e[3];
	T e_norm[3];
	point3d e_normalized[3];
	T R[3];
	point3d d[3];
	T d_norm[3];
	T C[3];
	point3d J[3];
	T C_tilde[3];
	point3d signed_solid_angle_gradient(0, 0, 0);
	for (unsigned int v = 0; v < 3; ++v)
		e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v)
		e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		e_normalized[v] = e[v] / e_norm[v];

	T signed_omega = 0.5;

	for (unsigned int v = 0; v < 3; ++v)
		R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v)
		d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v)
		C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	point3d Pt(-signed_omega * Nt);
	for (unsigned int v = 0; v < 3; ++v)
		Pt += point3d::cross(Nt, C[v] * d[v]);
	for (unsigned int v = 0; v < 3; ++v)
		J[v] = point3d::cross(e[(v + 2) % 3], e[(v + 1) % 3]);

	psi = 0;
	for (unsigned int v = 0; v < 3; ++v)
		psi -= C[v] * point3d::dot(J[v], Nt);
	psi_grad = -Pt;

	for (unsigned int v = 0; v < 3; ++v) {
		phi[v] = point3d::dot(Pt, J[v]) / (2.0 * At);
		phi_grad[v] = point3d::cross(Pt, d[v]) / (2.0 * At);
	}

	for (unsigned int v = 0; v < 3; ++v)
		C_tilde[v] = 1.0 / (2.0 * M_PI * (R[v] + d_norm[v]) * (R[v] - d_norm[v]));
	for (unsigned int v = 0; v < 3; ++v)
		signed_solid_angle_gradient += (C_tilde[v] * (1.0 / e_norm[(v + 1) % 3] + 1.0 / e_norm[(v + 2) % 3])) * J[v];

	mat33d const &Jacobian1 = mat33d::getFromCols(
									  C_tilde[0] * (e_normalized[1] + e_normalized[2]),
									  C_tilde[1] * (e_normalized[0] + e_normalized[2]),
									  C_tilde[2] * (e_normalized[1] + e_normalized[0])) *
					mat33d::getFromRows(d[0], d[1], d[2]) *
					mat33d::vectorial(Nt).getTranspose() +
			mat33d::tensor(signed_solid_angle_gradient, Nt);

	psi_H = -1.0 * Jacobian1;

	for (unsigned int v = 0; v < 3; ++v)
		phi_H[v] = (-0.5 / At) * mat33d::vectorial(d[v]) * Jacobian1;
}

// BIHARMONIC PART :
double bh_psi_subroutine_without_solid_angle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &n, double d) {
	point3d const &u_e = (v1 - v0).direction();
	double zeta0 = point3d::dot(eta - v0, u_e); // l0C0
	double zeta1 = point3d::dot(eta - v1, u_e); // l1C1
	double a_times_sign_ti = point3d::dot(point3d::cross(n, u_e), eta - v0);
	double D_squared = ((mat33d::Identity() - mat33d::tensor(u_e, u_e)) * (eta - v0)).sqrnorm();
	double le0 = (v0 - eta).norm();
	double le1 = (v1 - eta).norm();

	double alignment = (2 * d * d + D_squared);
	if (alignment < 0.00000000001)
		return 0.0;

	return (
				   0.5 * a_times_sign_ti * ((2 * d * d + D_squared) * log((le1 - zeta1) / (le0 - zeta0)) + le0 * zeta0 - le1 * zeta1)) /
			3;
}

double bh_psi(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();
	double d = point3d::dot(eta - v0, n); // Convention for the signed distance to the triangle: positive above the triangle, simpler to differentiate (no extra minus sign involved)
	double omega = get_signed_solid_angle(v2 - eta, v0 - eta, v1 - eta);

	return (1.0 / (8 * M_PI)) * (bh_psi_subroutine_without_solid_angle(eta, v0, v1, n, d) + bh_psi_subroutine_without_solid_angle(eta, v1, v2, n, d) + bh_psi_subroutine_without_solid_angle(eta, v2, v0, n, d) + d * d * d * omega / 3);
}

// HANDLE THIS AS A SPECIAL CASE TO MAKE SURE YOU HAVE ROBUST COMPUTATIONS THERE:
double bh_psi_on_triangle_corner(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2, unsigned int cornerIt) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();
	double d = point3d::dot(eta - v0, n); // Convention for the signed distance to the triangle: positive above the triangle, simpler to differentiate (no extra minus sign involved)

	if (cornerIt == 0)
		return (1.0 / (8 * M_PI)) * (bh_psi_subroutine_without_solid_angle(eta, v1, v2, n, d));
	if (cornerIt == 1)
		return (1.0 / (8 * M_PI)) * (bh_psi_subroutine_without_solid_angle(eta, v2, v0, n, d));
	if (cornerIt == 2)
		return (1.0 / (8 * M_PI)) * (bh_psi_subroutine_without_solid_angle(eta, v0, v1, n, d));

	// default (this works as well):
	return (1.0 / (8 * M_PI)) * (bh_psi_subroutine_without_solid_angle(eta, v0, v1, n, d) + bh_psi_subroutine_without_solid_angle(eta, v1, v2, n, d) + bh_psi_subroutine_without_solid_angle(eta, v2, v0, n, d));
}

double bh_psi_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();
	double d = point3d::dot(eta - v0, n); // Convention for the signed distance to the triangle: positive above the triangle, simpler to differentiate (no extra minus sign involved)

	return (1.0 / (8 * M_PI)) * (bh_psi_subroutine_without_solid_angle(eta, v0, v1, n, d) + bh_psi_subroutine_without_solid_angle(eta, v1, v2, n, d) + bh_psi_subroutine_without_solid_angle(eta, v2, v0, n, d));
}

point3d bh_psi_gradient_subroutine(point3d const &eta, point3d const &v0, point3d const &v1, double d, point3d const &n) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();

	double a_e = point3d::dot(r_e, eta - v0);

	point3d const &D_eVec = (mat33d::Identity() - mat33d::tensor(u_e, u_e)) * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	if (D_e_squared < 0.0000000000001)
		return (-l1 * zeta1 + l0 * zeta0) * r_e / 6.0;

	return ((2 * d * d + D_e_squared) * log((l1 - zeta1) / (l0 - zeta0)) - l1 * zeta1 + l0 * zeta0) * r_e / 6.0 +
			(a_e / 6.0) * log((l1 - zeta1) / (l0 - zeta0)) * (4 * d * n + 2 * D_eVec) +
			(a_e / 6.0) * (2 * d * d + D_e_squared) * ((u1 - u_e) / (l1 - zeta1) - (u0 - u_e) / (l0 - zeta0)) +
			(a_e / 6.0) * (zeta0 * u0 - zeta1 * u1 + (l0 - l1) * u_e);
}

point3d bh_psi_gradient(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();
	double d = point3d::dot(eta - v0, n); // Convention for the signed distance to the triangle: positive above the triangle, simpler to differentiate (no extra minus sign involved)
	double omega = get_signed_solid_angle(v2 - eta, v0 - eta, v1 - eta);
	point3d omega_gradient = get_signed_solid_angle_gradient(eta, v0, v1, v2);

	return (1.0 / (8 * M_PI)) * ((d * d * omega) * n + (d * d * d / 3.0) * omega_gradient + bh_psi_gradient_subroutine(eta, v0, v1, d, n) + bh_psi_gradient_subroutine(eta, v1, v2, d, n) + bh_psi_gradient_subroutine(eta, v2, v0, d, n));
}

point3d bh_psi_gradient_inside_triangle_subroutine(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &n) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();

	double a_e = point3d::dot(r_e, eta - v0);

	point3d const &D_eVec = (mat33d::Identity() - mat33d::tensor(u_e, u_e)) * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	if (D_e_squared < 0.0000000000001)
		return (-l1 * zeta1 + l0 * zeta0) * r_e / 6.0;

	return ((D_e_squared)*log((l1 - zeta1) / (l0 - zeta0)) - l1 * zeta1 + l0 * zeta0) * r_e / 6.0 +
			(a_e / 6.0) * log((l1 - zeta1) / (l0 - zeta0)) * (2 * D_eVec) +
			(a_e / 6.0) * (D_e_squared) * ((u1 - u_e) / (l1 - zeta1) - (u0 - u_e) / (l0 - zeta0)) +
			(a_e / 6.0) * (zeta0 * u0 - zeta1 * u1 + (l0 - l1) * u_e);
}

point3d bh_psi_gradient_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();

	return (1.0 / (8 * M_PI)) * (bh_psi_gradient_inside_triangle_subroutine(eta, v0, v1, n) + bh_psi_gradient_inside_triangle_subroutine(eta, v1, v2, n) + bh_psi_gradient_inside_triangle_subroutine(eta, v2, v0, n));
}

// HANDLE THIS AS A SPECIAL CASE TO MAKE SURE YOU HAVE ROBUST COMPUTATIONS THERE:
point3d bh_psi_gradient_subroutine_for_normal_derivatives_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &n) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();

	double a_e = point3d::dot(r_e, eta - v0);

	point3d const &D_eVec = (mat33d::Identity() - mat33d::tensor(u_e, u_e)) * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	return ((D_e_squared)*log((l1 - zeta1) / (l0 - zeta0)) - l1 * zeta1 + l0 * zeta0) * r_e / 6.0 +
			(a_e / 6.0) * log((l1 - zeta1) / (l0 - zeta0)) * (2 * D_eVec) +
			(a_e / 6.0) * (D_e_squared) * ((u1 - u_e) / (l1 - zeta1) - (u0 - u_e) / (l0 - zeta0)) +
			(a_e / 6.0) * (zeta0 * u0 - zeta1 * u1 + (l0 - l1) * u_e);
}
void bh_coordinates_normal_derivatives_inside_triangle(
		point3d const &eta,
		point3d const *tri,
		double &psi_normal_derivative,
		double *phi_normal_derivative) {
	psi_normal_derivative = 0; // fiou ! that was easy for once.

	point3d const &v0 = tri[0];
	point3d const &v1 = tri[1];
	point3d const &v2 = tri[2];
	point3d const &n = point3d::cross(tri[1] - tri[0], tri[2] - tri[0]).direction();
	double At = point3d::cross(tri[1] - tri[0], tri[2] - tri[0]).norm() / 2;

	point3d const &bh_psi_grad_special_case = -(1.0 / (8 * M_PI)) * (bh_psi_gradient_subroutine_for_normal_derivatives_inside_triangle(eta, tri[0], tri[1], n) + bh_psi_gradient_subroutine_for_normal_derivatives_inside_triangle(eta, tri[1], tri[2], n) + bh_psi_gradient_subroutine_for_normal_derivatives_inside_triangle(eta, tri[2], tri[0], n));

	point3d g0 = point3d::cross(n, v2 - v1);
	g0 /= point3d::dot(g0, v0 - v1);
	point3d g1 = point3d::cross(n, v0 - v2);
	g1 /= point3d::dot(g1, v1 - v2);
	point3d g2 = point3d::cross(n, v1 - v0);
	g2 /= point3d::dot(g2, v2 - v0);

	double gamma0 = point3d::dot(eta - v1, g0);
	double gamma1 = point3d::dot(eta - v2, g1);
	double gamma2 = point3d::dot(eta - v0, g2);

	phi_normal_derivative[0] = point3d::dot(bh_psi_grad_special_case, g0);
	phi_normal_derivative[1] = point3d::dot(bh_psi_grad_special_case, g1);
	phi_normal_derivative[2] = point3d::dot(bh_psi_grad_special_case, g2);

	{
		point3d e[3];
		double e_norm[3];
		point3d e_normalized[3];
		double R[3];
		point3d d[3];
		double d_norm[3];
		double C[3];
		for (unsigned int v = 0; v < 3; ++v)
			e[v] = tri[v] - eta;
		for (unsigned int v = 0; v < 3; ++v)
			e_norm[v] = e[v].norm();
		for (unsigned int v = 0; v < 3; ++v) {
			e_normalized[v] = e[v] / e_norm[v];
		}

		for (unsigned int v = 0; v < 3; ++v)
			R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
		for (unsigned int v = 0; v < 3; ++v)
			d[v] = tri[(v + 1) % 3] - tri[(v + 2) % 3];
		for (unsigned int v = 0; v < 3; ++v)
			d_norm[v] = d[v].norm();
		for (unsigned int v = 0; v < 3; ++v)
			C[v] = log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

		double psi = 2 * At * (C[0] * gamma0 + C[1] * gamma1 + C[2] * gamma2);

		phi_normal_derivative[0] += psi * gamma0 / 2;
		phi_normal_derivative[1] += psi * gamma1 / 2;
		phi_normal_derivative[2] += psi * gamma2 / 2;
	}
}

mat33d bh_psi_Hessian_subroutine(point3d const &eta, point3d const &v0, point3d const &v1, double d, point3d const &n) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();

	double a_e = point3d::dot(r_e, eta - v0);

	mat33d const &P_e = (mat33d::Identity() - mat33d::tensor(u_e, u_e));
	point3d const &D_eVec = P_e * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	mat33d U0 = (mat33d::Identity() - mat33d::tensor(u0, u0)) / l0;
	mat33d U1 = (mat33d::Identity() - mat33d::tensor(u1, u1)) / l1;

	point3d m0 = (u0 - u_e) / (l0 - zeta0);
	point3d m1 = (u1 - u_e) / (l1 - zeta1);
	mat33d M0 = U0 / (l0 - zeta0) - mat33d::tensor(m0, m0);
	mat33d M1 = U1 / (l1 - zeta1) - mat33d::tensor(m1, m1);

	mat33d contrib =
			(log((l1 - zeta1) / (l0 - zeta0)) / 6.0) * MT(r_e, 4 * d * n + 2 * D_eVec) + MT(r_e / 6, zeta0 * u0 - zeta1 * u1 + (l0 - l1) * u_e) + (a_e * log((l1 - zeta1) / (l0 - zeta0)) / 6.0) * (2 * MT(n, n) + 2 * P_e) + (1.0 / 6.0) * MT(6 * a_e * d * n + 3 * D_e_squared * r_e, m1 - m0) + (a_e / 6) * (zeta0 * U0 - zeta1 * U1 + MT(u0 - u1, u_e)) + (a_e * (2 * d * d + D_e_squared) / 6.0) * (M1 - M0);

	return contrib;
}

mat33d bh_psi_Hessian_subroutine_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &n) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();

	double a_e = point3d::dot(r_e, eta - v0);

	mat33d const &P_e = (mat33d::Identity() - mat33d::tensor(u_e, u_e));
	point3d const &D_eVec = P_e * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	mat33d U0 = (mat33d::Identity() - mat33d::tensor(u0, u0)) / l0;
	mat33d U1 = (mat33d::Identity() - mat33d::tensor(u1, u1)) / l1;

	point3d m0 = (u0 - u_e) / (l0 - zeta0);
	point3d m1 = (u1 - u_e) / (l1 - zeta1);
	mat33d M0 = U0 / (l0 - zeta0) - mat33d::tensor(m0, m0);
	mat33d M1 = U1 / (l1 - zeta1) - mat33d::tensor(m1, m1);

	mat33d contrib =
			(log((l1 - zeta1) / (l0 - zeta0)) / 6.0) * MT(r_e, 2 * D_eVec) + MT(r_e / 6, zeta0 * u0 - zeta1 * u1 + (l0 - l1) * u_e) + (a_e * log((l1 - zeta1) / (l0 - zeta0)) / 6.0) * (2 * MT(n, n) + 2 * P_e) + (1.0 / 6.0) * MT(3 * D_e_squared * r_e, m1 - m0) + (a_e / 6) * (zeta0 * U0 - zeta1 * U1 + MT(u0 - u1, u_e)) + (a_e * (D_e_squared) / 6.0) * (M1 - M0);

	return contrib;
}

mat33d bh_psi_Hessian_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();

	return (1.0 / (8 * M_PI)) * (bh_psi_Hessian_subroutine_inside_triangle(eta, v0, v1, n) + bh_psi_Hessian_subroutine_inside_triangle(eta, v1, v2, n) + bh_psi_Hessian_subroutine_inside_triangle(eta, v2, v0, n));
}

point3d bh_psi_gradient_in_tangent_plane_subroutine(point3d const &eta, point3d const &v0, point3d const &v1, double d, point3d const &n) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	double e = (v1 - v0).norm();
	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();

	double a_e = point3d::dot(r_e, eta - v0);

	point3d const &D_eVec = (mat33d::Identity() - mat33d::tensor(u_e, u_e)) * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	return -(2.0 * e * (l0 + l1)) * (d * d * d * d) * r_e / (3 * (l0 + l1 + e) * (l0 + l1 - e) * l0 * l1) +
			(3 * D_e_squared * log((l1 - zeta1) / (l0 - zeta0)) - l1 * zeta1 + l0 * zeta0) * r_e / 6.0 +
			(a_e / 6.0) * (2 * d * d + D_e_squared) * ((u1 - u_e) / (l1 - zeta1) - (u0 - u_e) / (l0 - zeta0)) +
			(a_e / 6.0) * (zeta0 * u0 - zeta1 * u1 + (l0 - l1) * u_e);
}

point3d bh_psi_gradient_in_tangent_plane(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();
	double d = point3d::dot(eta - v0, n); // Convention for the signed distance to the triangle: positive above the triangle, simpler to differentiate (no extra minus sign involved)

	return (1.0 / (8 * M_PI)) * (bh_psi_gradient_in_tangent_plane_subroutine(eta, v0, v1, d, n) + bh_psi_gradient_in_tangent_plane_subroutine(eta, v1, v2, d, n) + bh_psi_gradient_in_tangent_plane_subroutine(eta, v2, v0, d, n));
}

void bh_phi(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2, double &phi0, double &phi1, double &phi2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();
	double d = point3d::dot(eta - v0, n); // Convention for the signed distance to the triangle: positive above the triangle, simpler to differentiate (no extra minus sign involved)

	point3d r0 = point3d::cross(n, v2 - v1);
	r0 /= point3d::dot(r0, v0 - v1);
	point3d r1 = point3d::cross(n, v0 - v2);
	r1 /= point3d::dot(r1, v1 - v2);
	point3d r2 = point3d::cross(n, v1 - v0);
	r2 /= point3d::dot(r2, v2 - v0);

	double gamma0 = point3d::dot(eta - v1, r0);
	double gamma1 = point3d::dot(eta - v2, r1);
	double gamma2 = point3d::dot(eta - v0, r2);

	if (fabs(d) < 0.00000001) {
		// Test if eta is inside the triangle or outside. outside -> phi_k^t = 0
		if (gamma0 < 0 || gamma0 > 1 || gamma1 < 0 || gamma1 > 1 || gamma2 < 0 || gamma2 > 1) {
			phi0 = phi1 = phi2 = 0;
			return;
		}
	}

	point3d t[3] = { v0, v1, v2 };
	double psi_SingleHarmonic = h_psi(eta, &t[0]);
	//	point3d psi_gradient = bh_psi_gradient(eta, v0, v1, v2); // works just as well
	point3d psi_gradient = bh_psi_gradient_in_tangent_plane(eta, v0, v1, v2); // works just as well

	phi0 = -d * point3d::dot(psi_gradient, r0) + 0.5 * d * gamma0 * psi_SingleHarmonic;
	phi1 = -d * point3d::dot(psi_gradient, r1) + 0.5 * d * gamma1 * psi_SingleHarmonic;
	phi2 = -d * point3d::dot(psi_gradient, r2) + 0.5 * d * gamma2 * psi_SingleHarmonic;
}

void h_psi_linear(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2, double &psi0, double &psi1, double &psi2) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();

	point3d t[3] = { v0, v1, v2 };
	double psi_SingleHarmonic = h_psi(eta, &t[0]);
	//	point3d psi_gradient = bh_psi_gradient(eta, v0, v1, v2); // works just as well
	point3d psi_gradient = bh_psi_gradient_in_tangent_plane(eta, v0, v1, v2); // works just as well

	point3d r0 = point3d::cross(n, v2 - v1);
	r0 /= point3d::dot(r0, v0 - v1);
	point3d r1 = point3d::cross(n, v0 - v2);
	r1 /= point3d::dot(r1, v1 - v2);
	point3d r2 = point3d::cross(n, v1 - v0);
	r2 /= point3d::dot(r2, v2 - v0);

	double gamma0 = point3d::dot(eta - v1, r0);
	double gamma1 = point3d::dot(eta - v2, r1);
	double gamma2 = point3d::dot(eta - v0, r2);

	psi0 = -2 * point3d::dot(psi_gradient, r0) + gamma0 * psi_SingleHarmonic;
	psi1 = -2 * point3d::dot(psi_gradient, r1) + gamma1 * psi_SingleHarmonic;
	psi2 = -2 * point3d::dot(psi_gradient, r2) + gamma2 * psi_SingleHarmonic;
}

void bh_phi_grad_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2,
		point3d &phi0_grad, point3d &phi1_grad, point3d &phi2_grad) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();

	point3d t[3] = { v0, v1, v2 };
	double h_psi;
	double h_phi[3];
	h_coordinates_inside_triangle(eta, &t[0], h_psi, h_phi);

	point3d bh_psi_grad = bh_psi_gradient_inside_triangle(eta, v0, v1, v2);

	point3d r0 = point3d::cross(n, v2 - v1);
	r0 /= point3d::dot(r0, v0 - v1);
	point3d r1 = point3d::cross(n, v0 - v2);
	r1 /= point3d::dot(r1, v1 - v2);
	point3d r2 = point3d::cross(n, v1 - v0);
	r2 /= point3d::dot(r2, v2 - v0);

	double gamma0 = point3d::dot(eta - v1, r0);
	double gamma1 = point3d::dot(eta - v2, r1);
	double gamma2 = point3d::dot(eta - v0, r2);

	{
		phi0_grad = -point3d::dot(bh_psi_grad, r0) * n + 0.5 * gamma0 * h_psi * n;
		phi1_grad = -point3d::dot(bh_psi_grad, r1) * n + 0.5 * gamma1 * h_psi * n;
		phi2_grad = -point3d::dot(bh_psi_grad, r2) * n + 0.5 * gamma2 * h_psi * n;
		return;
	}
}

mat33d bh_phi_Hessian_subroutine(point3d const &eta, point3d const &v0, point3d const &v1, double d, point3d const &n, point3d const &r) {
	point3d const &u_e = (v1 - v0).direction();
	point3d const &r_e = point3d::cross(n, u_e);

	double zeta0 = point3d::dot(eta - v0, u_e), zeta1 = point3d::dot(eta - v1, u_e);
	double l0 = (v0 - eta).norm(), l1 = (v1 - eta).norm();
	double e = (v1 - v0).norm();

	point3d u0 = (eta - v0).direction(), u1 = (eta - v1).direction();
	mat33d U0 = (mat33d::Identity() - mat33d::tensor(u0, u0)) / l0;
	mat33d U1 = (mat33d::Identity() - mat33d::tensor(u1, u1)) / l1;

	double a_e = point3d::dot(r_e, eta - v0);

	mat33d const &P_e = (mat33d::Identity() - mat33d::tensor(u_e, u_e));
	point3d const &D_eVec = P_e * (eta - v0);
	double D_e_squared = D_eVec.sqrnorm();

	point3d m0 = (u0 - u_e) / (l0 - zeta0);
	point3d m1 = (u1 - u_e) / (l1 - zeta1);
	mat33d M0 = U0 / (l0 - zeta0) - mat33d::tensor(m0, m0);
	mat33d M1 = U1 / (l1 - zeta1) - mat33d::tensor(m1, m1);

	mat33d MU0r = -MT(r, u0) / (l0 * l0) + 2 * point3d::dot(r, u0) * mat33d::tensor(u0, u0) / (l0 * l0) - point3d::dot(r, u0) * U0 / l0;
	mat33d MU1r = -MT(r, u1) / (l1 * l1) + 2 * point3d::dot(r, u1) * mat33d::tensor(u1, u1) / (l1 * l1) - point3d::dot(r, u1) * U1 / l1;

	mat33d MM0r =
			MU0r / (l0 - zeta0) - MT(mat33d::tensor(u0 - u_e, r), U0) / ((l0 - zeta0) * (l0 - zeta0)) - point3d::dot(u0 - u_e, r) * U0 / ((l0 - zeta0) * (l0 - zeta0)) + 2 * point3d::dot(u0 - u_e, r) * mat33d::tensor(u0 - u_e, u0 - u_e) / ((l0 - zeta0) * (l0 - zeta0) * (l0 - zeta0));
	mat33d MM1r = MU1r / (l1 - zeta1) - MT(mat33d::tensor(u1 - u_e, r), U1) / ((l1 - zeta1) * (l1 - zeta1)) - point3d::dot(u1 - u_e, r) * U1 / ((l1 - zeta1) * (l1 - zeta1)) + 2 * point3d::dot(u1 - u_e, r) * mat33d::tensor(u1 - u_e, u1 - u_e) / ((l1 - zeta1) * (l1 - zeta1) * (l1 - zeta1));

	double rTre = point3d::dot(r_e, r);
	double twoErTre = (2.0 * e * rTre);
	double A1den = (square_t(l0 + l1) - e * e) * l0 * l1;
	double A1 = -twoErTre * (l0 + l1) / A1den;
	point3d A12brac = 2 * (l0 + l1) * l0 * l1 * (u0 + u1) + (square_t(l0 + l1) - e * e) * (l1 * u0 + l0 * u1);
	point3d A1grad = (-twoErTre / A1den) * (u0 + u1) + (twoErTre * (l0 + l1) / (A1den * A1den)) * A12brac;

	mat33d A_H =
			(4 * d * d * d / 3.0) * MT(A1grad, n) + 4 * d * d * A1 * mat33d::tensor(n, n) + (d * d * d * d / 3.0) * ((-twoErTre / A1den) * (U0 + U1) + (twoErTre / (A1den * A1den)) * MT(u0 + u1, A12brac) - (2 * twoErTre * (l0 + l1) / (A1den * A1den * A1den)) * mat33d::tensor(A12brac, A12brac) + (twoErTre * (l0 + l1) / (A1den * A1den)) * (2 * l0 * l1 * mat33d::tensor(u0 + u1, u0 + u1) + 2 * (l0 + l1) * MT(u0 + u1, l1 * u0 + l0 * u1) + (square_t(l0 + l1) - e * e) * (l1 * U0 + l0 * U1 + MT(u0, u1)) + 2 * (l0 + l1) * l0 * l1 * (U0 + U1)));

	mat33d reste = (1.0 / 6.0) * MT(mat33d::tensor(r_e, r), MT(u0 - u1, u_e) + zeta0 * U0 - zeta1 * U1) // A
			+ 0.5 * D_e_squared * MT(mat33d::tensor(r_e, r), M1 - M0) // B
			+ a_e * d * MT(mat33d::tensor(n, r), M1 - M0) // C
			+ MT(mat33d::tensor(D_eVec, r), mat33d::tensor(m1 - m0, r_e)) // D
			+ point3d::dot(m1 - m0, r) * d * MT(n, r_e) // E
			+ (a_e / 6.0) * MT(mat33d::tensor(u_e, r), U0 - U1) // F

			+ (point3d::dot(m1 - m0, r) * a_e) * (mat33d::tensor(n, n) + mat33d::tensor(r_e, r_e)) // reste 1
			+ log((l1 - zeta1) / (l0 - zeta0)) * rTre * P_e // reste 2
			+ (0.5 * D_e_squared * rTre) * (M1 - M0) // reste 3
			+ (rTre / 6.0) * (MT(u0 - u1, u_e) + zeta0 * U0 - zeta1 * U1) // reste 4
			+ ((a_e / 6.0) * (2 * d * d + D_e_squared)) * (MM1r - MM0r) // reste 5
			+ (a_e * point3d::dot(u_e, r) / 6.0) * (U0 - U1) // reste 6
			+ (a_e / 6.0) * (zeta0 * MU0r - zeta1 * MU1r) // reste 7
			;

	return A_H + reste;
}

void bh_phi_Hessian_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2,
		mat33d &phi0_H, mat33d &phi1_H, mat33d &phi2_H) {
	point3d const &n = point3d::cross(v1 - v0, v2 - v0).direction();

	point3d t[3] = { v0, v1, v2 };
	double h_psi;
	point3d h_psi_grad;
	mat33d h_psi_H;
	h_psi_with_derivatives_inside_triangle(eta, &t[0], h_psi, h_psi_grad, h_psi_H);

	mat33d bh_psi_H = bh_psi_Hessian_inside_triangle(eta, v0, v1, v2);

	point3d r0 = point3d::cross(n, v2 - v1);
	r0 /= point3d::dot(r0, v0 - v1);
	point3d r1 = point3d::cross(n, v0 - v2);
	r1 /= point3d::dot(r1, v1 - v2);
	point3d r2 = point3d::cross(n, v1 - v0);
	r2 /= point3d::dot(r2, v2 - v0);

	double Gamma0 = point3d::dot(eta - v1, r0);
	double Gamma1 = point3d::dot(eta - v2, r1);
	double Gamma2 = point3d::dot(eta - v0, r2);

	phi0_H =
			-MT(mat33d::tensor(n, r0), bh_psi_H) + 0.5 * h_psi * MT(n, r0) + 0.5 * Gamma0 * MT(h_psi_grad, n);

	phi1_H =
			-MT(mat33d::tensor(n, r1), bh_psi_H) + 0.5 * h_psi * MT(n, r1) + 0.5 * Gamma1 * MT(h_psi_grad, n);

	phi2_H =
			-MT(mat33d::tensor(n, r2), bh_psi_H) + 0.5 * h_psi * MT(n, r2) + 0.5 * Gamma2 * MT(h_psi_grad, n);
}

void compute_h_bh_coordinates(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2,
		double &h_psi_, double *h_phi_,
		double &bh_psi_, double *bh_phi_) {
	point3d t[3] = { v0, v1, v2 };
	h_coordinates(eta, &(t[0]), h_psi_, h_phi_);

	// Biharmonic part:
	bh_psi_ = bh_psi(eta, v0, v1, v2);
	bh_phi(eta, v0, v1, v2, bh_phi_[0], bh_phi_[1], bh_phi_[2]);
}

void compute_h_bh_coordinates_with_derivatives_inside_triangle(point3d const &eta, point3d const &v0, point3d const &v1, point3d const &v2,
		double &h_psi_, double *h_phi_,
		double &bh_psi_, double *bh_phi_,
		point3d &h_psi_grad_, point3d *h_phi_grad_,
		point3d &bh_psi_grad_, point3d *bh_phi_grad_,
		mat33d &h_psi_H_, mat33d *h_phi_H_,
		mat33d &bh_psi_H_, mat33d *bh_phi_H_) {
	point3d t[3] = { v0, v1, v2 };
	h_coordinates_with_derivatives_inside_triangle(eta, &(t[0]), h_psi_, h_phi_, h_psi_grad_, h_phi_grad_, h_psi_H_, h_phi_H_);

	// Biharmonic part:
	bh_psi_ = bh_psi_inside_triangle(eta, v0, v1, v2);
	bh_phi_[0] = bh_phi_[1] = bh_phi_[2] = 0;
	bh_psi_grad_ = bh_psi_gradient_inside_triangle(eta, v0, v1, v2);
	bh_phi_grad_inside_triangle(eta, v0, v1, v2, bh_phi_grad_[0], bh_phi_grad_[1], bh_phi_grad_[2]);
	bh_psi_H_ = bh_psi_Hessian_inside_triangle(eta, v0, v1, v2);
	bh_phi_Hessian_inside_triangle(eta, v0, v1, v2, bh_phi_H_[0], bh_phi_H_[1], bh_phi_H_[2]);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------------------------------------------------- //

template <class int_t>
void computeCoordinates(
		point3d const &eta,
		std::vector<std::vector<int_t>> const &cage_triangles, std::vector<point3d> const &cage_vertices,
		std::vector<double> &_h_phi, std::vector<double> &_h_psi,
		std::vector<double> &_bh_phi, std::vector<double> &_bh_psi) {
	_h_phi.clear();
	_h_phi.resize(cage_vertices.size(), 0.0);
	_h_psi.clear();
	_h_psi.resize(cage_triangles.size(), 0.0);

	_bh_phi.clear();
	_bh_phi.resize(cage_vertices.size(), 0.0);
	_bh_psi.clear();
	_bh_psi.resize(cage_triangles.size(), 0.0);

	// iterate over the triangles:
	for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
		double h_psi_;
		double h_phi_[3];

		double bh_psi_;
		double bh_phi_[3];

		compute_h_bh_coordinates(eta, cage_vertices[cage_triangles[t][0]], cage_vertices[cage_triangles[t][1]], cage_vertices[cage_triangles[t][2]],
				h_psi_, &(h_phi_[0]), bh_psi_, &(bh_phi_[0]));

		_h_psi[t] = h_psi_;
		for (unsigned int c = 0; c < 3; ++c)
			_h_phi[cage_triangles[t][c]] += h_phi_[c];

		_bh_psi[t] = bh_psi_;
		for (unsigned int c = 0; c < 3; ++c)
			_bh_phi[cage_triangles[t][c]] += bh_phi_[c];
	}
}

template <class int_t>
void computeCoordinatesOnCageVertices(
		std::vector<std::vector<int_t>> const &cage_triangles, std::vector<point3d> const &cage_vertices,
		std::vector<std::vector<double>> &_h_phi, std::vector<std::vector<double>> &_h_psi,
		std::vector<std::vector<double>> &_bh_phi, std::vector<std::vector<double>> &_bh_psi) {
	_h_phi.clear();
	_h_phi.resize(cage_vertices.size());
	_h_psi.clear();
	_h_psi.resize(cage_vertices.size());

	_bh_phi.clear();
	_bh_phi.resize(cage_vertices.size());
	_bh_psi.clear();
	_bh_psi.resize(cage_vertices.size());

	std::vector<point3d> cage_vertex_normals(cage_vertices.size(), point3d(0, 0, 0));
	for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
		point3d t0 = cage_vertices[cage_triangles[t][0]];
		point3d t1 = cage_vertices[cage_triangles[t][1]];
		point3d t2 = cage_vertices[cage_triangles[t][2]];
		point3d tN = point3d::cross(t1 - t0, t2 - t0);
		for (unsigned int c = 0; c < 3; ++c)
			cage_vertex_normals[cage_triangles[t][c]] += tN;
	}
	//#pragma omp parallel for
	for (auto &n : cage_vertex_normals)
		n.normalize();

	//#pragma omp parallel for
	for (int vPIt = 0; vPIt < cage_vertices.size(); ++vPIt) {
		point3d const &eta = cage_vertices[vPIt];

		_h_phi[vPIt].clear();
		_h_phi[vPIt].resize(cage_vertices.size());
		_bh_phi[vPIt].clear();
		_bh_phi[vPIt].resize(cage_vertices.size());
		_h_psi[vPIt].clear();
		_h_psi[vPIt].resize(cage_triangles.size());
		_bh_psi[vPIt].clear();
		_bh_psi[vPIt].resize(cage_triangles.size());

		// iterate over the triangles:
		for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
			double h_psi_ = 0;
			double h_phi_[3] = { 0, 0, 0 };
			double bh_psi_ = 0;
			double bh_phi_[3] = { 0, 0, 0 };

			point3d t0 = cage_vertices[cage_triangles[t][0]];
			point3d t1 = cage_vertices[cage_triangles[t][1]];
			point3d t2 = cage_vertices[cage_triangles[t][2]];

			if (cage_triangles[t][0] == vPIt || cage_triangles[t][1] == vPIt || cage_triangles[t][2] == vPIt) {
				int cornerIt = 0;
				if (cage_triangles[t][1] == vPIt)
					cornerIt = 1;
				if (cage_triangles[t][2] == vPIt)
					cornerIt = 2;
				point3d tri[3] = { t0, t1, t2 };

				h_coordinates_on_triangle_corner(eta, tri,
						h_psi_, h_phi_, cornerIt, cage_vertex_normals[vPIt]);

				bh_psi_ = bh_psi_on_triangle_corner(eta, tri[0], tri[1], tri[2], cornerIt);
				bh_phi_[0] = bh_phi_[1] = bh_phi_[2] = 0.0;
			} else {
				compute_h_bh_coordinates(eta, t0, t1, t2,
						h_psi_, &(h_phi_[0]), bh_psi_, &(bh_phi_[0]));
			}

			_h_psi[vPIt][t] = h_psi_;
			for (unsigned int c = 0; c < 3; ++c)
				_h_phi[vPIt][cage_triangles[t][c]] += h_phi_[c];

			_bh_psi[vPIt][t] = bh_psi_;
			for (unsigned int c = 0; c < 3; ++c)
				_bh_phi[vPIt][cage_triangles[t][c]] += bh_phi_[c];
		}
	}
}

struct SampleOnTri {
	int tri;
	point3d gamma;
	point3d eta;
};

template <class int_t>
void fill_one_sample_per_triangle(
		std::vector<std::vector<int_t>> const &cage_triangles, std::vector<point3d> const &cage_vertices, std::vector<SampleOnTri> &samples_on_cage_triangles) {
	samples_on_cage_triangles.resize(cage_triangles.size());
	for (int tPIt = 0; tPIt < cage_triangles.size(); ++tPIt) {
		point3d triIt[3] = { cage_vertices[cage_triangles[tPIt][0]], cage_vertices[cage_triangles[tPIt][1]], cage_vertices[cage_triangles[tPIt][2]] };
		samples_on_cage_triangles[tPIt].eta = (triIt[0] + triIt[1] + triIt[2]) / 3;
		samples_on_cage_triangles[tPIt].gamma[0] = samples_on_cage_triangles[tPIt].gamma[1] = samples_on_cage_triangles[tPIt].gamma[2] = 1.0 / 3.0;
		samples_on_cage_triangles[tPIt].tri = tPIt;
	}
}

template <class int_t>
void computeHarmonicCoordinatesOnCageTriangles(
		std::vector<std::vector<int_t>> const &cage_triangles, std::vector<point3d> const &cage_vertices,
		std::vector<SampleOnTri> &samples_on_cage_triangles,
		std::vector<std::vector<double>> &_h_phi, std::vector<std::vector<double>> &_h_psi) {
	_h_phi.clear();
	_h_phi.resize(samples_on_cage_triangles.size());
	_h_psi.clear();
	_h_psi.resize(samples_on_cage_triangles.size());

#pragma omp parallel for
	for (int sIt = 0; sIt < samples_on_cage_triangles.size(); ++sIt) {
		auto &sample = samples_on_cage_triangles[sIt];
		point3d const &eta = sample.eta;

		_h_phi[sIt].clear();
		_h_phi[sIt].resize(cage_vertices.size());
		_h_psi[sIt].clear();
		_h_psi[sIt].resize(cage_triangles.size());

		// iterate over the triangles:
		for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
			double h_psi_ = 0;
			double h_phi_[3] = { 0, 0, 0 };
			point3d tri[3] = { cage_vertices[cage_triangles[t][0]], cage_vertices[cage_triangles[t][1]], cage_vertices[cage_triangles[t][2]] };

			point3d t0 = cage_vertices[cage_triangles[t][0]];
			point3d t1 = cage_vertices[cage_triangles[t][1]];
			point3d t2 = cage_vertices[cage_triangles[t][2]];

			if (t == sample.tri) {
				h_coordinates_inside_triangle(eta, tri,
						h_psi_, h_phi_);
			} else {
				h_coordinates(eta, tri,
						h_psi_, h_phi_);
			}

			_h_psi[sIt][t] = h_psi_;
			for (unsigned int c = 0; c < 3; ++c)
				_h_phi[sIt][cage_triangles[t][c]] += h_phi_[c];
		}
	}
}

template <class int_t>
void computeCoordinatesOnCageTriangles(
		std::vector<std::vector<int_t>> const &cage_triangles, std::vector<point3d> const &cage_vertices,
		std::vector<SampleOnTri> &samples_on_cage_triangles,
		std::vector<std::vector<double>> &_h_phi, std::vector<std::vector<double>> &_h_psi,
		std::vector<std::vector<double>> &_bh_phi, std::vector<std::vector<double>> &_bh_psi) {
	_h_phi.clear();
	_h_phi.resize(samples_on_cage_triangles.size());
	_h_psi.clear();
	_h_psi.resize(samples_on_cage_triangles.size());
	_bh_phi.clear();
	_bh_phi.resize(samples_on_cage_triangles.size());
	_bh_psi.clear();
	_bh_psi.resize(samples_on_cage_triangles.size());

	//#pragma omp parallel for
	for (unsigned int sIt = 0; sIt < samples_on_cage_triangles.size(); ++sIt) {
		auto &sample = samples_on_cage_triangles[sIt];
		point3d const &eta = sample.eta;

		_h_phi[sIt].clear();
		_h_phi[sIt].resize(cage_vertices.size());
		_h_psi[sIt].clear();
		_h_psi[sIt].resize(cage_triangles.size());
		_bh_phi[sIt].clear();
		_bh_phi[sIt].resize(cage_vertices.size());
		_bh_psi[sIt].clear();
		_bh_psi[sIt].resize(cage_triangles.size());

		// iterate over the triangles:
		for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
			double h_psi_ = 0;
			double h_phi_[3] = { 0, 0, 0 };
			double bh_psi_ = 0;
			double bh_phi_[3] = { 0, 0, 0 };
			point3d tri[3] = { cage_vertices[cage_triangles[t][0]], cage_vertices[cage_triangles[t][1]], cage_vertices[cage_triangles[t][2]] };

			point3d t0 = cage_vertices[cage_triangles[t][0]];
			point3d t1 = cage_vertices[cage_triangles[t][1]];
			point3d t2 = cage_vertices[cage_triangles[t][2]];

			if (t == sample.tri) {
				h_coordinates_inside_triangle(eta, tri, h_psi_, h_phi_);
				bh_psi_ = bh_psi_inside_triangle(eta, tri[0], tri[1], tri[2]);
				// bh_phi_[3] = { 0,0,0 }; // is 0 on the triangle.
			} else {
				compute_h_bh_coordinates(eta, tri[0], tri[1], tri[2], h_psi_, h_phi_, bh_psi_, bh_phi_);
			}

			_h_psi[sIt][t] = h_psi_;
			_bh_psi[sIt][t] = bh_psi_;
			for (unsigned int c = 0; c < 3; ++c)
				_h_phi[sIt][cage_triangles[t][c]] += h_phi_[c];
			for (unsigned int c = 0; c < 3; ++c)
				_bh_phi[sIt][cage_triangles[t][c]] += bh_phi_[c];
		}
	}
}

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//

template <class int_t>
void computeConstrainedBiharmonicMatrices_13(
		std::vector<std::vector<int_t>> const &cage_triangles, std::vector<point3d> const &cage_vertices,
		Eigen::MatrixXd &C_11, Eigen::MatrixXd &C_12,
		Eigen::MatrixXd &C_21, Eigen::MatrixXd &C_22,
		double gamma_D = 1,
		int triangle_Laplacian_subdiv_count = 3, // default parameter should be posivite : many constraints at the center of the triangle. If 0: a single constraint
		int triangle_Dirichlet_subdiv_count = 3, // default parameter should be posivite : many constraints at the center of the triangle. If 0: a single constraint
		bool sample_Dirichlet_at_cage_vertices = false // default parameter should be false : we don't like sampling constraints at the cage vertices. coordinates are ill-defined here.
) {
	gamma_D = std::max<double>(gamma_D, 0.0); // making sure it's positive, we don't want to do weird things here...

	std::vector<std::vector<double>> h_phi_L, h_psi_L;
	std::vector<SampleOnTri> samples_on_triangles;
	{
		if (triangle_Laplacian_subdiv_count == 0)
			fill_one_sample_per_triangle(cage_triangles, cage_vertices, samples_on_triangles);
		else {
			struct BaryTri {
				point3d v0, v1, v2;
				BaryTri(point3d i0, point3d i1, point3d i2) :
						v0(i0), v1(i1), v2(i2) {}
			};
			std::vector<BaryTri> bary_tris;
			bary_tris.push_back(BaryTri(point3d(1, 0, 0), point3d(0, 1, 0), point3d(0, 0, 1)));

			for (unsigned int subdivIt = 0; subdivIt < triangle_Laplacian_subdiv_count; ++subdivIt) {
				std::vector<BaryTri> bary_tris_old = bary_tris;
				bary_tris.clear();
				for (auto &bt : bary_tris_old) {
					point3d v01 = (bt.v0 + bt.v1) / 2;
					point3d v02 = (bt.v0 + bt.v2) / 2;
					point3d v21 = (bt.v2 + bt.v1) / 2;
					bary_tris.push_back(BaryTri(bt.v0, v01, v02));
					bary_tris.push_back(BaryTri(bt.v1, v21, v01));
					bary_tris.push_back(BaryTri(bt.v2, v02, v21));
					bary_tris.push_back(BaryTri(v01, v21, v02));
				}
			}

			for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
				for (auto &bt : bary_tris) {
					SampleOnTri sample;
					sample.gamma = (bt.v0 + bt.v1 + bt.v2) / 3.0;
					sample.tri = t;
					auto &tric = cage_triangles[t];
					sample.eta = sample.gamma[0] * cage_vertices[tric[0]] + sample.gamma[1] * cage_vertices[tric[1]] + sample.gamma[2] * cage_vertices[tric[2]];
					samples_on_triangles.push_back(sample);
				}
			}
		}
	}

	computeHarmonicCoordinatesOnCageTriangles(cage_triangles, cage_vertices, samples_on_triangles, h_phi_L, h_psi_L);

	Eigen::MatrixXd Phi_L(samples_on_triangles.size(), cage_vertices.size()), Psi_L(samples_on_triangles.size(), cage_triangles.size());
	Eigen::MatrixXd Mass_L = Eigen::MatrixXd::Zero(samples_on_triangles.size(), cage_vertices.size());
	for (unsigned int sIt = 0; sIt < samples_on_triangles.size(); ++sIt) {
		auto &s = samples_on_triangles[sIt];
		for (unsigned int cornerIt = 0; cornerIt < 3; ++cornerIt)
			Mass_L(sIt, cage_triangles[s.tri][cornerIt]) = s.gamma[cornerIt];
		for (unsigned int t2 = 0; t2 < cage_triangles.size(); ++t2) {
			Psi_L(sIt, t2) = h_psi_L[sIt][t2];
		}
		for (unsigned int v = 0; v < cage_vertices.size(); ++v) {
			Phi_L(sIt, v) = h_phi_L[sIt][v];
		}
	}

	std::vector<std::vector<double>> h_phi_V, h_psi_V, bh_phi_V, bh_psi_V;
	if (sample_Dirichlet_at_cage_vertices)
		computeCoordinatesOnCageVertices(cage_triangles, cage_vertices, h_phi_V, h_psi_V, bh_phi_V, bh_psi_V);

	std::vector<SampleOnTri> samples_on_triangles_Dirichlet;
	{
		if (triangle_Dirichlet_subdiv_count == 0)
			fill_one_sample_per_triangle(cage_triangles, cage_vertices, samples_on_triangles_Dirichlet);
		else {
			struct BaryTri {
				point3d v0, v1, v2;
				BaryTri(point3d i0, point3d i1, point3d i2) :
						v0(i0), v1(i1), v2(i2) {}
			};
			std::vector<BaryTri> bary_tris;
			bary_tris.push_back(BaryTri(point3d(1, 0, 0), point3d(0, 1, 0), point3d(0, 0, 1)));

			for (unsigned int subdivIt = 0; subdivIt < triangle_Laplacian_subdiv_count; ++subdivIt) {
				std::vector<BaryTri> bary_tris_old = bary_tris;
				bary_tris.clear();
				for (auto &bt : bary_tris_old) {
					point3d v01 = (bt.v0 + bt.v1) / 2;
					point3d v02 = (bt.v0 + bt.v2) / 2;
					point3d v21 = (bt.v2 + bt.v1) / 2;
					bary_tris.push_back(BaryTri(bt.v0, v01, v02));
					bary_tris.push_back(BaryTri(bt.v1, v21, v01));
					bary_tris.push_back(BaryTri(bt.v2, v02, v21));
					bary_tris.push_back(BaryTri(v01, v21, v02));
				}
			}

			for (unsigned int t = 0; t < cage_triangles.size(); ++t) {
				for (auto &bt : bary_tris) {
					SampleOnTri sample;
					sample.gamma = (bt.v0 + bt.v1 + bt.v2) / 3.0;
					sample.tri = t;
					auto &tric = cage_triangles[t];
					sample.eta = sample.gamma[0] * cage_vertices[tric[0]] + sample.gamma[1] * cage_vertices[tric[1]] + sample.gamma[2] * cage_vertices[tric[2]];
					samples_on_triangles_Dirichlet.push_back(sample);
				}
			}
		}
	}

	unsigned int n_Dirichlet_constraints_on_cage_vertices = (sample_Dirichlet_at_cage_vertices ? cage_vertices.size() : 0);
	unsigned int n_Dirichlet_constraints = n_Dirichlet_constraints_on_cage_vertices + samples_on_triangles_Dirichlet.size();

	std::vector<std::vector<double>> h_phi_V_onTris, h_psi_V_onTris, bh_phi_V_onTris, bh_psi_V_onTris;

	computeCoordinatesOnCageTriangles(cage_triangles, cage_vertices, samples_on_triangles_Dirichlet, h_phi_V_onTris, h_psi_V_onTris, bh_phi_V_onTris, bh_psi_V_onTris);

	Eigen::MatrixXd H_Phi_V(n_Dirichlet_constraints, cage_vertices.size());
	Eigen::MatrixXd BH_Phi_V(n_Dirichlet_constraints, cage_vertices.size());
	Eigen::MatrixXd H_Psi_V(n_Dirichlet_constraints, cage_triangles.size());
	Eigen::MatrixXd BH_Psi_V(n_Dirichlet_constraints, cage_triangles.size());

	Eigen::MatrixXd Mass_V(n_Dirichlet_constraints, cage_vertices.size());
#pragma omp parallel for
	for (int i = 0; i < Mass_V.rows(); ++i) {
		for (int j = 0; j < Mass_V.cols(); ++j) {
			Mass_V(i, j) = 0.0;
		}
	}

	for (int v1 = 0; v1 < n_Dirichlet_constraints_on_cage_vertices; ++v1) {
		Mass_V(v1, v1) = 1.0;
		for (int v2 = 0; v2 < cage_vertices.size(); ++v2) {
			H_Phi_V(v1, v2) = h_phi_V[v1][v2];
			BH_Phi_V(v1, v2) = bh_phi_V[v1][v2];
		}
	}
	for (int v1 = 0; v1 < n_Dirichlet_constraints_on_cage_vertices; ++v1) {
		for (int t2 = 0; t2 < cage_triangles.size(); ++t2) {
			H_Psi_V(v1, t2) = h_psi_V[v1][t2];
			BH_Psi_V(v1, t2) = bh_psi_V[v1][t2];
		}
	}

	for (int v1 = 0; v1 < samples_on_triangles_Dirichlet.size(); ++v1) {
		auto &s = samples_on_triangles_Dirichlet[v1];
		for (int c = 0; c < 3; ++c)
			Mass_V(v1 + n_Dirichlet_constraints_on_cage_vertices, cage_triangles[s.tri][c]) = s.gamma[c];
		for (int v2 = 0; v2 < cage_vertices.size(); ++v2) {
			H_Phi_V(v1 + n_Dirichlet_constraints_on_cage_vertices, v2) = h_phi_V_onTris[v1][v2];
			BH_Phi_V(v1 + n_Dirichlet_constraints_on_cage_vertices, v2) = bh_phi_V_onTris[v1][v2];
		}
	}
	for (int v1 = 0; v1 < samples_on_triangles_Dirichlet.size(); ++v1) {
		for (int t2 = 0; t2 < cage_triangles.size(); ++t2) {
			H_Psi_V(v1 + n_Dirichlet_constraints_on_cage_vertices, t2) = h_psi_V_onTris[v1][t2];
			BH_Psi_V(v1 + n_Dirichlet_constraints_on_cage_vertices, t2) = bh_psi_V_onTris[v1][t2];
		}
	}

	// To compute:
	// A , B , C , C_L , C_D , C11 , C12 , C21 , C22:

	assert(BH_Phi_V.rows() == BH_Psi_V.rows());
	Eigen::MatrixXd A(BH_Phi_V.rows(), BH_Phi_V.cols() + BH_Psi_V.cols()); // A = ( bh_phi_V , bh_psi_V )
#pragma omp parallel for
	for (int r = 0; r < A.rows(); ++r) {
		for (unsigned int c = 0; c < BH_Phi_V.cols(); ++c) {
			A(r, c) = BH_Phi_V(r, c);
		}
		for (unsigned int c = 0; c < BH_Psi_V.cols(); ++c) {
			A(r, BH_Phi_V.cols() + c) = BH_Psi_V(r, c);
		}
	}

	//assert(H_Phi_L.rows() == H_Psi_L.rows());
	//assert(H_Phi_L.rows() == Mass_L.rows());
	//assert(H_Phi_L.cols() == Mass_L.cols());
	Eigen::MatrixXd B(Phi_L.rows(), Phi_L.cols() + Psi_L.cols()); // B = ( phi_L - M_L ; psi_L )
#pragma omp parallel for
	for (int r = 0; r < B.rows(); ++r) {
		for (unsigned int c = 0; c < Phi_L.cols(); ++c) {
			B(r, c) = Phi_L(r, c) - Mass_L(r, c);
		}
		for (unsigned int c = 0; c < Psi_L.cols(); ++c) {
			B(r, Phi_L.cols() + c) = Psi_L(r, c);
		}
	}

	Eigen::MatrixXd C;
	if (gamma_D < 1.0)
		C = gamma_D * (gamma_D * A.transpose() * A + B.transpose() * B).inverse() * A.transpose();
	else
		C = (A.transpose() * A + (1.0 / gamma_D) * B.transpose() * B).inverse() * A.transpose();

	Eigen::MatrixXd C_L(cage_vertices.size(), C.cols());
	Eigen::MatrixXd C_D(cage_triangles.size(), C.cols());
	assert(C_L.rows() + C_D.rows() == C.rows()); // C = (C_L  \\  C_D)

#pragma omp parallel for
	for (int r = 0; r < C_L.rows(); ++r) {
		for (unsigned int c = 0; c < C_L.cols(); ++c) {
			C_L(r, c) = C(r, c);
		}
	}
#pragma omp parallel for
	for (int r = 0; r < C_D.rows(); ++r) {
		for (unsigned int c = 0; c < C_D.cols(); ++c) {
			C_D(r, c) = C(r + C_L.rows(), c);
		}
	}

	C_11 = C_L * (Mass_V - H_Phi_V);
	C_12 = C_D * (Mass_V - H_Phi_V);
	C_21 = C_L * (-H_Psi_V);
	C_22 = C_D * (-H_Psi_V);
}

// this works for coordinates, gradients, Hessians, whatever:
template <class type_t>
void compute_13_blending_from_unconstrained_biharmonics(
		std::vector<type_t> const &h_phi_part, std::vector<type_t> const &h_psi_part,
		std::vector<type_t> const &bh_phi_part, std::vector<type_t> const &bh_psi_part,
		Eigen::MatrixXd &C_11, Eigen::MatrixXd &C_12,
		Eigen::MatrixXd &C_21, Eigen::MatrixXd &C_22,
		std::vector<type_t> &bh_13_alpha_part, std::vector<type_t> &bh_13_beta_part) {
	bh_13_alpha_part.resize(h_phi_part.size()); // no need to initialize to 0
	bh_13_beta_part.resize(h_psi_part.size()); // no need to initialize to 0

#pragma omp parallel for
	for (int cv = 0; cv < bh_13_alpha_part.size(); ++cv) {
		bh_13_alpha_part[cv] = h_phi_part[cv];
		// accumulating the biharmonic part:

		for (int cv2 = 0; cv2 < C_11.rows(); ++cv2) {
			bh_13_alpha_part[cv] += C_11(cv2, cv) * bh_phi_part[cv2];
		}
		for (int ct2 = 0; ct2 < C_12.rows(); ++ct2) {
			bh_13_alpha_part[cv] += C_12(ct2, cv) * bh_psi_part[ct2];
		}
	}

#pragma omp parallel for
	for (int ct = 0; ct < bh_13_beta_part.size(); ++ct) {
		bh_13_beta_part[ct] = h_psi_part[ct];
		// accumulating the biharmonic part:

		for (int cv2 = 0; cv2 < C_21.rows(); ++cv2) {
			bh_13_beta_part[ct] += C_21(cv2, ct) * bh_phi_part[cv2];
		}
		for (int ct2 = 0; ct2 < C_22.rows(); ++ct2) {
			bh_13_beta_part[ct] += C_22(ct2, ct) * bh_psi_part[ct2];
		}
	}
}

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
} //namespace BiharmonicCoordinates3D

#endif // BHC_H
