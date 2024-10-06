@tool
class_name VertexHandles extends Node3D

@onready var mesh = get_parent().mesh

@export var points := [] : set = _set_points

func _ready() -> void:
	assert(get_parent() is MeshInstance3D)
	get_parent().mesh = _to_array_mesh(get_parent().mesh)
	mesh = get_parent().mesh
	points = []
	for i in range(0,mesh.get_surface_count()):
		var arrays = mesh.surface_get_arrays(i)
		points.push_back(arrays[Mesh.ARRAY_VERTEX])

func _update_mesh():
	if is_node_ready() and not points.is_empty():
		mesh = get_parent().mesh
		var surface_arrays := []
		for i in range(0,mesh.get_surface_count()):
			var arrays = mesh.surface_get_arrays(i)
			surface_arrays.push_back(arrays)
		mesh.clear_surfaces()
		for i in range(0,surface_arrays.size()):
			surface_arrays[i][Mesh.ARRAY_VERTEX] = PackedVector3Array( points[i] )
			mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, surface_arrays[i])
		get_node("../../MeshMorph3D").apply_deformation_to_children()

func _set_points(value):
	points = value

func _to_array_mesh(_mesh:Mesh) -> ArrayMesh:
	var surface_arrays := []
	for i in range(0,_mesh.get_surface_count()):
		var arrays = _mesh.surface_get_arrays(i)
		surface_arrays.push_back(arrays)
	_mesh = ArrayMesh.new()
	for i in range(0,surface_arrays.size()):
		_mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, surface_arrays[i])
	return _mesh
