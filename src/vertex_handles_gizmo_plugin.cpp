#include "vertex_handles_gizmo_plugin.h"

#include "mesh_morph_3d.h"

#include <godot_cpp/classes/standard_material3d.hpp>

#include <godot_cpp/classes/geometry3d.hpp>

#include <godot_cpp/classes/mesh_instance3d.hpp>

#include <godot_cpp/variant/typed_array.hpp>

VertexHandlesGizmoPlugin::VertexHandlesGizmoPlugin() {
    create_material("main", Color(0, 1, 1), false, true);
    create_handle_material("handles", false);
}

bool VertexHandlesGizmoPlugin::_has_gizmo(Node3D *p_node) const {
    return cast_to<MeshMorph3D>(p_node);
}

String VertexHandlesGizmoPlugin::_get_gizmo_name() const {
    return "VertexHandlesGizmo";
}

String VertexHandlesGizmoPlugin::_get_handle_name(const Ref<EditorNode3DGizmo> &p_gizmo, int32_t p_handle_id, bool p_secondary) const {
    switch (p_handle_id) {
        case 0: return "Radius";
        case 1: return "Width";
        default: return String::num(p_handle_id);
    }
}

Variant VertexHandlesGizmoPlugin::_get_handle_value(const Ref<EditorNode3DGizmo> &p_gizmo, int32_t p_handle_id, bool p_secondary) const {
    Node3D *node3d = p_gizmo->get_node_3d();
    switch (p_handle_id) {
        case 0: return node3d->get("radius");
        case 1: return node3d->get("width");
        default: return Variant();
    }
}

void VertexHandlesGizmoPlugin::_set_handle(const Ref<EditorNode3DGizmo> &p_gizmo, int32_t p_handle_id, bool p_secondary, Camera3D *p_camera, const Vector2 &p_screen_pos) {
    Node3D *node3d = p_gizmo->get_node_3d();
    Transform3D gt = node3d->get_global_transform();
    Transform3D gi = gt.affine_inverse();
    Vector3 ray_from = p_camera->project_ray_origin(p_screen_pos);
    Vector3 ray_dir = p_camera->project_ray_normal(p_screen_pos);
    ray_from = gi.xform(ray_from);
    ray_dir = gi.basis.xform(ray_dir);
    Array points = node3d->get("points");
    PackedVector3Array point_array = points[0];
    Plane plane(p_camera->get_camera_transform().basis.get_column(2), point_array[p_handle_id]);
    Vector3 intersection_point;
    if (!plane.intersects_ray(ray_from, ray_dir, &intersection_point)) {
        return;   
    }
    point_array[p_handle_id] = intersection_point;
    points[0] = point_array;
    node3d->set("points", points);
    p_gizmo->_redraw();
    node3d->call("update_mesh");
}




void VertexHandlesGizmoPlugin::_redraw(const Ref<EditorNode3DGizmo> &p_gizmo) {
    p_gizmo->clear();
    Node3D *node3d = p_gizmo->get_node_3d();
    PackedVector3Array lines;
    lines.push_back(Vector3(1, 1, -1));
    lines.push_back(Vector3(0, 0, 0));
    PackedVector3Array handles;
    Array points = node3d->get("points");
    for (int i = 0; i < points.size(); ++i) {
        Array point_array = points[i];
        for (int j = 0; j < point_array.size(); ++j) {
            handles.push_back(point_array[j]);
        }
    }
    bool billboard = false;
    p_gizmo->add_handles(handles, get_material("handles", p_gizmo), Array(), billboard);
}