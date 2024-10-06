#ifndef VERTEX_HANDLES_GIZMO_PLUGIN_H
#define VERTEX_HANDLES_GIZMO_PLUGIN_H

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/camera3d.hpp>
#include <godot_cpp/classes/editor_node3d_gizmo_plugin.hpp>

using namespace godot;

class VertexHandlesGizmoPlugin : public EditorNode3DGizmoPlugin {
    GDCLASS(VertexHandlesGizmoPlugin, EditorNode3DGizmoPlugin);
protected:
    static void _bind_methods() {}
    
public:
    VertexHandlesGizmoPlugin();
    virtual bool _has_gizmo(Node3D *p_node) const override;
    virtual String _get_gizmo_name() const override;
    virtual String _get_handle_name(const Ref<EditorNode3DGizmo> &p_gizmo, int32_t p_handle_id, bool p_secondary) const override;
    virtual Variant _get_handle_value(const Ref<EditorNode3DGizmo> &p_gizmo, int32_t p_handle_id, bool p_secondary) const override;
    virtual void _set_handle(const Ref<EditorNode3DGizmo> &p_gizmo, int32_t p_handle_id, bool p_secondary, Camera3D *p_camera, const Vector2 &p_screen_pos) override;
    virtual void _redraw(const Ref<EditorNode3DGizmo> &p_gizmo) override;
};

#endif // VERTEX_HANDLES_GIZMO_PLUGIN_H
