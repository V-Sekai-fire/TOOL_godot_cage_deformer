#ifndef VERTEX_HANDLES_EDITOR_PLUGIN_H
#define VERTEX_HANDLES_EDITOR_PLUGIN_H

#include <godot_cpp/classes/editor_plugin.hpp>
#include <godot_cpp/classes/input_event.hpp>

#include "vertex_handles_gizmo_plugin.h"

using namespace godot;

class VertexHandlesEditorPlugin : public EditorPlugin {
    GDCLASS(VertexHandlesEditorPlugin, EditorPlugin);

    Ref<VertexHandlesGizmoPlugin> gizmo_plugin;
protected:
    static void _bind_methods();

public:
    virtual bool _forward_canvas_gui_input(const Ref<InputEvent> &p_event) override;
    virtual void _enter_tree() override {
        gizmo_plugin.instantiate();
        add_node_3d_gizmo_plugin(gizmo_plugin);
    }
    virtual void _exit_tree() override {
        remove_node_3d_gizmo_plugin(gizmo_plugin);
        gizmo_plugin.unref();
    }
};

#endif // VERTEX_HANDLES_EDITOR_PLUGIN_H
