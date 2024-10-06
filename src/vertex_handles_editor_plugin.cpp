#include "vertex_handles_editor_plugin.h"
#include <godot_cpp/core/class_db.hpp>

void VertexHandlesEditorPlugin::_bind_methods() {
}

bool VertexHandlesEditorPlugin::_forward_canvas_gui_input(const Ref<InputEvent> &p_event) {
    if (p_event->is_class("InputEventMouseMotion")) {
        update_overlays();
        return true;
    }
    return false;
}
