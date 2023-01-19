#include "scene.hpp"
#include "path_info.hpp"

#include <iostream>

using namespace cgp;


// Compute a 3D position of a 2D position given by its screen coordinates for orthographic projection
vec3 unproject(camera_projection_orthographic const& P, mat4 const& camera_view_inverse, vec2 const& p_screen)
{
    // Simple un-project assuming that the viewpoint is an orthogonal projection
    vec4 const p_proj = camera_view_inverse * P.matrix_inverse() * vec4(p_screen, 0.5f, 1.0f);
    return p_proj.xyz();
}

// Compute a 3D position of a 2D position given by its screen coordinates for perspective projection
vec3 unproject(camera_projection_perspective const& P, mat4 const& camera_view_inverse, vec2 const& p_screen)
{
    // Simple un-project assuming that the viewpoint is an orthogonal projection
    vec4 const p_proj = camera_view_inverse * P.matrix_inverse() *  vec4(p_screen, 0.5f, 1.0f);
    return p_proj.xyz() / p_proj.w;
}


void scene_structure::initialize()
{
    // Uncomment this line if using an orthographic projection
    //camera_projection = camera_projection_orthographic{ -1.1f, 1.1f, -1.1f, 1.1f, -10, 10, window.aspect_ratio() };

    camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
    camera_control.set_rotation_axis_y();
    camera_control.look_at({ 0.0f, 0.0f, 3.0f }, { 0,0,0 }, { 0,1,0 });
    global_frame.initialize_data_on_gpu(mesh_primitive_frame());

    std::string mesh_filename = path_info::assets + "fish.obj";
    mesh mesh_model = mesh_load_file_obj(mesh_filename);
    fish.initialize_data_on_gpu(mesh_model);
    fish.model.scaling = 0.2f;
    fish.model.rotation = rotation_transform::from_axis_angle({ 0,1,0 }, -cgp::Pi / 2.0f);
    fish.texture.load_and_initialize_texture_2d_on_gpu(path_info::assets + "fish.png");

    numarray<vec3> grid_points;
    int const grid_resolution = 5.0f;
    for (float i = 1; i < grid_resolution; ++i) {
        std::cout << "i in grid: " << i << std::endl;
        float p = i * 2.0f / grid_resolution - 1.0f;
        std::cout << "Grid index " << p << std::endl;
        // grid_points.push_back({p, -1.f, 0.f});
        // grid_points.push_back({p, 1.f, 0.f});
        // grid_points.push_back({-1.f, p, 0.f});
        // grid_points.push_back({1.f, p, 0.f});
        
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(p, -1.f)));
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(p, 1.f)));
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(-1.f, p)));
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(1.f, p)));
    }
    grid.initialize_data_on_gpu(grid_points);
    grid.display_type = curve_drawable_display_type::Segments;
    
    normals = curve_drawable_dynamic_extend();
    normals.initialize_data_on_gpu();
    normals.display_type = curve_drawable_display_type::Segments;
}

void scene_structure::display_frame()
{
    // Set the light to the current position of the camera
    environment.light = camera_control.camera_model.position();
    
    draw(fish, environment);

    if (gui.show_grid) {
        draw(grid, environment);
    }

    draw(normals, environment);
    
    for (int k = 0; k < sketch_drawable.size(); ++k) {
        draw(sketch_drawable[k], environment);
    }

    if (!gui.sketch_mode) {
        draw(global_frame, environment);
    }
}

void scene_structure::display_gui()
{
    bool cancel = ImGui::Button("Cancel last stroke");
    if (cancel)
    {
        // remove last stroke
        int const N_stroke = sketch_drawable.size();
        if (N_stroke > 0) {
            sketch_drawable[N_stroke - 1].clear();
            sketch_drawable.resize(N_stroke - 1);
            sketches.resize(N_stroke - 1);
        }
    }

    ImGui::Checkbox("Sketch Mode", &gui.sketch_mode);
    ImGui::Checkbox("Show grid", &gui.show_grid);
}

void scene_structure::mouse_move_event()
{
    if (gui.sketch_mode) {
        if (inputs.mouse.click.left) {
            // Add the new clicked position
            int k_sketch = sketch_drawable.size() - 1;
            vec3 const p = unproject(camera_projection, camera_control.camera_model.matrix_frame(), inputs.mouse.position.current);
            sketches[k_sketch].push_back(p);
            sketch_drawable[k_sketch].push_back(p);
            sketch_drawable[k_sketch].color = { 0.5f, 0.5f, 1.0f };
        }
    }
    else {
        if (!inputs.keyboard.shift)
            camera_control.action_mouse_move(environment.camera_view);
    }
}
void scene_structure::mouse_click_event()
{
    if (!inputs.mouse.on_gui) {
        if (gui.sketch_mode && inputs.mouse.click.last_action == last_mouse_cursor_action::click_left)
        {
            // Create new stroke (curve_dynamic_drawable)
            int k_sketch = sketch_drawable.size();
            sketch_drawable.push_back(curve_drawable_dynamic_extend());
            sketch_drawable[k_sketch].initialize_data_on_gpu();

            // Add the new clicked position
            std::cout << inputs.mouse.position.current << std::endl;
            vec3 const p = unproject(camera_projection, camera_control.camera_model.matrix_frame(), inputs.mouse.position.current);
            sketch_drawable[k_sketch].push_back(p);
            
            numarray<vec3> const new_curve = { p };
            sketches.push_back(new_curve);
        }
    }
}

void scene_structure::compute_tangents() {
    float resolution = 7;

    for (int i_sketch = 0; i_sketch < sketches.size(); ++i_sketch) {
        int n = sketches[i_sketch].size();
        for (float i_point = 1; i_point < resolution; ++i_point) {
            int i_cur = i_point / resolution * n;

            std::cout << i_sketch << " " << i_point << " " << i_cur << std::endl;
            std::cout << sketches[i_sketch][i_cur] << std::endl;

            vec3 normal = sketches[i_sketch][i_cur] - sketches[i_sketch][i_cur - 1];
            // float const L = norm(normal);
			// if(L > 1e-6f) {
			// 	   normal /= L;
            //     std::cout << "normalization of normal" << std::endl;
            // }

            // to draw normal
            normals.push_back(sketches[i_sketch][i_cur]);
            normals.push_back(13 * normal + sketches[i_sketch][i_cur]);
        }
    }
}

void scene_structure::keyboard_event()
{
    camera_control.action_keyboard(environment.camera_view);

    if (inputs.keyboard.is_pressed(GLFW_KEY_T)) {
        std::cout << "T is pressed, let's show tangents" << std::endl;
        compute_tangents();
    }
}

void scene_structure::idle_frame()
{
    camera_control.idle_frame(environment.camera_view);
}

