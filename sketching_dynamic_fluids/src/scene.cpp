#include "scene.hpp"
#include "path_info.hpp"

#include <iostream>
#include <fstream>

using namespace cgp;


// Compute a 3D position of a 2D position given by its screen coordinates for orthographic projection
vec3 unproject(camera_projection_orthographic const& P, mat4 const& camera_view_inverse, vec2 const& p_screen)
{
    // Simple un-project assuming that the viewpoint is an orthogonal projection
    vec4 const p_proj = camera_view_inverse * P.matrix_inverse() * vec4(p_screen, 0.0f, 1.0f);
    return p_proj.xyz();
}

// Compute a 3D position of a 2D position given by its screen coordinates for perspective projection
vec3 unproject(camera_projection_perspective const& P, mat4 const& camera_view_inverse, vec2 const& p_screen)
{
    // Simple un-project assuming that the viewpoint is an orthogonal projection
    vec4 const p_proj = camera_view_inverse * P.matrix_inverse() *  vec4(p_screen, 0.5f, 1.0f);
    vec3 res = p_proj.xyz() / p_proj.w;
    // std::cout << "x = " << p_screen.x << ", y = " << p_screen.y << "\n"
    //           << "x = " << p_proj.x << ", y = " << p_proj.y << "\n"
    //           << p_proj.w << "; " << res[2] << std::endl;
    // std::cout << "P:\n" << P.matrix_inverse() << std::endl;
    // std::cout << "V:\n" << camera_view_inverse << std::endl;
    return res;
}

// vec3 unproject_test(rotation_transform const& orientation, camera_projection_perspective const& P, vec2 const& p_screen) {
//     // return P.matrix_inverse() * vec4(orientation * vec3(p_screen, 0.0f), 0.0f);
//     return vec3(p_screen, 0.0f);
// }


void scene_structure::initialize()
{
    // Uncomment this line if using an orthographic projection
    // camera_projection = camera_projection_orthographic{ -1.1f, 1.1f, -1.1f, 1.1f, -10, 10, window.aspect_ratio() };

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
    float const grid_resolution = 20.0f;
    for (float i = 0; i <= grid_resolution; ++i) {
        float p = i * 2.0f / grid_resolution - 1.0f;

        // std::cout << "cube" << std::endl;
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(p, -1.f)));
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(p, 1.f)));
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(-1.f, p)));
        grid_points.push_back(unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(1.f, p)));
    }
    grid.initialize_data_on_gpu(grid_points);
    grid.display_type = curve_drawable_display_type::Segments;
    
    velocities_drawable = curve_drawable_dynamic_extend();
    velocities_drawable.initialize_data_on_gpu();
    velocities_drawable.display_type = curve_drawable_display_type::Segments;
    
    accelerations_drawable = curve_drawable_dynamic_extend();
    accelerations_drawable.initialize_data_on_gpu();
    accelerations_drawable.display_type = curve_drawable_display_type::Segments;
    accelerations_drawable.color = {0.0f, 1.0f, 0.0f};
}

void scene_structure::display_frame()
{
    // Set the light to the current position of the camera
    environment.light = camera_control.camera_model.position();
    
    draw(fish, environment);

    if (gui.show_grid) {
        draw(grid, environment);
    }

    draw(velocities_drawable, environment);
    draw(accelerations_drawable, environment);
    
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
            // how to remove velocities_drawable?
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
            // std::cout << "curve" << std::endl;
            vec3 const p = unproject(camera_projection, camera_control.camera_model.matrix_frame(), inputs.mouse.position.current);
            // vec3 const p = unproject_test(camera_control.camera_model.orientation(), camera_projection, inputs.mouse.position.current);
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
            // std::cout << "curve" << std::endl;
            vec3 const p = unproject(camera_projection, camera_control.camera_model.matrix_frame(), inputs.mouse.position.current);
            // vec3 const p = unproject_test(camera_control.camera_model.orientation(), camera_projection, inputs.mouse.position.current);
            sketch_drawable[k_sketch].push_back(p);
            
            numarray<vec3> const new_curve = { p };
            sketches.push_back(new_curve);
        }
    }
}

void scene_structure::compute_velocities(numarray<numarray<vec3>>& points, numarray<numarray<vec3>>& vels, curve_drawable_dynamic_extend& vel_drawable, float coef, bool print) {
    // for python optimization
    // if (print) {
        std::ofstream logFile("accelerations");
        logFile << sketches.size() << std::endl;
    // }
    
    for (int i_sketch = 0; i_sketch < sketches.size(); ++i_sketch) {
        // float resolution = 7;
        numarray<vec3> curve = sketches[i_sketch];
        numarray<vec3> point_set = points[i_sketch];
        int n = curve.size();

        // if (print) {
            logFile << n << std::endl;
        // }
        vels.push_back(numarray<vec3>());

        // for (float i_point = 1; i_point < resolution; ++i_point) {
        //     int i_cur = i_point / resolution * n;

        for (int i_cur = 0; i_cur < n; ++i_cur) {

            float const sigma = 0.001;

            vec3 vel(0.0f, 0.0f, 0.0f);

            float sum_w = 0.0;
            int j = i_cur;
            while (j > 0 && j + 1 < n && dist(curve[i_cur], curve[j]) <= 3 * sigma) {
                float d = dist(curve[i_cur], curve[j]);
                float w = exp(-d * d / (sigma * sigma));
                vel += w * (point_set[j] - point_set[j - 1]);
                sum_w += w;
                ++j;
            }

            j = i_cur - 1;
            while (j > 0 && j + 1 < n && dist(curve[i_cur], curve[j - 1]) <= 3 * sigma) {
                float d = dist(curve[i_cur], curve[j - 1]);
                float w = exp(-d * d / (sigma * sigma));
                vel += w * (point_set[j] - point_set[j - 1]);
                sum_w += 1;
                --j;
            }

            vel /= sum_w;
            
            // if (print) {
                logFile << curve[i_cur] << ";" << vel << std::endl;
            // }

            vels[i_sketch].push_back(vel);
            vel_drawable.push_back(curve[i_cur]);
            vel_drawable.push_back(vel * coef + curve[i_cur]);
        }
    }

    logFile.close();
}

void scene_structure::keyboard_event()
{
    camera_control.action_keyboard(environment.camera_view);

    if (inputs.keyboard.is_pressed(GLFW_KEY_T)) {
        std::cout << "T is pressed, let's show tangents" << std::endl;
        compute_velocities(sketches, velocities, velocities_drawable, 20.);
        compute_velocities(velocities, accelerations, accelerations_drawable, 140., true);
    }

    // if (inputs.keyboard.is_pressed(GLFW_KEY_R)) {
    //     velocities_drawable.clear();
    //     accelerations_drawable.clear();
    //     for (int i = 0; i < sketch_drawable.size(); ++i) {
    //         sketch_drawable[i].clear();
    //     }
    //     sketch_drawable.clear();
    //     sketches.clear();
    // }
}

void scene_structure::idle_frame()
{
    camera_control.idle_frame(environment.camera_view);
}

float scene_structure::dist(vec3 a, vec3 b) {
    return norm(a - b);
}