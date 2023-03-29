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
    vec4 const p_proj = camera_view_inverse * P.matrix_inverse() * vec4(p_screen, 0.5f, 1.0f);
    vec3 res = p_proj.xyz() / p_proj.w;
    // std::cout << "x = " << p_screen.x << ", y = "  < p_screen.y << "\n"
    //           << "x = " << p_proj.x << ", y = " << p_proj.y << "\n"
    //           << p_proj.w << "; " << res[2] << std::endl;
    // std::cout << "P:\n" << P.matrix_inverse() << std::endl;
    return res;
}

vec2 scene_structure::get_cell(vec3 p) {
    int x = (p.x - left_bot.x) / grid_step;
    int y = (p.y - left_bot.y) / grid_step;
    std::cout << "Point " << p << " is in the cell: " << x << ", " << y << std::endl;
    return {x, y};
}

void scene_structure::create_square_object(int h, vec3 p) {
    solid_boundary.clear();
    solid_inner.clear();

    vec2 p_cell = get_cell(p);
    int x = p_cell.x;
    int y = p_cell.y;

    for (int dx = -h; dx <= h; ++dx) {
        for (int dy = -h; dy <= h; ++dy) {
            vec2 body_point = vec2(x + dx, y + dy);
            bool boundary = false;
            vec2 normal = vec2(0.f, 0.f);

            if (x + dx == x - h) {
                boundary = true;
                normal.x = -1.;
            }
            if (x + dx == x + h) {
                boundary = true;
                normal.x = 1.;
            }
            if (p.y + dy == p.y - h) {
                boundary = true;
                normal.y = -1.;
            }
            if (p.y + dy == p.y + h) {
                boundary = true;
                normal.y = 1.;
            }
            
            if (boundary) {
                solid_boundary.push_back({body_point, normal});
                std::cout << "normal of " << body_point << " is: " << normal << std::endl;
            } else {
                solid_inner.push_back(body_point);
            }
        }
    }
}

numarray<vec3> scene_structure::solid_boundary_3d() {
    numarray<vec3> ret;
    for (auto& p : solid_boundary) {
        ret.push_back(vec3(p.first, 0.0f));
    }
    return ret;
}

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

    left_bot = unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(-1.f, -1.f));
    left_top = unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(-1.f, 1.f));
    right_top = unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(1.f, 1.f));
    right_bot = unproject(camera_projection, camera_control.camera_model.matrix_frame(), vec2(1.f, -1.f));

    grid_step = (right_top[0] - left_top[0]) / grid_resolution;

    // grid_cell_size = (grid_points[1] - grid_points[0]).y;
	// points_on_grid.resize(grid_resolution, grid_resolution);
    // points_on_grid.fill(std::set<int>());
    
    velocities_drawable = curve_drawable_dynamic_extend();
    velocities_drawable.initialize_data_on_gpu();
    velocities_drawable.display_type = curve_drawable_display_type::Segments;
    
    accelerations_drawable = curve_drawable_dynamic_extend();
    accelerations_drawable.initialize_data_on_gpu();
    accelerations_drawable.display_type = curve_drawable_display_type::Segments;
    accelerations_drawable.color = {0.0f, 1.0f, 0.0f};

    // create_square_object(1, {0, 0, 0});
    // solid_drawable.initialize_data_on_gpu(solid_boundary_3d());
}

void scene_structure::display_frame()
{
    // Set the light to the current position of the camera
    environment.light = camera_control.camera_model.position();
    
    draw(fish, environment);
    // draw(solid_drawable, environment);

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

            dtimes[k_sketch].push_back(timer.update());
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

            dtimes.push_back({numarray<float>()});
            dtimes[k_sketch].push_back(0.0);
            timer.start();
        }
    }
}

/* for each stroke compute the speed of changing the vectors in points,
   store it in vels and add to vel_drawable do draw */
void scene_structure::compute_velocities(numarray<numarray<vec3>>& points,
                                         numarray<numarray<vec3>>& vels,
                                         curve_drawable_dynamic_extend& vel_drawable,
                                         float coef)
{
    for (int i_sketch = 0; i_sketch < sketches.size(); ++i_sketch) {
        // float resolution = 7;
        numarray<vec3> curve = sketches[i_sketch];
        numarray<vec3> point_set = points[i_sketch];
        int n = curve.size();
        vels.push_back(numarray<vec3>());

        // for (float i_point = 1; i_point < resolution; ++i_point) {
        //     int i_point = i_point / resolution * n;

        for (int i_point = 0; i_point < n; ++i_point) {
            float const sigma = 0.001;

            vec3 vel(0.0f, 0.0f, 0.0f);

            float sum_w = 0.0;
            int j = i_point;
            while (j > 0 && j < n && dist(curve[i_point], curve[j]) <= 3 * sigma) {
                float d = dist(curve[i_point], curve[j]);
                float w = exp(-d * d / (sigma * sigma));
                vel += w * (point_set[j] - point_set[j - 1]) / dtimes[i_sketch][j];
                sum_w += w;
                ++j;
            }

            j = i_point - 1;
            while (j > 0 && j < n && dist(curve[i_point], curve[j - 1]) <= 3 * sigma) {
                float d = dist(curve[i_point], curve[j - 1]);
                float w = exp(-d * d / (sigma * sigma));
                vel += w * (point_set[j] - point_set[j - 1]) / dtimes[i_sketch][j];
                sum_w += w;
                --j;
            }

            vel /= sum_w;

            vels[i_sketch].push_back(vel);
            vel_drawable.push_back(curve[i_point]);
            vel_drawable.push_back(vel * coef + curve[i_point]);
        }
    }
}

void scene_structure::keyboard_event()
{
    camera_control.action_keyboard(environment.camera_view);

    if (inputs.keyboard.is_pressed(GLFW_KEY_T)) {
        std::cout << "T is pressed, let's show tangents" << std::endl;
        // computes all the velocities every time, not only last
        // not implemented yet: smooth_curves();
        compute_velocities(sketches, velocities, velocities_drawable, 0.1);
        compute_velocities(velocities, accelerations, accelerations_drawable, 0.005);
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

    if (inputs.keyboard.is_pressed(GLFW_KEY_S)) {
        std::cout << "S is pressed, let's print info about the solid" << std::endl;

        vec3 p = vec3(0, 0, 0);
        std::cout << "At position:\n" << p << std::endl;
        create_square_object(1, p);
        std::cout << "Boundary cells:\n";
        for (auto& v : solid_boundary) {
            std::cout << v.first << std::endl;
        }
        std::cout << "Inner cells:\n";
        for (vec2 v : solid_inner) {
            std::cout << v << std::endl;
        }
    }

    /* print info taken from sketches to file
       the format is:

       #num of body boundary cells
       #num of body inner cells
       #num of curves
       list of normals on boundary cells
       * then for each curve: *
        #num of points in a curve
        list of cells with body boundaries
        list of cells with body inner points
        body velocity
        body acceleration
    
     */
    if (inputs.keyboard.is_pressed(GLFW_KEY_P)) {
        std::cout << "P is pressed, let's print all points to the file" << std::endl;
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "../../accelerations_%H_%M_%S_%d_%m_%Y");
        std::ofstream logFile(oss.str(), std::ios::out);

        logFile << solid_boundary.size() << std::endl;
        logFile << solid_inner.size() << std::endl;
        logFile << sketches.size() << std::endl;

        for (auto& v : solid_boundary) {
            logFile << v.second << std::endl;
        }

        for (int i_sketch = 0; i_sketch < sketches.size(); ++i_sketch) {
            logFile << sketches[i_sketch].size() << std::endl;

            for (int i_point = 0; i_point < sketches[i_sketch].size(); ++i_point) {
                if (std::isnan(accelerations[i_sketch][i_point][0])) {
                    continue;
                }

                create_square_object(1, sketches[i_sketch][i_point]);

                for (auto& v : solid_boundary) {
                    logFile << v.first << std::endl;
                }

                logFile << solid_inner.size() << std::endl;
                for (vec2 v : solid_inner) {
                    logFile << v << std::endl;
                }

                // logFile << sketches[i_sketch][i_point] << ";" << velocities[i_sketch][i_point] << ";" << accelerations[i_sketch][i_point] << std::endl;
                logFile << velocities[i_sketch][i_point] << std::endl;
                logFile << accelerations[i_sketch][i_point] << std::endl;
            }
        }
        logFile.close();
    }
}

void scene_structure::idle_frame()
{
    camera_control.idle_frame(environment.camera_view);
}

float scene_structure::dist(vec3 a, vec3 b) {
    return norm(a - b);
}