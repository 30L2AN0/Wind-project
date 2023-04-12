#pragma once


#include "cgp/cgp.hpp"
#include "environment.hpp"


using cgp::mesh_drawable;


struct gui_parameters {
	bool sketch_mode = true;
	bool show_grid = false;
};

// The structure of the custom scene
struct scene_structure : cgp::scene_inputs_generic {
	
	// ****************************** //
	// Elements and shapes of the interface
	// ****************************** //
	camera_controller_orbit_euler camera_control;

	camera_projection_perspective camera_projection; // uncomment this line for projective perspective
	// camera_projection_orthographic camera_projection; // uncomment this line for orthographic projection

	window_structure window;

	mesh_drawable global_frame;          // The standard global frame
	environment_structure environment;   // Standard environment controler
	input_devices inputs;                // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                  // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	cgp::timer_basic timer;

	// Store the curve sketched on screen. 
	// Each new stroke (continuous click+motion of the mouse) is a new element of the buffer
	cgp::numarray<cgp::curve_drawable_dynamic_extend> sketch_drawable;
	// save the curve elements
	numarray<numarray<vec3>> sketches;

	// solid positions and normals
	numarray<std::pair<vec2, vec2>> solid_boundary;
	numarray<vec2> solid_inner;

	// grid info
	float grid_cell_size;
	curve_drawable grid;
	// grid_2D<std::set<int>> points_on_grid;
	vec3 left_top, left_bot, right_top, right_bot;
	double grid_step;
    float grid_resolution;

	numarray<numarray<float>> dtimes;
	numarray<numarray<vec3>> velocities;
	numarray<numarray<vec3>> accelerations;
	curve_drawable_dynamic_extend velocities_drawable;
	curve_drawable_dynamic_extend accelerations_drawable;
	curve_drawable solid_drawable;

	mesh_drawable fish;

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();    // Standard initialization to be called before the animation loop
	void display_frame(); // The frame display to be called within the animation loop
	void display_gui();   // The display of the GUI, also called within the animation loop

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

	void compute_velocities(numarray<numarray<vec3>>&, numarray<numarray<vec3>>&, curve_drawable_dynamic_extend&, float);
	// void smooth_the_curve();
	float dist(vec3 a, vec3 b);
	void create_square_object(int h, vec3 p);

	vec2 get_cell(vec3 p);
	numarray<vec3> solid_boundary_3d();
};
