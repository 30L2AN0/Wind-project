# Wind-project
From 2D strokes to a 3D fluid representation

The code cosists of two parts in directories `sketching_dynamic_fluids` and `optimization`. The first one is GUI where user can draw a curve, build velocities/accelerations from it by pressing "T" (for "tangents") and print the results to file by pressing "P" (for "print"). The second part is computing that builds the flow according to the curve using data from the file received by the first program. The flow is constructed frame by frame due to its change in time.

One can find the example results in the directory `examples`.
