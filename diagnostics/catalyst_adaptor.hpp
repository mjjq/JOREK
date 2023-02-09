#ifndef CATALYST_ADAPTOR_H_
#define CATALYST_ADAPTOR_H_

#if USE_CATALYST

extern "C" {
void catalyst_adaptor();

void catalyst_adaptor_initialise(char *a_catalyst_script);

void catalyst_adaptor_execute(int *a_step_index, double *a_time);

void catalyst_adaptor_finalise();

void catalyst_get_grid_params(int *a_nsub, int *a_n_elements);

void catalyst_interp_grid(int *a_nnos, int *a_nel, float *a_coords_RZst,
                          int *a_cell_points);
}

#endif /* USE_CATALYST */
#endif