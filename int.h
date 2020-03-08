void getdims(int *nx, int* ny, int *nz);
void getarr(double *arr);
void load_primordial_field(const char *fname);
void run_class(char *fname_ini, char *fname_pre);

/* Data structure for the renderer, generating perturbation theory grids */
struct renderer {

  /*! Array containing the primordial fluctuations on a grid */
  double *primordial_grid;
  double *primordial_dims;
  size_t primordial_grid_N;
};
