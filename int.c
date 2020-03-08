#include <stdlib.h>
#include <stdio.h>

/* Our header */
#include "int.h"
/* We use the HDF5 format for input/output of large arrays */
#include "hdf5.h"
/* We use CLASS to calculate cosmological perturbations */
#include "class.h"

#define NX 3
#define NY 4
#define NZ 9

void getdims(int *nx, int *ny, int *nz) {
     *nx = NX;
     *ny = NY;
     *nz = NZ;
}

void getarr(double *arr) {
     int nx = NX;
     int ny = NY;
     int nz = NZ;

     for (int x=0; x<nx; x++) {
         for (int y=0; y<ny; y++) {
             for (int z=0; z<nz; z++) {
                 // printf("%f\n", arr[z*NY*NX + y*NX + x]);
                 // printf("%f\n", arr[x*NY*NZ + y*NZ + z]);
                 arr[x*NY*NZ + y*NZ + z] = x*y*z;
             }
         }
     }
}


void load_primordial_field(const char *fname) {
  // Open the file containing the primordial fluctuation field
  const hid_t field_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (field_file < 0) {
    printf("Error opening the primordial field file.\n");
  }

  // Open the header group
  hid_t h_grp = H5Gopen(field_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) {
    printf("Error while opening file header\n");
  }

  // Read the physical dimensions of the box
  const hid_t hid_bsz = H5Aexists(h_grp, "BoxSize");
  if (hid_bsz < 0) {
    printf("Error while testing existance of 'BoxSize' attribute\n");
  }

  // double field_dims[3];
  // io_read_attribute(h_grp, "BoxSize", DOUBLE, field_dims);
  // rend->primordial_dims = (double *)malloc(3 * sizeof(double));
  // for (int i = 0; i < 3; i++) {
  //   rend->primordial_dims[i] = field_dims[i];
  // }

  // Now load the actual grid
  h_grp = H5Gopen(field_file, "/Field", H5P_DEFAULT);
  if (h_grp < 0) {
    printf("Error while opening field group\n");
  }

  hid_t h_data = H5Dopen(h_grp, "GaussianRandomField", H5P_DEFAULT);
  hid_t space = H5Dget_space(h_data);

  // The number of dimensions in the dataset (expected 3)
  const int rank = H5Sget_simple_extent_ndims(space);
  if (rank != 3) {
    printf("Incorrect dataset dimensions for primordial field.\n");
  }

  // Find the extent of each dimension (the grid size; not physical size)
  hsize_t grid_dims[rank];
  H5Sget_simple_extent_dims(space, grid_dims, NULL);

  // The grid must be cubic
  if (grid_dims[0] != grid_dims[1] || grid_dims[0] != grid_dims[2]) {
    printf("Primordial grid is not cubic.\n");
  }

  size_t N = grid_dims[0];
  // rend->primordial_grid_N = N;

  // Create a temporary array to read the data
  float grf[N][N][N];
  H5Dread(h_data, H5T_NATIVE_FLOAT, space, space, H5P_DEFAULT, grf);
  H5Dclose(h_data);

  // Allocate memory in the main program
  // rend->primordial_grid = malloc(N * N * N * sizeof(double));

  // Transfer the data to rend->primordial_grid
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      for (size_t k = 0; k < N; k++) {
        // rend->primordial_grid[k + j * N + i * N * N] = grf[i][j][k];
      }
    }
  }

  // printf("The primordial field has dimensions %fx%fx%f U_L^3\n", field_dims[0],
  //         field_dims[1], field_dims[2]);
  printf("The primordial grid has dimensions %zu x %zu x %zu\n",
          (size_t)grid_dims[0], (size_t)grid_dims[1], (size_t)grid_dims[2]);

  // Close the file
  H5Fclose(field_file);
}



void run_class(char *fname_ini, char *fname_pre) {
  /* Define the CLASS structures */
  struct precision pr;  /* for precision parameters */
  struct background ba; /* for cosmological background */
  struct thermo th;     /* for thermodynamics */
  struct perturbs pt;   /* for source functions */
  struct transfers tr;  /* for transfer functions */
  struct primordial pm; /* for primordial spectra */
  struct spectra sp;    /* for output spectra */
  struct nonlinear nl;  /* for non-linear spectra */
  struct lensing le;    /* for lensed spectra */
  struct output op;     /* for output files */
  ErrorMsg errmsg;      /* for CLASS-specific error messages */

  /* If no class .ini file was specified, infer parameters from the cosmology */
  if (fname_ini[0] == '\0') {
    printf("Error: no .ini file for CLASS given.\n");
    exit(0);
  } else {
    /* Otherwise, initialize CLASS with the parameter files */
    int class_argc = 2;
    char *class_argv[] = {"", fname_ini, fname_pre};

    printf("Reading CLASS parameters from '%s'.\n", fname_ini);
    printf("Reading CLASS precision parameters from '%s'.\n", fname_pre);

    if (input_init_from_arguments(class_argc, class_argv, &pr, &ba, &th, &pt,
                                  &tr, &pm, &sp, &nl, &le, &op,
                                  errmsg) == _FAILURE_) {
      printf("Error running input_init_from_arguments \n=>%s\n", errmsg);
    }
  }

  printf("Running CLASS.\n");

  if (background_init(&pr, &ba) == _FAILURE_) {
    printf("Error running background_init \n%s\n", ba.error_message);
  }

  if (thermodynamics_init(&pr, &ba, &th) == _FAILURE_) {
    printf("Error in thermodynamics_init \n%s\n", th.error_message);
  }

  if (perturb_init(&pr, &ba, &th, &pt) == _FAILURE_) {
    printf("Error in perturb_init \n%s\n", pt.error_message);
  }

  // /* Try getting a source */
  // int index_md = pt.index_md_scalars;      // scalar mode
  // int index_ic = 0;                        // index of the initial condition
  // int index_tp = pt.index_tp_delta_ncdm1;  // type of source function
  //
  // /* Size of the perturbations */
  // size_t k_size = pt.k_size[index_md];
  // size_t tau_size = pt.tau_size;
  //
  // /* The number of transfer functions to be read */
  // size_t n_functions = 1;
  //
  // /* Little h, which CLASS uses but Swift doesn't */
  // const double h = ba.h;
  //
  // /* Vector of the wavenumbers */
  // rend->transfer.k_size = k_size;
  // rend->transfer.k = (double *)calloc(k_size, sizeof(double));
  //
  // /* Vector of the conformal times at which the perturbation is sampled */
  // rend->transfer.tau_size = tau_size;
  // rend->transfer.log_tau = (double *)calloc(tau_size, sizeof(double));
  //
  // /* The number of transfer functions to be read */
  // rend->transfer.n_functions = n_functions;
  //
  // /* Vector with the transfer functions T(tau, k) */
  // rend->transfer.delta =
  //     (double *)calloc(n_functions * k_size * tau_size, sizeof(double));
  //
  // /* Read out the perturbation */
  // for (size_t index_tau = 0; index_tau < tau_size; index_tau++) {
  //   for (size_t index_k = 0; index_k < k_size; index_k++) {
  //     /* Convert k from h/Mpc to 1/U_L */
  //     double k = pt.k[index_md][index_k] * h / unit_length_factor;
  //     double p = pt.sources[index_md][index_ic * pt.tp_size[index_md] +
  //                                     index_tp][index_tau * k_size + index_k];
  //
  //     /* Convert transfer functions from CLASS format to CAMB/HeWon/icgen/
  //      *  Eisenstein-Hu format by multiplying by -1/k^2.
  //      */
  //     double T = -p / k / k;
  //
  //     rend->transfer.k[index_k] = k;
  //     rend->transfer.delta[index_tau * k_size + index_k] = T;
  //   }
  //   /* Convert tau from Mpc to U_T */
  //   double tau = pt.tau_sampling[index_tau] * unit_time_factor;
  //   rend->transfer.log_tau[index_tau] = log(tau);
  // }
  //
  // message("The neutrino density is sampled at %zu * %zu points.", k_size,
  //         tau_size);

  /* Pre-empt segfault in CLASS if there is no interacting dark radiation */
  if (ba.has_idr == _FALSE_) {
    pt.alpha_idm_dr = (double *)malloc(0);
    pt.beta_idr = (double *)malloc(0);
  }


  /* Close CLASS again */
  if (perturb_free(&pt) == _FAILURE_) {
    printf("Error in freeing class memory \n%s\n", pt.error_message);
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("Error in thermodynamics_free \n%s\n", th.error_message);
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("Error in background_free \n%s\n", ba.error_message);
  }

  printf("Done with CLASS.\n");
}
