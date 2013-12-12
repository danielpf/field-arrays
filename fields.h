/***************************************************************************
 *                         SOL: fields.c
 *--------------------------------------------------------------------------
 *  CREATED BY: danielf
 *  BEGIN     : Feb 12, 2013
 **************************************************************************/
/* field = one of the arrays associated with a variable.
 *
 * The functions in this directory accomplish three things:
 * - memory management of fields and arrays.
 * - operations on fields such as setting to a constant value or function,
 * sum, filtering, statistics (averages, norms), etc.
 * - initialization
 *
 * Follow this order for function arguments:
 * model_t* mdl, VarLoc_e l, int kvar, int korig, int kdest, double C
 *
 * Explanation:
 * mdl = pointer to model
 * l = location of the variable
 * kvar = index of variable in the variable's enum
 * kfrom = index of input field (eg FLD_VAR, FLD_AN, etc) with values
 * used in the operation.
 * kto = index of output field (whose values we wish to change)
 * C = input constant: becomes an argument when using a field function.
 *
 * Note: Don't use these functions inside iterations because they may
 * not be very fast.
 */
#ifndef FIELDS_H
#define FIELDS_H
/***************************************************************************/
/* HEADER FILES                                                            */
/***************************************************************************/
#include "MATH/math.h"

/***************************************************************************/
/* TYPE DEFINITIONS                                                        */
/***************************************************************************/
/* FIELD STORAGE INDEX */
enum {
  /* equation fields */
  FLD_EQSTART = 0,              /* of start equation fields */
  FLD_VAR = FLD_EQSTART,        /* main variable values */
  /* transient term */
  FLD_OLD,                  /* old values */
  FLD_OLD_OLD,              /* two time-levels before the time level of
                                   var: var_old_old->var_old->var */
  FLD_TMP,                  /* temporary array for the variable,
                                   which some time-stepping schemes need */
  FLD_EQEND = FLD_TMP, /* end of equation fields */
  /* iteration fields */
  FLD_RES,              /* system residuals in each cell's equation */
  FLD_ITDELTA,          /* difference in FLD_VAR from iter. to iter. */
  FLD_TD,               /* time derivative */
  /* errors & solution fields */
  FLD_AN,             /* analytic values */
  FLD_REF,            /* reference field */
  FLD_ERR,            /* errors of FLD_VAR */
  NFIELDS
};

typedef enum {
  CELL_C = 0,    // Cell centered
  FACE_C,        // Face centered
  VERT_C,        // Vertex centered
  NLOCS
} VarLoc_e;

/***************************************************************************/
/* GLOBAL FUNCTIONS                                                        */
/***************************************************************************/

/* Pointer to field function */
typedef double (*FLDfunc_t)(double x[],
    double A, double *ct, double *sig, double grd[]);

/* functions.c */
double constantfunc(double x[], double, double *, double *,  double g[]);
double linearfunc  (double x[], double, double *, double *,  double g[]);
double quadraticfunc  (double x[], double, double *, double *,  double g[]);
double quadratic1Dfunc(double x[], double, double *, double *,  double g[]);
double cubicfunc  (double x[], double, double *, double *,  double g[]);
double normalfunc (double x[], double, double *, double *,  double g[]);
double sinfunc    (double x[], double, double *, double *,  double g[]);
double multsinfunc(double x[], double, double *, double *,  double g[]);

/* alloc.c */
void alloc_initial_mdl_fields(model_t *);
void alloc_non_VAR_FLDs(model_t *mdl, int loc, int FLD, bool FULL);
void dealloc_non_VAR_FLDs(model_t *mdl, int loc, int FLD);
void empty_model_fields(model_t *);
void alloc_and_compute_error_fields(
    model_t *mdl,
    bool do_cells,      /* compute cell centered error fields */
    bool do_faces,
    bool do_verts,
    bool do_grads,      /* force grad error computation */
    bool do_all_vars);  /* compute errors for all available variables */
void dealloc_error_fields(model_t *mdl);

/* impose.c */
void init_model_fields(model_t *mdl);
void put_analytical_velocities(model_t *mdl);
void put_analytical_pressure(model_t *mdl);
void put_analytical_IncNavierStokes(model_t *mdl);
void put_analytical_FLD(model_t *mdl, int loc, int FLD);

/* fields.c */

/* Allocs new arrays for the specified variables in the model_t */
double *
FLD_alloc(model_t *mdl, VarLoc_e l, int kvar, int kfld);
/* Dealloc an array */
double *
FLD_dealloc(model_t*mdl, VarLoc_e l, int kvar, int kfld);
/* Dealloc all fields not set to NULL; set all to NULL. */
void
FLD_dealloc_all_fields(model_t *mdl);
/* print information about fields stored through this module; */
size_t
FLD_print_storage(model_t *mdl);
/* check field for errors */
int
FLD_check_for_not_finite(model_t*mdl, VarLoc_e l, int kvar, int kfld);


/* Operator is a field: ------------------------ */

/* copy the 'orig' field to the 'dest' field. */
double *
FLD_cpy_flds(model_t* mdl, VarLoc_e l, int kvar, int korig, int kdest);
/* swap two fields by exchanging pointers */
double *
FLD_swp_flds(model_t* mdl, VarLoc_e l, int kvar, int korig, int kdest);
/* sets the 'dest' field to it's sum with 'src' */
double *
FLD_sum_flds(model_t* mdl, VarLoc_e l,int kvar,int korig, int kdest, double C);
/* sets 'dest' to it's multiplication with 'src'*/
double *
FLD_mult_flds(model_t* mdl, VarLoc_e l,int kvar,int korig, int kdest, double C);
/* sets 'dest' to it's division with 'src'*/
double *
FLD_div_flds(model_t*mdl, VarLoc_e l, int kvar, int kfrom, int kto, double C);


/* Operator is a constant value: ------------------------ */

/* sets all the cell/face/vert values to constant C. */
double *
FLD_set_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C);
/* sets all the cell/face/vert values to constant C. */
double *
FLD_sum_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C);
/* sets all the cell/face/vert values to constant C. */
double *
FLD_mult_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C);
/* sets all the cell/face/vert values to constant C. */
double *
FLD_clip_below_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C);
/* sets all the cell/face/vert values to constant C. */
double *
FLD_clip_above_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C);


/* Operator is a function: ------------------------ */

/* sets the field to the analytical solution. */
double *
FLD_set_AnalyticalSol(model_t*mdl,VarLoc_e l,int kvar,int kfld);
/* sets the field to the analytical source. */
double *
FLD_set_AnalyticalSrcTerm(model_t*mdl,VarLoc_e l,int kvar,int kfld);
/* Add random field centered in zero with amplitude C */
double *
FLD_add_Random(model_t*mdl,VarLoc_e l,int kvar,int kfld, double C);
/* sets the field to the value returned by function pointer. */
double *
FLD_set_funcval(model_t*mdl, VarLoc_e l, int kvar, int kfld,
    FLDfunc_t func, double A, double *cent, double *sig);
/* sets the field to the gradient computed by function pointer. */
double *
FLD_set_funcgrad(model_t*mdl, VarLoc_e l, int kvar, int kfld,
    FLDfunc_t func, double A, double *cent, double *sig);


/* Field statistics: ------------------------ */

/* calculates the field's statistics (moments and norms). */
void
FLD_compute_stats(model_t*mdl, VarLoc_e l, int kvar, int kfld, int ivpl);

#endif /* FIELDS_H */
/***************************************************************************/
/*- END - END - END - END - END - END - END - END - END - END - END - END -*/
/***************************************************************************/
