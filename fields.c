/***************************************************************************
 *                         SOL: fields.c                                  
 *--------------------------------------------------------------------------
 *  CREATED BY: danielf
 *  BEGIN     : Feb 12, 2013                                                
 **************************************************************************/
/* This file contains functions for operations on fields */

/***************************************************************************/
/* INCLUDES                                                                */
/***************************************************************************/
#include "SYS/system.h"
#include "UTILS/utils.h"
#include "MODELS/models.h"
#include "MODELS/variables.h"
#include "MODELS/fields/fields.h"
#include "ANALYTICAL/analytical.h"

/***************************************************************************/
/* PRE-PROCESSOR                                                           */
/***************************************************************************/
/* todo: usar um #ifdef DEBUG ? */
#define _check_kvar_kfld(_kvar, _kfld) \
        if(_kvar<0)        fatal_error("fields: kvar is negative");\
        if(_kvar>=V_MAX)   fatal_error("fields: kvar>=V_MAX");\
        if(_kfld<0)        fatal_error("fields: kfld is negative");\
        if(_kfld>=NFIELDS) fatal_error("fields: kfld>=NFIELDS");\

/***************************************************************************/
/* TYPE DEFINITIONS                                                        */
/***************************************************************************/
/* IMPORTANT: Field Operation Structure
 * only visible inside of this file. */
struct op_t {
  int nk;             /* nº of elements - eg ntc */
  int nbulk;
  int *globalk;
  int nvpl;
  int kvar;              /* variable's index (VARPOS) */
  VarLoc_e loc;          /* cell/face or vert*/
  double *xyz;           /* centroid coordinates (access with 3*k) */
  double *to_vals;     /* pointer to 'destination' array */
  double **to_ptr;     /* pointer to pointer to dest */
  double *from_vals;     /* pointer to 'source' array, if needed */

  /* mathematical operations */
  enum {FLD_SUM, FLD_MULT, FLD_DIV, FLD_ULIMIT, FLD_LLIMIT, FLD_SET} math_op;

  /* input */
  enum {
    FLD_NONE,
    FLD_CONST,
    FLD_FLD,
    FLD_FUNCPTR_VAL,   /* function pointer */
    FLD_FUNCPTR_GRAD,
    FLD_FUNC_ANSOL,     /* analytical solution */
    FLD_FUNC_ANSRC,
    FLD_FUNC_FILTER,    /* spacial filter */
    FLD_FUNC_WALLDIST,  /* wall distance */
    FLD_FUNC_RANDOM     /* random values */
  } source_of_values;

  FLDfunc_t f;
  double *f_cent;
  double *f_sig;
  double constant;

  equation_t *eq;       /* todo: maybe it's not needed... */
  SampleMoments_t *mom;
  SampleNorms_t *norms;
  SampleNorms_t *norms_w;
};

/* \brief 'Types' of wall distance measuring    */
typedef enum {
  WDistance_TOP,        /* Distance to y=L */
  WDistance_BOTTOM,     /* Distance to y=0 */
  WDistance_LEFT,       /* Distance to x=0 */
  WDistance_RIGHT,      /* Distance to x=L */
  WDistance_NEAREST
} walldist_e;

/***************************************************************************/
/* FUNCTIONS                                                               */
/***************************************************************************/
void get_field_parameters(struct op_t *o,
        model_t* mdl, VarLoc_e l, int kvar, int kfrom, int kfld);
double *do_field_operation(struct op_t *o,
    model_t*mdl, VarLoc_e l, int kvar, int kfrom, int kto);


/***************************************************************************/
/*-------------------------------------------------------------------------*/
/* MEMORY MANAGEMENT - allocs & deallocs                                   */
/*-------------------------------------------------------------------------*/
/* CREATED BY : Daniel Fazeres                                             */
/* BEGIN      : Feb-2013                                                   */
/* CHANGES    :                                                            */
/***************************************************************************/

/* allocs new arrays for the specified variables in the
 * model_t. Also sets the corresponding equation_t var's pointers */
double *
FLD_alloc(model_t *mdl, VarLoc_e l, int kvar, int kfld)
{
  struct op_t o;

  o.source_of_values = FLD_NONE;
  get_field_parameters(&o, mdl, l, kvar, 0, kfld);

  /* error */
  if (*(o.to_ptr) != NULL) {
    printf("\nWARNING: cannot alloc; this is not a NULL pointer.\n");
    return *(o.to_ptr);
  }

  //**********************************
  /* alloc a new array */
  ALLOC(*(o.to_ptr), double, (o.nk*o.nvpl));
  //**********************************

  return *(o.to_ptr);   /* pointer to newly alloc'd array */
}

/* dealloc respective field */
double *
FLD_dealloc(model_t*mdl, VarLoc_e l, int kvar, int kfld)
{
  struct op_t o;

  o.source_of_values = FLD_NONE;
  get_field_parameters(&o, mdl, l, kvar, 0, kfld);

  DEALLOC(*(o.to_ptr));

  return NULL;
}

/* deallocs all fields not set to NULL */
void
FLD_dealloc_all_fields(model_t *mdl)
{
  int kvar, kfld, kloc;

  for (kvar = 0; kvar < V_MAX; kvar++)
    for (kfld = 0; kfld < NFIELDS; kfld++)
      for (kloc = 0; kloc < NLOCS; kloc++)
        if (mdl->vars[kvar].fld[kloc][kfld] != NULL)
        {
          FLD_dealloc(mdl, kloc, kvar, kfld);
        }

  return;
}

// sets all field pointers to NULL
void
init_FLD_pointers(model_t *mdl)
{
  int kvar, kfld, kloc;

  for (kvar = 0; kvar < V_MAX; kvar++)
    for (kfld = 0; kfld < NFIELDS; kfld++)
      for (kloc = 0; kloc < NLOCS; kloc++)
        mdl->vars[kvar].fld[kloc][kfld] = NULL;

  return;
}

/***************************************************************************/
/*-------------------------------------------------------------------------*/
/* INFORMATION ABOUT FIELDS - statistics, size, etc.                       */
/*-------------------------------------------------------------------------*/
/* CREATED BY : Daniel Fazeres                                             */
/* BEGIN      : 2013                                                       */
/* CHANGES    :                                                            */
/***************************************************************************/

// gets field usage information
size_t
FLD_print_storage(model_t *mdl)
{
  int kvar, kfld, kloc;
  size_t mem = 0;
  int total_nvars = 0;
  int nvars = 0;

  struct op_t o;
  o.source_of_values = FLD_NONE;

  printf("[FIELDS  V_MAX=%d, NFIELDS=%d]\n* ", V_MAX, NFIELDS);
  for (kfld = 0; kfld < NFIELDS; kfld++) {
    nvars = 0;
    for (kvar = 0; kvar < V_MAX; kvar++) {
      for (kloc = 0; kloc < NLOCS; kloc++) {
        if (mdl->vars[kvar].fld[kloc][kfld] != NULL) {
          get_field_parameters(&o, mdl, CELL_C, kvar, kfld, 0);
          mem += sizeof(double) * (size_t)(o.nk*o.nvpl);
          nvars++;
        }
      }
    }
    if (nvars > 0) { printf("#%d->%d,  ", kfld, nvars); }
    total_nvars += nvars;
  }
  /* todo: verificar se é mesmo 'bytes' */
  printf("\n* there are %d field arrays in block %d using %u \"bytes\"\n",
      total_nvars, mdl->myblock, (unsigned int) (mem));

  return mem;
}

/* static function used by other functions to compute norms and moments */
static bool get_field_stats(double *vals, int n,
    SampleMoments_t *mom, SampleNorms_t *norms, double *w)
{
  bool no_stats_were_computed = TRUE;

  /* non-weigthed norms and moments */
  if (mom) {
    compute_sample_moments(vals, n, mom, FALSE, NULL); /* todo weighted? */
    no_stats_were_computed = FALSE;
  }
  if (norms) {
    norms->which_i = compute_sample_norms(vals, n, w,
        &(norms->L1), &(norms->L2), &(norms->Linf), NULL);

    /* averaged */
    norms->L1_avg = norms->L1 / n;
    norms->L2_avg = norms->L2 / sqrt(n);
    no_stats_were_computed = FALSE;
  }
  return no_stats_were_computed;
}

/* calculates the field's statistics (SampleMoments_t and SampleNorms_t)
 * (one for each cartesian component) . */
void FLD_compute_stats(model_t*mdl,
    VarLoc_e l, int kvar, int kfld, int ivpl)
{
  struct op_t o;
  double *temp, *w, sum;

  o.source_of_values = FLD_NONE;
  get_field_parameters(&o, mdl, l, kvar, 0, kfld);

  /* Alloc temporary array for nvpl and face arrays
   * Using the only the bulk faces is needed for the statistics to
   * be correct. */
  if (o.nvpl != 1 || o.globalk != NULL) {
    int globalk, new_size;

    new_size = o.globalk==NULL ? o.nk : o.nbulk;
    ALLOC(temp, double, new_size);

    for (int k = 0; k<new_size; k++) {
      globalk = o.globalk==NULL? k : o.globalk[k];
      temp[k] = o.to_vals[o.nvpl*globalk];
    }
    o.nk = new_size;
  } else {
    if (ivpl > 1) {
      printf("\n* ERROR: requested component of single-valued field.\n\n");
      return;
    }
    temp = o.to_vals; /* array not needed */
  }

  /*----------*/

  /* not weighted */
  get_field_stats(temp, o.nk, o.mom+ivpl, o.norms+ivpl, NULL);

  /* weighted */
  if (l == CELL_C) {
    get_field_stats(temp, o.nk, NULL, o.norms_w+ivpl, kc_vol);
  }
  else if (l == FACE_C) {
    ALLOC(w, double, o.nk);
    for (int k = 0; k<o.nk; k++) { w[k] = normvN(face_nrml+3*k); }
    get_field_stats(temp, o.nk, NULL, o.norms_w+ivpl, w);
    DEALLOC(w);
  }
  else if (l == VERT_C) {
    /* what's the weight?*/
  }

  /* total amount of 'var' in domain - midpoint integration */
  if (l == CELL_C && o.nvpl == 1) {
    sum = 0.0;
    if (kfld == FLD_RES) {
      for (int k = 0; k<o.nk; k++) { sum += o.to_vals[k]; }
    } else {
      for (int k = 0; k<o.nk; k++) { sum += o.to_vals[k] * kc_vol[k]; }
    }
    mdl->vars[kvar].vol_sum[kfld] = sum;
  }

  /* dealloc temporary array */
  if (o.nvpl != 1) { DEALLOC(temp); }

  return;
}

/* returns the number of NaNs or INF values in the field */
int
FLD_check_for_not_finite(model_t*mdl, VarLoc_e l, int kvar, int kfld)
{
  int k, n, error_count;
  struct op_t o;
  o.source_of_values = FLD_NONE;
  get_field_parameters(&o, mdl, l, kvar, 0, kfld);

  n = o.nk*o.nvpl;
  for (error_count = 0, k = 0; k < n; k++)
    if (!is_finite(o.to_vals[k]))
      error_count++;

  return error_count;
}

/***************************************************************************/
/*-------------------------------------------------------------------------*/
/* FIELD OPERATIONS WITH OTHER FIELDS - copying, adding, etc.              */
/*-------------------------------------------------------------------------*/
/* CREATED BY : Daniel Fazeres                                             */
/* BEGIN      : Feb-2013                                                   */
/* CHANGES    :                                                            */
/***************************************************************************/

/* copy the 'from' field to the 'to' field. */
double *
FLD_cpy_flds(model_t*mdl, VarLoc_e l, int kvar, int kfrom, int kto)
{
  struct op_t o;

  // get parameters
  o.math_op = FLD_SET;
  o.source_of_values = FLD_FLD;
  get_field_parameters(&o, mdl, l, kvar, kfrom, kto);

  return memcpy(o.to_vals, o.from_vals, o.nk*o.nvpl*sizeof(double));
}
double *
FLD_swp_flds(model_t*mdl, VarLoc_e l, int kvar, int kfrom, int kto)
{
  struct op_t o;

  // get parameters
  o.math_op = FLD_SET;
  o.source_of_values = FLD_FLD;
  get_field_parameters(&o, mdl, l, kvar, kfrom, kto);

  SWAP(o.from_vals, o.to_vals, double*);

  return o.to_vals;
}
/* sets the 'dest' field to it's sum with 'src' */
double *
FLD_sum_flds(model_t*mdl, VarLoc_e l, int kvar, int kfrom,
    int kto, double C)
{
  struct op_t o;

  o.math_op = FLD_SUM;
  o.source_of_values = FLD_FLD;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, kfrom, kto);
}
/* sets 'dest' to it's multiplication with 'src'*/
double *
FLD_mult_flds(model_t*mdl, VarLoc_e l, int kvar, int kfrom,
    int kto, double C)
{
  struct op_t o;

  o.math_op = FLD_MULT;
  o.source_of_values = FLD_FLD;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, kfrom, kto);
}
/* sets 'dest' to it's division with 'src'*/
double *
FLD_div_flds(model_t*mdl, VarLoc_e l, int kvar, int kfrom,
    int kto, double C)
{
  struct op_t o;

  o.math_op = FLD_DIV;
  o.source_of_values = FLD_FLD;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, kfrom, kto);
}


/***************************************************************************/
/*-------------------------------------------------------------------------*/
/* FIELD OPERATIONS WITH CONSTANTS                                         */
/*-------------------------------------------------------------------------*/
/* CREATED BY : Daniel Fazeres                                             */
/* BEGIN      : Feb-2013                                                   */
/* CHANGES    :                                                            */
/***************************************************************************/

/* sets all the cell/face/vert values to constant C. */
double *
FLD_set_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C)
{
  struct op_t o;

  o.math_op = FLD_SET;
  o.source_of_values = FLD_CONST;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets all the cell/face/vert values to constant C. */
double *
FLD_sum_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C)
{
  struct op_t o;

  o.math_op = FLD_SUM;
  o.source_of_values = FLD_CONST;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets all the cell/face/vert values to constant C. */
double *
FLD_mult_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C)
{
  struct op_t o;

  o.math_op = FLD_MULT;
  o.source_of_values = FLD_CONST;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets all the cell/face/vert values to constant C. */
double *
FLD_clip_below_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C)
{
  struct op_t o;

  o.math_op = FLD_LLIMIT;
  o.source_of_values = FLD_CONST;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets all the cell/face/vert values to constant C. */
double *
FLD_clip_above_val(model_t*mdl, VarLoc_e l, int kvar, int kfld, double C)
{
  struct op_t o;

  o.math_op = FLD_ULIMIT;
  o.source_of_values = FLD_CONST;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}


/***************************************************************************/
/*-------------------------------------------------------------------------*/
/* FIELD OPERATIONS WITH ANALYTICAL FUNCTIONS - check functions.c          */
/*-------------------------------------------------------------------------*/
/* CREATED BY : Daniel Fazeres                                             */
/* BEGIN      : Feb-2013                                                   */
/* CHANGES    :                                                            */
/***************************************************************************/

/*  sets the field to the analytical solution. */
double *
FLD_set_AnalyticalSol(model_t*mdl, VarLoc_e l, int kvar, int kfld)
{
  struct op_t o;

  o.math_op = FLD_SET;
  o.source_of_values = FLD_FUNC_ANSOL;
  o.constant = 1.0;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets the field to the analytical source. */
double *
FLD_set_AnalyticalSrcTerm(model_t*mdl,VarLoc_e l,int kvar,int kfld)
{
  struct op_t o;

  o.math_op = FLD_SET;
  o.source_of_values = FLD_FUNC_ANSRC;
  o.constant = 1.0;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets the field to the value returned by function pointer. */
double *
FLD_set_funcval(model_t*mdl, VarLoc_e l, int kvar, int kfld,
    FLDfunc_t func, double A, double *cent, double *sig)
{
  struct op_t o;

  o.math_op = FLD_SET;
  o.source_of_values = FLD_FUNCPTR_VAL;
  o.f = func;
  o.f_cent = cent;
  o.f_sig = sig;
  o.constant = A;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* sets the field to the gradient computed by function pointer. */
double *
FLD_set_funcgrad(model_t*mdl, VarLoc_e l, int kvar, int kfld,
    FLDfunc_t func, double A, double *cent, double *sig)
{
  struct op_t o;

  o.math_op = FLD_SET;
  o.source_of_values = FLD_FUNCPTR_GRAD;
  o.f = func;
  o.f_cent = cent;
  o.f_sig = sig;
  o.constant = A;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}
/* Adds random noise */
double *
FLD_add_Random(model_t*mdl,VarLoc_e l,int kvar,int kfld, double C)
{
  struct op_t o;

  o.math_op = FLD_SUM;
  o.source_of_values = FLD_FUNC_RANDOM;
  o.constant = C;
  return do_field_operation(&o, mdl, l, kvar, 0, kfld);
}


/***************************************************************************/
/*-------------------------------------------------------------------------*/
/* INTERNAL FUNCTIONS        - only visible inside of this file            */
/*-------------------------------------------------------------------------*/
/* CREATED BY : Daniel Fazeres                                             */
/* BEGIN      : Feb-2013                                                   */
/* CHANGES    :                                                            */
/***************************************************************************/
void calc(int math_op, double from, double *to);

/* get_field_parameters:
 * - gets the properties of the storage array associated with the target val.
 * - inserts all the parameters related to this field in the 'o' struct.
 */
void
get_field_parameters(struct op_t *o, model_t* mdl, VarLoc_e l,
    int kvar, int kfrom, int kfld)
{
  _check_kvar_kfld(kvar, kfld);
  Var_t *var = mdl->vars+kvar;

  o->loc = l;
  o->kvar = kvar;
  o->eq = var->eq;

  /* todo: bulkfaces */

  /* centroids */
  switch(l)  {
  case CELL_C:
    o->nk = mdl->ntc;
    o->nbulk = o->nk;
    o->globalk = NULL;
    o->xyz = &(kc3_cent[mdl->fc]);
    break;
  case FACE_C:
    o->nk = mdl->nf;
    if (bulkfaces != NULL) {
      o->nbulk = nfbulk;
      o->globalk = bulkfaces;
    } else {
      o->nbulk = o->nk;
      o->globalk = NULL;
    }
    o->xyz = &(face_cent[mdl->ff]);
    break;
  case VERT_C:
    o->nk = mdl->nv;
    o->nbulk = o->nk;
    o->globalk = NULL;
    o->xyz = &(vert_xyz[mdl->fv]);
    /* todo bulk verts? */
    break;
  default: fatal_error("fields: invalid var location"); break;
  }
  o->nvpl = mdl->vars[kvar].nvpl;

  o->to_vals = var->fld[l][kfld];
  o->to_ptr  = &(var->fld[l][kfld]);
  o->mom       = var->mom[l][kfld];
  o->norms     = var->norms[l][kfld];
  o->norms_w   = var->norms_w[l][kfld];

  if (o->source_of_values == FLD_FLD) {
    _check_kvar_kfld(kvar, kfrom);
    o->from_vals = var->fld[l][kfrom];
  }

  return;
}

/* does all the necessary operations to the target field
 * and returns a pointer to it's array.
 *
 * Daniel Fazeres 2013
 * */
double *
do_field_operation(struct op_t *o, model_t*mdl, VarLoc_e l,
    int kvar, int kfrom, int kto)
{
  int k, i;
  double val[V_MAX];
  double C = o->constant;

  // get parameters
  get_field_parameters(o, mdl, l, kvar, kfrom, kto);
  if (o->to_vals == NULL) { fatal_error("fields: this array is NULL\n"); }

  // do operation
  switch (o->source_of_values) {
  case FLD_FLD:
    for (k = 0; k < (o->nk*o->nvpl); k++) {
      calc(o->math_op, C * o->from_vals[k], o->to_vals+k);
    }
    break;

  case FLD_CONST:
    if (C == 0.0 && (o->math_op == FLD_SET || o->math_op == FLD_MULT)) {
      memset(o->to_vals, 0, sizeof(double)*o->nk*o->nvpl); /* faster */
    } else {
      for (k = 0; k < (o->nk*o->nvpl); k++) {
        calc(o->math_op, C, o->to_vals+k);
      }
    }
    break;

  case FLD_FUNC_ANSOL:
    for (k = 0; k < o->nk; k++) {
      AnalyticalSolution(o->xyz+3*k, V_MAX, val);
      if (o->loc == FACE_C && kvar == V_vel) {
        /* COMPUTE Un (must have velocity components in 'val')
         * nvpl is 1;  */
        val[kvar] = dotprodvN(&(face_nrml[3*k]), val + V_vel);
      }
      // Calculate vector elements individually
      for (i = 0; i < o->nvpl; i++) {
        calc(o->math_op, val[kvar+i], o->to_vals + o->nvpl*k + i);
      }
    }
    break;

  case FLD_FUNC_ANSRC:
    for (k = 0; k < o->nk; k++) {
      AnalyticalSourceTerm(&(o->xyz[3*k]), V_MAX, val-V_SRC);
      for (i = 0; i < o->nvpl; i++) {
        calc(o->math_op, val[kvar+i], o->to_vals + o->nvpl*k + i);
      }
    }
    break;

  case FLD_FUNC_RANDOM:
    for (k = 0; k < o->nk; k++) {
      /* random function centered in zero. */
      val[kvar] = C * ((double)rand()/(double)RAND_MAX - 0.5);
      for (i = 0; i < o->nvpl; i++) {
        calc(o->math_op, val[kvar], o->to_vals + o->nvpl*k + i);
      }
    }
    break;

  case FLD_FUNCPTR_VAL:
    /* Use the value returned by the function-pointer. */
    for (k = 0; k < o->nk; k++) {
      val[0] = (o->f)(o->xyz+3*k, o->constant, o->f_cent, o->f_sig, val+1);
      for (i = 0; i < o->nvpl; i++) {
        calc(o->math_op, val[0], o->to_vals + o->nvpl*k + i);
      }
    }
    break;

  case FLD_FUNCPTR_GRAD:
    /* Use the function-pointer input-vector */
    for (k = 0; k < o->nk; k++) {
      o->f(o->xyz+3*k, o->constant, o->f_cent, o->f_sig, val);
      for (i = 0; i < o->nvpl; i++) {
        calc(o->math_op, val[i], o->to_vals + o->nvpl*k + i);
      }
    }
    break;

  case FLD_FUNC_WALLDIST:
    for (k = 0; k < o->nk*o->nvpl; k++) {
      calc(o->math_op, o->xyz[1], o->to_vals);
    }
    break;

  default: fatal_error("unrecognized option"); break;
  }

  return o->to_vals;
}

/* sum, multiply, divide, etc */
void calc(int math_op, double orig, double *dest)
{
  /* put value in dest */
  switch (math_op) {
  case FLD_SUM:  *dest += orig; break;
  case FLD_MULT: *dest *= orig; break;
  case FLD_DIV:
    if(is_equal(orig,0.0)) { printf("ERROR:fields: don't divide by zero"); }
    *dest = 1E50;
    break;
  case FLD_SET:  *dest = orig; break;
  case FLD_ULIMIT: if (is_lower(*dest, orig)) { *dest = orig; } break;
  case FLD_LLIMIT: if (is_greater(*dest, orig)){ *dest = orig; } break;
  default: fatal_error("fields: invalid mathematical calculation"); break;
  }

  return;
}
/***************************************************************************/
/*- END - END - END - END - END - END - END - END - END - END - END - END -*/
/***************************************************************************/
