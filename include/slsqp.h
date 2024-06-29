#ifndef SLSQP_H
#define SLSQP_H

#include "nlopt.h"
#include "nlopt-util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


nlopt_result nlopt_slsqp(unsigned n, nlopt_func f, void *f_data,
			 unsigned m, nlopt_constraint *fc,
			 unsigned p, nlopt_constraint *h,
			 const double *lb, const double *ub,
			 double *x, double *minf,
			 nlopt_stopping *stop);

void nnls_(double *a, int *mda, int *m, int *
           n, double *b, double *x, double *rnorm, double *w,
           double *z__, int *indx, int *mode);

void ldp_(double *g, int *mg, int *m, int *n,
          double *h__, double *x, double *xnorm, double *w,
          int *indx, int *mode);


void lsi_(double *e, double *f, double *g,
          double *h__, int *le, int *me, int *lg, int *mg,
          int *n, double *x, double *xnorm, double *w, int *
          jw, int *mode);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
