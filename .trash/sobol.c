#include <stdio.h>
#include <gsl/gsl_qrng.h>

int
main (void)
{
  int i;
  gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, 1);

  for (i = 0; i < 10; i++)
    {
      double v;
      gsl_qrng_get(q, &v);
      printf ("%.5f\n", v);
    }

  gsl_qrng_free (q);
  return 0;
}
