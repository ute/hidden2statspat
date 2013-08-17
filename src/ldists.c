// lsdist. c: locally scaled and Euclidean distance 
#include <R.h>
#include <Rinternals.h>

SEXP lspairdist(SEXP(x), SEXP(y), SEXP(la))
{
  const char *outnames[] = {"lsd", "eud", ""};
  int n = length(x);
  double *px, *py, *prola, *plsd, *peud;
  double xi, yi, rolai, dx, dy, dd, ldd;

  SEXP rola = PROTECT(duplicate(la));
  SEXP out = PROTECT (mkNamed(VECSXP, outnames));
  SEXP lsd = SET_VECTOR_ELT(out, 0, allocMatrix(REALSXP, n, n));
  SEXP eud = SET_VECTOR_ELT(out, 1, allocMatrix(REALSXP, n, n));
 
  px = REAL(x);
  py = REAL(y);
  prola = REAL(rola);
  plsd = REAL(lsd);
  peud = REAL(eud);
  
  for(int i = 0; i < n; i++) { prola[i] = sqrt(prola[i])/2; }
  for(int i = 0; i < n; i++) {
    xi = px[i]; yi = py[i]; rolai = prola[i];
    plsd[i * (n + 1)] = 0;
    peud[i * (n + 1)] = 0;
    for(int j = i+1; j < n; j++){
      dx = px[j] - xi;
      dy = py[j] - yi;
      dd = sqrt(dx * dx + dy * dy); 
      peud[i + n * j] = dd;
      peud[j + n * i] = dd;   
      ldd = dd * (prola[j] + rolai);
      plsd[i + n * j] = ldd;
      plsd[j + n * i] = ldd;   
    }
  }

  UNPROTECT(2);

  return(out);
}

