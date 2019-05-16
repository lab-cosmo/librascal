/**
 * file   cephes.h
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief defines bridge between the cephes C# library and
 * librascal
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef CEPHES_H
#define CEPHES_H


#ifdef __cplusplus
extern "C" {
  namespace rascal {
    namespace cephes {

#endif


  int airy(double x, double *ai, double *aip, double *bi, double *bip);

  double bdtrc(int k, int n, double p);
  double bdtr(int k, int n, double p);
  double bdtri(int k, int n, double y);

  double beta(double a, double b);
  double lbeta(double a, double b);

  double btdtr(double a, double b, double x);

  double cbrt(double x);
  double chbevl(double x, void *, int n);
  double chdtrc(double df, double x);
  double chdtr(double df, double x);
  double chdtri(double df, double y);
  double dawsn(double xx);

  double ellie(double phi, double m);
  double ellik(double phi, double m);
  double ellpe(double x);

  int ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph);
  double ellpk(double x);
  double exp10(double x);
  double exp1m(double x);
  double exp2(double x);

  double expn(int n, double x);

  double fdtrc(int a, int b, double x);
  double fdtr(int a, int b, double x);
  double fdtri(int a, int b, double y);

  int fresnl(double xxa, double *ssa, double *cca);
  double gamma(double x);
  double lgam(double x);
  double lgam_sgn(double x, int *sign);

  double gdtr(double a, double b, double x);
  double gdtrc(double a, double b, double x);
  double gdtri(double a, double b, double y);

  double hyp2f1(double a, double b, double c, double x);
  double hyperg(double a, double b, double x);
  double hyp2f0(double a, double b, double x, int type, double *err);
  double onef2(double a, double b, double c, double x, double *err);
  double threef0(double a, double b, double c, double x, double *err);

  double powi(double x, int nn);

  double i0(double x);
  double i0e(double x);
  double i1(double x);
  double i1e(double x);
  double igamc(double a, double x);
  double igam(double a, double x);
  double igam_fac(double a, double x);
  double igamci(double a, double q);
  double igami(double a, double p);

  double incbet(double aa, double bb, double xx);
  double incbi(double aa, double bb, double yy0);

  double iv(double v, double x);
  double j0(double x);
  double y0(double x);
  double j1(double x);
  double y1(double x);

  double jn(int n, double x);
  double jv(double n, double x);
  double k0(double x);
  double k0e(double x);
  double k1(double x);
  double k1e(double x);
  double kn(int nn, double x);

  double nbdtrc(int k, int n, double p);
  double nbdtr(int k, int n, double p);
  double nbdtri(int k, int n, double p);

  double ndtr(double a);
  double log_ndtr(double a);
  double erfc(double a);
  double erf(double x);
  double ndtri(double y0);

  double pdtrc(int k, double m);
  double pdtr(int k, double m);
  double pdtri(int k, double y);

  double psi(double x);

  double rgamma(double x);
  double round(double x);

  int shichi(double x, double *si, double *ci);
  int sici(double x, double *si, double *ci);

  double radian(double d, double m, double s);
  double sindg(double x);
  double cosdg(double x);
//   double sin(double x);
//   double cos(double x);

  double spence(double x);

  double stdtr(int k, double t);
  double stdtri(int k, double p);

  double yv(double v, double x);

  double tandg(double x);
  double cotdg(double x);

  double log1p(double x);
  double log1pmx(double x);
  double expm1(double x);
  double cosm1(double x);
  double lgam1p(double x);

  double yn(int n, double x);
  double zeta(double x, double q);
  double zetac(double x);

  double smirnov(int n, double e);
  double smirnovi(int n, double p);
  double kolmogorov(double x);
  double kolmogi(double p);

  double lanczos_sum_expg_scaled(double x);

  double owens_t(double h, double a);


#ifdef __cplusplus
    } //cephes
  }// rascal
} // extern "C"
#endif

#endif /* CEPHES_H */
