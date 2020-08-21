#include "mex.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

typedef enum { LEFT=0, RIGHT=1, UP=2, DOWN=3 } dir_t;

// Work out extents from a single mask shape
int get_mask_extents(int *M, int l, int *istart, int *istop, int *jstart, int *jstop)
{
  int i, j, f, p;

  p = 0;
  f = (l - 1) / 2;
  for (i=0; i < l; i++) {
    j = 0;
    istart[i] = j - f;
    while ((j < l) && M[j*l + i] <= 0.0) {
      j++;
      istart[i] = j - f;
    }
    istop[i] = istart[i] - 1;
    while ((j < l) && M[j*l + i] > 0.0) {
      istop[i] = j - f;
      j++;
    }
    j = 0;
    jstart[i] = j - f;
    while ((j < l) && M[i*l + j] <= 0.0) {
      j++;
      jstart[i] = j - f;
    }
    jstop[i] = jstart[i] - 1;
    while ((j < l) && M[i*l + j] > 0.0) {
      jstop[i] = j - f;
      j++;
    }
    p += jstop[i] - jstart[i] + 1;
  }

  return p;
}

// Efficient opening and/or closing using a fixed mask
void rankopen_2d(int *iarray, int w, int h, int l, int *M, int mthresh, bool dual, int ig1, int *oarray1, int ig2, int *oarray2)
{
  int i, j, k, ii, jj, f, p, m, e, *istart, *istop, *jstart, *jstop, maxval, minval;
  int *marray, *oldv, *newv, old_m1, old_m2, imin, imax, jmin, jmax, kmin, kmax;
  int loop, igr1, igr2, ig, *ia, *oa1, *oa2, *tarray1 = NULL, *tarray2 = NULL, *parray1 = NULL, *parray2 = NULL;
  bool forward, right_first;
  int mlim, minc, mtest, mmax, mmin, mk;
  bool threshold_mask = (mthresh > 0);
  dir_t dir;

  // return if silly values
  f = (l - 1) / 2;
  if ((f > (h / 2)) || (f > (w / 2))) return;

  // work out starting and ending points based on mask
  istart = new int[l];
  istop = new int[l];
  jstart = new int[l];
  jstop = new int[l];
  oldv = new int[l];
  newv = new int[l];
  p = get_mask_extents(M, l, istart, istop, jstart, jstop);
  imin = 0; while ((istart[imin] > istop[imin]) && (imin < (l - 1))) imin++;
  imax = l - 1; while ((istart[imax] > istop[imax]) && (imax > 0)) imax--;
  jmin = 0; while ((jstart[jmin] > jstop[jmin]) && (jmin < (l - 1))) jmin++;
  jmax = l - 1; while ((jstart[jmax] > jstop[jmax]) && (jmax > 0)) jmax--;
  if ((jmax - jmin) < (imax - imin)) right_first = true;
  else right_first = false;

  // Sort out centiles - zero-offset in C
  igr1 = p + 1 - ig1;
  igr2 = p + 1 - ig2;
  ig1 -= 1;
  ig2 -= 1;
  igr1 -= 1;
  igr2 -= 1;

  // Possibly additional arrays for thresholding
  if ((oarray1 != NULL) && threshold_mask) {
    parray1 = new int[w*h];
    for (i=0; i < w*h; i++) parray1[i] = 0;
  }
  if ((oarray2 != NULL) && threshold_mask) {
    parray2 = new int[w*h];
    for (i=0; i < w*h; i++) parray2[i] = 0;
  }

  // Possibly several passes
  for (loop=0; loop < 3; loop++) {
    if (loop == 0) {
      ia = iarray;
      forward = true;
      if (!dual) {
        oa1 = oarray1;
        oa2 = oarray2;
      } else {
        // will need a second pass to do each inverse
        // create temporary arrays
        if (oarray1 != NULL) tarray1 = new int[w*h];
        if (oarray2 != NULL) tarray2 = new int[w*h];
        oa1 = tarray1;
        oa2 = tarray2;
      }
    } else if (loop == 1) {
      if (oa1 == oarray1) break; // stop if this has already been done
      forward = false;
      ia = tarray1;
      oa1 = oarray1;
      oa2 = NULL;
    } else if (loop == 2) {
      if (oarray2 == NULL) break; // no further array to process
      ia = tarray2;
      oa1 = NULL;
      oa2 = oarray2;
    }

    // initialise at first pixel
    if (forward || !threshold_mask) {
      minval = ia[0];
      maxval = minval;
      for (i=1; i < w*h; i++) {
        if (ia[i] > maxval) maxval = ia[i];
        if (ia[i] < minval) minval = ia[i];
      }
    } else if (loop == 1) {
      minval -= 1; // to account for blanked values from forward pass
    }

    // create space for histogram, and fill with centre at (0, -1)
    e = maxval - minval + 1;
    if (loop > 0) delete[] marray;
    marray = new int[e];
    for (m=0; m < e; m++) marray[m] = 0;
    for (i=jmin; i <= jmax; i++) {
      for (j=jstart[i]; j <= jstop[i]; j++) {
        if ((i < (1 + f)) && (j < 0)) marray[ia[0] - minval]++;
        else if (i < (1 + f)) marray[ia[j] - minval]++;
        else if (j < 0) marray[ia[(i - 1 - f)*w] - minval]++;
        else marray[ia[(i - 1 - f)*w + j] - minval]++;
      }
    }

    // start walking through image using forwards then backwards approach to make the addition
    // of extra values to the sorted array more efficient
    old_m1 = 0;
    old_m2 = 0;
    i = 0;
    j = 0;
    dir = DOWN;
    kmin = imin;
    kmax = imax;
    while ((i < h) && (j < w)) {

      // Work out the list of old and new values we need to correct the sorted array with
      switch (dir) {

        case DOWN:
          for (k=kmin; k <= kmax; k++) {
            if ((j + k - f) < 0) jj=0;
            else if ((j + k - f) >= w) jj=w - 1;
            else jj=j + k - f;
            if ((istop[k] + i) >= h) newv[k] = ia[(h - 1)*w + jj];
            else if ((istop[k] + i) <= 0) newv[k] = ia[0 + jj];
            else newv[k] = ia[(i + istop[k])*w + jj];
            if ((istart[k] + i) <= 0) oldv[k] = ia[0 + jj];
            else if ((istart[k] + i) >= h) oldv[k] = ia[(h - 1)*w + jj];
            else oldv[k] = ia[(i + istart[k] - 1)*w + jj];
          }
          break;

        case UP:
          for (k=kmin; k <= kmax; k++) {
            if ((j + k - f) < 0) jj=0;
            else if ((j + k - f) >= w) jj=w - 1;
            else jj=j + k - f;
            if ((istart[k] + i) <= -1) newv[k] = ia[0 + jj];
            else if ((istart[k] + i) >= (h - 1)) newv[k] = ia[(h - 1)*w + jj];
            else newv[k] = ia[(i + istart[k])*w + jj];
            if ((istop[k] + i) >= (h - 1)) oldv[k] = ia[(h - 1)*w + jj];
            else if ((istop[k] + i) <= -1) oldv[k] = ia[0 + jj];
            else oldv[k] = ia[(i + istop[k] + 1)*w + jj];
          }
          break;

        case RIGHT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            if ((jstop[k] + j) >= w) newv[k] = ia[ii*w + w - 1];
            else if ((jstop[k] + j) <= 0) newv[k] = ia[ii*w + 0];
            else newv[k] = ia[ii*w + j + jstop[k]];
            if ((jstart[k] + j) <= 0) oldv[k] = ia[ii*w + 0];
            else if ((jstart[k] + j) >= w) oldv[k] = ia[ii*w + w - 1];
            else oldv[k] = ia[ii*w + j + jstart[k] - 1];
          }
          break;

        case LEFT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            if ((jstart[k] + j) <= 0) newv[k] = ia[ii*w + 0];
            else if ((jstart[k] + j) >= (w - 1)) newv[k] = ia[ii*w + w - 1];
            else newv[k] = ia[ii*w + j + jstart[k]];
            if ((jstop[k] + j) >= (w - 1)) oldv[k] = ia[ii*w + w - 1];
            else if ((jstop[k] + j) <= -1) oldv[k] = ia[ii*w + 0];
            else oldv[k] = ia[ii*w + j + jstop[k] + 1];
          }
          break;

      }

      // Now adjust for these new and old values

      // Remove old values and add new values to histogram
      for (k=kmin; k <= kmax; k++) {
        marray[oldv[k] - minval]--;
        marray[newv[k] - minval]++;
      }

      // Find centile from histogram
      if (oa1 != NULL) {

        if (forward) {
          ig = ig1;
        } else {
          ig = igr1;
        }

        if (!threshold_mask) {

          // Forward or inverse in simple case
          if (old_m1 < (e / 2)) {
            k=0;
            for (m=0; m < e; m++) {
              k += marray[m];
              if (k > ig) {
                old_m1 = m;
                oa1[i*w + j] = (m + minval);
                break;
              }
            }
          } else {
            k=p - 1;
            for (m=(e - 1); m >= 0; m--) {
              k -= marray[m];
              if (k < ig) {
                old_m1 = m;
                oa1[i*w + j] = (m + minval);
                break;
              }
            }
          }

        } else if (forward) {

          // Forward operation when using threshold_mask option
          // which has to be in the right direction according to ig, and stops at mtest
          mtest = iarray[i*w + j];
          mmax = (mtest + mthresh) - minval;
          if (mmax > (e - 1)) mmax = e - 1;
          mmin = (mtest - mthresh) - minval;
          if (mmin < 0) mmin = 0;
          mk = 0;
          if (ig <= (p / 2)) {
            k = 0;
            for (m=0; m < mmin; m++) mk += marray[m];
            for (; m <= mmax; m++) {
              k += marray[m];
              if (k > ig) {
                k = -2;
                break;
              }
            }
          } else {
            k = p - 1;
            mk = 0;
            for (m=(e - 1); m > mmax; m--) mk += marray[m];
            for (; m >= mmin; m--) {
              k -= marray[m];
              if (k < ig) {
                k = -2;
                break;
              }
            }
          }
          if (k == -2) {
            if (mk == 0) {
              oa1[i*w + j] = (m + minval);
            } else {
              parray1[i*w + j] = (m + minval);
              oa1[i*w + j] = minval - 1;
            }
          } else {
            parray1[i*w + j] = mtest;
            oa1[i*w + j] = minval - 1;
          }

        } else {

          // Inverse operation when using threshold_mask option
          // which has to be in the right direction according to ig, and stops at mtest
          mtest = ia[i*w + j] - minval;
          if (mtest == 0) mtest = parray1[i*w + j] - minval;
          if (ig < (p / 2)) {
            k=0;
            for (m=1; m < mtest; m++) {
              k += marray[m];
              if (k > ig) break;
            }
          } else {
            k=p - 1;
            for (m=(e - 1); m > mtest; m--) {
              k -= marray[m];
              if (k < ig) break;
            }
          }
          old_m1 = m;
          oa1[i*w + j] = (m + minval);

        }
      }

      if (oa2 != NULL) {

        if (forward) {
          ig = ig2;
        } else {
          ig = igr2;
        }

        if (!threshold_mask) {

          // Forward or inverse in simple case
          if (old_m2 < (e / 2)) {
            k=0;
            for (m=0; m < e; m++) {
              k += marray[m];
              if (k > ig) {
                old_m2 = m;
                oa2[i*w + j] = (m + minval);
                break;
              }
            }
          } else {
            k=p - 1;
            for (m=(e - 1); m >= 0; m--) {
              k -= marray[m];
              if (k < ig) {
                old_m2 = m;
                oa2[i*w + j] = (m + minval);
                break;
              }
            }
          }

        } else if (forward) {

          // Forward operation when using threshold_mask option
          // which has to be in the right direction according to ig, and stops at mtest
          mtest = iarray[i*w + j];
          mmax = (mtest + mthresh) - minval;
          if (mmax > (e - 1)) mmax = e - 1;
          mmin = (mtest - mthresh) - minval;
          if (mmin < 0) mmin = 0;
          mk = 0;
          if (ig <= (p / 2)) {
            k = 0;
            for (m=0; m < mmin; m++) mk += marray[m];
            for (; m <= mmax; m++) {
              k += marray[m];
              if (k > ig) {
                k = -2;
                break;
              }
            }
          } else {
            k = p - 1;
            mk = 0;
            for (m=(e - 1); m > mmax; m--) mk += marray[m];
            for (; m >= mmin; m--) {
              k -= marray[m];
              if (k < ig) {
                k = -2;
                break;
              }
            }
          }
          if (k == -2) {
            if (mk == 0) {
              oa2[i*w + j] = (m + minval);
            } else {
              parray2[i*w + j] = (m + minval);
              oa2[i*w + j] = minval - 1;
            }
          } else {
            parray2[i*w + j] = mtest;
            oa2[i*w + j] = minval - 1;
          }

        } else {

          // Inverse operation when using threshold_mask option
          // which has to be in the right direction according to ig, and stops at mtest
          mtest = ia[i*w + j] - minval;
          if (mtest == 0) mtest = parray2[i*w + j] - minval;
          if (ig < (p / 2)) {
            k=0;
            for (m=1; m < mtest; m++) {
              k += marray[m];
              if (k > ig) break;
            }
          } else {
            k=p - 1;
            for (m=(e - 1); m > mtest; m--) {
              k -= marray[m];
              if (k < ig) break;
            }
          }
          old_m2 = m;
          oa2[i*w + j] = (m + minval);

        }
      }

      // increment index for next time
      if (right_first) {
        switch (dir) {
          case DOWN:
            if (j == (w - 1)) {
              dir = LEFT;
              kmin = jmin;
              kmax = jmax;
              j--;
            } else {
              dir = RIGHT;
              kmin = jmin;
              kmax = jmax;
              j++;
            }
            break;
          case RIGHT:
            if (j == (w - 1)) {
              dir = DOWN;
              kmin = imin;
              kmax = imax;
              i++;
            } else {
              j++;
            }
            break;
          case LEFT:
            if (j == 0) {
              dir = DOWN;
              kmin = imin;
              kmax = imax;
              i++;
            } else {
              j--;
            }
            break;
        }
      } else {
        switch (dir) {
          case RIGHT:
            if (i == (h - 1)) {
              dir = UP;
              kmin = imin;
              kmax = imax;
              i--;
            } else {
              dir = DOWN;
              kmin = imin;
              kmax = imax;
              i++;
            }
            break;
          case DOWN:
            if (i == (h - 1)) {
              dir = RIGHT;
              kmin = jmin;
              kmax = jmax;
              j++;
            } else {
              i++;
            }
            break;
          case UP:
            if (i == 0) {
              dir = RIGHT;
              kmin = jmin;
              kmax = jmax;
              j++;
            } else {
              i--;
            }
            break;
        }
      }

    }
  }

  // free variables
  delete[] oldv;
  delete[] newv;
  delete[] marray;
  delete[] istart;
  delete[] istop;
  delete[] jstart;
  delete[] jstop;
  if (tarray1 != NULL) delete[] tarray1;
  if (tarray2 != NULL) delete[] tarray2;
  if (parray1 != NULL) delete[] parray1;
  if (parray2 != NULL) delete[] parray2;
}

// Efficient opening and/or closing using a fixed mask
void rankopen_2d(double *iarray, int w, int h, int l, int *M, double mthresh, bool dual, int ig1, double *oarray1, int ig2, double *oarray2)
{
  int i, j, k, ii, jj, f, p, m, *istart, *istop, *jstart, *jstop, minc;
  int left, right, old_m1, old_m2, imin, imax, jmin, jmax, kmin, kmax, loop, igr1, igr2, ig;
  double *marray, mtest, mlim, *oldv, *newv, *ia, *oa1, *oa2, *tarray1 = NULL, *tarray2 = NULL, *parray1 = NULL, *parray2 = NULL;
  bool forward, right_first;
  bool threshold_mask = (mthresh > 0);
  dir_t dir;

  // return if silly values
  f = (l - 1) / 2;
  if ((f > (h / 2)) || (f > (w / 2))) return;

  // work out starting and ending points based on mask
  istart = new int[l];
  istop = new int[l];
  jstart = new int[l];
  jstop = new int[l];
  oldv = new double[l];
  newv = new double[l];
  p = get_mask_extents(M, l, istart, istop, jstart, jstop);
  imin = 0; while ((istart[imin] > istop[imin]) && (imin < (l - 1))) imin++;
  imax = l - 1; while ((istart[imax] > istop[imax]) && (imax > 0)) imax--;
  jmin = 0; while ((jstart[jmin] > jstop[jmin]) && (jmin < (l - 1))) jmin++;
  jmax = l - 1; while ((jstart[jmax] > jstop[jmax]) && (jmax > 0)) jmax--;
  if ((jmax - jmin) < (imax - imin)) right_first = true;
  else right_first = false;

  // Sort out centiles - zero-offset in C
  igr1 = p + 1 - ig1;
  igr2 = p + 1 - ig2;
  ig1 -= 1;
  ig2 -= 1;
  igr1 -= 1;
  igr2 -= 1;

  // Possibly additional arrays for thresholding
  if ((oarray1 != NULL) && threshold_mask) {
    parray1 = new double[w*h];
    for (i=0; i < w*h; i++) parray1[i] = 0;
  }
  if ((oarray2 != NULL) && threshold_mask) {
    parray2 = new double[w*h];
    for (i=0; i < w*h; i++) parray2[i] = 0;
  }

  // Possibly several passes
  for (loop=0; loop < 3; loop++) {
    if (loop == 0) {
      ia = iarray;
      forward = true;
      if (!dual) {
        oa1 = oarray1;
        oa2 = oarray2;
      } else {
        // will need a second pass to do each inverse
        // create temporary arrays
        if (oarray1 != NULL) tarray1 = new double[w*h];
        if (oarray2 != NULL) tarray2 = new double[w*h];
        oa1 = tarray1;
        oa2 = tarray2;
      }
    } else if (loop == 1) {
      if (oa1 == oarray1) break; // stop if this has already been done
      forward = false;
      ia = tarray1;
      oa1 = oarray1;
      oa2 = NULL;
    } else if (loop == 2) {
      if (oarray2 == NULL) break; // no further array to process
      ia = tarray2;
      oa1 = NULL;
      oa2 = oarray2;
    }

    // initialise at first pixel
    // create space for the sorted list of length p, and fill it up with centre at (0,-1)
    if (loop > 0) delete[] marray;
    marray = new double[p];
    m = 0;
    for (i=jmin; i <= jmax; i++) {
      for (j=jstart[i]; j <= jstop[i]; j++) {
        if ((i < (1 + f)) && (j < 0)) marray[m++] = ia[0];
        else if (i < (1 + f)) marray[m++] = ia[j];
        else if (j < 0) marray[m++] = ia[(i - 1 - f)*w];
        else marray[m++] = ia[(i - 1 - f)*w + j];
      }
    }
    sort(marray, marray + p);

    // start walking through image using forwards then backwards approach to make the addition
    // of extra values to the sorted array more efficient
    old_m1 = 0;
    old_m2 = 0;
    i = 0;
    j = 0;
    dir = DOWN;
    kmin = imin;
    kmax = imax;
    while ((i < h) && (j < w)) {

      // Work out the list of old and new values we need to correct the sorted array with
      switch (dir) {

        case DOWN:
          for (k=kmin; k <= kmax; k++) {
            if ((j + k - f) < 0) jj=0;
            else if ((j + k - f) >= w) jj=w - 1;
            else jj=j + k - f;
            if ((istop[k] + i) >= h) newv[k] = ia[(h - 1)*w + jj];
            else if ((istop[k] + i) <= 0) newv[k] = ia[0 + jj];
            else newv[k] = ia[(i + istop[k])*w + jj];
            if ((istart[k] + i) <= 0) oldv[k] = ia[0 + jj];
            else if ((istart[k] + i) >= h) oldv[k] = ia[(h - 1)*w + jj];
            else oldv[k] = ia[(i + istart[k] - 1)*w + jj];
          }
          break;

        case UP:
          for (k=kmin; k <= kmax; k++) {
            if ((j + k - f) < 0) jj=0;
            else if ((j + k - f) >= w) jj=w - 1;
            else jj=j + k - f;
            if ((istart[k] + i) <= -1) newv[k] = ia[0 + jj];
            else if ((istart[k] + i) >= (h - 1)) newv[k] = ia[(h - 1)*w + jj];
            else newv[k] = ia[(i + istart[k])*w + jj];
            if ((istop[k] + i) >= (h - 1)) oldv[k] = ia[(h - 1)*w + jj];
            else if ((istop[k] + i) <= -1) oldv[k] = ia[0 + jj];
            else oldv[k] = ia[(i + istop[k] + 1)*w + jj];
          }
          break;

        case RIGHT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            if ((jstop[k] + j) >= w) newv[k] = ia[ii*w + w - 1];
            else if ((jstop[k] + j) <= 0) newv[k] = ia[ii*w + 0];
            else newv[k] = ia[ii*w + j + jstop[k]];
            if ((jstart[k] + j) <= 0) oldv[k] = ia[ii*w + 0];
            else if ((jstart[k] + j) >= w) oldv[k] = ia[ii*w + w - 1];
            else oldv[k] = ia[ii*w + j + jstart[k] - 1];
          }
          break;

        case LEFT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            if ((jstart[k] + j) <= 0) newv[k] = ia[ii*w + 0];
            else if ((jstart[k] + j) >= (w - 1)) newv[k] = ia[ii*w + w - 1];
            else newv[k] = ia[ii*w + j + jstart[k]];
            if ((jstop[k] + j) >= (w - 1)) oldv[k] = ia[ii*w + w - 1];
            else if ((jstop[k] + j) <= -1) oldv[k] = ia[ii*w + 0];
            else oldv[k] = ia[ii*w + j + jstop[k] + 1];
          }
          break;

      }

      // Now adjust for these new and old values

      // Go through sorted array, removing old values and adding new values
      for (k=kmin; k <= kmax; k++) {
        left = 0;
        right = p - 1;
        m = (left + right) / 2;
        while (oldv[k] != marray[m]) {
          if (oldv[k] > marray[m]) left = m + 1;
          else if (oldv[k] < marray[m]) right = m - 1;
          m = (left + right) / 2;
        }
        if (newv[k] > oldv[k]) {
          while ((m < (p - 1)) && (newv[k] > marray[m + 1])) {
            marray[m] = marray[m + 1];
            m++;
          }
        } else {
          while ((m > 0) && (newv[k] < marray[m - 1])) {
            marray[m] = marray[m - 1];
            m--;
          }
        }
        marray[m] = newv[k];
      }

      // Find centile from sorted array
      if (oa1 != NULL) {

        if (forward) {
          ig = ig1;
        } else {
          ig = igr1;
        }

        if (!threshold_mask) {

          // Forward or inverse in simple case
          oa1[i*w + j] = marray[ig];

        } else if (forward) {

          // Forward operation when using threshold_mask option
          mtest = iarray[i*w + j];
          if (ig <= (p / 2)) {
            m = 0;
            while (fabs(marray[m] - mtest) > mthresh) m++;
            m += ig;
            if ((m >= p) || (fabs(marray[m] - mtest) > mthresh)) {
              mlim = mtest;
            } else {
              mlim = marray[m];
            }
          } else {
            m = (p - 1);
            while (fabs(marray[m] - mtest) > mthresh) m--;
            m -= ((p - 1) - ig);
            if ((m < 0) || (fabs(marray[m] - mtest) > mthresh)) {
              mlim = mtest;
            } else {
              mlim = marray[m];
            }
          }
          if (m != ig) {
            parray1[i*w + j] = mlim;
            oa1[i*w + j] = -1e10; // blank value for next pass
          } else {
            oa1[i*w + j] = mlim;
          }

        } else {

          // Inverse operation when using threshold_mask option
          mtest = ia[i*w + j];
          if (mtest == -1e10) mtest = parray1[i*w + j];
          if (ig < (p / 2)) {
            for (m=0; m < p; m++) if (marray[m] != -1e10) break;
            m += ig;
            if (m >= p) {
              mlim = marray[p - 1];
              if (mlim < mtest) mlim = mtest;
            } else {
              mlim = marray[m];
            }
            if ((mlim == -1e10) || (mlim > mtest)) {
              oa1[i*w + j] = mtest;
            } else {
              oa1[i*w + j] = mlim;
            }
          } else {
            m = ig;
            for (; m < (p - 1); m++) if (marray[m] != -1e10) break;
            mlim = marray[m];
            if (m > ig) {
              if (mlim > mtest) mlim = mtest;
            }
            if ((mlim == -1e10) || (mlim < mtest)) {
              oa1[i*w + j] = mtest;
            } else {
              oa1[i*w + j] = mlim;
            }
          }

        }
      }

      if (oa2 != NULL) {

        if (forward) {
          ig = ig2;
        } else {
          ig = igr2;
        }

        if (!threshold_mask) {

          // Forward or inverse in simple case
          oa2[i*w + j] = marray[ig];

        } else if (forward) {

          // Forward operation when using threshold_mask option
          mtest = iarray[i*w + j];
          if (ig <= (p / 2)) {
            m = 0;
            while (fabs(marray[m] - mtest) > mthresh) m++;
            m += ig;
            if ((m >= p) || (fabs(marray[m] - mtest) > mthresh)) {
              mlim = mtest;
            } else {
              mlim = marray[m];
            }
          } else {
            m = (p - 1);
            while (fabs(marray[m] - mtest) > mthresh) m--;
            m -= (p - 1) - ig;
            if ((m < 0) || (fabs(marray[m] - mtest) > mthresh)) {
              mlim = mtest;
            } else {
              mlim = marray[m];
            }
          }
          if (m != ig) {
            parray2[i*w + j] = mlim;
            oa2[i*w + j] = -1e10; // blank value for next pass
          } else {
            oa2[i*w + j] = mlim;
          }

        } else {

          // Inverse operation when using threshold_mask option
          mtest = ia[i*w + j];
          if (mtest == -1e10) mtest = parray2[i*w + j];
          if (ig < (p / 2)) {
            for (m=0; m < p; m++) if (marray[m] != -1e10) break;
            m += ig;
            if (m >= p) {
              mlim = marray[p - 1];
              if (mlim < mtest) mlim = mtest;
            } else {
              mlim = marray[m];
            }
            if ((mlim == -1e10) || (mlim > mtest)) {
              oa2[i*w + j] = mtest;
            } else {
              oa2[i*w + j] = mlim;
            }
          } else {
            m = ig;
            for (; m < (p - 1); m++) if (marray[m] != -1e10) break;
            mlim = marray[m];
            if (m > ig) {
              if (mlim > mtest) mlim = mtest;
            }
            if ((mlim == -1e10) || (mlim < mtest)) {
              oa2[i*w + j] = mtest;
            } else {
              oa2[i*w + j] = mlim;
            }
          }

        }
      }

      // increment index for next time
      if (right_first) {
        switch (dir) {
          case DOWN:
            if (j == (w - 1)) {
              dir = LEFT;
              kmin = jmin;
              kmax = jmax;
              j--;
            } else {
              dir = RIGHT;
              kmin = jmin;
              kmax = jmax;
              j++;
            }
            break;
          case RIGHT:
            if (j == (w - 1)) {
              dir = DOWN;
              kmin = imin;
              kmax = imax;
              i++;
            } else {
              j++;
            }
            break;
          case LEFT:
            if (j == 0) {
              dir = DOWN;
              kmin = imin;
              kmax = imax;
              i++;
            } else {
              j--;
            }
            break;
        }
      } else {
        switch (dir) {
          case RIGHT:
            if (i == (h - 1)) {
              dir = UP;
              kmin = imin;
              kmax = imax;
              i--;
            } else {
              dir = DOWN;
              kmin = imin;
              kmax = imax;
              i++;
            }
            break;
          case DOWN:
            if (i == (h - 1)) {
              dir = RIGHT;
              kmin = jmin;
              kmax = jmax;
              j++;
            } else {
              i++;
            }
            break;
          case UP:
            if (i == 0) {
              dir = RIGHT;
              kmin = jmin;
              kmax = jmax;
              j++;
            } else {
              i--;
            }
            break;
        }
      }

    }
  }

  // free variables
  delete[] oldv;
  delete[] newv;
  delete[] marray;
  delete[] istart;
  delete[] istop;
  delete[] jstart;
  delete[] jstop;
  if (tarray1 != NULL) delete[] tarray1;
  if (tarray2 != NULL) delete[] tarray2;
  if (parray1 != NULL) delete[] parray1;
  if (parray2 != NULL) delete[] parray2;
}


// Call this function using [X X2] = rankopen2(A, M, R, T, DUAL)
#define A_IN prhs[0]
#define M_IN prhs[1]
#define R_IN prhs[2]
#define T_IN prhs[3]
#define DUAL_IN prhs[4]
#define X_OUT plhs[0]
#define X2_OUT plhs[1]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int rank, rank2;
  bool dual;
  int w, h, dw, dh, l, i;
  double *Ad, *Xd, *X2d, thresh;
  int *Au, *Xu, *X2u, *Tu, *D;
  unsigned char *Ac, *Xc, *X2c;
  mxArray *Darray;

  // check we have consistent input and output arguments
  if (nrhs < 3) {
    mexErrMsgTxt("A, M and R must all be defined.");
  } else if (nrhs > 5) {
    mexErrMsgTxt("No more than 5 input arguments.");
  } else {
    if (nlhs < 1) {
      mexErrMsgTxt("Must be at least one output argument.");
    } else if (nlhs > 2) {
      mexErrMsgTxt("No more than two output arguments.");
    }
  }

  // check input data is of the right type
  if (mxIsComplex(A_IN) || mxIsSparse(A_IN) || ((mxGetClassID(A_IN) != mxUINT8_CLASS) && (mxGetClassID(A_IN) != mxDOUBLE_CLASS) && (mxGetClassID(A_IN) != mxINT32_CLASS)) || (mxGetNumberOfDimensions(A_IN) != 2)) {
    mexErrMsgTxt("Input data must be a matrix of double, uint8 or int32 real values.");
  }
  if (mxIsComplex(M_IN) || mxIsSparse(M_IN) || (mxGetNumberOfDimensions(M_IN) != 2)) {
    mexErrMsgTxt("M must be an array real values.");
  }
  if (mxIsComplex(R_IN) || (mxGetNumberOfElements(R_IN) != 1)) {
    mexErrMsgTxt("R must be a scalar.");
  } else {
    rank = (int)mxGetScalar(R_IN);
  }
  if (nrhs > 3) {
    if (mxIsComplex(T_IN) || (mxGetNumberOfElements(T_IN) != 1)) {
      mexErrMsgTxt("T must be a scalar.");
    } else {
      thresh = (int)mxGetScalar(T_IN);
    }
  } else {
    thresh = 0;
  }
  if (nrhs > 4) {
    if ((!mxIsLogical(DUAL_IN) && !mxIsDouble(DUAL_IN)) || (mxGetNumberOfElements(DUAL_IN) != 1)) {
      mexErrMsgTxt("DUAL must be a scalar.");
    } else {
      dual = (mxGetScalar(DUAL_IN) != 0);
    }
  } else {
    dual = true;
  }

  // Get input data and domain
  dw = (int)mxGetM(M_IN);
  dh = (int)mxGetN(M_IN);
  w = (int)mxGetM(A_IN);
  h = (int)mxGetN(A_IN);

  // Check for obvious inconsistencies in filter size
  l = dw;
  if ((dw != dh) || (l % 2 != 1) || (l > h) || (l > w)) {
    mexErrMsgTxt("Filter length must be odd and less than input data dimensions.");
  }
  if (mxGetClassID(M_IN) == mxINT32_CLASS) {
    D = (int *)mxGetData(M_IN);
  } else {
    Darray = mxCreateNumericArray(mxGetNumberOfDimensions(M_IN), mxGetDimensions(M_IN), mxINT32_CLASS, mxREAL);
    D = (int *)mxGetData(Darray);
    switch (mxGetClassID(M_IN)) {
      case mxDOUBLE_CLASS:
        for (i=0; i < (dw*dh); i++) if (((double *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
        break;
      case mxLOGICAL_CLASS:
        for (i=0; i < (dw*dh); i++) if (((bool *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
        break;
      case mxUINT8_CLASS:
        for (i=0; i < (dw*dh); i++) if (((unsigned char *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
        break;
      default:
        mexErrMsgTxt("Only double, logical, uint8 and int32 classes supported for M.");
        break;
    }
  }

  // Check for obvious inconsistencies in rank
  int p=0;
  for (i=0; i < (dw*dh); i++) p+= (int)D[i];
  rank2 = p + 1 - rank;
  if ((rank < 1) || (rank > p)) {
    mexErrMsgTxt("Lowest rank is 1, largest is number of non-zero elements in M.");
  }

  // Get input and output data and call appropriate function
  if (mxGetClassID(A_IN) == mxDOUBLE_CLASS) {
    Ad = (double *)mxGetData(A_IN);
    X_OUT = mxCreateDoubleMatrix(w, h, mxREAL);
    Xd = (double *)mxGetData(X_OUT);
    if (nlhs > 1) {
      X2_OUT = mxCreateDoubleMatrix(w, h, mxREAL);
      X2d = (double *)mxGetData(X2_OUT);
      rankopen_2d(Ad, w, h, l, D, thresh, dual, rank, Xd, rank2, X2d);
    } else {
      rankopen_2d(Ad, w, h, l, D, thresh, dual, rank, Xd, -1, NULL);
    }
  } else if (mxGetClassID(A_IN) == mxUINT8_CLASS) {
    Ac = (unsigned char *)mxGetData(A_IN);
    Au = new int[w*h];
    for (i=0; i < w*h; i++) Au[i] = (int)Ac[i];
    X_OUT = mxCreateNumericMatrix(w, h, mxUINT8_CLASS, mxREAL);
    Xc = (unsigned char *)mxGetData(X_OUT);
    Xu = new int[w*h];
    if (nlhs > 1) {
      X2_OUT = mxCreateNumericMatrix(w, h, mxUINT8_CLASS, mxREAL);
      X2c = (unsigned char *)mxGetData(X2_OUT);
      X2u = new int[w*h];
      rankopen_2d(Au, w, h, l, D, (int)thresh, dual, rank, Xu, rank2, X2u);
      for (i=0; i < w*h; i++) X2c[i] = (unsigned char)X2u[i];
      delete[] X2u;
    } else {
      rankopen_2d(Au, w, h, l, D, (int)thresh, dual, rank, Xu, -1, NULL);
    }
    for (i=0; i < w*h; i++) Xc[i] = (unsigned char)Xu[i];
    delete[] Xu;
    delete[] Au;
  } else {
    Au = (int *)mxGetData(A_IN);
    X_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
    Xu = (int *)mxGetData(X_OUT);
    if (nlhs > 1) {
      X2_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
      X2u = (int *)mxGetData(X2_OUT);
      rankopen_2d(Au, w, h, l, D, (int)thresh, dual, rank, Xu, rank2, X2u);
    } else {
      rankopen_2d(Au, w, h, l, D, (int)thresh, dual, rank, Xu, -1, NULL);
    }
  }

  // Might need to destroy temporary domain array
  if (mxGetClassID(M_IN) != mxINT32_CLASS) mxDestroyArray(Darray);

}
