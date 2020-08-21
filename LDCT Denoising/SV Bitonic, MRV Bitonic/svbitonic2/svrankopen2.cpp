#include "mex.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

typedef struct {
  int i; // index into local array
  int k; // index into image
} hsvmf_t;

typedef struct {
  int v; // value
  int i; // index into local array
  int k; // index into image
} svmf_t;

typedef struct {
  double v; // value
  int i; // index into local array
  int k; // index into image
} dsvmf_t;

int qsort_dsvmf_compare(const void *a, const void *b)
{
  if (((dsvmf_t *)a)->v > ((dsvmf_t *)b)->v) return 1;
  else if (((dsvmf_t *)a)->v < ((dsvmf_t *)b)->v) return -1;
  else if (((dsvmf_t *)a)->i >((dsvmf_t *)b)->i) return 1;
  else if (((dsvmf_t *)a)->i < ((dsvmf_t *)b)->i) return -1;
  else return 0;
}

typedef enum { LEFT=0, RIGHT=1, UP=2, DOWN=3 } dir_t;

// Work out extents from a single mask shape
int get_mask_extents(int *M, int l, int m, int *istart, int *istop, int *jstart, int *jstop)
{
  int i, j, f, p, off;

  p = 0;
  off = m * l * l;
  f = (l - 1) / 2;
  for (i=0; i < l; i++) {
    j = 0;
    istart[i] = j - f;
    while ((j < l) && M[j*l + i + off] <= 0.0) {
      j++;
      istart[i] = j - f;
    }
    istop[i] = istart[i] - 1;
    while ((j < l) && M[j*l + i + off] > 0.0) {
      istop[i] = j - f;
      j++;
    }
    j = 0;
    jstart[i] = j - f;
    while ((j < l) && M[i*l + j + off] <= 0.0) {
      j++;
      jstart[i] = j - f;
    }
    jstop[i] = jstart[i] - 1;
    while ((j < l) && M[i*l + j + off] > 0.0) {
      jstop[i] = j - f;
      j++;
    }
    p += jstop[i] - jstart[i] + 1;
  }

  return p;
}

void svrankopen_2d(int *iarray, int w, int h, int l, int masks, int *M, int mthresh, bool dual, int *ig1, int *oarray1, int *marray1=NULL, int *moarray1=NULL, int *ig2=NULL, int *oarray2=NULL, int *marray2=NULL, int *moarray2=NULL);

// 2d filter to do rank filtering
// this version copes with multiple shapes at once, and also possibly different shapes at each location
void svrankopen_2d(int *iarray, int w, int h, int l, int masks, int *M, int mthresh, bool dual, int *ig1, int *oarray1, int *marray1, int *moarray1, int *ig2, int *oarray2, int *marray2, int *moarray2)
{
  int i, j, k, ii, jj, f, l2, p, m, e, *istart, *istop, *jstart, *jstop, maxval, minval, mii, m1, m2, v1, v2;
  int imin, imax, jmin, jmax, kmin, kmax, *mtotal, ttotal, ttotal2, *mp, moff, *mistart, *mistop, loop, *mmask;
  int *ia, *oa1, *oa2, *ma1=NULL, *ma2=NULL, *tarray1=NULL, *tarray2=NULL, *parray1 = NULL, *parray2 = NULL;
  int mlim, minc, mstart, mine, maxe, mtest;
  std::vector<int *> *mask;
  std::vector<int *>::iterator mi;
  std::vector<hsvmf_t> *harray;
  std::vector<hsvmf_t>::iterator hi;
  std::vector<hsvmf_t>::reverse_iterator rhi;
  svmf_t mval, *oldv, *newv;
  hsvmf_t hval;
  dir_t dir;
  bool right_first, open1, open2, inverse, forward;
  bool store_index = false;
  bool variable_ellipses = (masks > 1);
  bool threshold_mask = (mthresh > 0);

  // return if silly values
  f = (l - 1) / 2;
  if ((f > (h / 2)) || (f > (w / 2))) return;
  l2 = l * l;

  // Allocate various arrays
  if ((marray1 == NULL) && (oarray1 != NULL)) ma1 = new int[w*h];
  else ma1 = marray1;
  if ((marray2 == NULL) && (oarray2 != NULL)) ma2 = new int[w*h];
  else ma2 = marray2;
  if ((oarray1 != NULL) && threshold_mask) {
    parray1 = new int[w*h];
    for (i=0; i < w*h; i++) parray1[i] = 0;
  }
  if ((oarray2 != NULL) && threshold_mask) {
    parray2 = new int[w*h];
    for (i=0; i < w*h; i++) parray2[i] = 0;
  }
  istart = new int[l];
  istop = new int[l];
  jstart = new int[l];
  jstop = new int[l];
  mtotal = new int[masks];
  mp = new int[masks + 1];
  mistart = new int[l*masks]; // add space for a null mask at the end
  mistop = new int[l*masks];
  mask = new std::vector<int *>[l2];
  mmask = new int[l2*(masks + 1)]; // add space for a null mask at the end

  // Get extents for masks, and also a total mask list
  for (m=0; m < masks; m++) {
    p = get_mask_extents(M, l, m, istart, istop, jstart, jstop);
    for (i=0; i < l2; i++) mmask[m*l2 + i] = 0;
    for (j=0; j < l; j++) {
      for (i=istart[j]; i <= istop[j]; i++) {
        mask[(i + f)*l + j].push_back(mtotal + m);
        if ((i*i + (j - f)*(j - f)) <= (f - 1)*(f - 1) / 4) {
          mmask[m*l2 + (i + f)*l + j] = 2;
        } else {
          mmask[m*l2 + (i + f)*l + j] = 1;
        }
      }
      mistart[m*l + j] = istart[j];
      mistop[m*l + j] = istop[j];
    }
    mp[m] = p;
  }

  // set up null mask, only used for inverse operation
  for (i=0; i < l2; i++) mmask[masks*l2 + i] = 0;
  mmask[masks*l2 + (l2 / 2)] = 1;
  if (oarray1 != NULL) ig1[masks] = 0;
  if (oarray2 != NULL) ig2[masks] = 0;
  mp[masks] = 1;

  // Create the overall mask which represents the union of these shapes
  // this must also be convex, so fill in holes if necessary
  int *masksum = new int[l2];
  for (i=0; i < l2; i++) {
    if (mask[i].size() > 0) masksum[i] = 1;
    else masksum[i] = 0;
  }
  p = 0;
  for (j=0; j < l; j++) {
    i = 0;
    while ((i < l) && masksum[i*l + j] == 0) i++;
    istart[j] = i;
    i = l - 1;
    while ((i >= 0) && masksum[i*l + j] == 0) i--;
    istop[j] = i;
    // fill vertical holes
    for (i=istart[j]; i <= istop[j]; i++) masksum[i*l + j] = 1;
  }
  for (i=0; i < l; i++) {
    j = 0;
    while ((j < l) && masksum[i*l + j] == 0) j++;
    jstart[i] = j;
    j = l - 1;
    while ((j >= 0) && masksum[i*l + j] == 0) j--;
    jstop[i] = j;
    // fill horizontal holes
    for (j=jstart[i]; j <= jstop[i]; j++) masksum[i*l + j] = 1;
  }
  // horizontal hole filling might have changed vertical edges
  for (j=0; j < l; j++) {
    i = 0;
    while ((i < l) && masksum[i*l + j] == 0) i++;
    istart[j] = i;
    i = l - 1;
    while ((i >= 0) && masksum[i*l + j] == 0) i--;
    istop[j] = i;
    if (istop[j] >= istart[j]) p += (istop[j] - istart[j] + 1);
  }
  delete[] masksum;
  for (j=0; j < l; j++) {
    istart[j] -= f;
    istop[j] -= f;
    jstart[j] -= f;
    jstop[j] -= f;
  }
  imin = 0; while ((istart[imin] > istop[imin]) && (imin < (l - 1))) imin++;
  imax = l - 1; while ((istart[imax] > istop[imax]) && (imax > 0)) imax--;
  jmin = 0; while ((jstart[jmin] > jstop[jmin]) && (jmin < (l - 1))) jmin++;
  jmax = l - 1; while ((jstart[jmax] > jstop[jmax]) && (jmax > 0)) jmax--;
  if ((jmax - jmin) < (imax - imin)) right_first = true;
  else right_first = false;
  if (oarray1 != NULL) open1 = (ig1[0] < (mp[0] / 2));
  if (oarray2 != NULL) open2 = (ig2[0] < (mp[0] / 2));

  // create space for old and new values
  oldv = new svmf_t[l];
  newv = new svmf_t[l];

  // Might need to run the following code several times:
  for (loop=0; loop < 3; loop++) {
    if (loop == 0) {
      ia = iarray;
      forward = true;
      store_index = false;
      if (!dual) {
        // special case when want to do forward pass only
        inverse = false;
        oa1 = oarray1;
        oa2 = oarray2;
      } else {
        // will need a second pass to do each inverse
        inverse = false;
        // create temporary arrays
        if (oarray1 != NULL) tarray1 = new int[w*h];
        if (oarray2 != NULL) tarray2 = new int[w*h];
        oa1 = tarray1;
        oa2 = tarray2;
      }
    } else if (loop == 1) {
      if (oa1 == oarray1) break; // stop if this has already been done
      forward = false;
      inverse = true;
      store_index = true;
      ia = tarray1;
      oa1 = oarray1;
      oa2 = NULL;
    } else if (loop == 2) {
      if (oarray2 == NULL) break; // no further array to process
      ia = tarray2;
      oa1 = NULL;
      oa2 = oarray2;
    }

    // check for array range
    minval = ia[0];
    maxval = minval;
    for (i=1; i < w*h; i++) {
      if (ia[i] > maxval) maxval = ia[i];
      if (ia[i] < minval) minval = ia[i];
    }

    // create space for histogram, and fill with centre at (0, -1)
    e = maxval - minval + 1;
    if (loop > 0) delete[] harray;
    harray = new std::vector<hsvmf_t>[e];
    for (i=jmin; i <= jmax; i++) {
      for (j=jstart[i]; j <= jstop[i]; j++) {
        hval.i = (i - 1 - f)*l + j;
        if ((i < (1 + f)) && (j < 0)) k = 0;
        else if (i < (1 + f)) k = j;
        else if (j < 0) k = (i - 1 - f)*w;
        else k = (i - 1 - f)*w + j;
        ii = ia[k] - minval;
        if (store_index) {
          if (k == (i - 1 - f)*w + j) hval.k = k;
          else hval.k = -1;
        }
        harray[ii].push_back(hval);
      }
    }
    mine = 0;
    maxe = e - 1;

    // start walking through image using forwards then backwards approach to make the addition
    // of extra values to the sorted array more efficient
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
            mval.i = (i + istop[k])*l + j + k - f;
            if ((istop[k] + i) >= h) ii = (h - 1);
            else if ((istop[k] + i) <= 0) ii = 0;
            else ii = (i + istop[k]);
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((jj == (j + k - f)) && (ii == (i + istop[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + istart[k] - 1)*l + j + k - f;
            if ((istart[k] + i) <= 0) mval.v = ia[0 + jj];
            else if ((istart[k] + i) >= h) mval.v = ia[(h - 1)*w + jj];
            else mval.v = ia[(i + istart[k] - 1)*w + jj];
            oldv[k] = mval;
          }
          break;

        case UP:
          for (k=kmin; k <= kmax; k++) {
            if ((j + k - f) < 0) jj=0;
            else if ((j + k - f) >= w) jj=w - 1;
            else jj=j + k - f;
            mval.i = (i + istart[k])*l + j + k - f;
            if ((istart[k] + i) <= -1) ii = 0;
            else if ((istart[k] + i) >= (h - 1)) ii =(h - 1);
            else ii = (i + istart[k]);
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((jj == (j + k - f)) && (ii == (i + istart[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + istop[k] + 1)*l + j + k - f;
            if ((istop[k] + i) >= (h - 1)) mval.v = ia[(h - 1)*w + jj];
            else if ((istop[k] + i) <= -1) mval.v = ia[0 + jj];
            else mval.v = ia[(i + istop[k] + 1)*w + jj];
            oldv[k] = mval;
          }
          break;

        case RIGHT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            mval.i = (i + k - f)*l + j + jstop[k];
            if ((jstop[k] + j) >= w) jj = w - 1;
            else if ((jstop[k] + j) <= 0) jj = 0;
            else jj = j + jstop[k];
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((ii == (i + k - f)) && (jj == (j + jstop[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + k - f)*l + j + jstart[k] - 1;
            if ((jstart[k] + j) <= 0) mval.v = ia[ii*w + 0];
            else if ((jstart[k] + j) >= w) mval.v = ia[ii*w + w - 1];
            else mval.v = ia[ii*w + j + jstart[k] - 1];
            oldv[k] = mval;
          }
          break;

        case LEFT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            mval.i = (i + k - f)*l + j + jstart[k];
            if ((jstart[k] + j) <= 0) jj = 0;
            else if ((jstart[k] + j) >= (w - 1)) jj = w - 1;
            else jj = j + jstart[k];
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((ii == (i + k - f)) && (jj == (j + jstart[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + k - f)*l + j + jstop[k] + 1;
            if ((jstop[k] + j) >= (w - 1)) mval.v = ia[ii*w + w - 1];
            else if ((jstop[k] + j) <= -1) mval.v = ia[ii*w + 0];
            else mval.v = ia[ii*w + j + jstop[k] + 1];
            oldv[k] = mval;
          }
          break;
      }

      // Now adjust for these new and old values
      moff = ((i - f)*l + (j - f));

      // Remove old values and add new values to histogram
      for (k=kmin; k <= kmax; k++) {
        ii = oldv[k].v - minval;
        if (oldv[k].i == harray[ii].back().i) harray[ii].pop_back();
        else {
          for (hi=harray[ii].begin(); hi != harray[ii].end(); hi++) {
            if (hi->i == oldv[k].i) {
              *hi = harray[ii].back();
              harray[ii].pop_back();
              break;
            }
          }
        }
      }
      for (k=kmin; k <= kmax; k++) {
        ii = newv[k].v - minval;
        hval.i = newv[k].i;
        hval.k = newv[k].k;
        harray[ii].push_back(hval);
        if (ii > maxe) maxe = ii;
        if (ii < mine) mine = ii;
      }
      while (harray[mine].size() == 0) mine++;
      while (harray[maxe].size() == 0) maxe--;

      // Find required centile and mask from histogram
      if (forward) {

        // forward operation
        if (oa1 != NULL) {
          if (marray1 != NULL) m1 = ma1[i*w + j];
          else m1 = -1;
          mtest = iarray[i*w + j] - minval;
          ttotal2 = 0;
          if (open1) {
            if (m1 > -1) {
              ttotal = ig1[m1] + 1;
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig1[m] + 1;
              ttotal = masks;
            }
            minc = 1;
            mstart = mine;
            mlim = maxe + 1;
          } else {
            if (m1 > -1) {
              ttotal = mp[m1] - ig1[m1];
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig1[m];
              ttotal = masks;
            }
            minc = -1;
            mstart = maxe;
            mlim = mine - 1;
          }
          for (m=mstart; m != mlim; m+=minc) {
            for (hi=harray[m].begin(); hi != harray[m].end(); hi++) {
              mii = hi->i - moff;
              if (m1 > -1) {
                if (mmask[m1*l2 + mii] > 0) {
                  if (threshold_mask) {
                    ttotal2--;
                    if (abs(m - mtest) <= mthresh) {
                      ttotal--;
                    }
                  } else {
                    ttotal--;
                  }
                  if (ttotal == 0) break;
                }
              } else {
                for (mi=mask[mii].begin(); mi != mask[mii].end(); mi++) {
                  (**mi)--;
                  if (**mi == 0) {
                    ttotal--;
                    if (ttotal == 0) {
                      m1 = (int)(*mi - mtotal);
                      break;
                    }
                  }
                }
              }
              if (ttotal == 0) break;
            }
            if (ttotal == 0) break;
          }
          if (ttotal == 0) {
            v1 = m;
          } else {
            v1 = mtest;
          }
          if (threshold_mask && (ttotal2 != ttotal)) parray1[i*w + j] = 1;
          if (threshold_mask && (marray1 == NULL)) {
            int stop_search=false;
            if ((m + minval) > iarray[i*w + j]) {
              mlim = mine - 1;
              minc = -1;
            } else {
              mlim = maxe + 1;
              minc = 1;
            }
            for (; m != mlim; m+=minc) {
              for (hi=harray[m].begin(); hi != harray[m].end(); hi++) {
                mii = hi->i - moff;
                if ((abs((m + minval) - iarray[i*w + j]) <= mthresh)) {
                  if (m1 > -1) {
                    if (mmask[m1*l2 + mii] == 0) continue;
                  }
                  stop_search = true;
                  if (v1 != m) {
                    parray1[i*w + j] = 1;
                    v1 = m;
                  }
                  break;
                }
              }
              if (stop_search) break;
            }
          }
          v1 = v1 + minval;
        }
        if (oa2 != NULL) {
          if (marray2 != NULL) m2 = ma2[i*w + j];
          else m2 = -1;
          mtest = iarray[i*w + j] - minval;
          ttotal2 = 0;
          if (open2) {
            if (m2 > -1) {
              ttotal = ig2[m2] + 1;
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig2[m] + 1;
              ttotal = masks;
            }
            minc = 1;
            mstart = mine;
            mlim = maxe + 1;
          } else {
            if (m2 > -1) {
              ttotal = mp[m2] - ig2[m2];
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig2[m];
              ttotal = masks;
            }
            minc = -1;
            mstart = maxe;
            mlim = mine - 1;
          }
          for (m=mstart; m != mlim; m+=minc) {
            for (hi=harray[m].begin(); hi != harray[m].end(); hi++) {
              mii = hi->i - moff;
              if (m2 > -1) {
                if (mmask[m2*l2 + mii] > 0) {
                  if (threshold_mask) {
                    ttotal2--;
                    if (abs(m - mtest) <= mthresh) {
                      ttotal--;
                    }
                  } else {
                    ttotal--;
                  }
                  if (ttotal == 0) break;
                }
              } else {
                for (mi=mask[mii].begin(); mi != mask[mii].end(); mi++) {
                  (**mi)--;
                  if (**mi == 0) {
                    ttotal--;
                    if (ttotal == 0) {
                      m2 = (int)(*mi - mtotal);
                      break;
                    }
                  }
                }
                if (ttotal == 0) break;
              }
            }
            if (ttotal == 0) break;
          }
          if (ttotal == 0) {
            v2 = m;
          } else {
            v2 = mtest;
          }
          if (threshold_mask && (ttotal2 != ttotal)) parray2[i*w + j] = 1;
          if (threshold_mask && (marray2 == NULL)) {
            int stop_search=false;
            if ((m + minval) > iarray[i*w + j]) {
              mlim = mine - 1;
              minc = -1;
            } else {
              mlim = maxe + 1;
              minc = 1;
            }
            for (; m != mlim; m+=minc) {
              for (hi=harray[m].begin(); hi != harray[m].end(); hi++) {
                mii = hi->i - moff;
                if ((abs((m + minval) - iarray[i*w + j]) <= mthresh)) {
                  if (m2 > -1) {
                    if (mmask[m2*l2 + mii] == 0) continue;
                  }
                  stop_search = true;
                  if (v2 != m) {
                    parray2[i*w + j] = 1;
                    v2 = m;
                  }
                  break;
                }
              }
              if (stop_search) break;
            }
          }
          v2 = v2 + minval;
        }

      } else {

        // inverse operation
        if (oa1 != NULL) {
          v1 = -1;
          m1 = ma1[i*w + j];
          if (open1) {
            if (!variable_ellipses) {
              ttotal = ig1[m1] + 1;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig1[m] + 1;
              ttotal = masks;
            }
            minc = -1;
            mstart = maxe;
            mlim = mine - 1;
          } else {
            if (!variable_ellipses) {
              ttotal = mp[m1] - ig1[m1];
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig1[m];
              ttotal = masks;
            }
            minc = 1;
            mstart = mine;
            mlim = maxe + 1;
          }
          mlim = ia[i*w + j] - minval + minc;
          for (m=mstart; m != mlim; m+=minc) {
            for (hi=harray[m].begin(); hi != harray[m].end(); hi++) {
              if (hi->k == -1) continue;
              if (threshold_mask && (parray1[hi->k] > 0)) continue;
              m1 = ma1[hi->k];
              mii = hi->i - moff;
              mii = l2 - 1 - mii;
              if (mmask[m1*l2 + mii] > 0) {
                if (variable_ellipses) {
                  mtotal[m1]--;
                  if (mtotal[m1] == 0) {
                    v1 = m;
                    ttotal = 0;
                    break;
                  }
                } else {
                  ttotal--;
                  v1 = m;
                  if (ttotal == 0) break;
                }
              }
            }
            if (ttotal == 0) break;
          }
          if (m == mlim) {
            m -= minc;
            v1 = m;
          }
          v1 = v1 + minval;
        }
        if (oa2 != NULL) {
          v2 = -1;
          m2 = ma2[i*w + j];
          if (open2) {
            if (!variable_ellipses) {
              ttotal = ig2[m2] + 1;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig2[m] + 1;
              ttotal = masks;
            }
            minc = -1;
            mstart = maxe;
            mlim = mine - 1;
          } else {
            if (!variable_ellipses) {
              ttotal = mp[m2] - ig2[m2];
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig2[m];
              ttotal = masks;
            }
            minc = 1;
            mstart = mine;
            mlim = maxe + 1;
          }
          mlim = ia[i*w + j] - minval + minc;
          for (m=mstart; m != mlim; m+=minc) {
            for (hi=harray[m].begin(); hi != harray[m].end(); hi++) {
              if (hi->k == -1) continue;
              if (threshold_mask && (parray2[hi->k] > 0)) continue;
              mii = hi->i - moff;
              mii = l2 - 1 - mii;
              m2 = ma2[hi->k];
              if (mmask[m2*l2 + mii] > 0) {
                if (variable_ellipses) {
                  mtotal[m2]--;
                  if (mtotal[m2] == 0) {
                    v2 = m;
                    ttotal = 0;
                    break;
                  }
                } else {
                  ttotal--;
                  v2 = m;
                  if (ttotal == 0) break;
                }
              }
            }
            if (ttotal == 0) break;
          }
          if (m == mlim) {
            m -= minc;
            v2 = m;
          }
          v2 = v2 + minval;
        }
      }

      // Perform inverse operation with selected shape, or just write value to output array
      if (oa1 != NULL) {
        if (forward) ma1[i*w + j] = m1;
        oa1[i*w + j] = v1;
        if (moarray1 != NULL) moarray1[i*w + j] = m1;
      }
      if (oa2 != NULL) {
        if (forward) ma2[i*w + j] = m2;
        oa2[i*w + j] = v2;
        if (moarray2 != NULL) moarray2[i*w + j] = m2;
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
  delete[] harray;
  if (tarray1 != NULL) delete[] tarray1;
  if (tarray2 != NULL) delete[] tarray2;
  if ((marray1 == NULL) && (oarray1 != NULL)) delete[] ma1;
  if ((marray2 == NULL) && (oarray2 != NULL)) delete[] ma2;
  if (parray1 != NULL) delete[] parray1;
  if (parray2 != NULL) delete[] parray2;
  delete[] istart;
  delete[] istop;
  delete[] jstart;
  delete[] jstop;
  delete[] mask;
  delete[] mtotal;
  delete[] mp;
  delete[] mistart;
  delete[] mistop;
  delete[] mmask;

  return;
}

void svrankopen_2d(double *iarray, int w, int h, int l, int masks, int *M, double mthresh, bool dual, int *ig1, double *oarray1, int *marray1=NULL, int *moarray1=NULL, int *ig2=NULL, double *oarray2=NULL, int *marray2=NULL, int *moarray2=NULL);

// 2d filter to do rank filtering
// this version copes with multiple shapes at once, and also possibly different shapes at each location
void svrankopen_2d(double *iarray, int w, int h, int l, int masks, int *M, double mthresh, bool dual, int *ig1, double *oarray1, int *marray1, int *moarray1, int *ig2, double *oarray2, int *marray2, int *moarray2)
{
  int i, j, k, ii, jj, f, l2, p, m, *istart, *istop, *jstart, *jstop, mii, m1, m2, mlim, minc, mstart;
  int left, right, imin, imax, jmin, jmax, kmin, kmax, *mtotal, ttotal, ttotal2, *mp, moff, *mistart, *mistop, loop, *mmask;
  int *ma1=NULL, *ma2=NULL, *parray1 = NULL, *parray2 = NULL;
  double v1, v2, mtest, *ia, *oa1, *oa2, *tarray1=NULL, *tarray2=NULL;
  std::vector<int *> *mask;
  std::vector<int *>::iterator mi;
  dsvmf_t *marray, mval, *oldv, *newv;
  dir_t dir;
  bool right_first, open1, open2, inverse, forward;
  bool store_index = false;
  bool variable_ellipses = (masks > 1);
  bool threshold_mask = (mthresh > 0.0);

  // return if silly values
  f = (l - 1) / 2;
  if ((f > (h / 2)) || (f > (w / 2))) return;
  l2 = l * l;

  // Allocate various arrays
  if ((marray1 == NULL) && (oarray1 != NULL)) ma1 = new int[w*h];
  else ma1 = marray1;
  if ((marray2 == NULL) && (oarray2 != NULL)) ma2 = new int[w*h];
  else ma2 = marray2;
  if ((oarray1 != NULL) && threshold_mask) {
    parray1 = new int[w*h];
    for (i=0; i < w*h; i++) parray1[i] = 0;
  }
  if ((oarray2 != NULL) && threshold_mask) {
    parray2 = new int[w*h];
    for (i=0; i < w*h; i++) parray2[i] = 0;
  }
  istart = new int[l];
  istop = new int[l];
  jstart = new int[l];
  jstop = new int[l];
  mtotal = new int[masks];
  mp = new int[masks + 1];
  mistart = new int[l*masks]; // add space for a null mask at the end
  mistop = new int[l*masks];
  mask = new std::vector<int *>[l2];
  mmask = new int[l2*(masks + 1)]; // add space for a null mask at the end

  // Get extents for masks, and also a total mask list
  for (m=0; m < masks; m++) {
    p = get_mask_extents(M, l, m, istart, istop, jstart, jstop);
    for (i=0; i < l2; i++) mmask[m*l2 + i] = 0;
    for (j=0; j < l; j++) {
      for (i=istart[j]; i <= istop[j]; i++) {
        mask[(i + f)*l + j].push_back(mtotal + m);
        if ((i*i + (j - f)*(j - f)) <= (f - 1)*(f - 1) / 4) {
          mmask[m*l2 + (i + f)*l + j] = 2;
        } else {
          mmask[m*l2 + (i + f)*l + j] = 1;
        }
      }
      mistart[m*l + j] = istart[j];
      mistop[m*l + j] = istop[j];
    }
    mp[m] = p;
  }

  // set up null mask, only used for inverse operation
  for (i=0; i < l2; i++) mmask[masks*l2 + i] = 0;
  mmask[masks*l2 + (l2 / 2)] = 1;
  if (oarray1 != NULL) ig1[masks] = 0;
  if (oarray2 != NULL) ig2[masks] = 0;
  mp[masks] = 1;

  // Create the overall mask which represents the union of these shapes
  // this must also be convex, so fill in holes if necessary
  int *masksum = new int[l2];
  for (i=0; i < l2; i++) {
    if (mask[i].size() > 0) masksum[i] = 1;
    else masksum[i] = 0;
  }
  p = 0;
  for (j=0; j < l; j++) {
    i = 0;
    while ((i < l) && masksum[i*l + j] == 0) i++;
    istart[j] = i;
    i = l - 1;
    while ((i >= 0) && masksum[i*l + j] == 0) i--;
    istop[j] = i;
    // fill vertical holes
    for (i=istart[j]; i <= istop[j]; i++) masksum[i*l + j] = 1;
  }
  for (i=0; i < l; i++) {
    j = 0;
    while ((j < l) && masksum[i*l + j] == 0) j++;
    jstart[i] = j;
    j = l - 1;
    while ((j >= 0) && masksum[i*l + j] == 0) j--;
    jstop[i] = j;
    // fill horizontal holes
    for (j=jstart[i]; j <= jstop[i]; j++) masksum[i*l + j] = 1;
  }
  // horizontal hole filling might have changed vertical edges
  for (j=0; j < l; j++) {
    i = 0;
    while ((i < l) && masksum[i*l + j] == 0) i++;
    istart[j] = i;
    i = l - 1;
    while ((i >= 0) && masksum[i*l + j] == 0) i--;
    istop[j] = i;
    if (istop[j] >= istart[j]) p += (istop[j] - istart[j] + 1);
  }
  delete[] masksum;
  for (j=0; j < l; j++) {
    istart[j] -= f;
    istop[j] -= f;
    jstart[j] -= f;
    jstop[j] -= f;
  }
  imin = 0; while ((istart[imin] > istop[imin]) && (imin < (l - 1))) imin++;
  imax = l - 1; while ((istart[imax] > istop[imax]) && (imax > 0)) imax--;
  jmin = 0; while ((jstart[jmin] > jstop[jmin]) && (jmin < (l - 1))) jmin++;
  jmax = l - 1; while ((jstart[jmax] > jstop[jmax]) && (jmax > 0)) jmax--;
  if ((jmax - jmin) < (imax - imin)) right_first = true;
  else right_first = false;
  if (oarray1 != NULL) open1 = (ig1[0] < (mp[0] / 2));
  if (oarray2 != NULL) open2 = (ig2[0] < (mp[0] / 2));

  // create space for old and new values
  oldv = new dsvmf_t[l];
  newv = new dsvmf_t[l];

  // Might need to run the following code several times:
  for (loop=0; loop < 3; loop++) {
    if (loop == 0) {
      ia = iarray;
      forward = true;
      store_index = false;
      if (!dual) {
        // special case when want to do forward pass only
        inverse = false;
        oa1 = oarray1;
        oa2 = oarray2;
      } else {
        // will need a second pass to do each inverse
        inverse = false;
        // create temporary arrays
        if (oarray1 != NULL) tarray1 = new double[w*h];
        if (oarray2 != NULL) tarray2 = new double[w*h];
        oa1 = tarray1;
        oa2 = tarray2;
      }
    } else if (loop == 1) {
      if (oa1 == oarray1) break; // stop if this has already been done
      forward = false;
      inverse = true;
      store_index = true;
      ia = tarray1;
      oa1 = oarray1;
      oa2 = NULL;
    } else if (loop == 2) {
      if (oarray2 == NULL) break; // no further array to process
      ia = tarray2;
      oa1 = NULL;
      oa2 = oarray2;
    }

    // create space for the sorted list of length p, and fill it up with centre at (0,-1)
    if (loop > 0) delete[] marray;
    marray = new dsvmf_t[p];
    m = 0;
    for (i=jmin; i <= jmax; i++) {
      for (j=jstart[i]; j <= jstop[i]; j++) {
        mval.i = (i - 1 - f)*l + j;
        if ((i < (1 + f)) && (j < 0)) k = 0;
        else if (i < (1 + f)) k = j;
        else if (j < 0) k = (i - 1 - f)*w;
        else k = (i - 1 - f)*w + j;
        mval.v = ia[k];
        if (store_index) {
          if (k == (i - 1 - f)*w + j) mval.k = k;
          else mval.k = -1;
        }
        marray[m++] = mval;
      }
    }
    qsort(marray, p, sizeof(dsvmf_t), qsort_dsvmf_compare);

    // start walking through image using forwards then backwards approach to make the addition
    // of extra values to the sorted array more efficient
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
            mval.i = (i + istop[k])*l + j + k - f;
            if ((istop[k] + i) >= h) ii = (h - 1);
            else if ((istop[k] + i) <= 0) ii = 0;
            else ii = (i + istop[k]);
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((jj == (j + k - f)) && (ii == (i + istop[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + istart[k] - 1)*l + j + k - f;
            if ((istart[k] + i) <= 0) mval.v = ia[0 + jj];
            else if ((istart[k] + i) >= h) mval.v = ia[(h - 1)*w + jj];
            else mval.v = ia[(i + istart[k] - 1)*w + jj];
            oldv[k] = mval;
          }
          break;

        case UP:
          for (k=kmin; k <= kmax; k++) {
            if ((j + k - f) < 0) jj=0;
            else if ((j + k - f) >= w) jj=w - 1;
            else jj=j + k - f;
            mval.i = (i + istart[k])*l + j + k - f;
            if ((istart[k] + i) <= -1) ii = 0;
            else if ((istart[k] + i) >= (h - 1)) ii =(h - 1);
            else ii = (i + istart[k]);
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((jj == (j + k - f)) && (ii == (i + istart[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + istop[k] + 1)*l + j + k - f;
            if ((istop[k] + i) >= (h - 1)) mval.v = ia[(h - 1)*w + jj];
            else if ((istop[k] + i) <= -1) mval.v = ia[0 + jj];
            else mval.v = ia[(i + istop[k] + 1)*w + jj];
            oldv[k] = mval;
          }
          break;

        case RIGHT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            mval.i = (i + k - f)*l + j + jstop[k];
            if ((jstop[k] + j) >= w) jj = w - 1;
            else if ((jstop[k] + j) <= 0) jj = 0;
            else jj = j + jstop[k];
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((ii == (i + k - f)) && (jj == (j + jstop[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + k - f)*l + j + jstart[k] - 1;
            if ((jstart[k] + j) <= 0) mval.v = ia[ii*w + 0];
            else if ((jstart[k] + j) >= w) mval.v = ia[ii*w + w - 1];
            else mval.v = ia[ii*w + j + jstart[k] - 1];
            oldv[k] = mval;
          }
          break;

        case LEFT:
          for (k=kmin; k <= kmax; k++) {
            if ((i + k - f) < 0) ii=0;
            else if ((i + k - f) >= h) ii=h - 1;
            else ii=i + k - f;
            mval.i = (i + k - f)*l + j + jstart[k];
            if ((jstart[k] + j) <= 0) jj = 0;
            else if ((jstart[k] + j) >= (w - 1)) jj = w - 1;
            else jj = j + jstart[k];
            mval.v = ia[ii*w + jj];
            if (store_index) {
              if ((ii == (i + k - f)) && (jj == (j + jstart[k]))) mval.k = ii * w + jj;
              else mval.k = -1;
            }
            newv[k] = mval;
            mval.i = (i + k - f)*l + j + jstop[k] + 1;
            if ((jstop[k] + j) >= (w - 1)) mval.v = ia[ii*w + w - 1];
            else if ((jstop[k] + j) <= -1) mval.v = ia[ii*w + 0];
            else mval.v = ia[ii*w + j + jstop[k] + 1];
            oldv[k] = mval;
          }
          break;
      }

      // Now adjust for these new and old values
      moff = ((i - f)*l + (j - f));

      // Go through sorted array, removing old values and adding new values
      for (k=kmin; k <= kmax; k++) {
        left = 0;
        right = p - 1;
        m = (left + right) / 2;
        while (1) {
          if (oldv[k].v > marray[m].v) left = m + 1;
          else if (oldv[k].v < marray[m].v) right = m - 1;
          else if (oldv[k].i > marray[m].i) left = m + 1;
          else if (oldv[k].i < marray[m].i) right = m - 1;
          else break;
          m = (left + right) / 2;
        }
        if ((newv[k].v > oldv[k].v) || ((newv[k].v == oldv[k].v) && (newv[k].i > oldv[k].i))) {
          while ((m < (p - 1)) && (newv[k].v > marray[m + 1].v)) {
            marray[m] = marray[m + 1];
            m++;
          }
          while ((m < (p - 1)) && (newv[k].v == marray[m + 1].v) && (newv[k].i > marray[m + 1].i)) {
            marray[m] = marray[m + 1];
            m++;
          }
        } else {
          while ((m > 0) && (newv[k].v < marray[m - 1].v)) {
            marray[m] = marray[m - 1];
            m--;
          }
          while ((m > 0) && (newv[k].v == marray[m - 1].v) && (newv[k].i < marray[m - 1].i)) {
            marray[m] = marray[m - 1];
            m--;
          }
        }
        marray[m] = newv[k];
      }

      // Find required centile and mask from sorted array
      if (forward) {

        // forward operation
        if (oa1 != NULL) {
          if (marray1 != NULL) m1 = ma1[i*w + j];
          else m1 = -1;
          mtest = iarray[i*w + j];
          ttotal2 = 0;
          if (open1) {
            if (m1 > -1) {
              ttotal = ig1[m1] + 1;
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig1[m] + 1;
              ttotal = masks;
            }
            mstart = 0;
            mlim = p;
            minc = 1;
          } else {
            if (m1 > -1) {
              ttotal = mp[m1] - ig1[m1];
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig1[m];
              ttotal = masks;
            }
            mstart = p - 1;
            mlim = -1;
            minc = -1;
          }
          for (m=mstart; m != mlim; m+=minc) {
            mii = marray[m].i - moff;
            if (m1 > -1) {
              if (mmask[m1*l2 + mii] > 0) {
                if (threshold_mask) {
                  ttotal2--;
                  if (fabs(marray[m].v - mtest) <= mthresh) {
                    ttotal--;
                  }
                } else {
                  ttotal--;
                }
                if (ttotal == 0) break;
              }
            } else {
              for (mi=mask[mii].begin(); mi != mask[mii].end(); mi++) {
                (**mi)--;
                if (**mi == 0) {
                  ttotal--;
                  if (ttotal == 0) {
                    m1 = (int)(*mi - mtotal);
                    break;
                  }
                }
              }
            }
            if (ttotal == 0) break;
          }
          if (ttotal == 0) {
            v1 = marray[m].v;
          } else {
            v1 = mtest;
          }
          if (threshold_mask && (ttotal2 != ttotal)) parray1[i*w + j] = 1;
          if (threshold_mask && (marray1 == NULL)) {
            if (marray[m].v > iarray[i*w + j]) {
              mlim = -1;
              minc = -1;
            } else {
              mlim = p;
              minc = 1;
            }
            for (; m != mlim; m+=minc) {
              mii = marray[m].i - moff;
              if ((abs(marray[m].v - iarray[i*w + j]) <= mthresh) && (mmask[m1*l2 + mii] > 0)) {
                if (v1 != marray[m].v) {
                  parray1[i*w + j] = 1;
                  v1 = marray[m].v;
                }
                break;
              }
            }
          }
        }
        if (oa2 != NULL) {
          if (marray2 != NULL) m2 = ma2[i*w + j];
          else m2 = -1;
          mtest = iarray[i*w + j];
          ttotal2 = 0;
          if (open2) {
            if (m2 > -1) {
              ttotal = ig2[m2] + 1;
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig2[m] + 1;
              ttotal = masks;
            }
            mstart = 0;
            mlim = p;
            minc = 1;
          } else {
            if (m2 > -1) {
              ttotal = mp[m2] - ig2[m2];
              ttotal2 = ttotal;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig2[m];
              ttotal = masks;
            }
            mstart = p - 1;
            mlim = -1;
            minc = -1;
          }
          for (m=mstart; m != mlim; m+=minc) {
            mii = marray[m].i - moff;
            if (m2 > -1) {
              if (mmask[m2*l2 + mii] > 0) {
                if (threshold_mask) {
                  ttotal2--;
                  if (fabs(marray[m].v - mtest) <= mthresh) {
                    ttotal--;
                  }
                } else {
                  ttotal--;
                }
                if (ttotal == 0) break;
              }
            } else {
              for (mi=mask[mii].begin(); mi != mask[mii].end(); mi++) {
                (**mi)--;
                if (**mi == 0) {
                  ttotal--;
                  if (ttotal == 0) {
                    m2 = (int)(*mi - mtotal);
                    break;
                  }
                }
              }
            }
            if (ttotal == 0) break;
          }
          if (ttotal == 0) {
            v2 = marray[m].v;
          } else {
            v2 = mtest;
          }
          if (threshold_mask && (ttotal2 != ttotal)) parray2[i*w + j] = 1;
          if (threshold_mask && (marray2 == NULL)) {
            if (marray[m].v > iarray[i*w + j]) {
              mlim = -1;
              minc = -1;
            } else {
              mlim = p;
              minc = 1;
            }
            for (; m != mlim; m+=minc) {
              mii = marray[m].i - moff;
              if ((abs(marray[m].v - iarray[i*w + j]) <= mthresh) && (mmask[m2*l2 + mii] > 0)) {
                if (v2 != marray[m].v) {
                  parray2[i*w + j] = 1;
                  v2 = marray[m].v;
                }
                break;
              }
            }
          }
        }

      } else {

        // inverse operation
        if (oa1 != NULL) {
          m1 = ma1[i*w + j];
          v1 = ia[i*w + j];
          if (open1) {
            if (!variable_ellipses) {
              ttotal = ig1[m1] + 1;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig1[m] + 1;
              ttotal = masks;
            }
            mstart = p - 1;
            mlim = -1;
            minc = -1;
          } else {
            if (!variable_ellipses) {
              ttotal = mp[m1] - ig1[m1];
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig1[m];
              ttotal = masks;
            }
            mstart = 0;
            mlim = p;
            minc = 1;
          }
          for (m=mstart; m != mlim; m+=minc) {
            if (marray[m].k == -1) continue;
            mii = marray[m].i - moff;
            if (mii == (l2 - 1) / 2) {
              v1 = marray[m].v;
              break;
            }
            if (threshold_mask && (parray1[marray[m].k] > 0)) continue;
            m1 = ma1[marray[m].k];
            mii = l2 - 1 - mii;
            if (mmask[m1*l2 + mii] > 0) {
              if (variable_ellipses) {
                v1 = marray[m].v;
                mtotal[m1]--;
                if (mtotal[m1] == 0) {
                  ttotal = 0;
                  break;
                }
              } else {
                ttotal--;
                v1 = marray[m].v;
                if (ttotal == 0) break;
              }
            }
          }
        }
        if (oa2 != NULL) {
          m2 = ma2[i*w + j];
          v2 = ia[i*w + j];
          if (open2) {
            if (!variable_ellipses) {
              ttotal = ig2[m2] + 1;
            } else {
              for (m=0; m < masks; m++) mtotal[m] = ig2[m] + 1;
              ttotal = masks;
            }
            mstart = p - 1;
            mlim = -1;
            minc = -1;
          } else {
            if (!variable_ellipses) {
              ttotal = mp[m2] - ig2[m2];
            } else {
              for (m=0; m < masks; m++) mtotal[m] = mp[m] - ig2[m];
              ttotal = masks;
            }
            mstart = 0;
            mlim = p;
            minc = 1;
          }
          for (m=mstart; m != mlim; m+=minc) {
            if (marray[m].k == -1) continue;
            mii = marray[m].i - moff;
            if (mii == (l2 - 1) / 2) {
              v2 = marray[m].v;
              break;
            }
            if (threshold_mask && (parray2[marray[m].k] > 0)) continue;
            m2 = ma2[marray[m].k];
            mii = l2 - 1 - mii;
            if (mmask[m2*l2 + mii] > 0) {
              if (variable_ellipses) {
                v2 = marray[m].v;
                mtotal[m2]--;
                if (mtotal[m2] == 0) {
                  ttotal = 0;
                  break;
                }
              } else {
                ttotal--;
                v2 = marray[m].v;
                if (ttotal == 0) break;
              }
            }
          }
        }
      }

      // Perform inverse operation with selected shape, or just write value to output array
      if (oa1 != NULL) {
        if (forward) ma1[i*w + j] = m1;
        oa1[i*w + j] = v1;
        if (moarray1 != NULL) moarray1[i*w + j] = m1;
      }
      if (oa2 != NULL) {
        if (forward) ma2[i*w + j] = m2;
        oa2[i*w + j] = v2;
        if (moarray2 != NULL) moarray2[i*w + j] = m2;
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
  if (tarray1 != NULL) delete[] tarray1;
  if (tarray2 != NULL) delete[] tarray2;
  if ((marray1 == NULL) && (oarray1 != NULL)) delete[] ma1;
  if ((marray2 == NULL) && (oarray2 != NULL)) delete[] ma2;
  if (parray1 != NULL) delete[] parray1;
  if (parray2 != NULL) delete[] parray2;
  delete[] istart;
  delete[] istop;
  delete[] jstart;
  delete[] jstop;
  delete[] mask;
  delete[] mtotal;
  delete[] mp;
  delete[] mistart;
  delete[] mistop;
  delete[] mmask;

  return;
}


// Call this function using [X Mout X2 M2out] = svrankopen2(A, M, R, T, DUAL, MASKS, Min, M2in)
#define A_IN prhs[0]
#define M_IN prhs[1]
#define R_IN prhs[2]
#define T_IN prhs[3]
#define DUAL_IN prhs[4]
#define MASKS_IN prhs[5]
#define Min_IN prhs[6]
#define M2in_IN prhs[7]
#define X_OUT plhs[0]
#define M_X2_OUT plhs[1]
#define X2_OUT plhs[2]
#define M2_OUT plhs[3]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  bool dual, mask_out, mask_in, twice;
  int w, h, dw, dh, dm, l;
  long i, j, p;
  double *Ad, *Xd, *X2d, thresh;
  int *Au, *Xu, *X2u, *Tu, *D, *R, *R2, *Min, *M2in, *Mout, *M2out;
  unsigned char *Ac, *Xc, *X2c;
  mxArray *Darray, *Rarray, *R2array, *Marray, *M2array;

  // check we have consistent input and output arguments
  if (nrhs < 3) {
    mexErrMsgTxt("A, M and R must all be defined.");
  } else if (nrhs > 8) {
    mexErrMsgTxt("No more than 8 input arguments.");
  } else {
    if (nlhs < 1) {
      mexErrMsgTxt("Must be at least one output argument.");
    } else if ((nlhs != 1) && (nlhs != 2) && (nlhs != 4)) {
      mexErrMsgTxt("Must be either 1, 2 or 4 output arguments.");
    }
  }

  // check input data is of the right type
  if (mxIsComplex(A_IN) || mxIsSparse(A_IN) || ((mxGetClassID(A_IN) != mxDOUBLE_CLASS) && (mxGetClassID(A_IN) != mxUINT8_CLASS) && (mxGetClassID(A_IN) != mxINT32_CLASS)) || (mxGetNumberOfDimensions(A_IN) != 2)) {
    mexErrMsgTxt("Input data must be a matrix of double, uint8 or int32 real values.");
  }
  w = (int)mxGetM(A_IN);
  h = (int)mxGetN(A_IN);

  // Check domain size 
  if (mxIsComplex(M_IN) || mxIsSparse(M_IN) || (mxGetNumberOfDimensions(M_IN) < 2) || (mxGetNumberOfDimensions(M_IN) > 3)) {
    mexErrMsgTxt("M must be an array of logical or double real values.");
  }
  if (mxGetNumberOfDimensions(M_IN) == 2) {
    dm = 1;
    dw = (int)mxGetM(M_IN);
    dh = (int)mxGetN(M_IN);
  } else {
    dw = (int)(mxGetDimensions(M_IN)[0]);
    dh = (int)(mxGetDimensions(M_IN)[1]);
    dm = (int)(mxGetDimensions(M_IN)[2]);
  }
  l = dw;
  if ((dw != dh) || (l % 2 != 1) || (l > h) || (l > w)) {
    mexErrMsgTxt("M must be square, with odd side length and less than input data dimensions.");
  }

  // Check correct number of rank variables
  if (mxIsComplex(R_IN) || mxIsSparse(R_IN) || (mxGetNumberOfDimensions(R_IN) != 2) || (mxGetNumberOfElements(R_IN) != dm)) {
    mexErrMsgTxt("R must be a scalar vector of the same length as the number of masks in DOMAIN.");
  }

  // Check for threshold and read in immediately
  if (nrhs > 3) {
    if (mxIsComplex(T_IN) || (mxGetNumberOfElements(T_IN) != 1)) {
      mexErrMsgTxt("T must be a scalar.");
    } else {
      thresh = (int)mxGetScalar(T_IN);
    }
  } else {
    thresh = 0;
  }

  // Check for dual and read in immediately
  if (nrhs > 4) {
    if ((!mxIsLogical(DUAL_IN) && !mxIsDouble(DUAL_IN)) || (mxGetNumberOfElements(DUAL_IN) != 1)) {
      mexErrMsgTxt("DUAL must be a scalar.");
    } else {
      dual = (mxGetScalar(DUAL_IN) != 0);
    }
  } else {
    dual = true;
  }

  // Check for masks and read in immediately
  if (nrhs > 5) {
    if ((!mxIsLogical(MASKS_IN) && !mxIsDouble(MASKS_IN)) || (mxGetNumberOfElements(MASKS_IN) != 1)) {
      mexErrMsgTxt("MASKS must be a scalar.");
    } else {
      mask_out = (mxGetScalar(MASKS_IN) != 0);
    }
  } else {
    mask_out = false;
  }
  if (mask_out) {
    if (nlhs == 1) {
      mexErrMsgTxt("MASKS is true but there is only one output.");
    } else if (nlhs == 2) {
      twice = false;
    } else {
      twice = true;
    }
  } else {
    if (nlhs == 1) {
      twice = false;
    } else if (nlhs == 2) {
      twice = true;
    } else {
      mexErrMsgTxt("MASKS is false but there are more than two outputs.");
    }
  }

  // Check for input mask arrays
  if (nrhs > 6) {
    if (mxIsComplex(Min_IN) || mxIsSparse(Min_IN) || (mxGetNumberOfDimensions(Min_IN) != 2)) {
      mexErrMsgTxt("Min must be a matrix of real values.");
    }
    if ((mxGetM(Min_IN) != w) || (mxGetN(Min_IN) != h)) {
      mexErrMsgTxt("Min must be the same size as A.");
    }
    mask_in = true;
  } else {
    mask_in = false;
  }
  if (twice && mask_in) {
    if (nrhs > 7) {
      if (mxIsComplex(M2in_IN) || mxIsSparse(M2in_IN) || (mxGetNumberOfDimensions(M2in_IN) != 2)) {
        mexErrMsgTxt("M2in must be a matrix of real values.");
      }
      if ((mxGetM(M2in_IN) != w) || (mxGetN(M2in_IN) != h)) {
        mexErrMsgTxt("M2in must be the same size as A.");
      }
    } else {
      mexErrMsgTxt("Must also define M2in if two outputs.");
    }
  }

  // Convert input domain data to int
  Darray = mxCreateNumericArray(mxGetNumberOfDimensions(M_IN), mxGetDimensions(M_IN), mxINT32_CLASS, mxREAL);
  D = (int *)mxGetData(Darray);
  switch (mxGetClassID(M_IN)) {
    case mxDOUBLE_CLASS:
      for (i=0; i < (dw*dh*dm); i++) if (((double *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
      break;
    case mxLOGICAL_CLASS:
      for (i=0; i < (dw*dh*dm); i++) if (((bool *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
      break;
    case mxUINT8_CLASS:
      for (i=0; i < (dw*dh*dm); i++) if (((unsigned char *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
      break;
    case mxINT32_CLASS:
      for (i=0; i < (dw*dh*dm); i++) if (((unsigned char *)mxGetData(M_IN))[i]) D[i] = 1; else D[i] = 0;
      break;
    default:
      mexErrMsgTxt("Only double, logical, uint8 and int32 classes supported for M.");
      break;
  }

  // Convert input rank data to int and zero-offset
  // Note we need to turn this array around
  mwSize *Rsize = new mwSize[1];
  Rsize[0] = dm + 1;
  Rarray = mxCreateNumericArray(1, Rsize, mxINT32_CLASS, mxREAL);
  R = (int *)mxGetData(Rarray);
  switch (mxGetClassID(R_IN)) {
    case mxDOUBLE_CLASS:
      for (i=0; i < dm; i++) R[i] = (int)(((double *)mxGetData(R_IN))[i]) - 1;
      break;
    case mxUINT8_CLASS:
      for (i=0; i < dm; i++) R[i] = (int)(((unsigned char *)mxGetData(R_IN))[i]) - 1;
      break;
    case mxINT32_CLASS:
      for (i=0; i < dm; i++) R[i] = (int)(((int *)mxGetData(R_IN))[i]) - 1;
      break;
    default:
      mexErrMsgTxt("Only double, uint8 and int32 classes supported for R.");
      break;
  }

  // Convert input mask data to int
  if (mask_in) {
    Marray = mxCreateUninitNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
    Min = (int *)mxGetData(Marray);
    switch (mxGetClassID(Min_IN)) {
      case mxDOUBLE_CLASS:
        for (i=0; i < (w*h); i++) Min[i] = (int)(((double *)mxGetData(Min_IN))[i]) - 1;
        break;
      case mxUINT8_CLASS:
        for (i=0; i < (w*h); i++) Min[i] = (int)(((unsigned char *)mxGetData(Min_IN))[i]) - 1;
        break;
      case mxINT32_CLASS:
        for (i=0; i < (w*h); i++) Min[i] = (int)(((int *)mxGetData(Min_IN))[i]) - 1;
        break;
      default:
        mexErrMsgTxt("Only double, uint8 and int32 classes supported for Min.");
        break;
    }
  } else {
    Min = NULL;
  }
  if (twice && mask_in) {
    M2array = mxCreateUninitNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
    M2in = (int *)mxGetData(M2array);
    switch (mxGetClassID(M2in_IN)) {
      case mxDOUBLE_CLASS:
        for (i=0; i < (w*h); i++) M2in[i] = (int)(((double *)mxGetData(M2in_IN))[i]) - 1;
        break;
      case mxUINT8_CLASS:
        for (i=0; i < (w*h); i++) M2in[i] = (int)(((unsigned char *)mxGetData(M2in_IN))[i]) - 1;
        break;
      case mxINT32_CLASS:
        for (i=0; i < (w*h); i++) M2in[i] = (int)(((int *)mxGetData(M2in_IN))[i]) - 1;
        break;
      default:
        mexErrMsgTxt("Only double, uint8 and int32 classes supported for M2in.");
        break;
    }
  } else {
    M2in = NULL;
  }

  // Check for obvious inconsistencies in rank and mask values
  R2array = mxCreateNumericArray(1, Rsize, mxINT32_CLASS, mxREAL);
  R2 = (int *)mxGetData(R2array);
  for (j=0; j<dm; j++) {
    p=0;
    for (i=0; i < (dw*dh); i++) p += D[i + j*dw*dh];
    R2[j] = p - 1 - R[j];
    if ((R[j] < 0) || (R[j] >= p)) {
      mexErrMsgTxt("Lowest R is 1, largest is number of non-zero elements in M.");
    }
  }
  if (mask_in) {
    for (i=0; i < (w*h); i++) {
      if ((Min[i] < 0) || (Min[i] >= dm)) {
        mexErrMsgTxt("Masks in Min must be between 1 and the number of masks in M.");
      }
    }
  }
  if (mask_in && twice) {
    for (i=0; i < (w*h); i++) {
      if ((M2in[i] < 0) || (M2in[i] >= dm)) {
        mexErrMsgTxt("Masks in M2in must be between 1 and the number of masks in M.");
      }
    }
  }

  // Get input and output data and call appropriate function
  if (mask_out) {
    M_X2_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
    Mout = (int *)mxGetData(M_X2_OUT);
    if (twice) {
      M2_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
      M2out = (int *)mxGetData(M2_OUT);
    } else {
      M2out = NULL;
    }
  } else {
    Mout = NULL;
    M2out = NULL;
  }
  if (mxGetClassID(A_IN) == mxDOUBLE_CLASS) {
    Ad = (double *)mxGetData(A_IN);
    X_OUT = mxCreateDoubleMatrix(w, h, mxREAL);
    Xd = (double *)mxGetData(X_OUT);
    if (twice) {
      if (mask_out) {
        X2_OUT = mxCreateDoubleMatrix(w, h, mxREAL);
        X2d = (double *)mxGetData(X2_OUT);
      } else {
        M_X2_OUT = mxCreateDoubleMatrix(w, h, mxREAL);
        X2d = (double *)mxGetData(M_X2_OUT);
      }
      svrankopen_2d(Ad, w, h, l, dm, D, thresh, dual, R, Xd, Min, Mout, R2, X2d, M2in, M2out);
    } else {
      svrankopen_2d(Ad, w, h, l, dm, D, thresh, dual, R, Xd, Min, Mout);
    }
  } else if (mxGetClassID(A_IN) == mxUINT8_CLASS) {
    Ac = (unsigned char *)mxGetData(A_IN);
    Au = new int[w*h];
    for (i=0; i < w*h; i++) Au[i] = (int)Ac[i];
    X_OUT = mxCreateNumericMatrix(w, h, mxUINT8_CLASS, mxREAL);
    Xc = (unsigned char *)mxGetData(X_OUT);
    Xu = new int[w*h];
    if (twice) {
      if (mask_out) {
        X2_OUT = mxCreateNumericMatrix(w, h, mxUINT8_CLASS, mxREAL);
        X2c = (unsigned char *)mxGetData(X2_OUT);
      } else {
        M_X2_OUT = mxCreateNumericMatrix(w, h, mxUINT8_CLASS, mxREAL);
        X2c = (unsigned char *)mxGetData(M_X2_OUT);
      }
      X2u = new int[w*h];
      svrankopen_2d(Au, w, h, l, dm, D, (int)thresh, dual, R, Xu, Min, Mout, R2, X2u, M2in, M2out);
      for (i=0; i < w*h; i++) X2c[i] = (unsigned char)X2u[i];
      delete[] X2u;
    } else {
      svrankopen_2d(Au, w, h, l, dm, D, (int)thresh, dual, R, Xu, Min, Mout);
    }
    for (i=0; i < w*h; i++) Xc[i] = (unsigned char)Xu[i];
    delete[] Xu;
    delete[] Au;
  } else {
    Au = (int *)mxGetData(A_IN);
    X_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
    Xu = (int *)mxGetData(X_OUT);
    if (twice) {
      if (mask_out) {
        X2_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
        X2u = (int *)mxGetData(X2_OUT);
      } else {
        M_X2_OUT = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
        X2u = (int *)mxGetData(M_X2_OUT);
      }
      svrankopen_2d(Au, w, h, l, dm, D, (int)thresh, dual, R, Xu, Min, Mout, R2, X2u, M2in, M2out);
    } else {
      svrankopen_2d(Au, w, h, l, dm, D, (int)thresh, dual, R, Xu, Min, Mout);
    }
  }

  // Make output mask array not zero-offset
  if (mask_out) {
    for (i=0; i < (w*h); i++) Mout[i] += 1;
    if (twice) for (i=0; i < (w*h); i++) M2out[i] += 1;
  }

  // Might need to destroy temporary arrays
  mxDestroyArray(Darray);
  mxDestroyArray(Rarray);
  mxDestroyArray(R2array);
  if (mask_in) mxDestroyArray(Marray);
  if (twice && mask_in) mxDestroyArray(M2array);

}
