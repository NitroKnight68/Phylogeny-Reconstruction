#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "suffixArray.h"

#define SAIS_MYMALLOC(_num, _type) ((_type *)malloc((_num) * sizeof(_type)))
#define SAIS_MYFREE(_ptr, _num, _type) free((_ptr))
#define chr(_a) (cs == sizeof(int) ? ((int *)T)[(_a)] : ((unsigned char *)T)[(_a)])

// Qsort int comparison function
int int_cmp(const void *a, const void *b) {
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;
  return *ia - *ib;
}

// Find the Start/End of each bucket
static void getCounts(const void *T, int *C, int n, int k, int cs) {
  int i;
  for(i = 0; i < k; ++i) { C[i] = 0; }
  for(i = 0; i < n; ++i) { ++C[chr(i)]; }
}

static void getBuckets(const int *C, int *B, int k, int end) {
  int i, sum = 0;
  if(end) { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum; } }
  else { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum - C[i]; } }
}

// Sort all type LMS suffixes
static void LMSsort1(const void *T, int *SA, int *C, int *B, int n, int k, int cs) {
  int bb, i, j;
  int c0, c1;

  // Computes SAl
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); // Finds start of buckets
  j = n - 1;
  bb = B[c1 = chr(j)];
  --j;
  SA[bb++] = (chr(j) < c1) ? ~j : j;
  for(i = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      assert(chr(j) >= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = bb; bb = B[c1 = c0]; }
      assert(i < bb);
      --j;
      SA[bb] = (chr(j) < c1) ? ~j : j;
      ++bb;
      SA[i] = 0;
    } else if(j < 0) {
      SA[i] = ~j;
    }
  }

  // Computes SAs
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, bb = B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      assert(chr(j) <= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = bb; bb = B[c1 = c0]; }
      assert((bb) <= i);
      --j;
      SA[--bb] = (chr(j) > c1) ? ~(j + 1) : j;
      SA[i] = 0;
    }
  }
}

static int LMSpostproc1(const void *T, int *SA, int n, int m, int cs) {
  int i, j, p, q, plen, qlen, name;
  int c0, c1;
  int diff;

  // Compact all the sorted substrings into the first m items of SA
  assert(0 < n);
  for(i = 0; (p = SA[i]) < 0; ++i) { SA[i] = ~p; assert((i + 1) < n); }
  if(i < m) {
    for(j = i, ++i;; ++i) {
      assert(i < n);
      if((p = SA[i]) < 0) {
        SA[j++] = ~p; SA[i] = 0;
        if(j == m) { break; }
      }
    }
  }

  // Store the length of all substrings
  i = n - 1; j = n - 1; c0 = chr(n - 1);
  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
  for(; 0 <= i;) {
    do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
    if(0 <= i) {
      SA[m + ((i + 1) >> 1)] = j - i; j = i + 1;
      do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
    }
  }

  // Find the lexicographic names of all substrings
  for(i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
    p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
    if((plen == qlen) && ((q + plen) < n)) {
      for(j = 0; (j < plen) && (chr(p + j) == chr(q + j)); ++j) { }
      if(j == plen) { diff = 0; }
    }
    if(diff != 0) { ++name, q = p, qlen = plen; }
    SA[m + (p >> 1)] = name;
  }

  return name;
}

static void LMSsort2(const void *T, int *SA, int *C, int *B, int *D, int n, int k, int cs) {
  int *b, i, j, t, d;
  int c0, c1;
  assert(C != B);

  // Computes SAl
  getBuckets(C, B, k, 0); // Finds start of buckets
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  --j;
  t = (chr(j) < c1);
  j += n;
  *b++ = (t & 1) ? ~j : j;
  for(i = 0, d = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      if(n <= j) { d += 1; j -= n; }
      assert(chr(j) >= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      --j;
      t = c0; t = (t << 1) | (chr(j) < c1);
      if(D[t] != d) { j += n; D[t] = d; }
      *b++ = (t & 1) ? ~j : j;
      SA[i] = 0;
    } else if(j < 0) {
      SA[i] = ~j;
    }
  }
  for(i = n - 1; 0 <= i; --i) {
    if(0 < SA[i]) {
      if(SA[i] < n) {
        SA[i] += n;
        for(j = i - 1; SA[j] < n; --j) { }
        SA[j] -= n;
        i = j;
      }
    }
  }

  // Computes SAs
  getBuckets(C, B, k, 1); // Finds endsof buckets
  for(i = n - 1, d += 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      if(n <= j) { d += 1; j -= n; }
      assert(chr(j) <= chr(j + 1));
      if((c0 = chr(j)) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      --j;
      t = c0; t = (t << 1) | (chr(j) > c1);
      if(D[t] != d) { j += n; D[t] = d; }
      *--b = (t & 1) ? ~(j + 1) : j;
      SA[i] = 0;
    }
  }
}

static int LMSpostproc2(int *SA, int n, int m) {
  int i, j, d, name;

  // Compact all the sorted LMS substrings into the first m items of SA 
  assert(0 < n);
  for(i = 0, name = 0; (j = SA[i]) < 0; ++i) {
    j = ~j;
    if(n <= j) { name += 1; }
    SA[i] = j;
    assert((i + 1) < n);
  }
  if(i < m) {
    for(d = i, ++i;; ++i) {
      assert(i < n);
      if((j = SA[i]) < 0) {
        j = ~j;
        if(n <= j) { name += 1; }
        SA[d++] = j; SA[i] = 0;
        if(d == m) { break; }
      }
    }
  }
  if(name < m) {
    // Store the lexicographic names
    for(i = m - 1, d = name + 1; 0 <= i; --i) {
      if(n <= (j = SA[i])) { j -= n; --d; }
      SA[m + (j >> 1)] = d;
    }
  } else {
    for(i = 0; i < m; ++i) {
      if(n <= (j = SA[i])) { j -= n; SA[i] = j; }
    }
  }

  return name;
}

// Finding the Suffix Array
static int sais_main(const void *T, int *SA, int *LCP, int fs, int n, int k, int cs, int isbwt, int level0) {
  int *C, *B, *D, *RA, *PLCP, *PHI, *DELTA, *b;
  int i, j, m, p, q, t, name, pidx = 0, newfs; // m -> number of S*-suffixes
  int c0, c1;
  unsigned int flags;

  assert((T != NULL) && (SA != NULL));
  assert((0 <= fs) && (0 < n) && (1 <= k));

  if(k <= 256) {
    if((C = SAIS_MYMALLOC(k, int)) == NULL) { return -2; }
    if(k <= fs) {
      B = SA + (n + fs - k);
      flags = 1;
    } else {
      if((B = SAIS_MYMALLOC(k, int)) == NULL) { SAIS_MYFREE(C, k, int); return -2; }
      flags = 3;
    }
  } else if(k <= fs) {
    C = SA + (n + fs - k);
    if(k <= (fs - k)) {
      B = C - k;
      flags = 0;
    } else if(k <= (256 * 4)) {
      if((B = SAIS_MYMALLOC(k, int)) == NULL) { return -2; }
      flags = 2;
    } else {
      B = C;
      flags = 8;
    }
  } else {
    if((C = B = SAIS_MYMALLOC(k, int)) == NULL) { return -2; }
    flags = 4 | 8;
  }
  if((n <= 0x3fffffff) && (2 <= (n / k))) {
    if(flags & 1) { flags |= ((k * 2) <= (fs - k)) ? 32 : 16; }
    else if((flags == 0) && ((k * 2) <= (fs - k * 2))) { flags |= 32; }
  }

  // Sorts all the LMS Substrings
  getCounts(T, C, n, k, cs); getBuckets(C, B, k, 1); // Finds ends of buckets

  for(i = 0; i < n; ++i) { SA[i] = 0; }
  b = &t; i = n - 1; j = n; m = 0; c0 = chr(n - 1);
  do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
  for(; 0 <= i;) {
    do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
    if(0 <= i) {
      *b = j;
      b = SA + --B[c1]; j = i; ++m;
      do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
    }
  }

  if(1 < m) {
    if(flags & (16 | 32)) {
      if(flags & 16) {
        if((D = SAIS_MYMALLOC(k * 2, int)) == NULL) {
          if(flags & (1 | 4)) { SAIS_MYFREE(C, k, int); }
          if(flags & 2) { SAIS_MYFREE(B, k, int); }
          return -2;
        }
      } else {
        D = B - k * 2;
      }
      assert((j + 1) < n);
      ++B[chr(j + 1)];
      for(i = 0, j = 0; i < k; ++i) {
        j += C[i];
        if(B[i] != j) { assert(SA[B[i]] != 0); SA[B[i]] += n; }
        D[i] = D[i + k] = 0;
      }
      LMSsort2(T, SA, C, B, D, n, k, cs);
      name = LMSpostproc2(SA, n, m);
      if(flags & 16) { SAIS_MYFREE(D, k * 2, int); }
    } else {
      LMSsort1(T, SA, C, B, n, k, cs);
      name = LMSpostproc1(T, SA, n, m, cs);
    }
  } else if (m == 1) { // only one S*-suffix => set immediately
    *b = j + 1; // set entry in SA
    if (level0) { LCP[b-SA] = -1; } // mark first (=only) S*-suffix in bucket
    name = 1;
  } else {
    name = 0;
  }

  // Recurse if "name" is not Unique
  if(name < m) {
    if(flags & 4) { SAIS_MYFREE(C, k, int); }
    if(flags & 2) { SAIS_MYFREE(B, k, int); }
    newfs = (n + fs) - (m * 2);
    if((flags & (1 | 4 | 8)) == 0) {
      if((k + name) <= newfs) { newfs -= k; }
      else { flags |= 8; }
    }
    assert((n >> 1) <= (newfs + m));
    RA = SA + m + newfs;
    for(i = m + (n >> 1) - 1, j = m - 1; m <= i; --i) {
      if(SA[i] != 0) {
        RA[j--] = SA[i] - 1;
      }
    }
    
    if(sais_main(RA, SA, NULL, newfs, m, name, sizeof(int), 0, 0) != 0) {
      if(flags & 1) { SAIS_MYFREE(C, k, int); }
      return -2;
    }

    // Compute starting index of S*-Suffixes
    i = n - 1; j = m - 1; c0 = chr(n - 1);
    do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
    for(; 0 <= i;) {
      do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) <= c1));
      if(0 <= i) {
        RA[j--] = i + 1;
        do { c1 = c0; } while((0 <= --i) && ((c0 = chr(i)) >= c1));
      }
    }

    // Constructing LCP for Suffixes
    if (level0) {
      if (m < n/3) {
        PHI = LCP+m;
        RA[m] = n;
        j = SA[0];
        PHI[j<<1] = n-1;
        PHI[(j<<1)+1] = 0; 
        for (i = 1; i < m; ++i) {
          q = SA[i]; 
          p = q<<1; 
          PHI[p]=RA[j]; 
          PHI[p+1]=RA[j+1]-RA[j]; 
          j = q;
        }

        PLCP = PHI; 
        p = 0;
        j = 0;
        for (i = 0; i < n; ++i) {
          if (i == RA[j]) {
            if (p < 0) p = 0;
            int twoj = j << 1;
            while (chr(i+p) == chr(PHI[twoj]+p)) ++p;
            t = PHI[twoj+1];
            q = RA[j+1]-RA[j];
            PLCP[twoj] = p;
            ++j;
            p -= (t > q) ? t : q;
          }
        }

        // Translate PLCP-values to SA-order
        for (j = 0; j < m; ++j) LCP[j] = PLCP[SA[j]<<1];
      }
      else {
        PHI = LCP; 
        DELTA = LCP+m; 
        RA[m] = n; 
        j = SA[0]; 
        PHI[j] = n-1; 
        DELTA[j] = 0;
        for (i = 1; i < m; ++i) {
          q = SA[i]; 
          PHI[q]=RA[j]; 
          DELTA[q]=RA[j+1]-RA[j]; 
          j = q; 
        }

        PLCP = DELTA;
        p = 0;
        j = 0;
        for (i = 0; i < n; ++i) {
          if (i == RA[j]) {
            if (p < 0) p = 0;
            while (chr(i+p) == chr(PHI[j]+p)) ++p;
            t = PLCP[j];
            q = RA[j+1]-RA[j];
            PLCP[j++] = p;
            p -= (t > q) ? t : q;
          }
        }

        // Translate PLCP-values to SA-order
        for (j = 0; j < m; ++j) LCP[j] = PLCP[SA[j]];
      }
    }

    // Translate index in RA to index in T
    for(i = 0; i < m; ++i) SA[i] = RA[SA[i]];

    if(flags & 4) {
      if((C = B = SAIS_MYMALLOC(k, int)) == NULL) { return -2; }
    }
    if(flags & 2) {
      if((B = SAIS_MYMALLOC(k, int)) == NULL) {
        if(flags & 1) { SAIS_MYFREE(C, k, int); }
        return -2;
      }
    }
  }
  else if (level0) { // Only for Small Inputs    
    j = SA[0];
    for (i = 1; i < m; ++i) {
      p = 0;
      while (chr(SA[i]+p) == chr(j+p)) p++;
      LCP[i] = p;
      j = SA[i];
    }
  }

  // Induce the Result
  if(flags & 8) { getCounts(T, C, n, k, cs); }
  // Put all S*-suffixes into their buckets
  if(1 < m) {
    getBuckets(C, B, k, 1);
    i = m - 1, j = n, p = SA[m - 1], c1 = chr(p);
    if (level0) {
      newfs = LCP[m-1];
      do {
        q = B[c0 = c1];
        while(q < j) {
          SA[--j] = 0; LCP[j] = -2; // Set remaining entries in old bucket to 0/-2
        }

        do {
          SA[--j] = p; LCP[j] = newfs;
          if(--i < 0) break;
          newfs = LCP[i]; p = SA[i];
        } while((c1 = chr(p)) == c0);
        LCP[j] = -1; // Mark first S*-suffix in every bucket
      } while(0 <= i);
      while(0 < j) {
        SA[--j] = 0; LCP[j] = -2; // Set remaining entries in smallest buckets to 0/-2
      }
    }
    else {
        do {
          q = B[c0 = c1];
          while(q < j) SA[--j] = 0; // Set remaining entries in old bucket to 0
          do {
          SA[--j] = p;
          if(--i < 0) break;
          p = SA[i];
          } while((c1 = chr(p)) == c0);
        } while(0 <= i);
        while(0 < j) SA[--j] = 0; // Set remaining entries in 1st bucket to 0
    }
  } 

  if(flags & (1 | 4)) { SAIS_MYFREE(C, k, int); }
  if(flags & 2) { SAIS_MYFREE(B, k, int); }

  return pidx;
}

// Calling the Main Function
int sais(const unsigned char *T, int *SA, int* LCP, int n) {
  if((T == NULL) || (SA == NULL) || (LCP == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; LCP[0] = 0; } return 0; }
  return sais_main(T, SA, LCP, 0, n, 256, sizeof(unsigned char), 0,1);
}