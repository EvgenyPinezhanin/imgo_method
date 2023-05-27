#include <map.h>

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

void node(int is, int& n1, int nexp, int& l, int& iq, int *iu, int *iv) {
    // calculate iu = u[s], iv = v[s], l = l[s] by is = s

    int n, j, k1, k2, iff;

    n = n1 + 1;
    if (is == 0) {
        l = n1;
        for (int i = 0; i < n; i++) {
            iu[i] = -1; 
            iv[i] = -1;
        }
    } else if (is == (nexp - 1)) {
        l = n1;
        iu[0] = 1;
        iv[0] = 1;
        for (int i = 1; i < n; i++) {
            iu[i] = -1; 
            iv[i] = -1;
        }
        iv[n1] = 1;
    } else {
        iff = nexp;
        k1 = -1;
        for (int i = 0; i < n; i++) {
            iff /= 2;
            if (is >= iff) {
                if ((is == iff) && (is != 1)) {
                    l = i;
                    iq = -1;
                }
                is = is - iff;
                k2 = 1;
            } else {
                k2 = -1;
                if ((is == (iff - 1)) && (is != 0)) {
                    l = i;
                    iq = 1;
                }
            }
            j = -k1 * k2;
            iv[i] = j;
            iu[i] = j;
            k1 = k2;
        }
        iv[l] = iv[l] * iq;
        iv[n1] = -iv[n1];
    }
}

void mapd(double x, int m, double* y, int n, int key) {
    // mapping y(x) : 1 - center, 2 - line, 3 - node

    double d, mne, dd, dr;
    float p, r;
    int iw[11];
    int it, is = 0, i, j, k;
    int n1, nexp, l, iq, iu[10] = { 0 }, iv[10] = { 0 };

    p = 0.0;
    n1 = n - 1;
    for (nexp = 1, i = 0; i < n; nexp *= 2, i++); // nexp = 2**n
    d = x;
    r = 0.5;
    it = 0;
    dr = nexp;
    for (mne = 1, i = 0; i < m; mne *= dr, i++); // mne = dr**m
    for (i = 0; i < n; i++) {
        iw[i] = 1;
        y[i] = 0.0;
    }
    if (key == 2) {
        d = d * (1.0 - 1.0 / mne);
        k = 0;
    } else if (key > 2) {
        dd = floor(d * (mne - floor(mne / nexp)));
        d = floor(dd + (dd - 1.0) / (nexp - 1.0)) * (1.0 / mne);
    }
    for (j = 0; j < m; j++) {
        iq = 0;
        if (x == 1.0) {
            is = nexp - 1; 
            d = 0.0;
        } else {
              d = d * nexp;
              is = d;
              d = d - is;
        }
        i = is;
        node(i, n1, nexp, l, iq, iu, iv);
        i = iu[0];
        iu[0] = iu[it];
        iu[it] = i;
        i = iv[0];
        iv[0] = iv[it];
        iv[it] = i;
        if (l == 0) {
            l = it;
        } else if (l == it) {
            l = 0;
        }
        if ((iq > 0) || ((iq == 0) && (is == 0))) {
            k = l;
        } else if (iq < 0) {
            k = (it == n1) ? 0 : n1;
        }
        r = r * 0.5;
        it = l;
        for (i = 0; i < n; i++) {
            iu[i] = iu[i] * iw[i];
            iw[i] = -iv[i] * iw[i];
            p = r * iu[i];
            p = p + y[i];
            y[i] = p;
        }
    }
    if (key == 2) {
        if (is == (nexp - 1)) {
            i = -1;
        } else {
            i = 1;
        }
        p = 2 * i * iu[k] * r * d;
        p = y[k] - p;
        y[k] = p;
    } else if (key == 3) {
        for (i = 0; i < n; i++) {
            p = r * iu[i];
            p = p + y[i];
            y[i] = p;
        }
    }
}

void numbr(int *iss, int& n1, int nexp, int& l, int* iu, int* iv) {
    // calculate s(u) = is,l(u) = l,v(u) = iv by u = iu

    int n, is, iff, k1, k2, l1 = 0;

    n = n1 + 1;
    iff = nexp;
    is = 0;
    k1 = -1;
    for (int i = 0; i < n; i++) {
        iff = iff / 2;
        k2 = -k1 * iu[i];
        iv[i] = iu[i];
        k1 = k2;
        if (k2 < 0) {
        	l1 = i;
        } else {
        	is += iff;
        	l = i;
        }
    }
    if (is == 0) {
        l = n1;
    } else {
        iv[n1] = -iv[n1];
        if (is == (nexp - 1) ) {
        	l = n1;
        } else if (l1 == n1) {
        	iv[l] = -iv[l];
        } else {
        	l = l1;
        }
    }
    *iss = is;
}

void xyd(double *xx, int m, float y[], int n) {
    // calculate preimage x  for nearest level  m center to y
    // (x - left boundary point of level m interval)

    double x, r1;
    float r;
    int iw[10];
    int i, j, it, is;
    int n1, nexp, l, iu[10] = { 0 }, iv[10];

    n1 = n - 1;
    for (nexp = 1, i = 0; i < n; i++) {
        nexp *= 2; 
        iw[i] = 1;
    }
    r = 0.5;
    r1 = 1.0;
    x = 0.0;
    it = 0;
    for (j = 0; j < m; j++) {
        r *= 0.5;
        for (i = 0; i < n; i++) {
           iu[i] = (y[i] < 0) ? -1 : 1;
           y[i] -= r * iu[i];
           iu[i] *= iw[i];
        }
        i = iu[0];
        iu[0] = iu[it];
        iu[it] = i;
        numbr(&is, n1, nexp, l, iu, iv);
        i = iv[0];
        iv[0] = iv[it];
        iv[it] = i;
        for (i = 0; i < n; i++)
          	iw[i] = -iw[i] * iv[i];
        if (l == 0) {
        	l = it;
        } else if (l == it) {
        	l = 0;
        }
        it = l;
        r1 /= nexp;
        x += r1 * is;
    }
    *xx = x;
}

void invmad(int m, double xp[], int kp, int *kxx, double p[], int n, int incr) {
    // calculate kx preimage p node
    //   node type mapping m level

    // preimages calculation
    // - m - map level (number of partioning)
    // - xp - preimages to be calculated
    // - kp - number of preimages that may be calculated (size of xp)
    // - kxx - number of preimages being calculated
    // - p - image for which preimages are calculated
    // - n - dimension of image (size of p)
    // - incr - minimum number of map nodes that must be between preimages

    double mne, d1, dd, x, dr;
    float r, d, u[10], y[10];
    int i, k, kx, nexp;
    double del;

    kx = 0;
    kp--;
    for (nexp = 1, i = 0; i < n; i++) {
        nexp *= 2; 
        u[i] = -1.0;
    }
    dr = nexp;
    for (mne = 1, r = 0.5, i = 0; i < m; i++) {
        mne *= dr;
        r *= 0.5;
    }
    dr = floor(mne / nexp);
    del = 1.0 / (mne - dr);
    d1 = del * (incr + 0.5);
    for (kx = -1; kx < kp;) {
        for (i = 0; i < n; i++) { // label 2
            d = p[i];
            y[i] = d - r * u[i];
        }
        for (i = 0; (i < n) && (fabs(y[i]) < 0.5) ; i++);
        if (i >= n) {
            xyd(&x, m, y, n);
            x = (floor(x * mne) - floor(floor(x * mne) / nexp)) * del;
            if (kx > kp) break;
            k = kx++; // label 9
            if (kx == 0) {
                xp[0] = x;
            } else {
                while (k >= 0) {
                    dr = fabs(x - xp[k]); // label 11
                    if (dr <= d1) {
                        for (kx--; k < kx; k++, xp[k] = xp[k + 1]);
                        goto m6;
                    } else if (x <= xp[k]) {
                        xp[k + 1] = xp[k];
                        k--;
                    } else break;
                }
                xp[k + 1] = x;
            }
        }
        m6 : for (i = n - 1; (i >= 0) && (u[i] = (u[i] <= 0.0) ? 1 : -1) < 0; i--);
        if (i < 0) break;
    }
    *kxx = ++kx;
}
