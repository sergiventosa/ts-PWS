#include "cdotx.h"
#include "stdio.h"

void cdotx_dd (double *y, double *x0, unsigned int N, double *w0, unsigned int L, int c, int inc) {
	double *x, *w, aux;
	unsigned int l, l0, n, n0;

	if (N<L) {
		x = x0; x0 = w0; w0 = x;
		l = N; N = L; L = l;
	}
	for (n=0; n<N; n+=inc) {
		n0 = (N + n - c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		x = x0 + n0;
		w = w0;
		aux = 0;
		for (l=0; l<l0>>1; l++) {
			aux += w[0]*x[0] + w[1]*x[1];
			w+=2; x+=2;
		}
		if (l0&1) aux += *w++*x[0];
		x = x0;
		for (l=0; l<(L-l0)>>1; l++) {
			aux += w[0]*x[0] + w[1]*x[1];
			w+=2; x+=2;
		}
		if ((L-l0)&1) aux += w[0]*x[0];
		*y++ = aux;
	}
}

void cdotx_dc (double complex *y0, double *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc) {
	double complex *w, *y, aux;
	double *x;
	unsigned int l, l0, n, n0;

	if (N<L) {
		cdotx_cd(y0, w0, L, x0, N, c, inc);
		for (n=0; n<N/inc; n++) y0[n] = conj(y0[n]);    /* Integer div. truncates towards zero. */
	} else {
		y = y0;
		for (n=0; n<N; n+=inc) {
			n0 = (N + n - c) % N;
			l0 = N - n0;
			if (L < l0) l0 = L;

			x = x0 + n0;
			w = w0;
			aux = 0;
			for (l=0; l<l0>>2; l++) {
				aux += x[0]*w[0] + x[1]*w[1] + x[2]*w[2] + x[3]*w[3];
				x+=4; w+=4;
			}
			for (l=0; l<(l0&3); l++)
				aux += x[l]*w[l];

			w += l;
			x = x0;
			for (l=0; l<(L-l0)>>1; l++) {
				aux += x[0]*w[0] + x[1]*w[1];
				x+=2; w+=2;
			}
			if ((L-l0)&1)
				aux += x[0]*w[0];
			
			*y++ = conj(aux);
		}
	}
}

void cdotx_cd (double complex *y0, double complex *x0, unsigned int N, double *w0, unsigned int L, int c, int inc) {
	double complex *x, *y, aux;
	double *w;
	unsigned int l, l0, n;
	int n0;

	if (N<L) {
		cdotx_dc(y0, w0, L, x0, N, c, inc);
		for (n=0; n<N/inc; n++) y0[n] = conj(y0[n]);    /* Integer div. truncates towards zero. */
	} else {
		y = y0;
		for (n=0; n<N; n+=inc) {
			n0 = (N + n - c) % N;
			l0 = N - n0;
			if (L < l0) l0 = L;

			x = x0 + n0;
			w = w0;
			aux = 0;
			for (l=0; l<l0>>1; l++) {
				aux += w[0]*x[0] + w[1]*x[1];
				w+=2; x+=2;
			}
			if (l0&1) {
				aux += w[0]*x[0];
				w++;
			}
			x = x0;
			for (l=0; l<(L-l0)>>1; l++) {
				aux += w[0]*x[0] + w[1]*x[1];
				w+=2; x+=2;
			}
			if ((L-l0)&1)
				aux += w[0]*x[0];
			
			*y++ = aux;
		}
	}
}

void cdotx_cc (double complex *y, double complex *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc) {
	double complex *x, *w, aux;
	unsigned int l, l0, n;
	int n0;

	if (N<L) {
		x = x0; x0 = w0; w0 = x;
		l = N; N = L; L = l;
	}
	for (n=0; n<N; n+=inc) {
		n0 = (N + n -c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		x = x0 + n0;
		w = w0;
		aux = 0;
		for (l=0; l<l0; l++) aux += conj(*w++) * *x++;
		x = x0;
		for (    ; l<L; l++) aux += conj(*w++) * *x++;
		*y++ = aux;
	}
}

#if 0
void re_cdotx_cc (double *y, double complex *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc) {
	double complex *x, *w;
	double aux;
	unsigned int l, l0, n;
	int n0;

	if (N<L) {
		x = x0; x0 = w0; w0 = x;
		l = N; N = L; L = l;
	}
	for (n=0; n<N; n+=inc) {
		n0 = (N + n -c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		x = x0 + n0;
		w = w0;
		aux = 0;
		for (l=0; l<l0>>1; l++) {
			aux += creal(w[0])*creal(x[0]) + cimag(w[0])*cimag(x[0]) + creal(w[1])*creal(x[1]) + cimag(w[1])*cimag(x[1]);
			w+=2; x+=2;
		}
		if (l0&1) {
			aux += creal(w[0])*creal(x[0]) + cimag(w[0])*cimag(x[0]);
			w++;
		}
		x = x0;
		for (l=0; l<(L-l0)>>1; l++) {
			aux += creal(w[0])*creal(x[0]) + cimag(w[0])*cimag(x[0]) + creal(w[1])*creal(x[1]) + cimag(w[1])*cimag(x[1]);
			w+=2; x+=2;
		}
		if ((L-l0)&1)
			aux += creal(w[0])*creal(x[0]) + cimag(w[0])*cimag(x[0]);
		*y++ = aux;
	}
}
#else
void re_cdotx_cc (double *y, double complex *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc) {
	double *x, *w;
	double aux;
	unsigned int l, l0, n;
	int n0;

	if (N<L) {
		x = (double *)x0; x0 = (double complex *)w0; w0 = (double complex *)x;
		l = N; N = L; L = l;
	}
	for (n=0; n<N; n+=inc) {
		n0 = (N + n -c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		x = (double *)(x0 + n0);
		w = (double *)w0;
		aux = 0;
		for (l=0; l<l0>>1; l++) {
			aux += w[0]*x[0] + w[1]*x[1] + w[2]*x[2] + w[3]*x[3];
			w+=4; x+=4;
		}
		if (l0&1) {
			aux += w[0]*x[0] + w[1]*x[1];
			w+=2;
		}
		x = (double *)x0;
		for (l=0; l<(L-l0)>>1; l++) {
			aux += w[0]*x[0] + w[1]*x[1] + w[2]*x[2] + w[3]*x[3];
			w+=4; x+=4;
		}
		if ((L-l0)&1)
			aux += w[0]*x[0] + w[1]*x[1];
		*y++ = aux;
	}
}
#endif

void cdotx_upsampling_dd (double *y, unsigned int N, double *x0, unsigned int SX, double *w, unsigned int L, int c, int D) {
	/* SX is not used */
	double *x, ad;
	unsigned int l, l0, n, ai;
	int n0;
	
	ai = (N/D + D)*D; /* Just a number bigger than N with no reminder. */

	for (n=0; n<N; n++) {
		n0 = (N + n -c) % N;
		l0 = N - n0;

		ad = 0;
		x = x0 + (n0+D-1)/D;
		if (L < l0)
			for (l=(ai-n0)%D; l<L; l+=D) ad += w[l] * *x++;
		else {
			for (l=(ai-n0)%D; l<l0; l+=D) ad += w[l] * *x++;
			x = x0;
			for (l=l0; l<L; l+=D) ad += w[l] * *x++;
		}
		y[n] = D * ad;
	}
}

void cdotx_upsampling_dc (double complex *y, unsigned int N, double *x0, unsigned int SX, double complex *w, unsigned int L, int c, int D) {
	double complex ad;
	double *x;
	unsigned int l, l0, n, ai;
	int n0;
	
	ai = (N/D + D)*D; /* Just a number bigger than N with no reminder. */

	for (n=0; n<N; n++) {
		n0 = (N + n -c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		ad = 0;
		x = x0 + (n0+D-1)/D;
		for (l=(ai-n0)%D; l<l0; l+=D) ad += conj(w[l]) * *x++;

		x = x0;
		for (l=l0; l<L; l+=D) ad += conj(w[l]) * *x++;
		y[n] = D * ad;
	}
}

void cdotx_upsampling_cd (double complex *y, unsigned int N, double complex *x0, unsigned int SX, double *w, unsigned int L, int c, int D) {
	double complex ad, *x;
	unsigned int l, l0, n, ai;
	int n0;
	
	ai = (N/D + D)*D; /* Just a number bigger than N with no reminder. */

	for (n=0; n<N; n++) {
		n0 = (N + n -c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		ad = 0;
		x = x0 + (n0+D-1)/D;
		for (l=(ai-n0)%D; l<l0; l+=D) ad += w[l] * *x++;

		x = x0;
		for (l=l0; l<L; l+=D) ad += w[l] * *x++;
		y[n] = D * ad;
	}
}
void cdotx_upsampling_cc (double complex *y, unsigned int N, double complex *x0, unsigned int SX, double complex *w, unsigned int L, int c, int D) {
	double complex ad, *x;
	unsigned int l, l0, n, ai;
	int n0;
	
	ai = (N/D + D)*D; /* Just a number bigger than N with no reminder. */

	for (n=0; n<N; n++) {
		n0 = (N + n - c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		ad = 0;
		x = x0 + (n0+D-1)/D;
		for (l=(ai-n0)%D; l<l0; l+=D) ad += conj(w[l]) * *x++;

		x = x0;
		for (l=l0; l<L; l+=D) ad += conj(w[l]) * *x++;
		y[n] = D * ad;
	}
}

void re_cdotx_upsampling_cc (double *y, const unsigned int N, double complex *x0, unsigned int SX, const double complex *w, const unsigned int L, const int c, const int D) {
	double complex *x;
	double ad;
	unsigned int l, l0, n, ai;
	int n0;
	
	ai = (N/D + D)*D; /* Just a number bigger than N with no reminder. */
	
	for (n=0; n<N; n++) {
		n0 = (N + n - c) % N;
		l0 = N - n0;
		
		ad = 0;
		x = x0 + (n0+D-1)/D;
		if (L < l0) {
			for (l=(ai-n0)%D; l<L; l+=D) {
				/* ad += creal(conj(w[l]) * *x++); */
				ad += creal(w[l])*creal(x[0]) + cimag(w[l])*cimag(x[0]);
				x++;
			}
		} else {
			for (l=(ai-n0)%D; l<l0; l+=D) {
				/* ad += creal(conj(w[l]) * *x++); */
				ad += creal(w[l])*creal(x[0]) + cimag(w[l])*cimag(x[0]);
				x++;
			}
			x = x0;
			for (l=l0; l<L; l+=D) {
				/* ad += creal(conj(w[l]) * *x++); */
				ad += creal(w[l])*creal(x[0]) + cimag(w[l])*cimag(x[0]);
				x++;
			}
		}
		y[n] = D * ad;
	}
}
