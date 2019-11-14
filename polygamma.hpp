// This code uses the ideas from "PROGRAMS FOR COMPUTING THE LOGARITHM OF THE
// GAMMA FUNCTION, AND THE DIGAMMA FUNCTION, FOR COMPLEX ARGUMENT", Kolbig 1972
// with the error bounds derived in "Calculation of the Gamma Function by Stirling's 
// Formula", Spira 1971

// It defines pGamma(m, z) = d^(m+1)/dy^(m+1) ln(Gamma(y)) |_y=z

// Calculation is done with the asymptotic series 
// ln(Gamma(z)) = (z-1/2)ln(z)-z+1/2ln(2pi)+sum_k=1^K B_2k/(2k*(2k-1))z^(1-2k)
// and it's derivatives.

// The error on the m-th derivative is bounded by |B_2K|/(2K-1) * |z|^1-2m-k/(2m+k-1)^k

// Given this, this implementation uses K=7, |z| > 14 in the application of the series,
// mapping other values onto this range using identities of the polygama
// pGamma(m, z) = pGamma(m, z+1) - (-1)^m m! z^(-m-1)
// pGamma(m, z) = (-1)^m pGamma(m, 1-z)

// The choice of parameters should ensure that the approximation-induced error
// is below 2^-53, not accounting for floating-point operation inaccuracies.

// This implementation is partly inspired by the julia implementation at 
// https://github.com/JuliaMath/SpecialFunctions.jl/blob/master/src/gamma.jl

#include <complex>
#include <vector>
#include <cmath>

// Fast multiplication with (-1)^m (e.g. without using multiplication)
template<class T>
std::complex<T> signflip(int m, std::complex<T> z) {
	if (m % 2 == 0) return z;
	return -z;
}

// Factorial, needed in derivatives
template<class T>
T fac(int m) {
	if (m < 0) throw "Negative m";
	if (m == 0) return T(1);
	T res = T(1);
	for (int i=1; i<=m; i++) {
		res *= T(i);
	}
	return res;
}

// M-th derivative of cot(z)
// since d/dz cot(z) = -1-cot^2(z), it follows that the m-th derivatives of
// cot(z) is of the form p_m(cot(z)). Then working out, we find
// p_0(z) = z
// p_{m+1}(z) = -(1+z^2)p'_m(z)
// This function calculates the coeficients of these polynomials p, and then
// uses that to calculate d^m/dz^m cot(z).
template<class T>
std::complex<T> mdCot(int m, std::complex<T> z) {
	static std::vector<std::vector<T>> p = {{T(0), T(1)}};
	
	while (m >= (int)p.size()) {
		int c = p.size()-1;
		p.push_back(std::vector<T>(c+3, T(0)));
		for (int j=0; j<=c+1; j++) {
			if (j != 0) p[c+1][j-1] += T(j)*p[c][j];
			p[c+1][j+1] += T(j)*p[c][j];
		}
	}
	
	std::complex<T> res = T(0);
	std::complex<T> x = T(1)/std::tan(z);
	std::complex<T> pow = T(1);
	for (int i=0; i<=m; i++) {
		res += p[m][i]*pow;
		pow *= x;
	}
	
	return signflip(m, res);
}

// diGamma, see above for notes
// implemented separately, since this requires 0-th derivative of ln(z), which is not of
// form z^-m
template<class T>
std::complex<T> diGamma(std::complex<T> z) {
	// Ensure z > 0
	if (z.real() < T(0)) {
		return diGamma(std::complex<T>(T(1),T(0))-z)-M_PI/std::tan(M_PI*z);
	}
	
	// Ensure |z| > 14
	std::complex<T> result = T(0);
	while (z.real() < T(14)) {
		result -= std::pow(z,-1);
		z += T(1);
	}
	
	// ln(z)
	result += std::log(z);
	
	// -1/2z
	result -= T(1)/(T(2)*z);
	
	// Series terms
	// Series terms
	T B2k[] = {
		T(1)/T(1),
		T(1)/T(6),
		T(-1)/T(30),
		T(1)/T(42),
		T(-1)/T(30),
		T(5)/T(66),
		T(-691)/T(2730),
		T(7)/T(6)
	};
	for (int i=1; i<=7; i++) {
		result -= B2k[i]*std::pow(z, -2*i)/T(2*i);
	}
	return result;
}

// polygamma, see above for notes
template<class T>
std::complex<T> pGamma(int m, std::complex<T> z) {
	if (m < 0) throw "Negative m";
	if (m == 0) return diGamma(z);

	// Ensure z > 0
	if (z.real() < 0) {
		return signflip(m, pGamma(m, std::complex<T>(T(1), T(0))-z)) - std::pow(std::complex<T>(M_PI), m+1) * mdCot(m, M_PI * z);
	}

	// Ensure |z| > 2*K+m+1
	std::complex<T> result=T(0);	
	while (z.real() < T(2*7+m+1)) {
		result -= signflip(m, fac<T>(m) * std::pow(z, -m-1));
		z += T(1);
	}
	
	// m-th derivative ln(z)
	result += signflip(m-1, fac<T>(m-1) * std::pow(z, -m));
	
	// m-th derivative -1/2z
	result -= signflip(m, fac<T>(m) * std::pow(z, -m-1)/T(2));
	
	// Series terms
	T B2k[] = {
		T(1)/T(1),
		T(1)/T(6),
		T(-1)/T(30),
		T(1)/T(42),
		T(-1)/T(30),
		T(5)/T(66),
		T(-691)/T(2730),
		T(7)/T(6)
	};
	for (int i=1; i<=7; i++) {
		// Calculate derivative factor
		T dfac = T(1);
		for (int j=1; j<m; j++) {
			dfac *= T(2*i+j);
		}
		result -= signflip(m, dfac*B2k[i]*std::pow(z, -2*i-m));
	}
	
	return result;
}
