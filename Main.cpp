#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

// Standard normal probability density function
double norm_pdf(const double x) {
    return (1.0 / (std::pow(2 * M_PI, 0.5))) * std::exp(-0.5 * x * x);
}

// An approximation to the cumulative distribution function for the standard normal distribution
double norm_cdf(const double x) {
    double k = 1.0 / (1.0 + 0.2316419 * std::abs(x));
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));
    
    double cdf = 1.0 - (1.0 / (std::pow(2 * M_PI, 0.5))) * std::exp(-0.5 * x * x) * k_sum;
    return x >= 0.0 ? cdf : 1.0 - cdf;
}

// This calculates d_j, for j in {1, 2}. This term appears in the closed form solution for the European call or put price
double d_j(const int j, const double S, const double K, const double r, const double v, const double T) {
    return (std::log(S / K) + (r + (std::pow(-1, j - 1)) * 0.5 * v * v) * T) / (v * (std::pow(T, 0.5)));
}

// Calculate the European vanilla call price based on underlying S, strike K, risk-free rate r, volatility of underlying sigma and time to maturity T
double bs_call_price(const double S, const double K, const double r, const double sigma, const double T) {
    return S * norm_cdf(d_j(1, S, K, r, sigma, T)) - K * std::exp(-r * T) * norm_cdf(d_j(2, S, K, r, sigma, T));
}

// Calculate the European vanilla put price based on underlying S, strike K, risk-free rate r, volatility of underlying sigma and time to maturity T
double bs_put_price(const double S, const double K, const double r, const double sigma, const double T) {
    return K * std::exp(-r * T) * norm_cdf(-d_j(2, S, K, r, sigma, T)) - S * norm_cdf(-d_j(1, S, K, r, sigma, T));
}

// Calculate the Merton jump-diffusion call price based on a finite sum approximation to the infinite series solution, making use of the BS call price.
double bs_jd_call_price(const double S, const double K, const double r, const double sigma, const double T, const int N, const double m, const double lambda, const double nu) {
    double price = 0.0;  // Stores the final call price
    double factorial = 1.0;

    // Pre-calculate as much as possible
    double lambda_p = lambda * m;
    double lambda_p_T = lambda_p * T;

    // Calculate the finite sum over N terms
    for (int n = 0; n < N; ++n) {
        double sigma_n = std::sqrt(sigma * sigma + n * nu * nu / T);
        double r_n = r - lambda * (m - 1) + n * std::log(m) / T;

        // Calculate n!
        if (n > 0) {
            factorial *= n;
        }

        // Refine the jump price over the loop
        price += ((std::exp(-lambda_p_T) * std::pow(lambda_p_T, n)) / factorial) * bs_call_price(S, K, r_n, sigma_n, T);
    }

    return price;
}

// Calculate the Merton jump-diffusion put price based on a finite sum approximation to the infinite series solution, making use of the BS put price.
double bs_jd_put_price(const double S, const double K, const double r, const double sigma, const double T, const int N, const double m, const double lambda, const double nu) {
    double price = 0.0;  // Stores the final put price
    double factorial = 1.0;

    // Pre-calculate as much as possible
    double lambda_p = lambda * m;
    double lambda_p_T = lambda_p * T;

    // Calculate the finite sum over N terms
    for (int n = 0; n < N; ++n) {
        double sigma_n = std::sqrt(sigma * sigma + n * nu * nu / T);
        double r_n = r - lambda * (m - 1) + n * std::log(m) / T;

        // Calculate n!
        if (n > 0) {
            factorial *= n;
        }

        // Refine the jump price over the loop
        price += ((std::exp(-lambda_p_T) * std::pow(lambda_p_T, n)) / factorial) * bs_put_price(S, K, r_n, sigma_n, T);
    }

    return price;
}

int main() {
    // First we create the parameter list
    double S = 100.0;     // Option price
    double K = 100.0;     // Strike price
    double r = 0.05;      // Risk-free rate (5%)
    double v = 0.2;       // Volatility of the underlying (20%)
    double T = 1.0;       // One year until expiry
    int N = 50;           // Terms in the finite sum approximation
    double m = 1.083287;  // Scale factor for J
    double lambda = 1.0;  // Intensity of jumps
    double nu = 0.4;      // Stdev of lognormal jump process

    // Then we calculate the call and put black-scholes values
    double call_ bs = bs_call_price(S, K, r, sigma, T);
    double put_ bs = bs_put_price(S, K, r, sigma, T)
        
    // Then we calculate the call and put jump-diffusion values
    double call_jd = bs_jd_call_price(S, K, r, v, T, N, m, lambda, nu);
    double put_jd = bs_jd_put_price(S, K, r, v, T, N, m, lambda, nu);

    std::cout << "Call Price under BS: " << call_bs << std::endl;
    std::cout << "Put Price under BS:  " << put_bs << std::endl;
    
    std::cout << "Call Price under JD: " << call_jd << std::endl;
    std::cout << "Put Price under JD:  " << put_jd << std::endl;

    return 0;
}
