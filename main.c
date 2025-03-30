#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/*
 * @brief Monte carlo pricer for european option
 * 
 * @param n monte carlo samples
 * @param S_0 Spot price today 
 * @param k strike price
 * @param sigma Annual volatility
 * @param r Interest rate
 */
double monte_carlo_pricer(int n, double S_0, double k, double sigma, double r){
	double t = 1; // expire time
	double mu_t = (r - 0.5 * sigma * sigma) * t;
	double sigma_sqrt_t = sigma * sqrt(t);
	double total_payoff = 0;
		
	int seed = 69420;
	gsl_rng_env_setup();
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);
	
	double* rvs = malloc(n * sizeof(double));
	for (int i=0; i<n; i++){
        rvs[i] = gsl_ran_gaussian(rng, 1);
    }
    gsl_rng_free(rng);
	
	for (int i=0; i<n; i++){
		//double rv = gsl_ran_gaussian(rng, 1);
		double S_T = S_0 * exp(mu_t + sigma_sqrt_t * rvs[i]);
		total_payoff += fmax(S_T - k, 0);
	}
	free(rvs);
	
	return exp(-r * t) * (total_payoff / n);
} 


double monte_carlo_pricer_antithetic(int n, double S_0, double k, double sigma, double r){
	double t = 1; // expire time
	double mu_t = (r - 0.5 * sigma * sigma) * t;
	double sigma_sqrt_t = sigma * sqrt(t);
	double total_payoff = 0;
		
	int seed = 69420;
	gsl_rng_env_setup();
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);
	
	double* rvs = malloc(n * sizeof(double));
	for (int i=0; i<n; i+=2){
        double rv = gsl_ran_gaussian(rng, 1);
        rvs[i]   = rv;
        rvs[i+1] = -rv;
    }

    gsl_rng_free(rng);
	
	for (int i=0; i<n; i++){
		double S_T = S_0 * exp(mu_t + sigma_sqrt_t * rvs[i]);
		total_payoff += fmax(S_T - k, 0);
	}
	free(rvs);
	
	return exp(-r * t) * (total_payoff / n);
} 


double monte_carlo_pricer_sobol(int n, double S_0, double k, double sigma, double r){
	double t = 1; // expire time
	double mu_t = (r - 0.5 * sigma * sigma) * t;
	double sigma_sqrt_t = sigma * sqrt(t);
	double total_payoff = 0;

	gsl_qrng* qrng = gsl_qrng_alloc(gsl_qrng_sobol, 1);
	
	double* rvs = malloc(n * sizeof(double));
	double rv;
	for (int i=0; i<n; i++){
		gsl_qrng_get(qrng, &rv);
        rvs[i] = gsl_cdf_ugaussian_Pinv(rv);
    }
    gsl_qrng_free(qrng);
	
	for (int i=0; i<n; i++){
		double S_T = S_0 * exp(mu_t + sigma_sqrt_t * rvs[i]);
		total_payoff += fmax(S_T - k, 0);
	}
	free(rvs);
	
	return exp(-r * t) * (total_payoff / n);
} 

int main(){
	int n = 100000;
	double S_0 = 10;
	double k = 2;
	double sigma = 2;
	double r = 3;	
	//int seed = 69420;

	double val;
    val = monte_carlo_pricer(n, S_0, k, sigma, r);
    printf("default = %lf\n", val);	
    val = monte_carlo_pricer_antithetic(n, S_0, k, sigma, r);
    printf("antithetic = %lf\n", val);	
    val = monte_carlo_pricer_sobol(n, S_0, k, sigma, r);
    printf("sobol = %lf\n", val);	

	return 0;
}
