#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

	double* rvs = malloc(n * sizeof(double));

	const gsl_rng_type* T;
	gsl_rng* rng;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);

	for (int i=0; i<n; i++){
		rvs[i] = gsl_ran_gaussian(rng, 1);	
	}
	gsl_rng_free(rng);

	#pragma omp parallel num_threads(2) reduction(+:total_payoff) 
	{
		int id = omp_get_thread_num();
		int nt = omp_get_num_threads();

		printf("id = %d\n", id);
		int N = n/nt;
		for (int i=id*N; i<(id+1)*N; i++){
			//double rv = gsl_ran_gaussian(rng, 1);
			double S_T = S_0 * exp(mu_t + sigma_sqrt_t * rvs[i]);
			total_payoff += fmax(S_T - k, 0);
		}
	}
	
	free(rvs);
	return exp(-r * t) * (total_payoff / n);
} 

	
int main(){
	int n = 100000000;
	double S_0 = 10;
	double k = 2;
	double sigma = 2;
	double r = 3;	
	//int seed = 69420;		
		

	double val;
	val = monte_carlo_pricer(n, S_0, k, sigma, r);
	
	printf("val = %lf", val);
	
	
	return 0;
}
