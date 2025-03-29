#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


/*
 * @brief Wrapper for rand() to make uniform prng
 * 
 * @param[out] Uniform random number in [0, 1)
 */
double rand_uniform(){
	 return (double)rand() / (RAND_MAX + 1.0);
}

/*
 * @brief Custom normal prng
 * 
 * @param[out] Normal random number
 */
double rand_normal(){
	for (;;){
		double u1, u2;
		u1 = rand_uniform();
		u2 = rand_uniform();

		double v1 = (2*u1 - 1);
		double v2 = (2*u2 - 1);
		double s = v1*v1 + v2*v2;
		if (s > 1)
			continue;

		double z1 = v1 * sqrt( (-2 * log(s))/s);
		//double z2 = v2 * sqrt( (-2 * log(s))/s);
	
		return z1;
	}
}


/*
 * @brief Monte carlo pricer for european option
 * 
 * @param n monte carlo samples
 * @param S_0 Spot price today 
 * @param k strike price
 * @param sigma Annual volatility
 * @param r Interest rate
 */
double monte_carlo_pricer(int n, double S_0, double k, double sigma, double r, int seed){
	double t = 1; // expire time
	double mu_t = (r - 0.5 * sigma * sigma) * t;
	double sigma_sqrt_t = sigma * sqrt(t);
	double total_payoff = 0;

	#pragma omp parallel num_threads(4)
	{
		int id = omp_get_thread_num();
		int nt = omp_get_num_threads();
		int N = n/nt;

		srand(seed);
		double* rvs = malloc(N * sizeof(double));
		for (int i=0; i<(id + 1)*N; i++){
			rand_normal();	
		}
		for (int i=0; i<N; i++){
			rvs[i] = rand_normal();	
		}

		double local_payoff = 0;
		double S_T;
		for (int i=0; i<N; i++){
			S_T = S_0 * exp(mu_t + sigma_sqrt_t * rvs[i]);
			local_payoff += fmax(S_T - k, 0);
		}
		
		free(rvs);
		
		#pragma omp atomic	
		total_payoff += local_payoff;
	}
	
	return exp(-r * t) * (total_payoff / n);
} 


int main(){
	int n = 1000000;
	double S_0 = 10;
	double k = 2;
	double sigma = 2;
	double r = 3;	
	int seed = 69420;		
	
	double val;
	val = monte_carlo_pricer(n, S_0, k, sigma, r, seed);
	
	printf("val = %lf", val);
	
	
	return 0;
}
