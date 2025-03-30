#include <stdio.h>
#include <math.h>
#include <mpi.h>
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
	int proc;
    int num_procs;
    //MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	double t = 1; // expire time
	double mu_t = (r - 0.5 * sigma * sigma) * t;
	double sigma_sqrt_t = sigma * sqrt(t);
	double total_payoff = 0;
	double local_payoff = 0;

	int N = n/num_procs;
	printf("N = %d\n", N);
	
	double* rvs = malloc(N * sizeof(double));

	const gsl_rng_type* T;
	gsl_rng* rng;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);

	for (int i=0; i<proc*N; i++){
		gsl_ran_gaussian(rng, 1);	
	}
	for (int i=0; i<N; i++){
		rvs[i] = gsl_ran_gaussian(rng, 1);	
	}

	gsl_rng_free(rng);
	
	printf("proc = %d\n", proc);
	for (int i=proc*N; i<(proc+1)*N; i++){
		double S_T = S_0 * exp(mu_t + sigma_sqrt_t * rvs[i]);
		local_payoff += fmax(S_T - k, 0);
	}
	free(rvs);
	
	MPI_Reduce(&local_payoff, &total_payoff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
	return exp(-r * t) * (total_payoff / n);
} 

	
int main(){
	int rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0){
		int n = 100;
		double S_0 = 10;
		double k = 2;
		double sigma = 2;
		double r = 3;	
		//int seed = 69420;		
			
		double val;
		val = monte_carlo_pricer(n, S_0, k, sigma, r);
		
		printf("val = %lf", val);
	}
	MPI_Finalize();
	return 0;
	
}
