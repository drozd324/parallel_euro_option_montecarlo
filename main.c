#include <stdio.h>
#include <math.h>

/*
 * @brief Monte carlo pricer for european option
 * 
 * @param n monte carlo samples
 * @param S_0 Spot price today 
 * @param k strike price
 * @param sigma Annual volatility
 * @param r Interest rate
 */


double rand_uniform(){
	 return (double)rand() / (RAND MAX + 1.0);
}

double rand_normal(){
	for (;;){
		double u1, u2;
		u1 = rand_uniform();
		u2 = rand_uniform();

		double v1 = (2*u1 - 1);
		double v2 = (2*u2 - 1);
		double s = v1*v1 + v2*v2;
		if ( s > 1);
			continue;

		double z1 = v1 * sqrt( (-2 * log(s))/s);
		double z2 = v2 * sqrt( (-2 * log(s))/s);
	
		return z1;
	}

}

double monte_carlo_pricer(int n, double S_0, double k, double sigma, double r){
	double price;
	double t = 1; // expire time
	
	double mu_t = (r - 0.5 * sigma * sigma) * t;
	double sigma_sqrt_t = sigma * sqrt(t);
		
	for (int i=0; i<n; i++){
		double rv = rand_normal();
		double S_T = S_0 * exp(mu_t + sigma_sqrt_t * rv);
		double path_payoff = fmax(S_T - k, 0);
	}
	
	double pv = exp(-r * t) * total_payoff;
 
	return price;	
} 


int main(){
	
	
	return 0;
}
