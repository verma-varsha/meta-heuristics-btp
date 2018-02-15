#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <bits/stdc++.h>

using namespace std;

double costFunc(vector<double> arr, int num_dimensions){
	double sum=0;
	for(int i=0; i<num_dimensions; i++){
		sum = sum + (arr[i]*arr[i]);
	}
	return sum;
}

struct Particle{
		vector<double> position;
		vector<double> velocity;
		vector<double> best_position_in;
		double best_cost_in;
		double cost_in;
	};

int main(){
	srand(time(NULL));
	int num_particles, num_dimensions, num_iter;
	double bound_lower, bound_upper;
	cout<<"Number of particles: ";
	cin>>num_particles;
	cout<<"Number of dimensions: ";
	cin>>num_dimensions;
	cout<<"Lower bound: ";
	cin>>bound_lower;
	cout<<"Upper bound: ";
	cin>>bound_upper;
	cout<<"Number of iterations: ";
	cin>>num_iter;
	
	//SWARM PARTICLES INITIALIZATION
	vector<Particle> swarm;
	for(int i=0; i<num_particles; i++){
		Particle p1;
		unsigned int seed=1;
		for(int j=0; j<num_dimensions; j++){
			p1.position.push_back(rand()%100+1);
			p1.velocity.push_back(rand()%100+1);
			p1.best_cost_in=-1;
			p1.cost_in=-1;
			seed++;
		}
		swarm.push_back(p1);	
	}

	vector<double> best_position_g;
	double best_cost_g=-1;

	double w=1; //CONSTANT INERTIA WEIGHT (HOW MUCH TO WEIGH THE PREVIOUS VELOCITY)
	double c1=2; //COGNITIVE CONSTANT
	double c2=2; //SOCIAL CONSTANT

	//PERFORMING REQUIRED ITERATIONS OF OPTIMIZATION LOOP
	for(int i=0; i<num_iter; i++){

		//CALCULATING THE GLOBAL BEST PARTICLE WALKING THROUGH & EVAUATING COST OF ALL PARTICLES
		for(int j=0; j<num_particles; j++){
			swarm[j].cost_in = costFunc(swarm[j].position, num_dimensions);
			//UPDATING THE CURRENT BEST
			if(swarm[j].cost_in<swarm[j].best_cost_in || swarm[j].best_cost_in==-1){
				swarm[j].best_cost_in=swarm[j].cost_in;
				swarm[j].best_position_in = swarm[j].position;
			}
			//UPDATING THE GLOBAL BEST
			if(swarm[j].cost_in<best_cost_g || best_cost_g==-1){
				//cout<<"Hola from global best."<<endl;
				best_cost_g=swarm[j].cost_in;
				best_position_g = swarm[j].position;
			}
		}

		//WALKING THROUGH THE SWARM AND UPDATING VELOCITIES AND POSITIONS
		for(int j=0; j<num_particles; j++){

			//UPDATE VELOCITY
			unsigned int seed = 1;
			for(int k=0; k<num_dimensions; k++){
				double r1=((double)rand()/(double)(2<<30-1));
				double r2=((double)rand()/(double)(2<<30-1));
				double vel_cognitive = r1*c1*(swarm[j].best_position_in[k]-swarm[j].position[k]);
				double vel_social = r2*c2*(best_position_g[k]-swarm[j].position[k]);
				swarm[j].velocity[k] = (w*swarm[j].velocity[k]) + vel_cognitive + vel_social;
				//cout<<vel_cognitive<<" "<<vel_social<<" "<<w*swarm[j].velocity[k]<<endl;
				seed++;
			}

			//UPDATE POSITION
			for(int k=0; k<num_dimensions; k++){
				swarm[j].position[k] += swarm[j].velocity[k];

				//CHECKING FOR EXCEED OF LOWER BOUND
				if(swarm[j].position[k] < bound_lower)
					swarm[j].position[k] = bound_lower;

				//CHECKING FOR EXCEED OF UPPER BOUND
				if(swarm[j].position[k] > bound_upper)
					swarm[j].position[k] = bound_upper;
			}
		}
	}

	//DISPLAYING THE RESULTS
	cout<<"BEST POSITION = ";
	for(int i=0; i<num_dimensions; i++){
		cout<<best_position_g[i]<<" ";
	}
	cout<<endl;

	cout<<"BEST COST = "<<best_cost_g<<endl;
}