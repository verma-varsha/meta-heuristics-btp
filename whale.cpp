#include <bits/stdc++.h>
using namespace std;

const int D = 2;

struct particle{
    double pos[D];
    double fitness;
};

double evaluate_fitness(double pos[]){
    double sum=0;
    for(int i=0; i<D; i++){
        //fitness function being used here is f = (x-4)^2 + (y-4)^2
        sum+= (pos[i]-4)*(pos[i]-4);
    }
    return sum;
}

bool compareByFitness(const particle &a, const particle &b)
{
    return a.fitness < b.fitness;
}

double cal_magnitude(double arr[]){
    double value = 0.0;
    for(int i=0; i<D; i++){
        value+= arr[i]*arr[i];
    }
    return value;
}

int main(){
    int num_iter, num_particles;
    cout<<"Enter the number of iterations:";
    cin>>num_iter;
    cout<<"Enter the number of particles:";
    cin>>num_particles;
    particle p_best;
    double a, l, p;
    double A[D];
    double C[D];
    double bound_lower, bound_upper;
    cout<<"Enter lower bound:";
    cin>>bound_lower;
    cout<<"Enter upper bound:";
    cin>>bound_upper;
    //POPULATION INITIALIZATION
    vector<particle> population;
    unsigned int seed=1;
    for(int i=0; i<num_particles; i++){
        particle p1;
        p1.fitness = -1;
        for(int j=0; j<D; j++){
            p1.pos[j]=(double)(rand_r(&seed))/(double)((2<<25)-1) -10;
            //cout<<p1.pos[j]<<" ";
            seed++;
        }
        population.push_back(p1);
    }
    //EVALUATING FITNESS FOR EACH SEARCH AGENT
    for(int i=0; i<num_particles; i++){
        population[i].fitness = evaluate_fitness(population[i].pos);
    }
    //BEST SEARCH AGENT
    sort(population.begin(), population.end(), compareByFitness);
    p_best.fitness = population[0].fitness;
    for(int i=0; i<D; i++){
        p_best.pos[i] = population[0].pos[i];
    }
    //OPTIMIZATION LOOP BEGINS..
    for(int k=1; k<=num_iter; k++){
        for(int i=0; i<num_particles; i++){
            //UPDATE a, A, C, l and p
            a =2.0 -  2.0*k/num_iter;
            for(int j=0; j<D; j++){
                A[j] = 2*a*((double)(rand_r(&seed))/(double)((2<<30)-1)) - a;
                seed++;
                C[j] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
                seed++;
            }
            l = (double)(rand_r(&seed))/(double)((2<<29)-1) - 1;
            seed++;
            p = (double)(rand_r(&seed))/(double)((2<<30)-1);
            seed++;

            if(p<0.5){
                double mag_A = cal_magnitude(A);
                if(mag_A<1){
                    for(int j=0; j<D; j++){
                        population[i].pos[j] = p_best.pos[j] - A[j]*(C[j]*p_best.pos[j] - population[i].pos[j]);
                    }
                }
                else{
                    int r_i = rand()%num_particles;
                    for(int j=0; j<D; j++){
                        population[i].pos[j] = population[r_i].pos[j] - A[j]*(C[j]*population[r_i].pos[j] - population[i].pos[j]);
                    }
                }
            }
            else{
                for(int j=0; j<D; j++){
                    population[i].pos[j] = (p_best.pos[j] - population[i].pos[j])*exp(0.206*l)*cos(2*3.141592*l) + p_best.pos[j];
                }
            }

            //CHECKING FOR EXCEED OF BOUNDS
            for(int j=0; j<D; j++){
                if(population[i].pos[j] < bound_lower){
                    population[i].pos[j] = bound_lower;
                }
                if(population[i].pos[j] > bound_upper){
                    population[i].pos[j] = bound_upper;
                }
            }
            population[i].fitness = evaluate_fitness(population[i].pos);
        }
        sort(population.begin(), population.end(), compareByFitness);
        if(population[0].fitness<p_best.fitness){
            p_best.fitness = population[0].fitness;
            for(int j=0; j<D; j++){
                p_best.pos[j] = population[0].pos[j];
            }
        }
        cout<<"p_best fitness at iter "<<k<<"= "<<p_best.fitness<<endl;
    }
    

    cout<<"BEST SOLUTION: ";
    for(int j=0; j<D; j++){
        cout<<p_best.pos[j]<<" ";
    }
    cout<<endl;
    return 0;
}