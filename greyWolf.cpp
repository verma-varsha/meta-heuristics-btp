#include <bits/stdc++.h>
using namespace std;
const int D = 5;

struct particle{
    double fitness;
    double pos[D];
};

double evaluate_fitness(double pos[]){
    double sum=0;
    for(int i=0; i<D; i++){
        //fitness function being used here is f = (x-4)^2 + (y-4)^2
        sum+= (pos[i]-(i+1))*(pos[i]-(i+1));
    }
    return sum;
}

bool compareByFitness(const particle &a, const particle &b)
{
    return a.fitness < b.fitness;
}

int main(){
    int num_iter, num_particles;
    cout<<"Enter the number of iterations:";
    cin>>num_iter;
    cout<<"Enter the number of particles:";
    cin>>num_particles;
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
            cout<<p1.pos[j]<<" ";
            seed++;
        }
        population.push_back(p1);
    }

    //INITIALIZE a, A, C
    double a = 2.0;
    double A[D];
    double C[D];
    for(int i=0; i<D; i++){
        A[i] = 2*a*((double)(rand_r(&seed))/(double)((2<<30)-1)) - a;
        seed++;
        C[i] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
        seed++;
    }

    //CALCULATE FITNESS FOR EACH SEARCH AGENT
    for(int i=0; i<num_particles; i++){
        population[i].fitness = evaluate_fitness(population[i].pos);
    }

    //IDENTIFY THE THREE BEST DEARCH AGENTS
    sort(population.begin(), population.end(), compareByFitness);
    particle alpha;
    alpha.fitness = population[0].fitness;
    for(int i=0; i<D; i++){
        alpha.pos[i] = population[0].pos[i];
    }
    particle beta;
    beta.fitness = population[1].fitness;
    for(int i=0; i<D; i++){
        beta.pos[i] = population[1].pos[i];
    }
    particle delta;
    delta.fitness = population[2].fitness;
    for(int i=0; i<D; i++){
        delta.pos[i] = population[2].pos[i];
    }

    //OPTIMIZATION LOOP BEGINS
    double x1[D];
    double x2[D];
    double x3[D];
    for(int k=1; k<=num_iter; k++){
        //UPDATING POSITION OF EACH SEARCH AGENT USING ALPHA, BETA, DELTA
        for(int i=0; i<num_particles; i++){
            for(int j=0; j<D; j++){
                x1[j] = alpha.pos[j] - A[j]*(C[j]*alpha.pos[j] - population[i].pos[j]);
                x2[j] = beta.pos[j] - A[j]*(C[j]*beta.pos[j] - population[i].pos[j]);
                x3[j] = delta.pos[j] - A[j]*(C[j]*delta.pos[j] - population[i].pos[j]);
                population[i].pos[j] = (x1[j]+x2[j]+x3[j])/3;
            }
        }
        //UPDATE a, A, C
        a = 2.0 - (2*k)/num_iter;
        for(int i=0; i<D; i++){
            A[i] = 2*a*((double)(rand_r(&seed))/(double)((2<<30)-1)) - a;
            seed++;
            C[i] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
            seed++;
        }
        //CALCULATE THE FITNESS OF ALL SEARCH AGENTS
        for(int i=0; i<num_particles; i++){
            population[i].fitness = evaluate_fitness(population[i].pos);
        }
        //UPDATE ALPHA, BETA, DELTA
        sort(population.begin(), population.end(), compareByFitness);
        if(alpha.fitness > population[0].fitness){
            alpha.fitness = population[0].fitness;
            for(int i=0; i<D; i++){
                alpha.pos[i] = population[0].pos[i];
            }
        }
        if(population[1].fitness < beta.fitness && population[1].fitness > alpha.fitness){
            beta.fitness = population[1].fitness;
            for(int i=0; i<D; i++){
                beta.pos[i] = population[1].pos[i];
            }
        }
        if(population[2].fitness > alpha.fitness && population[2].fitness > beta.fitness && population[2].fitness < delta.fitness){
            delta.fitness = population[2].fitness;
            for(int i=0; i<D; i++){
                delta.pos[i] = population[2].pos[i];
            }
        }
        cout<<"Fitness at iter "<<k<<" = "<<alpha.fitness<<endl;
    }
    cout<<"Predicted output:";
    for(int l=0; l<D; l++){
        cout<<alpha.pos[l]<<" ";
    }
    cout<<endl;
    return 0;
}