#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

const int D = 5;

struct particle{
    double pos[D];
    double fitness;
};

double evaluate_fitness(double pos[]){
    double sum=0;
    for(int i=0; i<D; i++){
        //fitness function being used here is f = (a-1)^2 + (b-2)^2 + (c-3)^2 + (d-4)^2 + (e-5)^2
        sum+= (pos[i]-(i+1))*(pos[i]-(i+1));
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
    /************************************************************************************************
        VARIABLES COMMON TO BOTH GREYWOLF AND WHALE
    *************************************************************************************************/
    int num_particles, num_iter;
    double bound_lower, bound_upper;
    cout<<"Enter the number of particles:"<<"\n";
    cin>>num_particles;
    num_iter = 100;
    bound_lower = -10;
    bound_upper = 20;
    unsigned int seed=1;

    /*************************************************************************************************
        VARIABLES FOR USAGE WITH GREYWOLF OPTIMIZATION
    **************************************************************************************************/
    //POPULATION INITIALIZATION for grey wolf
    vector<particle> populationGrey;
    for(int i=0; i<num_particles; i++){
        particle p1;
        p1.fitness = -1;
        for(int j=0; j<D; j++){
            p1.pos[j]=(double)(rand_r(&seed))/(double)((2<<25)-1) -10;
            cout<<p1.pos[j]<<" ";
            seed++;
        }
        populationGrey.push_back(p1);
    }
    //INITIALIZE a, A, C
    double aGrey = 2.0;
    double AGrey[D];
    double CGrey[D];
    for(int i=0; i<D; i++){
        AGrey[i] = 2*aGrey*((double)(rand_r(&seed))/(double)((2<<30)-1)) - aGrey;
        seed++;
        CGrey[i] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
        seed++;
    }
    double x1[D];
    double x2[D];
    double x3[D];

    /*************************************************************************************************
        VARIABLES FOR USAGE WITH WHALE OPTIMIZATION
    **************************************************************************************************/
    particle p_best_W;
    //POPULATION INITIALIZATION for whale
    vector<particle> populationWhale;
    for(int i=0; i<num_particles; i++){
        particle p1;
        p1.fitness = -1;
        for(int j=0; j<D; j++){
            p1.pos[j]=(double)(rand_r(&seed))/(double)((2<<25)-1) -10;
            //cout<<p1.pos[j]<<" ";
            seed++;
        }
        populationWhale.push_back(p1);
    }
    double aWhale, l, p;
    double AWhale[D];
    double CWhale[D];

    

    //CALCULATE FITNESS FOR EACH SEARCH AGENT IN populationGrey
    for(int i=0; i<num_particles; i++){
        populationGrey[i].fitness = evaluate_fitness(populationGrey[i].pos);
    }
    //CALCULATE FITNESS FOR EACH SEARCH AGENT IN populationWhale
    for(int i=0; i<num_particles; i++){
        populationWhale[i].fitness = evaluate_fitness(populationWhale[i].pos);
    }

    //IDENTIFY THE THREE BEST SEARCH AGENTS-->for grey wolf
    sort(populationGrey.begin(), populationGrey.end(), compareByFitness);
    particle alpha;
    alpha.fitness = populationGrey[0].fitness;
    for(int i=0; i<D; i++){
        alpha.pos[i] = populationGrey[0].pos[i];
    }
    particle beta;
    beta.fitness = populationGrey[1].fitness;
    for(int i=0; i<D; i++){
        beta.pos[i] = populationGrey[1].pos[i];
    }
    particle delta;
    delta.fitness = populationGrey[2].fitness;
    for(int i=0; i<D; i++){
        delta.pos[i] = populationGrey[2].pos[i];
    }

    //IDENTIFY THE BEST SEARCH AGENT-->for whale
    p_best_W.fitness = populationWhale[0].fitness;
    for(int i=0; i<D; i++){
        p_best_W.pos[i] = populationWhale[0].pos[i];
    }

    //OPTIMIZATION LOOP BEGINS
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                //SECTION CORRESPONDING TO GREY WOLF OPTIMIZATION
                for(int kg=1; kg<=num_iter/2; kg++){
                //UPDATING POSITION OF EACH SEARCH AGENT USING ALPHA, BETA, DELTA
                    for(int i=0; i<num_particles; i++){
                        for(int j=0; j<D; j++){
                            x1[j] = alpha.pos[j] - AGrey[j]*(CGrey[j]*alpha.pos[j] - populationGrey[i].pos[j]);
                            x2[j] = beta.pos[j] - AGrey[j]*(CGrey[j]*beta.pos[j] - populationGrey[i].pos[j]);
                            x3[j] = delta.pos[j] - AGrey[j]*(CGrey[j]*delta.pos[j] - populationGrey[i].pos[j]);
                            populationGrey[i].pos[j] = (x1[j]+x2[j]+x3[j])/3;
                        }
                    }
                    //UPDATE a, A, C
                    aGrey = 2.0 - (2*kg)/num_iter;
                    for(int i=0; i<D; i++){
                        AGrey[i] = 2*aGrey*((double)(rand_r(&seed))/(double)((2<<30)-1)) - aGrey;
                        seed++;
                        CGrey[i] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
                        seed++;
                    }
                    //CALCULATE THE FITNESS OF ALL SEARCH AGENTS
                    for(int i=0; i<num_particles; i++){
                        populationGrey[i].fitness = evaluate_fitness(populationGrey[i].pos);
                    }
                    //UPDATE ALPHA, BETA, DELTA
                    sort(populationGrey.begin(), populationGrey.end(), compareByFitness);
                    if(alpha.fitness > populationGrey[0].fitness){
                        alpha.fitness = populationGrey[0].fitness;
                        for(int i=0; i<D; i++){
                            alpha.pos[i] = populationGrey[0].pos[i];
                        }
                    }
                    if(populationGrey[1].fitness < beta.fitness && populationGrey[1].fitness > alpha.fitness){
                        beta.fitness = populationGrey[1].fitness;
                        for(int i=0; i<D; i++){
                            beta.pos[i] = populationGrey[1].pos[i];
                        }
                    }
                    if(populationGrey[2].fitness > alpha.fitness && populationGrey[2].fitness > beta.fitness && populationGrey[2].fitness < delta.fitness){
                        delta.fitness = populationGrey[2].fitness;
                        for(int i=0; i<D; i++){
                            delta.pos[i] = populationGrey[2].pos[i];
                        }
                    }
                    //cout<<"Fitness at iter "<<k<<" = "<<alpha.fitness<<endl;
                }
            }
            #pragma omp section
            {
                //SECTION CORRESPONDING TO WHALE OPTIMIZATION
                for(int kw=1; kw<=num_iter/2; kw++){
                    for(int i=0; i<num_particles; i++){
                        //UPDATE a, A, C, l and p
                        aWhale =2.0 -  2.0*kw/num_iter;
                        for(int j=0; j<D; j++){
                            AWhale[j] = 2*aWhale*((double)(rand_r(&seed))/(double)((2<<30)-1)) - aWhale;
                            seed++;
                            CWhale[j] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
                            seed++;
                        }
                        l = (double)(rand_r(&seed))/(double)((2<<29)-1) - 1;
                        seed++;
                        p = (double)(rand_r(&seed))/(double)((2<<30)-1);
                        seed++;

                        if(p<0.5){
                            double mag_A = cal_magnitude(AWhale);
                            if(mag_A<1){
                                for(int j=0; j<D; j++){
                                    populationWhale[i].pos[j] = p_best_W.pos[j] - AWhale[j]*(CWhale[j]*p_best_W.pos[j] - populationWhale[i].pos[j]);
                                }
                            }
                            else{
                                int r_i = rand()%num_particles;
                                for(int j=0; j<D; j++){
                                    populationWhale[i].pos[j] = populationWhale[r_i].pos[j] - AWhale[j]*(CWhale[j]*populationWhale[r_i].pos[j] - populationWhale[i].pos[j]);
                                }
                            }
                        }
                        else{
                            for(int j=0; j<D; j++){
                                populationWhale[i].pos[j] = (p_best_W.pos[j] - populationWhale[i].pos[j])*exp(0.206*l)*cos(2*3.141592*l) + p_best_W.pos[j];
                            }
                        }

                        //CHECKING FOR EXCEED OF BOUNDS
                        for(int j=0; j<D; j++){
                            if(populationWhale[i].pos[j] < bound_lower){
                                populationWhale[i].pos[j] = bound_lower;
                            }
                            if(populationWhale[i].pos[j] > bound_upper){
                                populationWhale[i].pos[j] = bound_upper;
                            }
                        }
                        populationWhale[i].fitness = evaluate_fitness(populationWhale[i].pos);
                    }
                    sort(populationWhale.begin(), populationWhale.end(), compareByFitness);
                    if(populationWhale[0].fitness<p_best_W.fitness){
                        p_best_W.fitness = populationWhale[0].fitness;
                        for(int j=0; j<D; j++){
                            p_best_W.pos[j] = populationWhale[0].pos[j];
                        }
                    }
                    //cout<<"p_best fitness at iter "<<k<<"= "<<p_best_W.fitness<<endl;
                }
            }
        }
        #pragma omp barrier
    }

    //OPTIMIZATION PHASE 2 BEGINS..
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                //SECTION CORRESPONDING TO GREY WOLF OPTIMIZATION
                for(int kg=(num_iter/2)+1; kg<=num_iter; kg++){
                //UPDATING POSITION OF EACH SEARCH AGENT USING ALPHA, BETA, DELTA
                    for(int i=0; i<num_particles; i++){
                        for(int j=0; j<D; j++){
                            x1[j] = alpha.pos[j] - AGrey[j]*(CGrey[j]*alpha.pos[j] - populationWhale[i].pos[j]);
                            x2[j] = beta.pos[j] - AGrey[j]*(CGrey[j]*beta.pos[j] - populationWhale[i].pos[j]);
                            x3[j] = delta.pos[j] - AGrey[j]*(CGrey[j]*delta.pos[j] - populationWhale[i].pos[j]);
                            populationWhale[i].pos[j] = (x1[j]+x2[j]+x3[j])/3;
                        }
                    }
                    //UPDATE a, A, C
                    aGrey = 2.0 - (2*kg)/num_iter;
                    for(int i=0; i<D; i++){
                        AGrey[i] = 2*aGrey*((double)(rand_r(&seed))/(double)((2<<30)-1)) - aGrey;
                        seed++;
                        CGrey[i] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
                        seed++;
                    }
                    //CALCULATE THE FITNESS OF ALL SEARCH AGENTS
                    for(int i=0; i<num_particles; i++){
                        populationWhale[i].fitness = evaluate_fitness(populationWhale[i].pos);
                    }
                    //UPDATE ALPHA, BETA, DELTA
                    sort(populationWhale.begin(), populationWhale.end(), compareByFitness);
                    if(alpha.fitness > populationWhale[0].fitness){
                        alpha.fitness = populationWhale[0].fitness;
                        for(int i=0; i<D; i++){
                            alpha.pos[i] = populationWhale[0].pos[i];
                        }
                    }
                    if(populationWhale[1].fitness < beta.fitness && populationWhale[1].fitness > alpha.fitness){
                        beta.fitness = populationWhale[1].fitness;
                        for(int i=0; i<D; i++){
                            beta.pos[i] = populationWhale[1].pos[i];
                        }
                    }
                    if(populationWhale[2].fitness > alpha.fitness && populationWhale[2].fitness > beta.fitness && populationWhale[2].fitness < delta.fitness){
                        delta.fitness = populationWhale[2].fitness;
                        for(int i=0; i<D; i++){
                            delta.pos[i] = populationWhale[2].pos[i];
                        }
                    }
                    //cout<<"Fitness at iter "<<k<<" = "<<alpha.fitness<<endl;
                }
            }
            #pragma omp section
            {
                //SECTION CORRESPONDING TO WHALE OPTIMIZATION
                for(int kw=(num_iter/2)+1; kw<=num_iter; kw++){
                    for(int i=0; i<num_particles; i++){
                        //UPDATE a, A, C, l and p
                        aWhale =2.0 -  2.0*kw/num_iter;
                        for(int j=0; j<D; j++){
                            AWhale[j] = 2*aWhale*((double)(rand_r(&seed))/(double)((2<<30)-1)) - aWhale;
                            seed++;
                            CWhale[j] = 2*((double)(rand_r(&seed))/(double)((2<<30)-1));
                            seed++;
                        }
                        l = (double)(rand_r(&seed))/(double)((2<<29)-1) - 1;
                        seed++;
                        p = (double)(rand_r(&seed))/(double)((2<<30)-1);
                        seed++;

                        if(p<0.5){
                            double mag_A = cal_magnitude(AWhale);
                            if(mag_A<1){
                                for(int j=0; j<D; j++){
                                    populationGrey[i].pos[j] = p_best_W.pos[j] - AGrey[j]*(CWhale[j]*p_best_W.pos[j] - populationGrey[i].pos[j]);
                                }
                            }
                            else{
                                int r_i = rand()%num_particles;
                                for(int j=0; j<D; j++){
                                    populationGrey[i].pos[j] = populationGrey[r_i].pos[j] - AWhale[j]*(CWhale[j]*populationGrey[r_i].pos[j] - populationGrey[i].pos[j]);
                                }
                            }
                        }
                        else{
                            for(int j=0; j<D; j++){
                                populationGrey[i].pos[j] = (p_best_W.pos[j] - populationGrey[i].pos[j])*exp(0.206*l)*cos(2*3.141592*l) + p_best_W.pos[j];
                            }
                        }

                        //CHECKING FOR EXCEED OF BOUNDS
                        for(int j=0; j<D; j++){
                            if(populationGrey[i].pos[j] < bound_lower){
                                populationGrey[i].pos[j] = bound_lower;
                            }
                            if(populationGrey[i].pos[j] > bound_upper){
                                populationGrey[i].pos[j] = bound_upper;
                            }
                        }
                        populationGrey[i].fitness = evaluate_fitness(populationGrey[i].pos);
                    }
                    sort(populationGrey.begin(), populationGrey.end(), compareByFitness);
                    if(populationGrey[0].fitness<p_best_W.fitness){
                        p_best_W.fitness = populationGrey[0].fitness;
                        for(int j=0; j<D; j++){
                            p_best_W.pos[j] = populationGrey[0].pos[j];
                        }
                    }
                    //cout<<"p_best fitness at iter "<<k<<"= "<<p_best_W.fitness<<endl;
                }
            }
        }
        #pragma omp barrier
    }
    
    //COMPARISON OF THE RESULTS OBTAINED FROM THE TWO OPTIMIZATION PARADIGMS
    if(alpha.fitness > p_best_W.fitness){
        cout<<"Best Solution:"<<"\n";
        cout<<"Fitness:"<<alpha.fitness<<"\n";
        for(int s=0; s<D; s++){
            cout<<alpha.pos[s]<<" ";
        }
        cout<<endl;
    }
    else{
        cout<<"Best Solution:"<<"\n";
        cout<<"Fitness:"<<p_best_W.fitness<<"\n";
        for(int s=0; s<D; s++){
            cout<<p_best_W.pos[s]<<" ";
        }
        cout<<endl;
    }


    return 0;
}