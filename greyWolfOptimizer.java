import java.io.*;
import java.util.*;
import java.text.*;
import java.math.*;
import java.util.regex.*;

class greyWolf implements Comparable<greyWolf>{
    public static final int dimensions = 2;
    double[] pos_i = new double[dimensions];
    double fitness = -1;

    greyWolf(){
        for(int j = 0; j<dimensions; j++){
            pos_i[j] = (20.0 * Math.random()) -10 ;
        }
    }

    public int compareTo(greyWolf gg){
        double fitgg = gg.fitness;
        return (int)(this.fitness-fitgg);
    }

    double evaluate_fitness(){
        double sum = 0;
        for(int i=0; i<dimensions; i++){
            sum = sum + (pos_i[i] - 4)*(pos_i[i] - 4);    
        }
        return sum;
    }
}

public class greyWolfOptimizer{

    static void sort_with_fitness(greyWolf[] population){
        Arrays.sort(population);
    }

    public static void main(String[] args){
        Scanner sc = new Scanner(System.in);
        int num_wolves, num_iters;
        System.out.println("Enter number of iterations:");
        num_iters = sc.nextInt();
        System.out.println("Enter wolf population size:");
        num_wolves = sc.nextInt();

        //INITIALIZE THE GREY WOLF POPULATION
        greyWolf[] population = new greyWolf[num_wolves];
        //List<greyWolf> population = new ArrayList<greyWolf>();
        for(int i=0; i<num_wolves; i++){
            //greyWolf x = new greyWolf();
            //population.add(x);
            population[i] = new greyWolf();       
        }
        //INITIALIZE a, A and C
        double a = 2.0;
        double[] A = new double[greyWolf.dimensions];
        double[] C = new double[greyWolf.dimensions];
        for(int i =0; i<greyWolf.dimensions; i++){
            A[i] = 2*a*Math.random() - a;
            C[i] = 2*Math.random();
        } 

        //CALCULATE FITNESS FOR EACH SEARCH AGENT
        for(int i=0; i<num_wolves; i++){
            population[i].fitness = population[i].evaluate_fitness();
        }

        //IDENTIFY THE THREE BEST SEARCH AGENTS
        sort_with_fitness(population);
        greyWolf alpha = new greyWolf();
        alpha.fitness = population[0].fitness;
        alpha.pos_i = population[0].pos_i.clone();
        greyWolf beta = new greyWolf();
        beta.fitness = population[1].fitness;
        beta.pos_i = population[1].pos_i.clone();
        greyWolf delta = new greyWolf();
        delta.fitness = population[2].fitness;
        delta.pos_i = population[2].pos_i.clone();

        //OPTIMIZATION LOOP BEGINS
        for(int iter=1; iter<=num_iters; iter++){
            //UPDATING POSITION OF EACH SEARCH AGENT USING ALPHA, BETA AND DELTA
            double[] d_alpha = new double[greyWolf.dimensions];
            double[] d_beta = new double[greyWolf.dimensions];
            double[] d_delta = new double[greyWolf.dimensions];
            double[] x1 = new double[greyWolf.dimensions];
            double[] x2 = new double[greyWolf.dimensions];
            double[] x3 = new double[greyWolf.dimensions];
            for(int i =0; i<num_wolves; i++){
                for(int j=0; j<greyWolf.dimensions; j++){
                    d_alpha[j] = C[j]*alpha.pos_i[j] - population[i].pos_i[j];
                    x1[j] = alpha.pos_i[j] - A[j]*d_alpha[j];
                    d_beta[j] = C[j]*beta.pos_i[j] - population[i].pos_i[j];
                    x2[j] = beta.pos_i[j] - A[j]*d_beta[j];
                    d_delta[j] = C[j]*delta.pos_i[j] - population[i].pos_i[j];
                    x3[j] = delta.pos_i[j] - A[j]*d_delta[j];
                    population[i].pos_i[j] = (x1[j]+x2[j]+x3[j])/3;
                }
            }
            //PRINTING VALUE OF A
            double mag_A = A[0]*A[0] + A[1]*A[1];
            //System.out.println("Magnitude of A at iter "+ iter+" = "+mag_A);
            //System.out.println("Magnitude of a at iter "+iter+" = "+a);
            //UPDATE a, A, C
            a = 2.0 - (2.0*iter)/num_iters;
            for(int i=0; i<greyWolf.dimensions; i++){
                A[i]=2*a*Math.random() - a;
                C[i]=2*Math.random();
            }
            //CALCULATE THE FITNESS OF ALL SEARCH AGENTS
            for(int i=0; i<num_wolves; i++){
                population[i].fitness = population[i].evaluate_fitness();
            }
            //UPDATE ALPHA, BETA AND DELTA
            sort_with_fitness(population);
            if(alpha.fitness > population[0].fitness){
                System.out.println("Updating alpha from "+alpha.fitness+" to"+population[0].fitness);
                alpha.fitness = population[0].fitness;
                alpha.pos_i = population[0].pos_i.clone();
            }
            if(population[1].fitness < beta.fitness && population[1].fitness > alpha.fitness){
                beta.fitness = population[1].fitness;
                beta.pos_i = population[1].pos_i.clone();
            }
            if(population[2].fitness > alpha.fitness && population[2].fitness > beta.fitness && population[2].fitness < delta.fitness){
                delta.fitness = population[2].fitness;
                delta.pos_i = population[2].pos_i.clone();
            }
              
            System.out.println("Fitness at iter "+iter+"="+alpha.fitness);
        }
        //PRINTING OUTPUT
        System.out.println("Predicted output:");
        for(int l=0; l<greyWolf.dimensions; l++){
            System.out.print(alpha.pos_i[l]+" ");
        }
        System.out.println();
        System.out.println("Fitness of the best solution = "+alpha.fitness);
    }
}