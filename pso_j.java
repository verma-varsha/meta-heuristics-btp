import java.io.*;
import java.util.*;
import java.text.*;
import java.math.*;
import java.util.regex.*;

class Particle{
	public static final int num_dimensions = 3;
	double[] position_i = new double[num_dimensions];
	double[] velocity_i = new double[num_dimensions];
	double[] pos_best_i = new double[num_dimensions];
	double err_best_i = -1;
	double err_i = -1;	

	Particle(){
		Random randNum = new Random();
		for(int i=0; i<num_dimensions; i++){
			velocity_i[i] = -2 + randNum.nextInt(2);
			position_i[i] = 5;
		}
	}
	
	void evaluate(){
		err_i = 0;
		for(int i = 0; i<num_dimensions; i++){
			err_i+=(position_i[i]*position_i[i]);
		}
		if(err_i < err_best_i || err_best_i == -1){
			for(int j = 0; j<num_dimensions; j++){
				pos_best_i[j] = position_i[j];			
			}
			err_best_i = err_i;
		}
	}

	void update_velocity(double[] pos_best_g){
		double w = 0.5;
		double c1 = 1;
		double c2 = 2;
		Random randNum = new Random();
		for(int i=0; i<num_dimensions; i++){
			double r1 = -2 + randNum.nextInt(2);
			double r2 = -2 + randNum.nextInt(2);

			double vel_cognitive = c1*r1*(pos_best_i[i] - position_i[i]);
			double vel_social = c2*r2*(pos_best_g[i] - position_i[i]);
			velocity_i[i] = w*velocity_i[i] + vel_cognitive + vel_social;
		}
	}
	
	void update_position(double bound_upper, double bound_lower){
		for(int i=0; i<num_dimensions; i++){
			position_i[i] = position_i[i] + velocity_i[i];
			if(position_i[i]>bound_upper)
				position_i[i] = bound_upper;
			if(position_i[i]<bound_lower)
				position_i[i] = bound_lower;
		}
	}
}

public class hello{

    public static void main(String[] args) {
	Scanner sc = new Scanner(System.in);
    int num_particles, num_iter;
	double bound_lower, bound_upper;
	System.out.println("Number of particles: ");
	num_particles = sc.nextInt();
	System.out.println("Lower bound: ");
	bound_lower = sc.nextDouble();
	System.out.println("Upper bound: ");
	bound_upper = sc.nextDouble();
	System.out.println("Number of iterations: ");
	num_iter = sc.nextInt();
	double err_best_g = -1;
	double[] pos_best_g = new double[Particle.num_dimensions];

	Particle[] swarm = new Particle[num_particles];
	for(int j = 0; j<num_particles; j++){
		Particle p1 = new Particle();
		swarm[j] = p1;
	}

	for(int i = 0; i<num_iter; i++){
		for(int j=0; j<num_particles; j++){
			swarm[j].evaluate();

			if(swarm[j].err_i < err_best_g || err_best_g == -1){
				for(int k = 0; k<Particle.num_dimensions; k++){
					pos_best_g[k] = swarm[j].position_i[k];
				} 
				err_best_g = swarm[j].err_i;
			}
		}

		for(int j =0; j<num_particles; j++){
			swarm[j].update_velocity(pos_best_g);
			swarm[j].update_position(bound_upper, bound_lower);
		}
		System.out.println("error "+ err_best_g+"\n");
		System.out.println("position :");
		for(int g =0; g<Particle.num_dimensions; g++){
			System.out.print(pos_best_g[g]+" ");
		}
		System.out.print("\n");
	}

	System.out.println("Best Position:");
	for(int i =0; i<Particle.num_dimensions; i++){
		System.out.print(pos_best_g[i]+" ");
	}
	System.out.println("\nBest Error:");
	System.out.println(err_best_g);

    }
}	
