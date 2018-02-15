import random

def costFunction(x):
	sum=0
	for i in range(len(x)):
		sum+= x[i]*x[i]
	return sum

class Particle:
	def __init__(self, x0):
		self.position=[]
		self.velocity=[]
		self.best_pos_in=[]
		self.best_cost_in=-1
		self.cost=-1

		for i in range(0, num_dimensions):
			self.velocity.append(random.uniform(-1, 1))
			self.position.append(x0[i])

	def calculate_cost(self, costFunc):
		self.cost=costFunc(self.position)
		if self.cost < self.best_cost_in or self.best_cost_in==-1:
			self.best_cost_in=self.cost
			self.best_pos_in=self.position

	def update_velocity(self, best_pos_g):
		w=0.5
		c1=2
		c2=2

		for i in range(0, num_dimensions):
			r1=random.random()
			r2=random.random()

			vel_cognitive=c1*r1*(self.best_pos_in[i]-self.position[i])
			vel_social= c2*r2*(best_pos_g[i]-self.position[i])
			self.velocity[i]= w*self.velocity[i]+vel_social+vel_cognitive;


	def update_position(self, bounds):
		for i in range(0, num_dimensions):
			self.position[i]+=self.velocity[i]

			if self.position[i]<bounds[i][0]:
				self.position[i]=bounds[i][0]

			if self.position[i]>bounds[i][1]:
				self.position[i]=bounds[i][1]


class PSO():
	def __init__(self, costFunc, x0, bounds, num_particles, num_iter):
		global num_dimensions
		num_dimensions=len(x0)

		best_cost_g=-1
		best_pos_g=[]

		swarm=[]
		for i in range(0, num_particles):
			swarm.append(Particle(x0))

		for i in range(num_iter):
			for j in range(0, num_particles):
				swarm[j].calculate_cost(costFunc)

				if swarm[j].cost< best_cost_g or best_cost_g==-1:
					best_cost_g=float(swarm[j].cost)
					best_pos_g=list(swarm[j].position)

			for j in range(0, num_particles):
				swarm[j].update_velocity(best_pos_g)
				swarm[j].update_position(bounds)

		print 'Best position : '
		print best_pos_g
		print 'Best cost : '
		print best_cost_g

x0=[5,5]
bounds=[(-10, 10), (-10, 10)]
PSO(costFunction, x0, bounds, num_particles=15, num_iter=100)

