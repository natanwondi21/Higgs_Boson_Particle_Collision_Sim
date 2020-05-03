#!/usr/bin/env python
# coding: utf-8

# In[9]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[1]:


import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations
get_ipython().run_line_magic('matplotlib', '')


# In[2]:


class Particle:
    """A class representing a two-dimensional particle."""

    def __init__(self, x, y, vx, vy, radius=0.01, styles=None,Name='P'):
        """Initialize the particle's position, velocity, and radius.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor.

        """
        self.name=Name
        self.t=0
        self.r = np.array((x, y))
        self.v = np.array((vx, vy))
        self.radius = radius
        self.styles = styles
        if not self.styles:
            # Default circle styles
            self.styles = {'edgecolor': 'b', 'fill': True}

    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value

    def overlaps(self, other):
        """Does the circle of this Particle overlap that of other?"""

        return np.hypot(*(self.r - other.r)) < self.radius + other.radius

    def draw(self, ax):
        """Add this Particle's Circle patch to the Matplotlib Axes ax."""

        circle = Circle(xy=self.r, radius=self.radius, **self.styles)
        ax.add_patch(circle)
        plt.draw()
        
        #label = ax.annotate(self.name, xy=self.r, fontsize=10, ha="center")
        return circle #label

    def advance(self, dt):
        """Advance the Particle's position forward in time by dt."""

        self.r += self.v * dt

class Simulation:
    """A class for a simple hard-circle molecular dynamics simulation.

    The simulation is carried out on a square domain: 0 <= x < 1, 0 <= y < 1.

    """

    def __init__(self, n,P,radius=0.01, styles=None):
        """Initialize the simulation with n Particles with radii radius.
        radius can be a single value or a sequence with n values.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor when drawing
        the Particles.

        """
        self.t=0
        self.name=[]
        self.n = n
        self.particles = []
        for i in range(len(P)):
            self.particles.append(P[i])
            
        #self.init_particles(n, radius, styles)

    #def init_particles(self, n, radius, styles=None):
            
    def handle_collisions(self):
        """Detect and handle any collisions between the Particles.

        When two Particles collide, they do so inelastically: their velocities
        change such that both energy and momentum are conserved.

        """

        def change_velocities(p1, p2):
            """
            Particles p1 and p2 have collided inelastically: update their
            velocities.

            """

            m1, m2 = p1.radius**2, p2.radius**2
            M = m1 + m2
            r1, r2 = p1.r, p2.r
            v1, v2 = p1.v, p2.v
            v=0.1*(m1*v1+m2*v2)/M
            p1.v=v*1
            p1.r=(r1+r2)/2

        # We're going to need a sequence of all of the pairs of particles when
        # we are detecting collisions. combinations generates pairs of indexes
        # into the self.particles list of Particles on the fly.    
        if self.n==2:
            pairs = combinations(range(self.n), 2)
            for i,j in pairs:
                if self.particles[i].overlaps(self.particles[j]):
                    change_velocities(self.particles[i], self.particles[j])
                    self.particles.pop()
                    self.particles[0].radius=0.05
                    self.particles[0].styles = {'edgecolor': 'r', 'fill':True,'color':'r'}
                    self.particles[0].name='H'
                    self.t=0
        self.n=len(self.particles)        
        #print(self.n)            
    def advance_animation(self, dt):
        """Advance the animation by dt, returning the updated Circles list."""
        # decay process  - generating the TAU's after 1.5 sec 
        
        if self.n==1 :
            self.t+=dt
            if self.t>=1.5:
                T1=self.particles[0]
                T2=self.particles[0]
                T1.r+=T1.radius
                T2.r-=T2.radius
                T1.radius=2/3*T1.radius
                T2.radius=2/3*T2.radius
                T1.v=np.array([0.01,0.01])
                T2.v=np.array([-0.01,-0.01])
                H0=self.particles[0]
                radii = np.array([H0.radius,T1.radius])
                P0=np.array([H0.r[0],H0.r[1],H0.v[0],H0.v[1]])
                P1=np.array([T1.r[0],T1.r[1],0.05,0.05])
                P2=np.array([T2.r[0],T2.r[1],-0.05,-0.05])
                #P1=np.array([0.02,0.02,0.1,0.1])
                #P2=np.array([0.02,0.98,0.1,-0.1])
                self.particles=[]
                p0 = Particle(x=P0[0], y=P0[1], vx=0.2*P0[2], vy=0*P0[3], radius=0*radii[0],styles = {'edgecolor': 'tab:orange', 'fill':True,'color':'tab:orange'})
                p1 = Particle(x=P1[0], y=P1[1], vx=P1[2], vy=P1[3], 
                              radius=radii[1],styles = {'edgecolor': 'y', 'fill':True,'color':'y'},Name='T')
                p2 = Particle(x=P2[0], y=P2[1], vx=P2[2], vy=P2[3], 
                              radius=radii[1],styles = {'edgecolor': 'y', 'fill':True,'color':'y'},Name='T-')
                self.particles.append(p0)
                self.particles.append(p1)
                self.particles.append(p2)
                self.n=len(self.particles) 
                self.t=0
        for i, p in enumerate(self.particles):
            p.advance(dt)
            self.circles[i].center = p.r
        if self.n==2:    
            self.handle_collisions()      
        if self.n>=3 :
            self.t+=dt
            

        return self.circles

    def advance(self, dt):
        """Advance the animation by dt."""
        for i, p in enumerate(self.particles):
            p.advance(dt)
        if self.n==2:    
            self.handle_collisions()
    

    def init(self):
        """Initialize the Matplotlib animation."""
        #print(self.t)
        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax))
        return self.circles

    def animate(self,i):
        """The function passed to Matplotlib's FuncAnimation routine."""
        #while (self.t<10 and self.n<=3):
        if (self.n>=3 and self.t>10):
            end()
        else:
            self.advance_animation(0.01)
            
            return self.circles
            

    def do_animation(self, save=True):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        fig, self.ax = plt.subplots(figsize=(10,7))
        plt.draw()
        #fig = plt.figure()
        self.ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
        line, = self.ax.plot([], [], lw=2)
        self.ax.xaxis.set_ticks([])
        self.ax.yaxis.set_ticks([])
        Mh = 125.35 * pow(10,9)
        Mtau = 1776.86 * pow(10,6)
        Q = Mh - (2*Mtau)
        Ptau = (np.sqrt(pow(((Q + (2*Mtau)) / 2),2) - pow(Mh / pow(10,9), 2)))
        C=299792458
        E = np.sqrt(pow(Ptau, 2) + pow(Mtau,2))
        KEtau = E - Mtau
        Vtau = np.sqrt(1 - pow(Mtau/(KEtau + Mtau),2))
        print(Vtau)
        #Velocity=np.sqrt(1+pow((Mtau)/(KE+Mtau)),2)
        textstr = '\n' " ".join((r'The Mass of Higgs Boson is $%.2f$ GeV' % (Mh/1e9, ),
                                 r'The Mass of Tauon is $%.2f$MeV' % (Mtau/1e6, ),
                                 r'The Q Value of the  interaction is $%.2f$ GeV' % (Q/1e9, ),
                                 r'The momentum of one of the tauons is $%.2f$GeV/c' % (Ptau/1e9, ),
                                 r'The Velocity of Tau is $%.6f$c' % (Vtau, ),
                                 r'The kinetic energy of Tau is $%.2f$ GeV'% (KEtau/1e9,),
                                 r'The energy of Tau is $%.02f$GeV/c^2'%(E/1e9,))) 
        
        textstr2 ='\n' " ".join((r'Legend',r'',
                               r'        -> H',r'', 
                               r'        -> P',r'',
                               r'        -> T'))
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)

        # place a text box in upper left in axes coords
        self.ax.text(0.05, 0.95, textstr, transform=self.ax.transAxes, fontsize=7,
            verticalalignment='top', bbox=props)
        self.ax.text(0.85, 0.25, textstr2, transform=self.ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
        circleP = Circle((0.87,0.12), 0.015   ,edgecolor= 'g', fill=True ,color='g')
        circleH = Circle((0.87,0.18), 0.02  ,edgecolor= 'r', fill=True ,color='r')
        circleT = Circle((0.87,0.06), 0.012  ,edgecolor= 'y', fill=True ,color='y')
        self.ax.add_patch(circleP)
        self.ax.add_patch(circleH)
        self.ax.add_patch(circleT)
        
        #plt.draw()
        
        
        
        anim = animation.FuncAnimation(fig, self.animate, init_func=self.init,
                               frames=10, interval=1, blit=True)
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=30, bitrate=1800)
            anim.save('collision.mp4', writer=writer)
            
        else:
            plt.show()


# In[5]:

plt.ion() 
get_ipython().run_line_magic('matplotlib', '')
if __name__ == '__main__':
    nparticles = 2
    P1=np.array([0.02,0.5,0.1,0])
    P2=np.array([0.98,0.5,-0.1,0])
    radii = np.array([0.03,0.03])
    particle1 = Particle(x=P1[0], y=P1[1], vx=P1[2], vy=P1[3], radius=radii[0],styles = {'edgecolor': 'g', 'fill':True,'color':'g'},Name='P')
    particle2 = Particle(x=P2[0], y=P2[1], vx=P2[2], vy=P2[3], radius=radii[1],styles = {'edgecolor': 'g', 'fill':True,'color':'g'},Name='P')
    P= [particle1, particle2]
    #radii = np.array([0.05,0.3])
    styles = {'edgecolor': 'C0', 'linewidth': 2, 'fill': None}
    sim = Simulation(nparticles,P, radii, styles)
    sim.do_animation(save=True)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




