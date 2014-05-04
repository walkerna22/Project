from __future__ import division
import numpy as np
from pylab import plot,xlabel,ylabel,show



n_steps = 1000

def init_arrays(n_points,start_time, end_time):
    """
    
    """
    time, delta_t = np.linspace(start_time, end_time, num=n_points, retstep=True)
    posx = np.zeros_like(time)
    velx = np.zeros_like(time)
    posy = np.zeros_like(time)
    vely = np.zeros_like(time)
    return time, delta_t, posx, velx, posy, vely
    #print end_time
def init_conditions(n_points):
    time, delta_t, posx, velx, posy, vely = init_arrays(n_points, 0.0, 5000)
    posx[0] = 4.0e12
    velx[0] = 0.0
    posy[0] = 0.0
    vely[0] = 500.0
      

    #print "init_conditions setting posy[0] to "+str(posx[0])
    return time, delta_t, posx, velx, posy, vely

def integrate_pendulum(n_steps, integration_method):
    
    time, delta_t, posx, velx, posy, vely = init_conditions(n_steps)
    time,pos_outx, pos_outy, vel_outx, vel_outy = integration_method(time,posx, posy, velx, vely)
    
    return time, pos_outx, vel_outx, pos_outy, vel_outy
#print g
def rk4(present_posx, present_velx, present_posy, present_vely, acceleration_func, delta_t):
    # Implement 4th Order RK
    pos1x = present_posx
    pos1y = present_posy
    vel1x = present_velx
    vel1y = present_vely
    accel1x, accel1y = acceleration_func(pos1x,pos1y)
    #print (present_velx, present_vely)
    #print (pos1x, pos1y)
    
    pos2x = present_posx + 0.5*vel1x*delta_t
    pos2y = present_posy + 0.5*vel1y*delta_t
    vel2x = present_velx + 0.5*accel1x*delta_t
    vel2y = present_vely + 0.5*accel1y*delta_t
    accel2x, accel2y = acceleration_func(pos2x,pos2y)
    
    
    pos3x = present_posx + 0.5*vel2x*delta_t
    pos3y = present_posy + 0.5*vel2y*delta_t
    vel3x = present_velx + 0.5*accel2x*delta_t
    vel3y = present_vely + 0.5*accel2y*delta_t
    accel3x,accel3y = acceleration_func(pos3x,pos3y)
    
    

    pos4x = present_posx + vel3x*delta_t
    pos4y = present_posy + vel3y*delta_t
    vel4x = present_velx + accel3x*delta_t
    vel4y = present_vely + accel3y*delta_t
    accel4x,accel4y = acceleration_func(pos4x,pos4y)
    
    
    
    new_posx = present_posx + (delta_t/6.0)*(vel1x + 2*vel2x + 2*vel3x + vel4x)
    new_posy = present_posy + (delta_t/6.0)*(vel1y + 2*vel2y + 2*vel3y + vel4y)
    new_velx = present_velx + (delta_t/6.0)*(accel1x + 2*accel2x + 2*accel3x + accel4x)
    new_vely = present_vely + (delta_t/6.0)*(accel1y + 2*accel2y + 2*accel3y + accel4y)

    

    return new_posx, new_velx, new_posy, new_vely

def pendulum_linear_runge_kutta4(time, pos_inputx, vel_inputx, pos_inputy, vel_inputy
                                 ):
   
    # Define a function for the acceleration of the pendulum (assuming
    # small angle linear approximation)
    def linear_accel(posx,posy,
            g = 6.67e-11, # N(m/kg)^2
            m = 1.99e30): #kilograms)
        return -(g*m*posx)/(np.sqrt(posx**2 + posy**2))**3, -(g*m*posy)/(np.sqrt(posx**2 + posy**2))**3
    
    # Copy the initial array to avoid clobbering it in case I need to re-run
    # this procedure
    posx = pos_inputx.copy()
    velx = vel_inputx.copy()
    posy = pos_inputy.copy()
    vely = vel_inputy.copy()
   

    # Compute position and velocities of pendulum using 2nd Order Runge-Kutta
    for i in np.arange(0, len(time)-1):
        delta_t = time[i+1] - time[i]
        posx[i+1], velx[i+1], posy[i+1], vely[i+1] = rk4(posx[i], velx[i], posy[i], vely[i], linear_accel, delta_t)
        
    print("y position:")
    print posy
    #print ("y position:")
    #print vely
    #print delta_t
    return time, posx, velx, posy, vely
    

