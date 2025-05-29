import numpy as np
from math import *
import matplotlib.pyplot as plt

#############################
# Variable Declarations     #
#############################

Fdt = 16.72
Nmc = 1000
m = 0.44
r = 0.11
rho = 1.2
g = 9.81
A = (r**2) * np.pi
k = (rho * A) / (2 * m)
Cx = 0.47
t_0 = 0.0
t_end = 5.0
h = 0.001
distance = 35

#############################
# Function Definitions      #
#############################

def Impact():
    '''
    Randomly generates impact parameters:
        Yc = y-coordinate of impact point as a percentage of radius
        Zc = z-coordinate of impact point as a percentage of radius
        alphaHit = shot angle in degrees (x-y plane)
        betaHit = shot angle in degrees (x-y-r plane)
    '''
    randomY = 0.1 * np.random.randn(1, 1)
    randomZ = 0.1 * np.random.randn(1, 1)
    randomalpha = 0.1 * np.random.randn(1, 1)
    randombeta = 0.1 * np.random.randn(1, 1)

    Yc = randomY[0][0]
    Zc = randomZ[0][0]
    alphaHit = randomalpha[0][0]
    betaHit = randombeta[0][0]

    return Yc, Zc, alphaHit, betaHit


def rotation(Ypercent, Zpercent, r, m, Hitalpha, Hitbeta, Fdt):
    '''
    Calculates the velocity components and angular rotation
    based on the impact point on the ball and impact angles.
    '''
    b = r * Ypercent
    c = r * Zpercent
    a = -np.sqrt(r**2 - c**2 - b**2)

    sinA = np.sin((Hitalpha * np.pi) / 180)
    sinB = np.sin((Hitbeta * np.pi) / 180)
    cosA = np.cos((Hitalpha * np.pi) / 180)
    cosB = np.cos((Hitbeta * np.pi) / 180)

    Fx = Fdt * cosB * cosA
    Fy = Fdt * cosB * sinA
    Fz = Fdt * sinB

    Wx = 5 / (2 * m * (r**2)) * (b * Fz - c * Fy)
    Wy = 5 / (2 * m * (r**2)) * (-a * Fz + c * Fx)
    Wz = 5 / (2 * m * (r**2)) * (a * Fy - b * Fx)

    vx = Fx / m
    vy = Fy / m
    vz = Fz / m

    alpha = (atan(vy / vx) * 180) / np.pi
    beta = (atan(vz / (np.sqrt(vx**2 + vy**2))) * 180) / np.pi

    return vx, vy, vz, alpha, beta, Wx, Wy, Wz


def runge_kutta_4th_order(f, y0, t_0, t_end, h, v, C_MC, C_MB, sinA, cosA, cosB, distance):
    '''
    Solves the system using 4th-order Runge-Kutta method.
    '''
    num_steps = int((t_end - t_0) / h) + 1
    t_values = np.linspace(t_0, t_end, num_steps)
    y_values = np.zeros((num_steps, len(y0)))
    y_values[0] = y0

    for i in range(1, num_steps):
        k1 = h * f(t_values[i - 1], y_values[i - 1], v, C_MC, C_MB, sinA, cosA, cosB)
        k2 = h * f(t_values[i - 1] + 0.5 * h, y_values[i - 1] + 0.5 * k1, v, C_MC, C_MB, sinA, cosA, cosB)
        k3 = h * f(t_values[i - 1] + 0.5 * h, y_values[i - 1] + 0.5 * k2, v, C_MC, C_MB, sinA, cosA, cosB)
        k4 = h * f(t_values[i - 1] + h, y_values[i - 1] + k3, v, C_MC, C_MB, sinA, cosA, cosB)

        y_values[i] = y_values[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Check if the ball is out of bounds
        if (y_values[i, 0] < 0) or (y_values[i, 2] < -10) or (y_values[i, 2] > 10) or y_values[i, 0] > distance:
            return y_values[:i+1]

        # Check if the ball hits the ground
        if y_values[i, 4] < 0:
            return y_values[:i+1]

        # Check if the ball scores a goal
        if y_values[i, 0] == distance and (y_values[i, 4] > 2.44) and (-3.66 < y_values[i, 2] < 3.66):
            return y_values[:i+1]
        
        # Check if it hits the wall
        if y_values[i, 0] == 9.1 and ((y_values[i, 4] < 2.00) and (-1.5 < y_values[i, 2] < 0)):
            return y_values[:i+1]

    return y_values


def f(t, Y, v, C_MC, C_MB, sinA, cosA, cosB):
    '''
    Defines the differential equations of motion.
    '''
    x, x1 = Y[0], Y[1]
    y, y1 = Y[2], Y[3]
    z, z1 = Y[4], Y[5]

    x_dot = x1
    y_dot = y1
    z_dot = z1

    x1_dot = -k * v * (Cx * x1 + C_MC * v * sinA - C_MB * z1 * cosA)
    y1_dot = -k * v * (Cx * y1 - C_MC * v * cosA - C_MB * z1 * sinA)
    z1_dot = -k * v * (Cx * z1 + C_MB * v * cosB) - g

    return np.array([x_dot, x1_dot, y_dot, y1_dot, z_dot, z1_dot])


def ballon_affiche(points, sol):
    '''
    Displays impact points on the ball surface.
    '''
    theta = np.linspace(0, 2*np.pi, 100)
    x1 = r * np.cos(theta)
    x2 = r * np.sin(theta)
    fig, ax = plt.subplots(1)
    ax.plot(x1, x2)
    ax.set_aspect(1)

    plt.xlim(-0.12, 0.12)
    plt.ylim(-0.12, 0.12)
    plt.grid(linestyle='--')

    for i in points:
        x, y = i
        if not ((x/r, y/r) in sol):
            plt.scatter(x, y, c='black', marker='x', s=20)
        else:
            plt.scatter(x, y, c='green', marker='o', s=30)

    plt.title('Ball Impact Points', fontsize=8)


def cherche_max(lst):
    '''
    Finds the maximum Y and Z values from a list of tuples.
    '''
    if not lst:
        return [], []

    max_Y = lst[0][0]
    max_Z = lst[0][1]

    for Y, Z in lst:
        if Y > max_Y:
            max_Y = Y
        if Z > max_Z:
            max_Z = Z

    return [max_Y], [max_Z]

#######################
# 3D Field Definition #
#######################

fig, cx = plt.subplots()
cx = plt.axes(projection='3d')

# Define the crossbar of the goal
goal_z_p = np.linspace(2.44, 2.44, 50)  # z-coordinates of the crossbar (fixed height at 2.44m)
goal_x_p = np.linspace(distance, distance, 50)  # x-coordinates of the crossbar (fixed distance from the kicker)
goal_y_p = np.linspace(-3.66, 3.66, 50)  # y-coordinates spanning the width of the goal (-3.66m to 3.66m)

# Plot the crossbar
cx.plot3D(goal_x_p, goal_y_p, goal_z_p, '-b', label='goal')

# Define the right goalpost
goal_z_pr = np.linspace(0, 2.44, 50)  # height from ground to crossbar
goal_x_pr = np.linspace(distance, distance, 50)  # fixed x-position at goal distance
goal_y_pr = np.linspace(-3.66, -3.66, 50)  # y-position of the right post

# Plot the right goalpost
cx.plot3D(goal_x_pr, goal_y_pr, goal_z_pr, '-b')

# Define the left goalpost
goal_z_pl = np.linspace(0, 2.44, 50)  # height from ground to crossbar
goal_x_pl = np.linspace(distance, distance, 50)  # fixed x-position at goal distance
goal_y_pl = np.linspace(3.66, 3.66, 50)  # y-position of the left post

# Plot the left goalpost
cx.plot3D(goal_x_pl, goal_y_pl, goal_z_pl, '-b')

# Define the top of the defensive wall
wall_z_top = np.linspace(2.0, 2.0, 50)  # height of the top of the wall (2 meters)
wall_x_top = np.linspace(9.1, 9.1, 50)  # fixed x-position of the wall
wall_y_top = np.linspace(0, -1.5, 50)  # span of the wall in y-direction

# Plot the top of the wall
cx.plot3D(wall_x_top, wall_y_top, wall_z_top, '-g', label='wall of players')

# Define the right side of the wall
wall_z_r = np.linspace(2.0, 0.0, 50)  # vertical line from 2m to ground
wall_x_r = np.linspace(9.1, 9.1, 50)  # fixed x-position
wall_y_r = np.linspace(-1.5, -1.5, 50)  # y-position of the right side

# Plot the right side of the wall
cx.plot3D(wall_x_r, wall_y_r, wall_z_r, '-g')

# Define the left side of the wall
wall_z_l = np.linspace(2.0, 0.0, 50)  # vertical line from 2m to ground
wall_x_l = np.linspace(9.1, 9.1, 50)  # fixed x-position
wall_y_l = np.linspace(0.0, 0.0, 50)  # y-position of the left side

# Plot the left side of the wall
cx.plot3D(wall_x_l, wall_y_l, wall_z_l, '-g')

# Label axes
cx.set_xlabel('x (m)')
cx.set_ylabel('y (m)')
cx.set_zlabel('z (m)')

# Set axis limits
cx.set_xlim(0, distance)
cx.set_ylim(-10, 10)
cx.set_zlim(0, 15)



####################
# MAIN LOOP        #
####################

def main():
    sol=[]
    ang_sol=[]
    points=[]
    for i in range(Nmc):
        Yc, Zc, alphaHit, betaHit = Impact()
        vx, vy, vz, alpha, beta, Wx, Wy, Wz = rotation(Yc, Zc, 0.11, 0.41, alphaHit, betaHit, Fdt)
        v = np.sqrt(vx**2 + vy**2 + vz**2)
        angVy = Wy
        angVz = Wz
        angV = np.sqrt(Wx**2 + Wy**2 + Wz**2)
        C_MC = 1 / (2 + (v / (r * angVz)))
        C_MB = 1 / (2 + (v / (r * angVy)))

        sinA = np.sin((alpha * np.pi) / 180) 
        sinB = np.sin((beta * np.pi) / 180)
        cosA = np.cos((alpha * np.pi) / 180)
        cosB = np.cos((beta * np.pi) / 180)

        y0 = np.array([0, vx, 0, vy, 0, vz])
        y_values = runge_kutta_4th_order(f, y0, t_0,t_end , h,v,C_MC,C_MB,sinA,cosA,cosB,distance)
        x = y_values[:, 0]
        y = y_values[:, 2]
        z = y_values[:, 4]

        if y_values[-1, 0]>=distance and 2<y_values[-1, 4] < 2.44 and (-3.66 < y_values[-1, 2] <3.66):
            cx.plot3D(x, y, z, color='green')
            sol.append((Yc,Zc))
            ang_sol.append((alphaHit,betaHit))
                           
        # elif y_values[-1, 0] ==9.1 and ((y_values[-1, 4] < 2.00) and (-1.5 < y_values[-1, 2] < 0)):
        #     cx.plot3D(x, y, z, color='red')

        # else:
        #     cx.plot3D(x, y, z, color='black')

        points.append((r*Yc,r*Zc))
        
    ballon_affiche(points,sol)
    return sol,ang_sol

if __name__ == "__main__":
    sol, ang_sol = main()

    if sol:
        L_Y, L_Z = cherche_max(sol)
        L_A, L_B = cherche_max(ang_sol)
        print("Max impact coordinates on the ball:", L_Y, L_Z)
        print("Max shot angles:", L_A, L_B)
    else:
        print("No trajectory resulted in a goal.")

plt.show()















