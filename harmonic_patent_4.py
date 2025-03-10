import numpy as np
import matplotlib.pyplot as plt

# ----------------
# INPUT PARAMETERS
# ----------------
m = 1           # module
Nr = 32        # Number of teeth (ring gear)

# ----------------
# CALCULATED PARAMETERS
# ----------------
n = 2           # Number of lobes on wave generator (strain inducer) - Fixed at 2
Ns = Nr - n     # Number of teeth for strain gear
ratio = Ns / n  # Reduction Ratio
P = 1 / m       # Diametrical Pitch
p = np.pi / P   # Circular Pitch
d = n / P       # deflection
a = 0.139 * n * p   # addendum
b = 0.179 * n * p   # dedendum
c = 0.04 * n * p    # clearance
t = 11/(8*P)    # tooth thickness at pitch line
root_diameter_delta = d * (((2/3)+(9/16))/2)    # Offsets for root diameter
tip_diameter_delta  = d * ((7/16)-(1/32))       # Offsets for root diameter

# ----------------
# For Mechanical Analysis
# ----------------
Ti = 0.15       # Input torque, Nm
E = 3600        # Elastic Modulus, MPa
S = 60          # Tensile Strength, MPa
l = 10          # width of ring
thick = 1       # thickness of ring
f1 = 0.05       # Coefficient of sliding friction, gear teeth
f2 = 0.0015     # Coefficient of friction, strain inducer (wave generator)
F = l           # Face width of gear

def circle_coords(x_center, y_center, diameter):
    radius = diameter / 2
    # Generate points along the circle
    theta = np.linspace(0, 2 * np.pi, 300)
    x = x_center + radius * np.cos(theta)
    y = y_center + radius * np.sin(theta)
    return x, y

def find_intersection(m, b, x_c, y_c, diameter):
    R = diameter / 2
    # Coefficients of the quadratic equation A*x^2 + B*x + C = 0
    A = 1 + m**2
    B = -2*x_c + 2*m*(b - y_c)
    C = x_c**2 + (b - y_c)**2 - R**2
    # Compute the discriminant
    discriminant = B**2 - 4*A*C
    # Check the number of intersections based on the discriminant
    if discriminant < 0:
        #print("No intersection between the line and the circle.")
        x_intersections = []
        y_intersections = []
    elif np.isclose(discriminant, 0):
        #print("The line is tangent to the circle (one intersection).")
        x_int = -B / (2*A)
        x_intersections = [x_int]
        y_intersections = [m * x_int + b]
    else:
        #print("The line intersects the circle at two points.")
        sqrt_disc = np.sqrt(discriminant)
        x_int1 = (-B + sqrt_disc) / (2*A)
        x_int2 = (-B - sqrt_disc) / (2*A)
        x_intersections = [x_int1, x_int2]
        y_intersections = [m * x_int1 + b, m * x_int2 + b]
    
    for x_int, y_int in zip(x_intersections, y_intersections):
        #print(f"({x_int}, {y_int})")
        if x_int > 0:
            return x_int, y_int

def rotate_points(points, angle):
    # Function to rotate a set of points by a given angle (in radians)
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, -s], [s, c]])
    return points.dot(R.T)
    
def plot_points(ax, pitch_angle, points_tooth, N, type):
  # Plot each tooth by rotating the base tooth shape
    for i in range(N):
        angle_rotation = i * pitch_angle
        tooth_points = rotate_points(points_tooth, angle_rotation)
        # Close the polygon by appending the first point at the end
        #tooth_points = np.vstack([tooth_points, tooth_points[0]])
        ax.plot(tooth_points[:, 0], tooth_points[:, 1], type, lw=1)

def main():
    # Ring Gear (Internal gear) parameters 
    Dr = Nr * m                          # Pitch Diameter for ring gear
    r_ringGear = Dr / 2                    # Radius of Pitch Circle
    theta_r = np.arctan(1.091/n)           # press angle ring gear (in rad)
    theta_r_deg = np.rad2deg(theta_r)      # press angle in degrees
    Root_Diameter_RingGear = Dr + 2 * root_diameter_delta
    Tip_Diameter_RingGear = Dr - 2 * tip_diameter_delta

    # Strain Gear (External gear) parameters
    Deffective = Dr - d                # Effective Pitch Diameter
    r_effective = Deffective / 2         # Effective Pitch Radius 
    Drelaxed = Deffective - 0.416*d      # Relaxed Pitch Diameter
    r_relaxed = Drelaxed / 2
    theta_s = theta_r + np.arctan(0.458*d*n / r_effective)
    theta_s_deg = np.rad2deg(theta_s)
    Root_Diameter_strainGear = Drelaxed - 2 * root_diameter_delta
    Tip_Diameter_strainGear  = Drelaxed + 2 * tip_diameter_delta

    # Output Gear (Internal Gear)
    Do = Dr             # Match Pitch Diameter of Ring Gear
    r_output = Do / 2
    No = Ns             # Match Teeth of Strain Gear
    theta_o = theta_r   # Match Press Angle of Ring Gear
    Root_Diameter_Output = Do + 2 * root_diameter_delta
    Tip_Diameter_Output = Do - 2 * tip_diameter_delta 

    # Determine the angular half-width of the tooth for Strain / Ring / Output Gear
    half_tooth_angle_strainGear = t / (2 * r_relaxed)  # in radians
    half_tooth_angle_RingGear = t / (2 * r_ringGear)  # in radians
    half_tooth_angle_output = t / (2 * r_output)

    # The pitch angle between consecutive teeth on the strain / Ring gear
    pitch_angle_StrainGear = 2* np.pi / Ns
    half_pitch_angle_StrainGear = pitch_angle_StrainGear / 2

    pitch_angle_RingGear = 2* np.pi / Nr

    pitch_angle_output = 2 * np.pi / No

    # Point on the Pitch Diameter
    x0_StrainGear = r_relaxed * np.cos(half_tooth_angle_strainGear)
    y0_StrainGear = r_relaxed * np.sin(half_tooth_angle_strainGear)

    x0_RingGear = r_ringGear * np.cos(half_tooth_angle_RingGear)
    y0_RingGear = r_ringGear * np.sin(half_tooth_angle_RingGear)

    x0_output = r_output * np.cos(half_tooth_angle_output)
    y0_output = r_output * np.sin(half_tooth_angle_output) 

    # Absolut angle of edge line
    angle_edge_StrainGear = half_pitch_angle_StrainGear - theta_s
    slope_StrainGear = np.tan(angle_edge_StrainGear)

    angle_edge_RingGear = theta_r
    slope_RingGear = np.tan(angle_edge_RingGear)

    angle_edge_output = theta_o
    slope_output = np.tan(angle_edge_output)

    # Find line, y = mx + b
    b_StrainGear = y0_StrainGear - slope_StrainGear * x0_StrainGear
    b_RingGear = y0_RingGear - slope_RingGear * x0_RingGear
    b_output = y0_output - slope_output * x0_output

    # Coordinates of root circles
    x_root_StrainGear, y_root_StrainGear = circle_coords(0,0,Root_Diameter_strainGear)

    x_root_RingGear, y_root_RingGear = circle_coords(0,0,Root_Diameter_RingGear)

    x_root_output, y_root_output = circle_coords(0,0,Root_Diameter_Output)

    # Find intersection between edge line and root diameter / tip diameter
    # These will form 4 vertices of tooth
    x_root_int_StrainGear, y_root_int_StrainGear = find_intersection(slope_StrainGear, b_StrainGear, 0, 0, Root_Diameter_strainGear)
    x_tip_int_StrainGear, y_tip_int_StrainGear = find_intersection(slope_StrainGear, b_StrainGear, 0, 0, Tip_Diameter_strainGear)

    x_root_int_RingGear, y_root_int_RingGear = find_intersection(slope_RingGear, b_RingGear, 0, 0, Root_Diameter_RingGear)
    x_tip_int_RingGear, y_tip_int_RingGear = find_intersection(slope_RingGear, b_RingGear, 0, 0, Tip_Diameter_RingGear)

    x_root_int_output, y_root_int_output = find_intersection(slope_output, b_output, 0, 0, Root_Diameter_Output)
    x_tip_int_output, y_tip_int_output = find_intersection(slope_output, b_output, 0, 0, Tip_Diameter_Output)

    ################################################
    # Find Tooth Point - trapezoidal shape, 4 vertices ( 4 points)
    ################################################
    points_tooth_StrainGear = np.array([
        [x_root_int_StrainGear, y_root_int_StrainGear],  
        [x_tip_int_StrainGear, y_tip_int_StrainGear],   
        [x_tip_int_StrainGear, -y_tip_int_StrainGear],    
        [x_root_int_StrainGear, -y_root_int_StrainGear]     
    ])

    points_tooth_RingGear = np.array([
        [x_root_int_RingGear, y_root_int_RingGear],  
        [x_tip_int_RingGear, y_tip_int_RingGear],   
        [x_tip_int_RingGear, -y_tip_int_RingGear],    
        [x_root_int_RingGear, -y_root_int_RingGear]     
    ])

    points_tooth_output = np.array([
        [x_root_int_output, y_root_int_output],  
        [x_tip_int_output, y_tip_int_output],   
        [x_tip_int_output, -y_tip_int_output],    
        [x_root_int_output, -y_root_int_output]     
    ])

    ################################################
    # Mechanical Analysis
    ################################################
    sigma_max = 1.03 * thick * d * E / r_relaxed**2     # max tensile stress in deflected ring at outside crest of wave
    sigma_1 = 0.59 * thick * d * E / r_relaxed**2       # tensile stress in deflected ring at wave base
    W1 = 0.56 * d * l * thick**3 * E / r_relaxed**3     # radial load required to deflec ring
    efficiency = (1/ratio) * ( (1-f1*np.tan(theta_s) ) / (np.tan(theta_s)+f1) ) * (( 1 - (f2 * 0.458 * d * n / r_relaxed)) / (((0.458*d*n)/r_relaxed) + f2 ) ) 
    To = Ti*ratio*efficiency
    sigma_s = 1.6 * To / (r_relaxed**2 * F)
    print(f'Ratio: {ratio}, Input Torque: {Ti} Nm, Output Torque: {round(To,1)} Nm, Efficiency: {round(efficiency*100,1)} %, Max stress: {round(sigma_max,1)} MPa')
    ################################################
    ################################################

    # Plot
    # Create plot
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot Strain Gear
    plot_points(ax, pitch_angle_StrainGear, points_tooth_StrainGear, Ns, 'r-')            # Teeth
    ax.plot(x_root_StrainGear, y_root_StrainGear, 'r-', lw=1, label="Strain Gear")    # Root Diameter
    
    # Plot Ring Gear
    ax.plot(x_root_RingGear, y_root_RingGear, 'b-', lw=1, label="Ring Gear")    # Root Diameter
    plot_points(ax, pitch_angle_RingGear, points_tooth_RingGear, Nr, 'b-')            # Teeth

    # Plot Output Gear
    ax.plot(x_root_output, y_root_output, 'g-', lw=1, label="Output Gear")    # Root Diameter
    plot_points(ax, pitch_angle_output, points_tooth_output, No, 'g-')            # Teeth

    #plt.savefig("strain_gear.svg")  # or "gear_profile.pdf"
    ax.legend()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Strain Gear (External Gear) Profile")
    ax.axis('equal')
    ax.grid(True)
    plt.show()

if __name__ == "__main__":
    main()


