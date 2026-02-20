import numpy as np
import matplotlib.pyplot as plt
import math

def reverseArray(array):
    newArray = []
    numElements = len(array)
    #print("number of array elements: " + str(numElements))
    #appends the reverse of the original array to an empty array
    for i in range((numElements-1),-1,-1):
        newArray.append(array[i])
        #print(newArray)
    return newArray

def mirrorArray(array):
    newArray = []
    numElements = len(array)
    #print("number of array elements: " + str(numElements))
    #appends the reverse of the original array to an empty array
    for i in range((numElements-1),-1,-1):
        newArray.append(array[i])
        #print(newArray)
    normal_Array = array
    for values in normal_Array:
        newArray.append(values)
    return newArray

def reflectArray(array):
    newArray = []
    numElements = len(array)
    #print("number of array elements: " + str(numElements))
    #appends the reverse of the original array to an empty array
    for i in range((numElements-1),0,-1):
        newArray.append(-array[i])
        #print(newArray) 
    normal_Array = array
    for values in normal_Array:
        #print(values)
        newArray.append(values)
    return newArray

def bendingStiffness(E_array,r_inner,r_outer,numLayers):
    t_p = (r_outer-r_inner)/numLayers
    t_p = 0.00625 

    Radii = [0, r_inner] #r_0 and r_1
    # Add 8 composite layers (r2 to r9)
    for i in range(numLayers):
        Radii.append(Radii[i+1] + t_p)

    H33ci = np.zeros(numLayers+1)

    for i in range(numLayers+1):
        H33ci[i] = E_array[i] * (math.pi / 4) * (Radii[i+1]**4 - Radii[i]**4)
    
    #total bending stiffness
    H33c_total = np.sum(H33ci)

    return (H33c_total)

def Axial_Stress_From_Bending(Mom_Array,Ex1x1_array,H33_c,r_inner,r_outer,numLayers):
    print("Ex1x1 array: " +str(Ex1x1_array))
    Ex1x1_array=mirrorArray(Ex1x1_array)
    print("Ex1x1 array (mirrored): " + str(Ex1x1_array))
    t_p = (r_outer-r_inner)/numLayers
    #t_p = 0.00625 
    print("t_p = "+ str(t_p))

    r_array = np.linspace(-r_outer,r_outer,len(Mom_Array))

    Radii = [0, r_inner] #r_0 and r_1
    # Add 8 composite layers (r2 to r9)
    for i in range(numLayers):
        Radii.append(Radii[i+1] + t_p)
    print("Radii: "+str(Radii))
    Radii = reflectArray(Radii)
    print("Radii (reflected): " + str(Radii))
    E_radial_array=[0]
    #Ex1x1_array.append(Ex1x1_array[len(Ex1x1_array)-1])

    for j in range(len(Radii)-1):
        for i in range(len(r_array)-1):
            #makes sure that we're in the correct range
            if (r_array[i] < Radii[j+1] and r_array[i] >= Radii[j]):
                #append an array to have the needed E value for that layer
                print(i)
                E_radial_array.append(Ex1x1_array[j])
            #print("r: "+str(r_array[i])+ "  E: "+str(E_radial_array[i]))
            #print("E_value"+" index: "+str(i))
    E_radial_array.append(Ex1x1_array[len(Ex1x1_array)-1])

    omega_1 = np.zeros(len(r_array))
    print("length of r_array: "+str(len(r_array))+"     length of E_radial_array: "+str(len(E_radial_array)))
    for i in range(len(r_array)-1):
        #print("r: "+str(r_array[i])+ "      E: "+str(E_radial_array[i])+ "      index: " + str(i))
        omega_1_section = -E_radial_array[i] * r_array[i]*Mom_Array[len(Mom_Array)-1]/H33_c
        omega_1[i] = omega_1_section
    
    #adds mirror 
    return (omega_1,r_array)




def symmetric_layup(layupArray):
    #I don't think this differentiates between the layups that have an
    #odd or even number of elements
    numElements = len(layupArray)
    #check if there's an even number of layers
    if ((numElements % 2)== 0):
        laminate = layupArray
        
        #add the reverse of the array to the end so it shows the whole layup
        layupReverse = reverseArray(layupArray)
        print("Layup array: " + str(layupArray))
        print("Reverse array: " + str(layupReverse))
        #need to go step-by-step instead of appending the whole array
        for i in range(numElements):
            laminate.append(layupReverse[i])
        print("Full laminate array: " + str(laminate))
    return(laminate)

def degToRad(degrees):
    radians = degrees * math.pi/180
    return radians
    
def orientation_Transform(degrees):
    #T-matrix found on page 4 of calculation_bendingstiffness_laminate_EAE135.pdf
    Beta = degToRad(degrees)
    T_1_1 = (math.cos(Beta))**2
    T_1_2 = (math.sin(Beta))**2
    T_1_3 = 2*math.sin(Beta)*math.cos(Beta)
    T_2_1 = T_1_2
    T_2_2 = T_1_1
    T_2_3 = -T_1_3
    T_3_1 = (1/2)*T_2_3
    T_3_2 = -T_3_1
    T_3_3 = T_1_1 -T_1_2

    T = np.array([[T_1_1,T_1_2,T_1_3],[T_2_1,T_2_2,T_2_3],[T_3_1,T_3_2,T_3_3]])

    #print("Transformation matrix for "+str(degrees)+" degrees: \n" + str(T))
    return T

    

class laminate:
    #initialize with a material class and layup Array
    def __init__(self,material,layupArray):
        self.material = material
        self.layupArray = layupArray
        self.fullLayupArray = symmetric_layup(layupArray)
    def Q_Matrix(self):
        denominator = (1-self.material.v_12*self.material.v_21)
        Q_1_1 = self.material.E1/denominator
        Q_1_2 = (self.material.v_21*self.material.E1)/denominator
        Q_1_3 = 0
        Q_2_1 = (self.material.v_12*self.material.E2)/denominator
        Q_2_2 = (self.material.E2)/denominator
        Q_2_3 = 0
        Q_3_1 = 0
        Q_3_2 = 0
        Q_3_3 = self.material.G_12

        Q = np.array([[Q_1_1,Q_1_2,Q_1_3],[Q_2_1,Q_2_2,Q_2_3],[Q_3_1,Q_3_2,Q_3_3]])
        #print("Q matrix (GPa): \n"+str(Q))
        return Q
    
    def Q_bar(self,degrees):
        #print("Start of Q_bar calculations")
        R = np.array([[1,0,0],[0,1,0],[0,0,2]])
        #print("R_matrix: \n"+str(R))
        T=orientation_Transform(degrees)
        Q=self.Q_Matrix()

        T_inv = np.linalg.inv(T)
        R_inv = np.linalg.inv(R)

        #Q_bar = T^-1 * Q * R * T * R^-1
        Q_bar = T_inv @ Q @ R @ T @ R_inv
        #print("T = \n" + str(T))
        #print("R = \n" + str(R))
        #print("Q_bar = \n" + str(Q_bar))
        return Q_bar
    
    def Q_bar_array(self):
        Q_bars = []
        #print("Q_bars initial: " + str(Q_bars))
        for i in range(len(layup)):
            #takes the sheet orientations from layup array
            #and calculates the Q_bar for each sheet sequentially
            #stores these matrices in Q_bars
            
            
            #print("Q_bar "+ str(i) + " :\n"+str(self.Q_bar(layup[i])))
            Q_bars.append(self.Q_bar(layup[i]))
            #print("printing Q_bars")
            #for Q_bar in Q_bars:
                #print("\n" + str(Q_bar))
        
        return Q_bars
      

class Material:
    def __init__(self,rho,E1,E2,v_12,v_23,
                 G_12,G_23,
                 omega_axial_T,omega_axial_C,
                 omega_transverse_T,
                 omega_transverse_C):
        
        #Initialize Material with given properties
            self.rho = rho
            self.E1 = E1
            self.E2 = E2
            self.v_12 = v_12
            self.v_23 = v_23
            self.v_21 = v_12*E2/E1
            self.G_12 = G_12
            self.G_23 = G_23
            self.omega_axial_T = omega_axial_T
            self.omega_axial_C = omega_axial_C
            self.omega_transverse_T = omega_transverse_T
            self.omega_transverse_C = omega_transverse_C
    
    

#Conversion Factors
ksi_to_MPa = 6.89476

#starting t+5 seconds
T_max_1 = 726 #kN
D_outer = 128 / 100 #cm -> m
D_inner = 118 / 100 #cm -> m

t_maxQ = 36 #seconds

#Table 1 Input Data
AOA = 20 #degree
l_rocket = 5*D_outer #m
t_layer = (1/8) * (50) / 1000 #m
alpha = degToRad(AOA)
l_cp_to_cg = D_inner #m
F_axial_O = 600 #kN
W_rocket = 8500 * 9.8 / 1000 #kg to #kN
Lift = 1400 #kN
safety_factor = 1.25

#Problem Statement
x_1 = np.linspace(0,l_rocket,1000)

#Material Properties
#AS4/epoxy
c_rho = 1.52 #g/cm^3
c_E1 = 148 #GPa
c_E2 = 10.50 #GPa
c_v_12 = 0.30
c_v_23 = 0.59
c_G_12 = 5.61 #GPa
c_G_23 = 3.17 #GPa
c_omega_axial_T = 2137 #MPa
c_omega_axial_C = -184 * ksi_to_MPa #MPa
c_omega_transverse_T = 53.4 #MPa
c_omega_transverse_C = -24.4 * ksi_to_MPa #MPa

E_HTPB = 9.036 * 10**(-3) # GPa

carbon_epoxy=Material(c_rho,c_E1,c_E2,c_v_12,c_v_23,c_G_12,c_G_23,
                c_omega_axial_T,c_omega_axial_C,
                c_omega_transverse_T,c_omega_transverse_C)

layup = [0,45,-45,90]

#symmetric_layup(layup)

black_Aluminum=laminate(carbon_epoxy,layup)
black_Aluminum.Q_Matrix()
#orientation_Transform(45)
black_Aluminum.Q_bar(45)

Q_bars=[]

Q_bars = black_Aluminum.Q_bar_array()
#print("Black Aluminum Q_bars: ")
S=[]

E_x1x1 = [E_HTPB] #GPa
#print("S matrix (should be inverse of Q): ")
for Q in Q_bars:
    S_bar = np.linalg.inv(Q)
    S.append(S_bar)
    print(S_bar)
    E_x1x1.append(1/(S_bar[0,0]))
#print("S_bar should be inverse of Q_bar")
print("Ex1x1: "+str(E_x1x1) + " GPa")

H_33c = bendingStiffness(E_x1x1,D_inner/2,D_outer/2,len(black_Aluminum.fullLayupArray)) #check units here
print("H33_c: "+ str(H_33c))


#Bending
#Need unit checks for the rest of this

#Integrate Moment Equation Once
print("L_rocket: " + str(l_rocket) + "\nL_rocket^2: " + str(l_rocket**2))
c_3 = -((1/2)*(Lift*math.cos(alpha)-W_rocket)*l_rocket**2 + W_rocket*D_outer*l_rocket)

#Integrate Moment Equation Twice
c_4= -((1/6)*(Lift*math.cos(alpha))*l_rocket**3 + (1/2)*W_rocket*D_outer*l_rocket**2 + c_3*l_rocket)

print("c_3 = " + str(c_3)+ "\nc_4 = "+str(c_4))

#plug into displacement
u2_x1 = (1/((10**3)*H_33c))*((1/6)*(Lift*math.cos(alpha)*(x_1)**3)+(1/2)*(W_rocket)*(D_outer)*(x_1)**2+c_3*x_1+c_4) #micrometers


Mom_Array = [] #Moments
for x in x_1:
    if x < D_outer:
        Mom_near_end = Lift*math.cos(degToRad(AOA)) * x
        Mom_Array.append(Mom_near_end)
        #print("X = "+ str(x) + "      Moment: "+str(Mom_near_end))
    else:
        Mom_near_free = Lift*math.cos(degToRad(AOA)) * x - W_rocket*(x-D_outer)
        Mom_Array.append(Mom_near_free)
        #print("X = "+ str(x) + "      Moment: "+str(Mom_near_free))
    



#Axial bending stress

(omega_1, r_array) = Axial_Stress_From_Bending(Mom_Array,E_x1x1,H_33c,D_inner/2,D_outer/2,8)

fig, ax1 = plt.subplots()

color = 'tab:blue'
#plt.figure(figsize=(16,9))
#ax1.title("Vertical Displacement of Beam")
ax1.set_xlabel("x_beam (m)")
ax1.set_ylabel("u_2 (um)")
ax1.plot(x_1,u2_x1,marker='o',markersize=10,color='blue')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:red'
ax2.set_ylabel("Moment (I dont know the units are)")
ax2.plot(x_1,Mom_Array,marker='o',markersize=5,color='red')
ax2.tick_params(axis='y', labelcolor=color)
plt.grid()
plt.show()

plt.figure(figsize=(16,9))
plt.title("Axial Stress versus radius at Critical Cross Section")
plt.xlabel("Axial Stress")
plt.ylabel("radius")
plt.plot(omega_1,r_array)
plt.grid()
plt.show()

