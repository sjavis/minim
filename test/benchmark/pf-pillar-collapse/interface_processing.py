import numpy as np
import matplotlib.pyplot as plt

nx = 100
ny = 100
nz = 100

def data_extract(filename):
    #loading the in the data obtained from the simulation
    
    
    data = np.loadtxt(filename).reshape((nx,ny,nz,3))
    #data = np.swapaxes(data, 0, 1) # Swap x-y axes for visualisation
    solid = data[:,:,:,0]
    liquid = data[:,:,:,1]
    gas = data[:,:,:,2]

    return solid,liquid,gas

def region_definition(solid, liquid, gas):
    #explicity define the regions using the following definition
    # solid = 1, gas = 0 and liquid = -1
    
    
    explicit_structure = np.zeros((nx,ny,nz))
    
    explicit_structure[solid > 0.5] = 1
    explicit_structure[np.logical_and(solid <= 0.5, liquid > gas)] = -1
    
    return explicit_structure
    
def solid_visualisation(solid, liquid, gas):
    #determines where the top of the solid surface is to visualise structure 
    
    
    solid_structure = np.zeros((nx,ny,nz))
    solid_structure[np.logical_and(solid > liquid, solid > gas)] = 1

    solid_offset1 = solid_structure[:,:-1,:]
    solid_offset2 = solid_structure[:,1:,:]

    solid_top = np.where(np.logical_and(solid_offset1 + solid_offset2 == 1, solid_offset2 == 0) == True)
    coords = np.array((solid_top[0],solid_top[2],solid_top[1])).T
    
    max_y_val_dict = {}
    
    for i in coords:
        x, z, y = i
        
        if (x,z) not in max_y_val_dict or y > max_y_val_dict[(x,z)]:
            max_y_val_dict[(x,z)] = y
            
    max_y_val = np.array([[x,z,y] for (x,z), y in max_y_val_dict.items()])   
            
    X = np.reshape(max_y_val[:,0],(nx,nz))
    Y = np.reshape(max_y_val[:,2],(nx,nz))
    Z = np.reshape(max_y_val[:,1],(nx,nz))

    plt.figure()
    ax = plt.axes(projection = '3d')
    #ax.plot_wireframe(X,Z,Y, color = 'black') 
    ax.scatter(X,Z,Y, cmap = 'viridis', linewidth = 0.5) 
    #ax.plot_surface(X,Z,Y, rstride = 1, cstride = 1, cmap = 'viridis', edgecolor = 'none') 
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_zlabel('y')    
    
def liquid_visualisation(solid,liquid,gas):
    #determines where the bottom of the liquid surface is to visualise structure
    
    
    liquid_structure = np.zeros((nx,ny,nz))
    liquid_structure[np.logical_and(liquid > solid, liquid > gas)] = 1
    
    liquid_offset1 = liquid_structure[:,:-1,:]
    liquid_offset2 = liquid_structure[:,1:,:]
    
    liquid_bottom = np.where(np.logical_and(liquid_offset1 + liquid_offset2 == 1, liquid_offset1 == 0) == True)
    coords = np.array((liquid_bottom[0],liquid_bottom[2],liquid_bottom[1])).T
    
    min_y_val_dict = {}
    
    for i in coords:
        x, z, y = i
        
        if (x,z) not in min_y_val_dict or y < min_y_val_dict[(x,z)]:
            min_y_val_dict[(x,z)] = y
            
    min_y_val = np.array([[x,z,y] for (x,z), y in min_y_val_dict.items()])   
            
    X = np.reshape(min_y_val[:,0],(nx,nz))
    Y = np.reshape(min_y_val[:,2],(nx,nz))
    Z = np.reshape(min_y_val[:,1],(nx,nz))

    plt.figure()
    ax = plt.axes(projection = '3d')  
    ax.plot_surface(X,Z,Y, rstride = 1, cstride = 1, cmap = 'viridis', edgecolor = 'none') 
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_zlabel('y')
    
def side_profile(explicit_structure,cut, dimensions = 3, figuresave = False, figfilename = ''):
    #a snapshot of the y-z plane for a specific x
    
    cut = int(np.floor(cut))
    plt.figure()
    ax = plt.axes()
    
    
    if dimensions == 3:
        region_of_interest = explicit_structure[:,:,cut]
        ax.contour(region_of_interest)
    else:
        region_of_interest = explicit_structure[:,:,0]
        ax.contour(region_of_interest)
        #if figuresave == True:
            #plt.savefig(figfilename[:-4]+'/.png', dpi = 500)
    

def sl_interface(explicit_structure):
    #determines the lowest point of contact between the solid and liquid phases
    #to work out 'h' for 2D cone geometries
    ex_structure = np.squeeze(explicit_structure).T
    
    es_offset1 = ex_structure[:-1,:]
    es_offset2 = ex_structure[1:,:]

    sl_interface = np.logical_and(es_offset1 + es_offset2 == 0, es_offset2 == -1)

    if np.sum(sl_interface) > 0:
        points_of_contact = np.sum(sl_interface, axis = 1)
        minimum_contact = np.min(np.where(np.sum(sl_interface, axis = 1) > 0)) + 1
        
        minimum_contact = ny - minimum_contact
        print(minimum_contact)
        return minimum_contact
    
def solid_frac(nx,spacing):
    phi = (nx**2-spacing**2)/nx**2
    return phi    

def critP(spacing):
    critical_pressure = -2*100*np.cos(112*np.pi/180)/spacing
    return critical_pressure

def max_volume(contact_angle,a):
    
    contact_angle = contact_angle/180*np.pi
    volume = np.pi/3*(a/np.sin(contact_angle))**3*(2+np.cos(contact_angle))*(1-np.cos(contact_angle))**2*1e9
    
    return volume

def apparent_ca(contact_angle, solid_frac):
    contact_angle = contact_angle/180*np.pi
    apparent_theta = np.arccos(-1+solid_frac*(np.cos(contact_angle)+1))
    apparent_theta = apparent_theta*180/np.pi
    
    return apparent_theta


#%%   

solid, liquid, gas = data_extract("res-1/st-100/s-14/ca-105/p-0.000000.txt")

explicit_structure = region_definition(solid, liquid, gas)
if nz > 1:
    #solid_visualisation(solid, liquid, gas) 
    #liquid_visualisation(solid,liquid,gas) 
    side_profile(explicit_structure,21)
else:
    min_contact = sl_interface(explicit_structure)
    side_profile(explicit_structure,0, dimensions = 2)


#%%    
solid_1d = np.ravel(solid)

filename_s = 'solid_output.txt'

outfile = open(filename_s,'w')

for number in range(len(solid_1d)):
    outfile.write(str(solid_1d[number])+',')
    
outfile.close()

liquid_1d = np.ravel(liquid)

filename_l = 'liquid_output.txt'

outfile = open(filename_l,'w')

for number in range(len(liquid_1d)):
    outfile.write(str(liquid_1d[number])+',')
    
outfile.close()
    
#%%

spacings = np.linspace(20,70,6)
phis = solid_frac(nx, spacings)

critP20 = np.array([1826.2,1611.3,1435.6,1240.2,1064.5,869.1])
critP30 = np.array([3642.6,3017.6,2568.4,2216.8,1923.8,1650.4])
critP60 = np.array([9824.2,7294.9,5771.5,4775.4,4072.3,3564.5])
critP90 = np.array([14355.5,10058.6,7753.9,6308.6,5302.7,4599.6])
critP120 = np.array([16894.5,11425.8,8613.3,6933.6,5791.0,4990.2])
critP130 = np.array([17441.4,11640.6,8808.6,7050.8,5888.7,5048.8])
critP150 = np.array([17949.2,11972.7,9003.9,7207.0,5996.1,5136.7])

reduced_critP = critP130/1000*100/100

#%%
plt.figure()
plt.plot(phis,critP20,label = '20°')
plt.plot(phis,critP30,label = '30°')
plt.plot(phis,critP60,label = '60°')
plt.plot(phis,critP90,label = '90°')
plt.plot(phis,critP120,label = '120°')
plt.plot(phis,critP130,label = '130°')
plt.plot(phis,critP150,label = '150°')
plt.legend()
plt.xlabel('Solid Fraction Φ')
plt.ylabel('Critical Pressure (Pa)')
plt.rcParams['figure.dpi'] = 500
#plt.savefig('testfig.png', dpi = 1000)
plt.show()


#%%
theta = np.array([15,30,45,60,90,95,100,105,110,115,120,130,140])

d = np.array([60,80,100])

critP_1 = np.array([0.65,1.92,3.08,4.03,5.22,5.34,5.44,5.52,5.60,5.65,5.69,5.79,5.85])

critP_2 = np.array([0.56,1.53,2.39,3.12,3.99,4.07,4.15,4.21,4.25,4.29,4.31,4.37,4.42])

critP_3 = np.array([0.46,1.28,1.96,2.53,3.23,3.31,3.35,3.39,3.43,3.45,3.47,3.51,3.54])

critP_4 = np.array([0.361,0.95,1.46,1.85,2.35,2.39,2.43,2.45,2.47,2.49,2.49,2.51,2.53])

test = [critP_1, critP_2, critP_3]

critP_all = np.zeros([len(d),len(theta)])

for i in range(len(d)):
    critP_all[i,:] = test[i]

x_adjusted = np.zeros_like(critP_all)

for i in range(len(theta)):
    x_adjusted[:,i] = 4*100*np.sin(theta[i]*np.pi/180)/d
    


#%%

plt.figure()
plt.plot(theta,critP_1*6/14, marker = 'o',label = 'P_1')
plt.plot(theta,critP_2*8/14, marker = 'o', label = 'P_2')
plt.plot(theta,critP_3*10/14, marker = 'o', label = 'P_3')
plt.plot(theta,critP_4, marker = 'o', label = 'P_4')
plt.xlabel('Contact Angle (°)')
plt.ylabel('Critical Pressure (Scaled to a = 140)')
plt.legend()

critP_ratio1 = 1 - critP_1*0.75/critP_2

critP_ratio2 = 1 - critP_1*0.6/critP_3

critP_ratio3 = 1 - critP_2*0.8/critP_3

critP_ratio4 = 1 - critP_3/critP_4*100/140


plt.figure()
plt.plot(theta,critP_ratio1, marker = 'o', label = 'a_1 v a_2')
plt.plot(theta,critP_ratio2, marker = 'o', label = 'a_1 v a_3')
plt.plot(theta,critP_ratio3, marker = 'o', label = 'a_2 v a_3')
plt.plot(theta,critP_ratio4, marker = 'o', label = 'a_3 v a_4')
plt.xlabel('Contact Angle (°)')
plt.ylabel('1 - Ratio of Pc x Ratio of a')
plt.legend()

#%%
theta_sin = np.sin(theta*np.pi/180)

plt.figure()
plt.plot(theta_sin,critP_1*0.6, marker = 'o')
plt.plot(theta_sin,critP_2*0.8, marker = 'o')    
plt.plot(theta_sin,critP_3, marker = 'o')  

#%%
solid_frac_inf = [15,30,50,70]
critp_mine = [185,78.1,41.2,21.3]  
critp_inf = [114,67.5,38,30]

solid_frac_theory = np.linspace(15,70,200)

critp_hybrid = (1+1/np.sqrt(2))*2*78.2e-3/((solid_frac_theory-2.8)*1e-6)*np.cos((90-42)*np.pi/180)/100
critp_diagonal = 4/np.sqrt(2)*78.2e-3/((solid_frac_theory-2.8)*1e-6)*np.cos((90-42)*np.pi/180)/100
critp_side = 4*78.2e-3/((solid_frac_theory-2.8)*1e-6)*np.cos((90-42)*np.pi/180)/100

plt.figure()
plt.plot(solid_frac_inf,critp_mine, marker = 'o', label = 'My Simulations')
plt.plot(solid_frac_inf,critp_inf, marker = 'o', label = 'Measurements')

plt.plot(solid_frac_theory,critp_hybrid, label = 'Hybrid')
plt.plot(solid_frac_theory,critp_diagonal, label = 'Diagonal')
plt.plot(solid_frac_theory,critp_side, label = 'Side')

plt.xlabel('δ (μm)')
plt.ylabel('ΔP (mbar)')
plt.title('θ = 30°')
plt.legend()

#%%
square_theta = np.array([10,15,20,30,40,45,50,60,70,80,90,95,100,105,110,115,120,130,140,150,160])
square_theta_complete = np.linspace(10,160,1000)
    
square_alpha = np.array([0.382958237,0.488178913,0.569995668,0.665,0.707854341,0.72266313,0.726459157,0.747668599,
     0.767272174, 0.796094464, 0.8225,0.839695294,0.863620333,0.887749325,0.919981684, 0.961593856,
     1.006321519,1.146800304,1.377593449,1.785,2.609495427])  

critP_4_2 = np.array([0.19,0.361,0.557,0.95,1.3,1.46,1.59,1.85,2.06,2.24,2.35,2.39,
                      2.43,2.45,2.47,2.49,2.49,2.51,2.53,2.55,2.55]) 

alpha_hybrid = np.zeros_like(square_theta_complete)
alpha_hybrid[:] = 2*(1+1/np.sqrt(2))/4

alpha_diagonal = np.zeros_like(square_theta_complete)
alpha_diagonal[:] = 4/np.sqrt(2)/4

alpha_side = np.zeros_like(square_theta_complete)
alpha_side[:] = 1 

plt.figure()
plt.plot(square_theta,square_alpha, marker = 'o', label = 'a = 140')
plt.plot(square_theta_complete, alpha_hybrid, label = 'hybrid')
plt.plot(square_theta_complete, alpha_diagonal, label = 'diagonal')
plt.plot(square_theta_complete, alpha_side, label = 'side')    

plt.xlabel('Contact Angle (°)')
plt.ylabel('Alpha')
plt.legend()

plt.figure()
plt.plot(square_theta,critP_4_2, marker = 'o')
plt.xlabel('Contact Angle (°)')
plt.ylabel('ΔP (A.U.)')

#%%
x = np.array([60,80,100,140,200])
y = np.array([0.98553691,0.995351864,1.001702717,1.006321519,1.009237138])

plt.figure()
plt.plot(x,y, marker = 'o')

plt.xlabel('δ - d')
plt.ylabel('Alpha')


#%%
x_real = np.array([12.2,27.2,47.2,67.2])
y_infineon = np.array([114,67,38,30])


x_real_detailed = np.linspace(12.2,70,100)
y_sim = 1.14258 * (72.8e-3/100)/(x_real_detailed*1e-6/240)/100
y_sim2 = 0.76*4*72.8e-3*np.sin(60*np.pi/180)/(x_real_detailed*1e-6)/100

y_theoretical_side = 4*72.8e-3/(x_real_detailed*1e-6)*np.sin(60*np.pi/180)/100
y_theoretical_diagonal = 0.71*4*72.8e-3/(x_real_detailed*1e-6)*np.sin(60*np.pi/180)/100
y_theoretical_hybrid = 0.85*4*72.8e-3/(x_real_detailed*1e-6)*np.sin(60*np.pi/180)/100

y_error = [9,6,2,1.5]
y_error = [y_error, y_error]


plt.figure()
plt.errorbar(x_real,y_infineon, marker = 'o', label = 'infineon values', yerr = y_error, linewidth = 2)
plt.plot(x_real_detailed,y_sim2, marker = ',', label = 'simulation results', linewidth = 2)
plt.plot(x_real_detailed, y_theoretical_hybrid,label = 'hybrid', linewidth = 0.5)
plt.plot(x_real_detailed, y_theoretical_diagonal,label = 'diagonal', linewidth = 0.5)
plt.plot(x_real_detailed, y_theoretical_side,label = 'side', linewidth = 0.5)

plt.xlabel('Hole Size (μm)')
plt.ylabel('Critical Pressure (mbar)')
plt.title('θ = 60°')
plt.legend()


#%%
theta_predict = np.array([30,45,60,90,110,120,130])
alpha_predict = np.array([0.71,0.8,0.8,0.87,0.95,1.01,1.16])
alpha_predict_90 = np.array([0.71,0.8,0.8,0.87,0.89,0.88,0.89])


plt.figure()
plt.plot(theta_predict,alpha_predict,marker='o', label = 'θ > 90°')
plt.plot(theta_predict,alpha_predict_90,marker = 'o', label = 'θ = 90°')

plt.plot(np.linspace(30,130,100),0.85*np.ones(100), label = 'hybrid')
plt.plot(np.linspace(30,130,100),0.71*np.ones(100), label = 'diagonal')
plt.plot(np.linspace(30,130,100),np.ones(100), label = 'side')

plt.xlabel('Contact Angle (°)')
plt.ylabel('Predicted Alpha')

plt.legend()

#%%

intrinsic_ca = np.array([10,15,20,30,40,45,50,60,75,90,100,110,120,130,140,150])

critp_15 = np.array([17.48,31.47,47.78,80.42,110.72,124.71,136.36,157.34,182.98,199.30,203.96,208.62,208.62,210.95,213.28,213.28])
critp_30 = np.array([7.84,14.11,21.43,36.07,49.66,55.93,61.16,70.57,82.07,89.39,91.48,93.57,93.57,94.62,95.66,95.66])
critp_50 = np.array([4.52,8.13,12.35,20.79,28.62,32.23,35.25,40.67,47.30,51.51,52.72,53.92,53.92,54.53,55.13,55.13])
critp_70 = np.array([3.17,5.71,8.68,14.60,20.10,22.64,24.76,28.56,33.22,36.18,37.03,37.87,37.87,38.30,38.72,38.72])

plt.figure()
plt.plot(intrinsic_ca,critp_15,marker = 'o', label = 'Hole size = 12.2 μm')
plt.plot(intrinsic_ca,critp_30,marker = 'o', label = 'Hole size = 27.2 μm')
plt.plot(intrinsic_ca,critp_50,marker = 'o', label = 'Hole size = 47.2 μm')
plt.plot(intrinsic_ca,critp_70,marker = 'o', label = 'Hole size = 67.2 μm')

plt.xlabel('Intrinsic contact angle (°)')
plt.ylabel('Water critical pressure (mbar)')

plt.legend()