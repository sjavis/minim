import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import os
import glob

nx = 120
ny = 120
nz = 1

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
        region_of_interest = explicit_structure[cut,:,:]
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


#%%   

solid, liquid, gas = data_extract("outputs/res-0.001/st-100/s-35/ca-105/p-0.000000.txt")

explicit_structure = region_definition(solid, liquid, gas)
if nz > 1:
    #solid_visualisation(solid, liquid, gas) 
    #liquid_visualisation(solid,liquid,gas) 
    side_profile(explicit_structure,50)
else:
    min_contact = sl_interface(explicit_structure)
    side_profile(explicit_structure,0, dimensions = 2)
    
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
all_critP = np.array([critP150,critP130,critP120,critP90,critP60,critP30,critP20])
phis = solid_frac(nx, spacings)

x,y = np.meshgrid([phis],[20,30,60,90,120,130,150])

all_critP = all_critP[:-1, :-1]

fig, ax = plt.subplots()
c = ax.pcolormesh(x,y,all_critP, cmap = 'viridis')
ax.axis([x.min(),x.max(),y.min(),y.max()])
fig.colorbar(c, ax=ax)

#%%

solid_1d = np.ravel(solid)

filename_s = "solid_output.txt"

outfile = open(filename_s,'w')

for number in range(len(solid_1d)):
    outfile.write(str(solid_1d[number])+ ",")

outfile.close()

liquid_1d = np.ravel(liquid)

filename_l = "liquid_output.txt"

outfile = open(filename_l,'w')

for number in range(len(liquid_1d)):
    outfile.write(str(liquid_1d[number])+ ",")

outfile.close()

#%%
spacing_theory = np.linspace(10,60,100)
spacing_sim = np.linspace(10,60,11)
critP_theory = -2*100*np.cos(112*np.pi/180)/(spacing_theory*0.001)
critP_sim = [8330,5283,3896,3096,2568,2178,1904,1689,1514,1377,1260]

ca_theory = np.linspace(60,120,100)
ca_sim = [60,85,100,110,120]
critP_ca_theory = -2*100*np.cos(ca_theory*np.pi/180)/(50*0.001)
critP_ca_sim = [-2021,-342,713,1377,2021]

#%%
plt.figure()
plt.plot(spacing_theory,critP_theory, label = 'Theory')
plt.plot(spacing_sim,critP_sim,marker = 'o',linestyle = '')
plt.legend()
plt.xlabel('Spacing (mm)')
plt.ylabel('Critical Pressure (Pa)')

plt.figure()
plt.plot(ca_theory,critP_ca_theory, label = 'Theory')
plt.plot(ca_sim,critP_ca_sim,marker = 'o',linestyle = '')
plt.legend()
plt.xlabel('Contact Angle (°)')
plt.ylabel('Critical Pressure (Pa)')

#%%
specific_path = "outputs/res-0.001/st-100/s-10/ca-105/"

for name in glob.glob(specific_path+"p-*[0-9].*"):
    print(name)
    solid, liquid, gas = data_extract(name)

    explicit_structure = region_definition(solid, liquid, gas)
    if nz > 1:
        #solid_visualisation(solid, liquid, gas) 
        #liquid_visualisation(solid,liquid,gas) 
        side_profile(explicit_structure,50)
    else:
        #min_contact = sl_interface(explicit_structure)
        side_profile(explicit_structure,0, dimensions = 2, figuresave = True, figfilename = name)
    
    
    
    
    
    