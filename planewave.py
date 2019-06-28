import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 24 # pixels/μm
time = 200

cell_size = mp.Vector3(14,10,0)

pml_layers = [mp.PML(thickness=2,direction=mp.X)]

# rotation angle (in degrees) of planewave, counter clockwise (CCW) around z-axis
rot_angle = np.radians(10)

fsrc = 2.5 # frequency of planewave (wavelength = 1/fsrc)

n = 1  # refractive index of homogeneous material
default_material = mp.Medium(index=n)

geometry = [mp.Block(size = mp.Vector3(5,10,0),
                     center = mp.Vector3(3.5,0,0),
                     material = mp.Medium(epsilon=1.5))]


#We will use p-polarized wave (Polarized in the z axis)
k_point = mp.Vector3(fsrc*n).rotate(mp.Vector3(z=1), rot_angle)

sources = [mp.EigenModeSource(src=mp.ContinuousSource(fsrc),
                              center=mp.Vector3(-1.5,0,0),
                              size=mp.Vector3(y=10),
                              direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION, 
                              #we need to set the direction to no direction when we want to rotate our plane wave
                              eig_kpoint=k_point)]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    boundary_layers=pml_layers,
                    geometry = geometry,
                    sources=sources,
                    k_point=k_point)


fcen = 2.5
df= 0.1
nfreq = 100           

          
refl_fr = mp.FluxRegion(center=mp.Vector3(0,0,0), size=mp.Vector3(0,10,0))      
refl = sim.add_flux(fcen, df, nfreq,refl_fr)


tran_fr = mp.FluxRegion(center=mp.Vector3(0,0,0), size=mp.Vector3(0,10,0))
tran = sim.add_flux(fcen, df, nfreq,tran_fr)



sim.run(until=time)


# for normalization run, save flux fields data for reflection plane
straight_refl_flux = mp.get_fluxes(refl)
#print(straight_refl_data)


# save incident power for transmission plane
straight_tran_flux = mp.get_fluxes(tran)
#print(straight_tran_flux)






nonpml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(10,10,0))
ez_data = sim.get_array(vol=nonpml_vol, component=mp.Ez)

plt.figure()
plt.imshow(np.flipud(np.transpose(np.real(ez_data))), interpolation='spline36', cmap='RdBu')
plt.axis('off')
plt.show()




#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
##################       Second Simulation      #############################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################



sim.reset_meep()


cell_size = mp.Vector3(14,10,0)

pml_layers = [mp.PML(thickness=2,direction=mp.X)]



fsrc = 2.5 # frequency of planewave (wavelength = 1/fsrc)

n = 1  # refractive index of homogeneous material
default_material = mp.Medium(index=n)

k_point = mp.Vector3(fsrc*n).rotate(mp.Vector3(z=1), rot_angle)

geometry = [mp.Block(size = mp.Vector3(5,10,0),
                     center = mp.Vector3(3.5,0,0),
                     material = mp.Medium(epsilon=1))]

sources = [mp.EigenModeSource(src=mp.ContinuousSource(fsrc),
                              center=mp.Vector3(-1.5,0,0),
                              size=mp.Vector3(y=10),
                              direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION, 
                              eig_kpoint=k_point)]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    boundary_layers=pml_layers,
                    geometry = geometry,
                    sources=sources,
                    k_point=k_point)


refl_fr2 = mp.FluxRegion(center=mp.Vector3(0,0,0), size=mp.Vector3(0,10,0))      
refl2 = sim.add_flux(fcen, df, nfreq,refl_fr2)


tran_fr2 = mp.FluxRegion(center=mp.Vector3(0,0,0), size=mp.Vector3(0,10,0))
tran2 = sim.add_flux(fcen, df, nfreq,tran_fr2)



sim.run(until=time)


# for normalization run, save flux fields data for reflection plane
bend_refl_flux = mp.get_fluxes(refl2)
#print(straight_refl_data)


# save incident power for transmission plane
bend_tran_flux = mp.get_fluxes(tran2)
#print(straight_tran_flux)


flux_freqs = mp.get_flux_freqs(refl)


nonpml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(10,10,0))
ez_data = sim.get_array(vol=nonpml_vol, component=mp.Ez)

plt.figure()
plt.imshow(np.flipud(np.transpose(np.real(ez_data))), interpolation='spline36', cmap='RdBu')
plt.axis('off')
plt.show()



wl = []
Rs = []
Ts = []
for i in range(1):
    wl = np.append(wl, 1/flux_freqs[i])
    Ts = np.append(Ts,abs(straight_tran_flux[i]/bend_tran_flux[i]))
    Rs = np.append(Rs,1-Ts[i])
       

if mp.am_master():
    plt.figure()
    plt.plot(wl,Rs,'bo-',label='reflectance')
    plt.plot(wl,Ts,'ro-',label='transmittance')
    plt.axis([0, 1, 0, 1])
    plt.xlabel("wavelength (μm)")
    plt.legend(loc="upper right")
    plt.show()
    print(wl,Ts)