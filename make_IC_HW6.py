################################################################################
###### This is an example script to generate HDF5-format ICs for GIZMO
######  The specific example below is obviously arbitrary, but could be generalized
######  to whatever IC you need. 
################################################################################
################################################################################

## load libraries we will use 
import numpy as np
import h5py as h5py

# the main routine. this specific example builds an N-dimensional box of gas plus 
#   a collisionless particle species, with a specified mass ratio. the initial 
#   gas particles are distributed in a uniform lattice; the initial collisionless 
#   particles laid down randomly according to a uniform probability distribution 
#   with a specified random velocity dispersion
#
def make_IC():
    '''
    This is an example subroutine provided to demonstrate how to make HDF5-format
    ICs for GIZMO. The specific example here is arbitrary, but can be generalized
    to whatever IC you need
    '''
    x_min = 0
    x_max = 2
    q0 = 4
    q1 = 0.5
    gamma = 1.4
    N = 200
    DIMS=1; # number of dimensions 
    
    fname='sodShockHW2.hdf5'; # output filename 

    Lbox = 2.0 # box side length
    dx = (x_max - x_min)/N
    
    #xs = np.linspace(x_min, x_max, N+1)
    
    center = np.linspace(x_min+dx/2, x_max-dx/2, N)
    rhos = np.array([1 if i <= 0.75 else 0.125 for i in center])
    Ps = np.array([1 if i<=0.75 else 0.1 for i in center])
    
    es = Ps*(gamma-1)**-1*rhos**-1
    cs = (gamma*Ps/rhos)**0.5
    qs = np.array([0 for i in rhos])
    
    # first we set up the gas properties (particle type 0)
    
    # make a regular 1D grid for particle locations (with N_1D elements and unit length)
    x0=np.arange(-0.5,0.5,1./N); x0+=0.5*(0.5-x0[-1]);
    
    xv_g=x0; yv_g = 0.0*xv_g; zv_g = 0.0*xv_g; 
        
    # the gas particle number is the lattice size: this should be the gas particle number
    Ngas = xv_g.size
    # flatten the vectors (since our ICs should be in vector, not matrix format): just want a
    #  simple list of the x,y,z positions here. Here we multiply the desired box size in
    xv_g=xv_g.flatten()*Lbox; yv_g=yv_g.flatten()*Lbox; zv_g=zv_g.flatten()*Lbox; 
    
    # set the initial velocity in x/y/z directions (here zero)
    vx_g=0.*xv_g; vy_g=0.*xv_g; vz_g=0.*xv_g;
    
    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=0.*xv_g; by_g=0.*xv_g; bz_g=0.*xv_g;
    
    # set the particle masses. Here we set it to be a list the same length, with all the same mass
    mv_g = (x_max-x_min)*rhos/N
    
    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    uv_g=Ps*(gamma-1)**-1*rhos**-1
    # set the gas IDs: here a simple integer list
    id_g=np.arange(1,N+1)

    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
    file = h5py.File(fname,'w') 

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    ##npart = np.array([Ngas,0,0,Ngrains,0,0]) # we have gas and particles we will set for type 3 here, zero for all others
    npart = np.array([Ngas,0,0,0,0,0])
    # now we make the Header - the formatting here is peculiar, for historical (GADGET-compatibility) reasons
    h = file.create_group("Header");
    # here we set all the basic numbers that go into the header
    # (most of these will be written over anyways if it's an IC file; the only thing we actually *need* to be 'correct' is "npart")
    h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs 
                                         #  if we make a multi-part IC file. with this simple script, we aren't equipped to do that.
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since 
                                        # we set masses explicitly by-particle this should be zero. that is more flexible anyways, as it 
                                        # allows for physics which can change particle masses 
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    ## no other parameters will be read from the header if we are parsing an HDF5 file for reading in
    
    # Now, the actual data!
    #   These blocks should all be written in the order of their particle type (0,1,2,3,4,5)
    #   If there are no particles of a given type, nothing is needed (no block at all)
    #   PartType0 is 'special' as gas. All other PartTypes take the same, more limited set of information in their ICs
    
    # start with particle type zero. first (assuming we have any gas particles) create the group 
    p = file.create_group("PartType0")
    # now combine the xyz positions into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=xv_g; q[:,1]=yv_g; q[:,2]=zv_g;
    # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q)
    # similarly, combine the xyz velocities into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=vx_g; q[:,1]=vy_g; q[:,2]=vz_g;
    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=q)
    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)
    # write particle masses to the Masses block
    p.create_dataset("Masses",data=mv_g)
    # write internal energies to the InternalEnergy block
    p.create_dataset("InternalEnergy",data=uv_g)
    # combine the xyz magnetic fields into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=bx_g; q[:,1]=by_g; q[:,2]=bz_g;
    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with 
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can 
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=q)

    # no PartType1 for this IC
    # no PartType2 for this IC
    # no PartType3
    # no PartType4 for this IC
    # no PartType5 for this IC

    # close the HDF5 file, which saves these outputs
    file.close()
    # all done!
