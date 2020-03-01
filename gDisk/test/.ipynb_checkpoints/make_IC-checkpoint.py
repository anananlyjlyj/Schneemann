import numpy as np
import h5py as h5py

def make_IC():
    DIMS=3;                     # number of dimensions
    
    def is_float(x):
        try:
            num = float(x)
        except ValueError:
            return False
        return True
    def globalize(command):
        exec(command, globals())
    for line in open('./my.params'):
        if not line.lstrip().startswith("%"):
            if line.strip():
                li = line.split()
                if is_float(li[1]):
                    globalize("%s =%.10f"%(li[0], float(li[1])))
    
    ########################### units set in parameter file #########################    
    G = GravityConstantInternal #39.4429 # G = 6.6743e-11*(UnitMass_in_g/1000)/(UnitLength_in_cm/100*(UnitVelocity_in_cm_per_s/100)**2)
    gamma_eos = 1.001#5./3.           # polytropic index of ideal equation of state the run will assume
    k_B = 1.38064852e-23/(((UnitVelocity_in_cm_per_s)/100)**2 * (UnitMass_in_g/1000))   # 1.38064852e-23 m2 kg s-2 K-1
    mu = 2.3;   proton_mass = 1.6726219e-27/(UnitMass_in_g/1000);   mean_molecular_weight = mu*proton_mass
    
    ########################### physical setup ######################################
    M = 1.0                     # mass of central star
    a = 50; b = 50; c = 50      # position of central star 
    M_disk = 0.1                # toaal mass of disk
    r_o = 20.0;   r_i = 5.0     # disc outer & inner radius;
    Ngas = 250000
    
    ########################### initial profile setup ###############################
    # in the following: cylindical coordinate, r**2 = x**2 + y**2
    # Surface density \propto (r/r0)**p =============================================
    p = -1;     Sigma_0 = M_disk*(p+2) /(2.0*np.pi)/((r_o/r0)**(p+2)-(r_i/r0)**(p+2))  # Sigma_0 = Sigma(r = r0) 
    # \int{\Sigma dA} = \int{\Sigma_0 * r^p * (2*pi*r*dr)} = M_disk -> determine \Sigma_0
    def Sigma(r):
        if r < r_i or r > r_o:
            return 0
        return Sigma_0 * (r/r0)**p
    vSigma = np.vectorize(Sigma)
    
    # Temperature \propto (r/r0)**q ==================================================
    Q_for_normalize = 2         # Q = c_s*Omega/(np.pi*G*Sigma), normalised s.t. Q(r=r_o) = Q_for_normalize 
    q = -0.5
    def Temp(r):
        if r < r_i or r > r_o:
            return 200
        res = (Q_for_normalize*np.pi*G*Sigma(r_o)/Omega_K(r_o))**2/(k_B*(r_o/r0)**q/mean_molecular_weight) * (r/r0)**q #/100
        return res
    vT = np.vectorize(Temp)
    
    # vertical structure of the disk h = r * AspectRatio as a Gaussian ==============
    def c_s(r):                 # speed of sound
        return (k_B * Temp(r)/mean_molecular_weight)**0.5     
    def Omega_K(r):             # Keplerian angular speed
        return (G*M/r**3)**0.5 
    def vertical_scale_h(r):
        return c_s(r) / Omega_K(r)
    def AspectRatio(r):         # for verification: representative value = 0.05
        return vertical_scale_h(r)/r
    vO_K = np.vectorize(Omega_K); vh = np.vectorize(vertical_scale_h); vAR = np.vectorize(AspectRatio)
    
    def rho(r, z):
        rho_0 = Sigma(r)/(vertical_scale_h(r)*(2*np.pi)**0.5)
        return rho_0 * np.exp(-0.5*(z/vertical_scale_h(r))**2)
    
    def Q(r):                   # for verification
        return c_s(r)*Omega_K(r)/(np.pi*G*Sigma(r))
    
    ########################### initial coordinate ##################################
    # --- cylindical shell method ---
    def Mshell_intg(r):
        return ((r/r0)**(p+2)) /(p+2) * Sigma_0 * 2 * np.pi *r0
    # radii of particles are inverse function of Mshell_intg
    r_g = r0 * np.exp(np.log((p+2.0)*(Mshell_intg(r_i)+(Mshell_intg(r_o)-Mshell_intg(r_i))*np.random.uniform(0, 1, Ngas))/(2*np.pi*r0*Sigma_0))/(p+2.0))
    z_g = vertical_scale_h(r0) * (r_g/r0)**((3+q)/2) * np.random.normal(0, 1, Ngas)
    phi_g =  np.random.uniform(-np.pi, np.pi, Ngas)   # randomrized by pdf = uniform distribution
    xv_g = r_g*np.cos(phi_g) + a;   yv_g = r_g*np.sin(phi_g) + b;   zv_g = z_g + c
    '''
    # --- acceptance method: ---
    # see https://stackoverflow.com/questions/51050658/
    def uniform_proposal(x, delta=2.0):
        return np.random.uniform(x - delta, x + delta)
    def metropolis_sampler(p, nsamples, proposal=uniform_proposal):
        r = r_i;  z = -r*0.5                            # start somewhere
        for i in range(nsamples):
            trial_r = proposal(r)                           # random neighbour from the proposal distribution
            acceptance_r = p(trial_r, z)/p(r, z)
            if np.random.uniform() < acceptance_r:          # accept the move conditionally
                r = trial_r
            trial_z = proposal(z)
            acceptance_z = p(r, trial_z)/p(r, z)
            if np.random.uniform() < acceptance_z:          # accept the move conditionally
                z = trial_z
            yield r, z
    result = list(metropolis_sampler(rho, Ngas));       result = np.array(result)
    
    r_g = result[:,0];     z_g = result[:,1]   
    phi_g =  np.random.uniform(-np.pi, np.pi, Ngas)   # randomrized by pdf = uniform distribution
    xv_g = r_g*np.cos(phi_g) + a;   yv_g = r_g*np.sin(phi_g) + b;   zv_g = z_g + c
    
    # ---meshgrid method: ----
    # (x, y, z)
    #create meshgrid:
    x0 = np.linspace(-r_o,r_o,200)
    [x_grid,y_grid,z_grid] = np.meshgrid(x0, x0, x0)
    xi,yi,zi = x_grid.ravel(),y_grid.ravel(),z_grid.ravel()
    
    #create normalized pdf:
    pdf = rho((x_grid**2 + y_grid**2)**0.5, z_grid);
    pdf = pdf/np.sum(pdf);
    
    #obtain indices of randomly selected points, as specified by pdf:
    randices = np.random.choice(np.arange(x_grid.ravel().shape[0]), Ngas, replace = False,p = pdf.ravel())
    
    #positions
    xv_g = xi[randices]+a;     yv_g = yi[randices]+b;     zv_g = zi[randices]+c
    r2_g = (xv_g-a)**2 + (yv_g-b)**2;     r_g = r2_g**0.5;   z_g = zv_g-c
    
    # (r, phi, z)
    #create meshgrid:
    r0 = np.linspace(r_i,r_o,200)
    phi0 = np.linspace(-np.pi, np.pi, 200, endpoint = False)
    z0 = np.linspace(-0.5*r_o, 0.5*r_o, 100)
    [r_grid,phi_grid,z_grid] = np.meshgrid(r0, phi0, z0)
    ri,phii,zi = r_grid.ravel(),phi_grid.ravel(),z_grid.ravel()
    
    #create normalized pdf:
    pdf = rho(r_grid, z_grid);
    pdf = pdf/np.sum(pdf);
    
    #obtain indices of randomly selected points, as specified by pdf:
    randices = np.random.choice(np.arange(z_grid.ravel().shape[0]), Ngas, replace = False,p = pdf.ravel())
    
    #positions
    r_g = ri[randices];     phi_g = phii[randices];     z_g = zi[randices]
    xv_g = r_g*np.cos(phi_g)+a;   yv_g = r_g*np.sin(phi_g)+b;        zv_g = z_g+c
    '''
    ############################# initial mass ########################################
    ##### uniform mass = M_disk / Ngas
    mv_g = M_disk/Ngas + 0.*xv_g
    
    ############################# initial velocity ####################################
    # initial velocity from equillibrium, c_s account for pressure, F_sg self-gravity
    #F_sg = 0
    #vx_g=(G*M/r_g+r_g*F_sg-c_s**2)**0.5*(xv_g-b)/r_g; vy_g=(G*M/r_g+r_g*F_sg-c_s**2)**0.5*(yv_g-b)/r_g; vz_g=0.*zv_g;
    # equation 24 from Armitage "Physical Processes in Protoplanetary disks"
    #v_phi = v_K*(1-(11/4)*AspectRatio(r_g)**2)**0.5
    #v_phi += (2*np.pi+G*Sigma(r_g)*r_g)**0.5
    # see https://arxiv.org/pdf/1102.0671.pdf
    Omega = vOK(r_g)*((1+q) - q*r_g/(r_g**2+z_g**2)**0.5 + (p-(3+q)/2+q) * vAR(r_g)**2)**0.5
    vx_g=-Omega*(yv_g-b);     vy_g=Omega*(xv_g-a);      vz_g=0.*zv_g;
    
    # set the initial velocity in x/y/z direction = (-Omega * r_y, Omega*r_y, 0)
    #vx_g=-Omega*(yv_g-a); vy_g=Omega*(xv_g-b) ; vz_g=0.*zv_g;
    
    ############################ initial internal energy ##############################
    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    #uv_g=P_desired/((gamma_eos-1.)*rho_desired)
    uv_g = k_B * vT(r_g)/(mean_molecular_weight*(gamma_eos-1))
    
    
    ############################ initial magnetic field ################################    
    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=0.*xv_g; by_g=0.*xv_g; bz_g=0.*xv_g;    
    
    # set the gas IDs: here a simple integer list
    id_g=np.arange(1,Ngas+1)
    
    fname='ics.hdf5'; # output filename 
    
    '''
    # particles
    # vgrainrms = 1.0 # rms velocity of collisionless particles
    # dust_to_gas_ratio = 0.01 # mass ratio of collisionless particles to gas
    
    # now we set the properties of the collisionless particles: we will assign these to particle type '3', 
    #   but (barring special compile-time flags being set) GIZMO will treat all collisionless particle types
    #   the same. so the setup would be identical for any of the particle types 1,2,3,4,5
    
    # set the desired number of particles (here to about twice as many as the gas particles, because we feel like it)
    Ngrains = int(np.round(2. * (1.*N_1D)**DIMS))
    # set the x/y/z positions: again a simple list for each: here to random numbers from a uniform distribution
    xv_d = (np.random.rand(Ngrains)-0.5)*Lbox
    yv_d = (np.random.rand(Ngrains)-0.5)*Lbox
    zv_d = (np.random.rand(Ngrains)-0.5)*Lbox
    # set the IDs: these must be unique, so we start from the maximum gas ID and go up
    id_d = np.arange(Ngas+1,Ngrains+Ngas+1)
    # set the velocities. again we will set to a random value, here a Gaussian-distributed one
    vx_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vy_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vz_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    # set the masses, again a list with all the same mass
    mv_d = dust_to_gas_ratio * (1.*Ngas)/(1.*Ngrains) * mv_g[0] + 0.*xv_d
    # set the types for grains. GrainType = 1: Epstein/Stokes; 2: Charged Epstein/Stokes; 3: Cosmic Rays
    type_d = (np.ones(Ngrains) * 3).astype("int")
    '''



    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
    file = h5py.File(fname,'w') 

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    #npart = np.array([Ngas,0,0,Ngrains,0,0]) # we have gas and particles we will set for type 3 here, zero for all others
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
    ## all of the parameters below will be overwritten by whatever is set in the run-time parameterfile if
    ##   this file is read in as an IC file, so their values are irrelevant. they are only important if you treat this as a snapshot
    ##   for restarting. Which you shouldn't - it requires many more fields be set. But we still need to set some values for the code to read
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = 1.0; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
    ## ok, that ends the block of 'useless' parameters
    
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

    '''
    # no PartType1, 2, 4, 5 for this IC
    
    # now assign the collisionless particles to PartType3. note that this block looks exactly like 
    #   what we had above for the gas. EXCEPT there are no "InternalEnergy" or "MagneticField" fields (for 
    #   obvious reasons). 
    p = file.create_group("PartType3")
    q=np.zeros((Ngrains,3)); q[:,0]=xv_d; q[:,1]=yv_d; q[:,2]=zv_d;
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((Ngrains,3)); q[:,0]=vx_d; q[:,1]=vy_d; q[:,2]=vz_d;
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_d)
    p.create_dataset("Masses",data=mv_d)
    p.create_dataset("GrainType",data=type_d)
    
    '''

    # close the HDF5 file, which saves these outputs
    file.close()
    # all done!

make_IC()
