import ipyvolume as ipv
import numpy as np
import matplotlib.cm as cm
import yt
import h5py

def aniplot(folder, frame):
    if frame is None:
         frame = len(fnmatch.filter(os.listdir('/u/yali/'+folder+'/test/output'), '*.hdf5'))
            
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    v = []
    
    for i in np.arange(frame):
        fname = str(format(i, '03d'))
        f = h5py.File("/u/yali/"+folder+"/test/output/snapshot_"+fname+".hdf5", 'r') # Read-only access to the file
        #intE = []
        #density = []
        x.append(f['PartType0']['Coordinates'][:,0])
        y.append(f['PartType0']['Coordinates'][:,1])
        z.append(f['PartType0']['Coordinates'][:,2])
        #intE.append(f['Partype0']['InternalEnergy'][:])
        vx.append(f['PartType0']['Velocities'][:,0])
        vy.append(f['PartType0']['Velocities'][:,1])
        vz.append(f['PartType0']['Velocities'][:,2])
        #density.append(f['PartType0']['Density'][:])
    
    # and normalize
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    vx = np.array(vx)
    vy = np.array(vy)
    vz = np.array(vz)
    v = np.array(np.sqrt(vx**2+vy**2+vz**2))
    v -= v.min()
    v /= v.max()
    
    # map the normalized values to rgba colors
    cmap = cm.Reds
    colors = np.array([cmap(k) for k in v])
    colors.shape
    
    # plot animation
    fig = ipv.figure()
    ipv.style.use('dark')
    # use just rgb colors, not rgba
    quiver = ipv.quiver(x, y, z, vx, vy, vz, size=5, color=colors[:,:,:3])
    # create the animation widgets/slider
    ipv.animation_control(quiver, interval=500)
    ipv.show()
    
