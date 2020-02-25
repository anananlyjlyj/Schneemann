import yt
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3
plt.rcParams["animation.html"] = "jshtml"

unit_base = {'UnitLength_in_cm'         : 1, #1.496e13,  
             'UnitMass_in_g'            : 1, #1.989e33,  
             'UnitVelocity_in_cm_per_s' : 1} #474318.3259}

bbox_lim = 80*1.496e13

bbox = [[0,bbox_lim],
        [0,bbox_lim],
        [0,bbox_lim]]

def PartProjAni(folder):
    yt.funcs.mylog.setLevel(50)
    ts = yt.load("/u/yali/"+folder+"/test/output/snapshot_*",unit_base=unit_base,bounding_box=bbox)
    
    #sorted(ds.field_list)
    #print(dir(ds.fields.PartType0))

    #plot = yt.ProjectionPlot(ts[0], 'x', 'density')
    plot = yt.ParticleProjectionPlot(ts[0], 'x', fields = "Masses")
    #plot.set_zlim('density', 8e-29, 3e-26)
    
    fig = plot.plots['Masses'].figure

    def animate(i):
        ds = ts[i]
        plot._switch_ds(ds)
        plot.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
        plot.annotate_scale(corner='upper_right')
    
    plot.set_buff_size((200, 200))
    ani = FuncAnimation(fig, animate, frames=len(ts))
    return ani

def pp(folder):
    ds = yt.load("/u/yali/"+folder+"/test/output/snapshot_001.hdf5",unit_base=unit_base,bounding_box=bbox)
    #plot = yt.ParticleProjectionPlot(ds, 'z', fields = "Masses")
    #plot = yt.ProjectionPlot(ds, "x", "density")
    
    # Modify the projection
    # The argument specifies the region along the line of sight for which particles will be gathered.
    # 1.0 signifies the entire domain in the line of sight
    # p.annotate_particles(1.0)
    # but in this case we only go 10 Mpc in depth
    
    #plot.annotate_particles((10, 'km'))
    plot = yt.ParticlePhasePlot(ds.all_data(), ('PartType3', 'particle_position_y'), ('PartType3', 'particle_position_z'),('PartType3','Masses'))
    #fig = plot.plots['PartType3','Density'].figure
    #fig.show()
    #ad = ds.all_data()
    #density = ad[("PartType0","density")]
    #plot = yt.ProjectionPlot(ds, 'x', ('gas', 'density'))
    
    plot.show()
    
    
def radial_qntty(what, file, quantity, rad_max, rad_num, cut = 0):
    # e.g. file = yt.load("/u/yali/gDisk1/test/output_softened_Potential/snapshot_000.hdf5",unit_base=unit_base, bounding_box=bbox)
    if what is "intg" or "avg":
        plot = yt.ProjectionPlot(file, 'z', ('gas', quantity), center = [50,50,50], width = [100,100,100])
    if what is "cut":
        plot = yt.SlicePlot(file, "z", ('gas', quantity), center=[50,50,50+cut], width = [100,100,100])        

    qntty_mtx = plot.frb[('gas', quantity)]
    qntty = qntty_mtx.flatten()

    center_pxl = np.array([plot.buff_size[0]/2, plot.buff_size[1]/2])   # buff_size = number of pixel in 1d
    dpx = (plot.bounds[1]-plot.bounds[0])/plot.buff_size[0]    
    dpy = (plot.bounds[3]-plot.bounds[2])/plot.buff_size[1]
    dist_x = np.array([np.abs(np.arange(plot.buff_size[0]) - center_pxl[0]) * dpx] * plot.buff_size[1])
    dist_y = np.array([np.abs(np.arange(plot.buff_size[1]) - center_pxl[1]) * dpy] * plot.buff_size[0]).transpose()
    dist_mtx = (dist_x**2 + dist_y**2)**0.5
    dist = dist_mtx.flatten()

    EPSILON = 0.005    # average quantity in reasonable range
    rad_list = np.linspace(0, rad_max, rad_num)
    qntty_list = []
    for r in rad_list:
        x = qntty[np.abs(dist - r) < EPSILON]
        qntty_list.append(x.sum()/np.size(x))
    res = np.asarray(qntty_list)
    
    if what is "avg":
        height = 1
        res /= height
    
    return  ts0.arr(rad_list, 'code_length'), res

def qntty_movie(what, files, qntty1, qntty2, rad_max, rad_num, cut = 0, frame = 1):
    fig = plt.figure(figsize=(7,7))
    ax = plt.axes(xlim=(0,rad_max)) #,ylim=(size1, size2))
    #dens_plot = ax.loglog(rp.x.value, rp["density"].value)
    res_plot = ax.scatter([],[])
    ax.set_xlabel(qntty1)#r"$\mathrm{r\ (AU)}$")
    ax.set_ylabel(qntty2)
    #ax.set_yscale("log")
    ax.set_ylim([0,100])
    
    def update(i):
        cur = files[i]
        #disk = cur.disk(cur.domain_center, [1,0,0], (50., "AU"), (20., "AU")) # Get a disk object
        # Bin up the data from the disk into a radial profile
        #rp = yt.create_profile(disk, 'radius', [('gas', 'density'), ('gas', 'temperature'), ('gas', 'angular_momentum_magnitude')],
        #                        n_bins = 128,
        #                        units = {'radius': 'AU'}, logs = {'radius': False})
        if qntty1 is "radius":
            rad, q2 = radial_qntty(what, cur, qntty2, rad_max, rad_num, cut)
            res_plot.set_offsets(np.transpose(np.asarray([rad, q2])))#rp.x.value, rp[qntty2].value)
            #ax.set_ylim([0,np.max(q2)])
        else: 
            if qntty2 is "radius":
                rad, q1 = radial_qntty(what, cur, qntty1, rad_max, rad_num, cut)
                res_plot.set_offsets(np.transpose(np.asarray([rad, q1])))
                #ax.set_ylim([0,np.max(q1)])
            else:
                rad, q1 = radial_qntty(what, cur, qntty1, rad_max, rad_num, cut)
                rad, q2 = radial_qntty(what, cur, qntty2, rad_max, rad_num, cut)
                res_plot.set_offsets(np.transpose(np.asarray([q1, q2])))   
                #ax.set_ylim([0,np.max(q2)])
                ax.set_xlim([0,np.max(q1)])
        return res_plot,
    
    plt.close()
    ani = FuncAnimation(fig, update, frame, blit=True)
    return ani
        