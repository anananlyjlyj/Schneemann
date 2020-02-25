import numpy as np
import h5py
import mpl_toolkits.mplot3d.axes3d as p3
from IPython.display import HTML
import fnmatch
import os
import io
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.widgets import Slider, Button, RadioButtons
plt.rcParams["animation.html"] = "jshtml"
import io
from plot_static import plot_single


def ani_top(folder, frame, start):
    fig = plt.figure(figsize=(10,10))

    def update(i):
        fname = str(format(start + i, "03d"))
        f = h5py.File("/u/yali/"+folder+"/test/output/snapshot_"+fname+".hdf5", 'r')
        coord = f['PartType0']['Coordinates'][:,:]
        points.set_data(np.asarray([coord[:,0], coord[:,1]]))
        ax.set_title('= %.5f'%(np.max(f['PartType0']['InternalEnergy'][:])))
        return points,

    ax = fig.add_subplot()
    points, = ax.plot([], [], linestyle="", marker='.', markersize = 0.3, alpha=0.8)       
    ax.set_xlim(30, 70);    ax.set_ylim(30, 70)
    plt.close()
    ani = FuncAnimation(fig, update, frame, blit=True)    
    return ani

def ani_edge(folder, frame):
    fig = plt.figure(figsize=(10,10))

    def update(i):
        fname = str(format(i, "03d"))
        f = h5py.File("/u/yali/"+folder+"/test/output/snapshot_"+fname+".hdf5", 'r')
        coord = f['PartType0']['Coordinates'][:,:]
        points.set_data(np.asarray([coord[:,1][np.abs(coord[:,0] - 50) < 6], coord[:,2][np.abs(coord[:,0] - 50) < 6]]))
        ax.set_title('= %.5f'%(np.max(f['PartType0']['InternalEnergy'][:])))
        return points,

    ax = fig.add_subplot()
    points, = ax.plot([], [], linestyle="", marker='.', markersize = 0.3, alpha=0.8)       
    ax.set_xlim(30, 70);    ax.set_ylim(30, 70)
    plt.close()
    ani = FuncAnimation(fig, update, frame, blit=True)    
    return ani


def ani_3d(folder, frame):
    fig = plt.figure(figsize=(10,10))

    def update(i):
        fname = str(format(i, "03d"))
        f = h5py.File("/u/yali/"+folder+"/test/output/snapshot_"+fname+".hdf5", 'r')
        coord = f['PartType0']['Coordinates'][:,:]
        points.set_data(np.asarray([coord[:,0], coord[:,1]]))  
        points.set_3d_properties(np.asarray(coord[:,2]))
        ax.set_title('= %.5f'%(np.max(f['PartType0']['InternalEnergy'][:])))
        return points,

    ax = fig.add_subplot(111, projection='3d')
    points, = ax.plot([], [], [], linestyle="", marker='.', markersize = 0.1, alpha=0.6)       
    ax.set_xlim(30, 70);    ax.set_ylim(30, 70);      ax.set_zlim(30, 70)
    ax.view_init(elev=35., azim=20.)
    plt.close()
    ani = FuncAnimation(fig, update, frame, blit=True)    
    return ani


class plots:
    def __init__(self, folder, output):
        self.folder = folder
        self.output = output
    def scatter_qty(self, face, qty = "None", frame = "None"):
        if frame is None:
            frame = len(fnmatch.filter(os.listdir('/u/yali/'+self.folder+'/test/output'+self.output), '*.hdf5'))
        if face is "3d":
            fig = plt.figure(figsize=(10,10))
        else:
            fig = plt.figure(figsize=(10,8))
        
        def update(i):
            cur = plot_single(self.folder, self.output, i)
            if face is "3d":
                points.set_data(np.asarray([cur.coordx, cur.coordy]))
                points.set_3d_properties(np.asarray(cur.coordz))
                if i == 0:
                    xmax = np.max(cur.coordx); xmin = np.min(cur.coordx); ymax = np.max(cur.coordy); ymin = np.min(cur.coordy); 
                    zmax = np.max(cur.coordz); zmin = np.min(cur.coordz);
                    ax.set_xlim3d(xmin-0.5*(xmax-xmin),xmax+0.5*(xmax-xmin));    ax.set_ylim3d(ymin-0.5*(ymax-ymin),ymax+0.5*(ymax-ymin));   
                    ax.set_zlim3d(zmin-0.5*(zmax-zmin),zmax+0.5*(zmax-zmin))
                #ax.set_title('%03d MAX '%(i) + qty + '= %.2f'%(np.max(getattr(cur, qty))))
                #line.set_color(cur.dens)
            else:
                if face is "xy":
                    res1 = cur.coordx; res2 = cur.coordy
                    ax.set_xlabel('x');    ax.set_ylabel('y');
                    #plt.text(cur.coordx[cur.ids == 0], cur.coordy[cur.ids == 0], "particle 0", fontsize=12)
                elif face is "yz":
                    res1 = cur.coordy; res2 = cur.coordz
                    ax.set_xlabel('y');    ax.set_ylabel('z');
                elif face is "xz":
                    res1 = cur.coordx; res2 = cur.coordz
                    ax.set_xlabel('x');    ax.set_ylabel('z');
                elif face is "rz":
                    res1 = cur.rad; res2 = cur.height
                    ax.set_xlabel('r');    ax.set_ylabel('z');
                points.set_offsets(np.transpose([res1, res2]))
                if i == 0:
                    xmax = np.max(res1); xmin = np.min(res1); ymax = np.max(res2); ymin = np.min(res2)
                    ax.set_xlim(xmin-0.5*(xmax-xmin),xmax+0.5*(xmax-xmin));        ax.set_ylim(ymin-0.5*(ymax-ymin),ymax+0.5*(ymax-ymin))
                points.set_array(getattr(cur, qty))
                ax.set_title('%03d MAX '%(i) + qty + '= %.5f'%(np.max(getattr(cur, qty))))
                            
            return points,
        
        if face is "3d": 
            ax = fig.add_subplot(111, projection='3d')
            #ax = p3.Axes3D(fig)
            points, = ax.plot([], [], [], linestyle="", marker='.', markersize = 1, alpha=0.6)
            ax.view_init(elev=35., azim=20.)
            ax.set_xlabel('x');    ax.set_ylabel('y');    ax.set_zlabel('z')
        else:
            ax = fig.add_subplot()
            points = ax.scatter([],[], c=[])
            norm = Normalize(0.0, 1.0)
            colormap = cm.ScalarMappable(norm, 'gnuplot')
            fig.colorbar(colormap)

        plt.close()
        ani = FuncAnimation(fig, update, frame, blit=True)
        print(frame)
    
        return ani
    
    def qty_qty(self, qty1, qty2, frame = "None", size1 = 0, size2 = 50):
        if frame is None:
            frame = len(fnmatch.filter(os.listdir('/u/yali/'+self.folder+'/test/output'), '*.hdf5'))
        fig = plt.figure(figsize=(10,10))
        
        def update(i):
            cur = plot_single(self.folder, self.output, i)
            res1 = np.transpose(getattr(cur, qty1));   res2 = np.transpose(getattr(cur, qty2))
            xmax = np.max(res1); xmin = np.min(res1); ymax = np.max(res2); ymin = np.min(res2)
            if qty1 is "dens":
                ax.set_xscale('log')
            else:
                ax.set_xlim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin));         ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
            points.set_offsets(np.transpose([res1, res2]))
            ax.set_title(qty2+" vs "+qty1+" of snapshot "+cur.fname)      
            return points,
        ax = fig.add_subplot()#xlim=(size1,size2))#,ylim=(size1, size2))
        points = ax.scatter([],[], c=[])
        plt.xlabel(qty1);          plt.ylabel(qty2)
        plt.close()
        ani = FuncAnimation(fig, update, frame, blit=True)
        print(frame)
    
        return ani
    
    def surfDens(self, frame):
        if frame is None:
            frame = len(fnmatch.filter(os.listdir('/u/yali/'+self.folder+'/test/output'+self.output), '*.hdf5'))
        fig = plt.figure(figsize=(10,10))
        
        n = 200; m = 200;       eps_rad = 21/n;  eps_phi = 2*np.pi/m
        grid_rad = np.arange(4, 25, eps_rad);  grid_phi = np.arange(-np.pi, np.pi, eps_phi)
        grid = (np.array(np.meshgrid(grid_rad, grid_phi)).T.reshape(-1,2))
        grid_x = np.multiply(grid[:,0], np.cos(grid[:,1])).ravel();  grid_y = np.multiply(grid[:,0], np.sin(grid[:,1])).ravel()
        def sigma(obj, a, b):       
            arr = obj.height[((obj.rad - a <= eps_rad/2)&(a - obj.rad< eps_rad/2)) & ((obj.phi - b <= eps_phi/2)&(b - obj.phi< eps_phi/2))]
            if len(arr) == 0:
                return 0
            return len(arr)*(0.1/250000)/(np.pi*((a + eps_rad/2)**2-(a - eps_rad/2)**2)*(eps_phi/(2*np.pi)))
        
        def update(i):
            cur = plot_single(self.folder, self.output, i)
            res = np.array([sigma(cur, i, j) for i, j in grid])
            points.set_array(res)
            ax.set_title('%03d MAX '%(i) + 'surface density = %.5f'%(np.max(res)))      
            return points,
        ax = fig.add_subplot()
        points = ax.scatter(grid_x,grid_y, c=[])
        plt.xlabel("x");          plt.ylabel("y")
        plt.close()
        ani = FuncAnimation(fig, update, frame, blit=True)
        print(frame)
    
        return ani
        
    

