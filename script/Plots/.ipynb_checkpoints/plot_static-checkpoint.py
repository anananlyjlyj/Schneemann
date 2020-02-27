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

class plot_single:
    
    def __init__(self, folder, output, i):      # read data
        def is_float(val):
            try:
                num = float(val)
            except ValueError:
                return False
            return True
        for line in open('/u/yali/'+folder+'/test/output'+output+'/parameters-usedvalues'):
            li = line.split()
            if is_float(li[1]):
                exec("self.%s = %.10f" % (li[0],float(li[1])))
                
        self.mu = 2.3;     # mean molecular mass
        self.gamma_eos = 1.001#(5./3.)           # polytropic index of ideal equation of state the run will assume

        # G = 6.6743e-11 m3 kg-1 s-2 *(UnitMass_in_g/1000)/(UnitLength_in_cm/100*(UnitVelocity_in_cm_per_s/100)**2)
        self.k_B = 1.38064852e-23/((self.UnitVelocity_in_cm_per_s/100)**2 * (self.UnitMass_in_g/1000))   # 1.38064852e-23 m2 kg s-2 K-1
        self.proton_mass = 1.6726219e-27/(self.UnitMass_in_g/1000);                 self.mean_molecular_weight = self.mu*self.proton_mass
        self.r_i = 5.;   self.r_o = 20.;  self.m_disk = 0.1
        self.star = [50, 50, 50];  self.m_star = 1
        self.Sigma_0 = self.m_disk/(2*np.pi*(self.r_o-self.r_i));      self.power = -1

        self.fname = str(format(i, '03d'));   self.time = i*self.TimeBetSnapshot * (self.UnitLength_in_cm / self.UnitVelocity_in_cm_per_s)
        f = h5py.File("/u/yali/"+folder+"/test/output"+output+"/snapshot_"+self.fname+".hdf5", 'r') # Read-only access to the file
        
        self.ids =  f['PartType0']['ParticleIDs'][:]
        
        self.coord = f['PartType0']['Coordinates'][:,:]
        self.coordx = self.coord[:,0];     self.coordy = self.coord[:,1];   self.coordz = self.coord[:,2];
        self.rad = ((self.coordx - self.star[0])**2 + (self.coordy - self.star[1])**2)**0.5   #  + self.coordz**2
        self.phi = np.arctan2((self.coordy - self.star[1]),(self.coordx - self.star[0]))
        self.height = self.coordz - self.star[2]
        
        self.vel = f['PartType0']['Velocities'][:,:]
        self.velx = self.vel[:,0];         self.vely = self.vel[:,1];       self.velz = self.vel[:,2];
        self.vel_mag = (self.velx**2 + self.vely**2 + self.velz**2)**0.5
        self.Omega_K = (self.GravityConstantInternal*self.m_star/self.rad**3)**0.5
        
        self.dens = f['PartType0']['Density'][:]
        
        self.intE = f['PartType0']['InternalEnergy'][:]
        self.temp = (self.mean_molecular_weight/self.k_B) * (self.gamma_eos-1) * self.intE  
        #self.temp_dimless = 
        self.c_s = (self.k_B /self.mean_molecular_weight * self.temp)**0.5
        self.vertical_scale_h = np.divide(self.c_s, self.Omega_K)
        self.AspectRatio = np.divide(self.vertical_scale_h, self.rad)
        
        self.specAngMom = self.velx * (self.coordx - self.star[0]) + self.vely * (self.coordy - self.star[1])  + self.velz + (self.coordz - self.star[2])
        
        
        self.Omega_K_10AU = (self.GravityConstantInternal*self.m_star/10**3)**0.5; self.vel_10AU = self.Omega_K_10AU * 10
        self.Omega_dimless = self.Omega_K/self.Omega_K_10AU;    self.vel_dimless = self.vel_mag/self.vel_10AU
        
    def scatter_qty(self, face, qty = None):
        fig = plt.figure(figsize=(10,8))
        ax = plt.axes()
        if qty is "None":
            clr = self.coordx * 0.0
        else:
            clr = getattr(self, qty)
        if face is "xy":
            res1 = self.coordx; res2 = self.coordy
            ax.set_xlabel('x');    ax.set_ylabel('y');
        elif face is "yz":
            res1 = self.coordy; res2 = self.coordz
            ax.set_xlabel('y');    ax.set_ylabel('z');
        elif face is "xz":
            res1 = self.coordx; res2 = self.coordz
            ax.set_xlabel('x');    ax.set_ylabel('z');
        elif face is "rz":
            res1 = self.rad; res2 = self.height
            ax.set_xlabel('r');    ax.set_ylabel('z');
        if face is "3d":
            ax = fig.add_subplot(111, projection='3d')
            #ax.view_init(elev=15., azim=15.)
            ax.scatter(self.coordx, self.coordy, self.coordz, c =  clr)
            plt.xlabel('x');          plt.ylabel("y");        plt.xlabel('z')
        else:
            ax.scatter(res1, res2, c = clr)
            xmax = np.max(res1); xmin = np.min(res1); ymax = np.max(res2); ymin = np.min(res2)
            ax.set_xlim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin));         ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
        plt.title(qty+" of t = %.1f"%(self.time)+" s (snapshot"+self.fname+") (Max = %.5f)"%(np.max(clr)))
        #ax.xaxis.set_minor_locator(MultipleLocator(5))
        norm = Normalize(0.0, 1.0);        colormap = cm.ScalarMappable(norm, 'gnuplot')
        fig.colorbar(colormap, label = qty)
        plt.show()
        #fig.savefig("./images/"+quantity+str(i)+'.png', dpi = 400)
        plt.close()

    def qty_qty(self, qty1, qty2):
        plt.figure(figsize=(10,10))
        plt.scatter(getattr(self, qty1), getattr(self, qty2))
        plt.title(qty2+" vs "+qty1+" of t = %.1f"%(self.time) + "s (snapshot "+self.fname+")")
        #ax.xaxis.set_minor_locator(MultipleLocator(5))
        plt.xlabel(qty1);          plt.ylabel(qty2)
        plt.show()
        plt.savefig("./images/"+qty1+"_"+qty2+"_"+self.fname+'.png', dpi = 400)
        plt.close()
        
    def grid_qty(self, qty, what):
        n = 200; m = 200;       eps_rad = 21/n;  eps_phi = 2*np.pi/m
        grid_rad = np.arange(0.1, 30, eps_rad);  grid_phi = np.arange(-np.pi, np.pi, eps_phi)
        grid = (np.array(np.meshgrid(grid_rad, grid_phi)).T.reshape(-1,2))
        grid_x = np.multiply(grid[:,0], np.cos(grid[:,1])).ravel();  grid_y = np.multiply(grid[:,0], np.sin(grid[:,1])).ravel()
        
        def Sigma(a, b):
            arr = self.height[((self.rad - a <= eps_rad/2)&(a - self.rad< eps_rad/2)) & ((self.phi-b <= eps_phi/2)&(b-self.phi< eps_phi/2))]
            if len(arr) == 0:
                return 0
            return len(arr)*(0.1/250000)/(np.pi*((a + eps_rad/2)**2-(a - eps_rad/2)**2)*(eps_phi/(2*np.pi)))
        
        def grid_qty(qty, a, b):   
            arr1 = getattr(self, qty)
            arr = arr1[((self.rad - a <= eps_rad/2)&(a - self.rad< eps_rad/2)) & ((self.phi-b <= eps_phi/2)&(b-self.phi< eps_phi/2))]
            if len(arr) == 0:
                return 0
            return np.average(arr)
              
        if qty is "surfDen":
            res = [Sigma(i, j) for i, j in grid]        
        elif qty is "Toomre":
            def Q_Toomre(a, b, what = None):
                if what is "exp":
                    x = self.Sigma_0/a
                else:
                    x = Sigma(a, b)
                if x == 0:
                    return 0
                Tmp = grid_qty("temp",a, b)
                c_s = (self.k_B/self.mean_molecular_weight * Tmp)**0.5
                Omega_K = (self.GravityConstantInternal*self.m_star/a**3)**0.5 
                return((self.c_s * self.Omega_K/(np.pi*self.GravityConstantInternal))/x)
            res = [Q_Toomre(i, j) for i, j in grid]
        else:
            res = [grid_qty(qty, i, j) for i, j in grid]

        def Sigma_init(r):
            res = self.Sigma_0 * r**self.power
            res[(r > self.r_o) | (r < self.r_i)] = 0
            return res

        plt.figure(figsize = (10, 8))
        if what is "radial":
            res_avg = [sum(res[i:i+m])/m for i in range(0, len(res), m)]
            plt.plot(grid_rad, res_avg, "rx", label = "averaged by radius")
            plt.scatter(grid[:,0], res, s = 5.0, label = "measured by cells")
            if qty is "surfDen":
                plt.plot(grid_rad, Sigma_init(grid_rad), 'y', linewidth = 3, label = "expected initial condition")
            xmax = 25; xmin = 4; ymax = np.max(res); ymin = np.min(res)
            plt.xlim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin));         plt.ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
            plt.legend()
            plt.xlabel("radius");          plt.ylabel("Surface Density (M_sun / AU^2)")
        if what is "xy":
            res = np.array(res)
            plt.scatter(grid_x, grid_y, s = 2.0, c = res)
            plt.colorbar()
        plt.title(qty+" of t = %.1f"%(self.time)+" s (snapshot"+self.fname+") ")
        #ax.xaxis.set_minor_locator(MultipleLocator(5))
        plt.show()
        #plt.savefig("./images/"+qty1+"_"+qty2+"_"+self.fname+'.png', dpi = 400)
        plt.close()

        
    def slice_qty(self, face, qty = None, cut = 0):
        fig = plt.figure(figsize=(10,8))
        ax = plt.axes()
        if qty is "None":
            clr = self.coordx * 0.0
        else:
            clr = getattr(self, qty)
        if face is "xy":
            res1 = self.coordx; res2 = self.coordy; res3 = self.coordz
            ax.set_xlabel('x');    ax.set_ylabel('y');
        elif face is "yz":
            res1 = self.coordy; res2 = self.coordz; res3 = self.coordx
            ax.set_xlabel('y');    ax.set_ylabel('z');
        elif face is "xz":
            res1 = self.coordx; res2 = self.coordz; res3 = self.coordy
            ax.set_xlabel('x');    ax.set_ylabel('z');
        elif face is "rz":
            res1 = self.rad; res2 = self.height; res3 = self.phi
            ax.set_xlabel('r');    ax.set_ylabel('z');
        eps = 1e-1
        plotx = res1[(np.abs(res3 - cut) < eps)]; ploty = res2[(np.abs(res3 - cut) < eps)]
        ax.scatter(plotx, ploty, c = clr[(np.abs(res3 - cut) < eps)])
        xmax = np.max(plotx); xmin = np.min(plotx); ymax = np.max(ploty); ymin = np.min(ploty)
        ax.set_xlim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin));         ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
        plt.title(qty+" of t = %.1f"%(self.time)+" s (snapshot"+self.fname+")  (Max = %.5f)"%(np.max(clr)))
        #ax.xaxis.set_minor_locator(MultipleLocator(5))
        norm = Normalize(0.0, 1.0);        colormap = cm.ScalarMappable(norm, 'gnuplot')
        fig.colorbar(colormap, label = qty)
        plt.show()
        #fig.savefig("./images/"+quantity+str(i)+'.png', dpi = 400)
        plt.close()

        #fig, ax = plt.subplots()
        #CS = ax.contourf(x, y, gridres, levels = 20, cmap=plt.cm.viridis)
        #ax.xaxis.set_minor_locator(MultipleLocator(5))   
        #plt.colorbar(CS, label = qty)
        #plt.title(qty+"of snapshot "+self.fname)
        #plt.legend()
        #plt.xlabel('x');       plt.ylabel('y')
        #plt.show()
        #print(plt.cm.cmap_d.keys())
        #matplotlib.colors.Normalize(vmin=-1.,vmax=1.)
        #fig.savefig("./images/"+quantity+str(i)+'.png', dpi = 400)
        #plt.close()
        

def make_video(folder, frame, qty1, qty2):
    for i in range(frame):
        current = plot_single(folder, i)
        current.qty_qty(qty1, qty2)
    os.system("ffmpeg -f image2 -r 1 -i ./images/*_%03d.png -vcodec mpeg4 -y test.mp4")
        

        
def slides(folder):
    frame = len(fnmatch.filter(os.listdir('/u/yali/'+folder+'/test/output'), '*.hdf5'))

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)
    
    f = h5py.File("/u/yali/"+folder+"/test/output/snapshot_000.hdf5", 'r') # Read-only access to the file
    coord0 = f['PartType0']['Coordinates'][:,:]
    l, = plt.plot(coord0[:,0], coord0[:,1], linestyle="", marker='.', alpha=0.6)
    #l, = plt.scatter(coordx0[0], coordy0[0], s = 0.5)
    
    axcolor = 'lightgoldenrodyellow'
    axt = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    
    st = Slider(axt, 't', 5, frame, valinit=0, valstep=1)
    
    def update(val):
        fname = str(format(i, '03d'))
        f2 = h5py.File("/u/yali/"+folder+"/test/output/snapshot_"+fname+".hdf5", 'r') # Read-only access to the file
        coord = f2['PartType0']['Coordinates'][:,:]
        l.set_data(coord[:,0],coord[:,1])
        fig.canvas.draw_idle()
    
    st.on_changed(update)
    
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    def reset(event):
        #print("Amplitude: ",  st.val)
        st.reset()
    button.on_clicked(reset)
    
    plt.show()        
        
        
        
        
        
        
