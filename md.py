
# coding: utf-8

# In[37]:

import math
import numpy as np
import random
import itertools
import scipy.spatial.distance as scipy_dist

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# In[38]:


Pi=3.1415926535897932384
core_r=0.500
density=0.140

period=4
clas="fcc"
place=[]
xi=xj=xk=0
yi=yj=yk=0
zi=zj=zk=0
xa=ya=za=0
lx=1.0
ly=1.0
lz=1.0

for i in range(period):
    xi=0
    yi=0
    zi=lz*i
    for j in range(period):
        xj=0
        yj=ly*j
        zj=0
        for k in range(period):
            xk=lx*k
            yk=0
            zk=0
            for b in range(2):
                if b==0:
                    xb=lx/2
                    yb=ly/2
                    zb=0.0
                elif b==1:
                    xb=lx/2
                    yb=0.0
                    zb=lz/2                 
                elif b==2:
                    xb=0.0
                    yb=ly/2
                    zb=lz/2
                xijk=xi+xj+xk+xb
                yijk=yi+yj+yk+yb
                zijk=zi+zj+zk+zb
                place.append([xijk,yijk,zijk])

Lx=xijk+lx/2
Ly=yijk+ly
Lz=zijk+lz/2
Lv=Lx*Ly*Lz
print "\t{0}\t{1}\t{2}\t{3} -> ".format(Lx,Ly,Lz,Lv)

place=np.array(place)
num = place.size / 3
v=core_r*core_r*core_r*4.0*math.pi*float(num)/(3.0*density)
ratio3=v / Lv
ratio=math.pow(ratio3,(1.0/3.0))

Lx*=ratio
Ly*=ratio
Lz*=ratio
Lv=Lx*Ly*Lz
print "\t{0}\t{1}\t{2}\t{3}".format(Lx,Ly,Lz,Lv)
place[:,0]*=ratio
place[:,1]*=ratio
place[:,2]*=ratio




# In[40]:

class md:
    def __init__(self):
        print "num : {0}\tdens : {1}".format(num,density)
        atom_num=num
        step=1000
        self.step=step
        self.atom_num=atom_num
        self.dt=0.01
        self.xyz=np.zeros((step,atom_num,3))
        self.vxyz=np.zeros((step,atom_num,3))
        self.fxyz=np.zeros((step,atom_num,3))
        self.dist=scipy_dist.euclidean
        
    def init_set(self):
        self.xyz[0,:,:]=place[:,:]
        self.vxyz[0,:,:]=2*np.random.rand(self.atom_num,3)-1
        K=0
        for i in xrange(num):
            K+=(self.vxyz[0,i,0]**2+self.vxyz[0,i,1]**2+self.vxyz[0,i,2]**2)/2
        Tb=K*num*(1/(1.38*10**(-23)))*2/3
        self.vxyz[0,:,:]=self.vxyz[0,:,:]*math.sqrt(1.0/Tb)
        K=0
        for i in xrange(num):
            K+=(self.vxyz[0,i,0]**2+self.vxyz[0,i,1]**2+self.vxyz[0,i,2]**2)/2
        Ta=K*num*(1/(1.38*10**(-23)))*2/3
        print "T = {0} -> {1}".format(Tb,Ta)
        
        for i in xrange(self.atom_num):
            self.fxyz[0,i,:]=self.force(self.xyz[0,:,:],i)
            
    def start(self):
        for i in xrange(0,self.step-1):
            self.verlet(i)


    def verlet(self,istep):
        for i in xrange(0,self.atom_num):
#             print "before\t{0}".format(self.xyz[istep,i,:])
            self.xyz[istep+1,i,:]=self.xyz[istep,i,:]+self.dt*self.vxyz[istep,i,:]+self.dt**2*self.fxyz[istep,i,:]/2
            self.pbc(istep)
            self.fxyz[istep+1,i,:]=self.fxyz[istep,i,:]+self.force(self.xyz[istep,:,:],i)
            self.vxyz[istep+1,i,:]=self.vxyz[istep,i,:]+(1/2)*self.dt*(self.fxyz[istep+1,i,:]-self.fxyz[istep,i,:])

    def pbc(self,step):
        for i in xrange(self.atom_num):
            if abs(self.xyz[step,i,0]) > Lx:
                self.xyz[step,i,0] = Lx - self.xyz[step,i,0]
            if abs(self.xyz[step,i,1]) > Ly:
                self.xyz[step,i,1] = Ly - self.xyz[step,i,1]
            if abs(self.xyz[step,i,2]) > Lz:
                self.xyz[step,i,2] = Lz - self.xyz[step,i,2]
            if abs(self.xyz[step,i,0]) < 0:
                self.xyz[step,i,0] = Lx + self.xyz[step,i,0]
            if abs(self.xyz[step,i,1]) < 0:
                self.xyz[step,i,1] = Ly + self.xyz[step,i,1]
            if abs(self.xyz[step,i,2]) < 0:
                self.xyz[step,i,2] = Lz + self.xyz[step,i,2]
        
    def force(self,atom_xyz,index):
        f=0
        tmpf=0
        sig=1.0
        A=12*sig**12
        B=6*sig**6
        eps=1.0
        for i,tmp_xyz in enumerate(atom_xyz):
            if index !=i:
#                 print "pair\t{0}\t{1}".format(index,i)
#                 print "xyz\t{0}\t{1}".format(atom_xyz[i],atom_xyz[index])
                dr=self.dist(atom_xyz[index],atom_xyz[i])
                dxyz=atom_xyz[index]-atom_xyz[i]
                if dr == 0.0:
                    dr=0.0001
                    dxyz=[0.0001,0.0001,0.0001]
                    dxyz=np.array(dxyz)
                k=4*eps*(A*(1/dr)**13-B*(1/dr)**7)
                tmpf=k*dxyz
                f=tmpf+f
#                 print k
#                 print tmpf
#                 dr=self.dist(atom_xyz[index],atom_xyz[i])
#                 u+=4*eps*(A*(1/dr)**12-B*(1/dr)**6)
#                 self.fxyz[i]=4*eps*(A*(1/dr)**13-B*(1/dr)**7)
        return f

    def show(self):
        for istep in xrange(self.step):
            for i in xrange(self.atom_num):
#                 print "R\t{0}\tF\t{1}".format(self.xyz[istep,i,:],self.fxyz[istep,i,:])
                print self.xyz[istep,:,:]
            print "\n"
    
    def plot(self):
        for istep in xrange(self.step):
            fig = plt.figure()
            ax = Axes3D(fig)
            box_size=2
            ax.set_aspect('equal')
            ax.set_xlim3d(-box_size, box_size)
            ax.set_ylim3d(-box_size, box_size)
            ax.set_zlim3d(-box_size, box_size)
            ax.view_init(9, 45)
            for i in xrange(self.atom_num):
                ax.scatter(self.xyz[istep,i,0], self.xyz[istep,i,1], self.xyz[istep,i,2], s=5.0,alpha=1.0)
                # c=clr[tag], edgecolor=clr[tag], 
            plt.savefig("./img/%05d.png" % istep)
        
    def plot2(self):
        fig = plt.figure()
        ax = Axes3D(fig)
        box_size=10

        trace1 = go.Scatter3d(
            x=self.xyz[:,0,0], 
            y=self.xyz[:,0,1],
            z=self.xyz[:,0,2],
            mode='markers',
            marker=dict(
                sizemode='diameter',
                colorscale = 'Portland',
                line=dict(color='rgb(255, 255, 255)'),
                opacity=0.9,
                size=5 
            )
        )

        data=[trace1]
        layout=dict(height=540, width=960, title='place')
        fig=dict(data=data, layout=layout)
        offline.iplot(fig, show_link=False)
        
    def save(self):
#         print self.xyz[self.step-1,:,:]
        for i in xrange(self.step-1):
            np.savetxt('./img/place{0}n{1:03d}.dat'.format(int(num),i),self.xyz[i,:,:],delimiter='\t',fmt='%.11f')
            



# In[41]:

if __name__=="__main__":
    spawn=md()
    spawn.init_set()
    spawn.start()
#     spawn.show()
#     spawn.plot()
    spawn.save()

