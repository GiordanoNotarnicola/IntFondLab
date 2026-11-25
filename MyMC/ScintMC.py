import numpy as np
import matplotlib.pyplot as plt
import time

class Detector:
    def __init__(self, xmin, xmax, ymin, ymax, z0):
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.z0=z0

    def hit(self, ray):
        x,y=ray.pos_at_z(self.z0)
        inx=(self.xmin<x) & (x<self.xmax)
        iny=(self.ymin<y) & (y<self.ymax)
        return inx & iny

class Cristal:
    def __init__(self, xmin, xmax, ymin, ymax, zup, zdown):
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.zup=zup
        self.zdown=zdown

    def intersect(self, ray):
        x1=(self.xmin-ray.x0)/ray.tx
        x2=(self.xmax-ray.x0)/ray.tx
        y1=(self.ymin-ray.y0)/ray.ty
        y2=(self.ymax-ray.y0)/ray.ty
        xin=np.maximum(x1, x2)
        xout=np.minimum(x1,x2)
        yin=np.maximum(y1,y2)
        yout=np.minimum(y1, y2)
        z_enter = np.minimum(np.minimum(xin, yin), self.zup)
        z_exit = np.maximum(np.maximum(xout, yout), self.zdown)
        return z_enter, z_exit  

    def hit(self, ray):
        z_enter, z_exit=self.intersect(ray)
        return z_enter>z_exit

    def length(self, ray):
        z_enter, z_exit=self.intersect(ray)
        dz=np.maximum(0., z_enter- z_exit)
        return dz/ray.costh
              
class Ray:

    def __init__(self, x0, y0, costh, phi):
        self.x0=x0
        self.y0=y0
        self.costh=costh
        self.phi=phi
        tanth=np.sqrt(1-self.costh**2)/self.costh
        self.tx=tanth*np.cos(self.phi)
        self.ty=tanth*np.sin(self.phi)

    def pos_at_z(self, z):
        x=self.x0+z*self.tx
        y=self.y0+z*self.ty
        return x,y
    
#sezione geometria detector da modificare a piacimento
expo=200000
sizex=0.53
sizey=0.27
detcen=Detector(0, 0.51, -0.09, 0.31, -0.08)
cri1=Cristal(0.09, 0.47, 0.01, 0.115, -0.105, -0.21)
cri2=Cristal(0.09, 0.47, 0.115, 0.22, -0.105, -0.21)
cri3=Cristal(0.09, 0.47, 0.01, 0.115, -0.21, -0.315)
cri4=Cristal(0.09, 0.47, 0.115, 0.22, -0.21, -0.315)
cri5=Cristal(0.09, 0.47, 0.01, 0.115, -0.315, -0.42)
cri6=Cristal(0.09, 0.47, 0.115, 0.22, -0.315, -0.42)
detinf=Detector(0, 0.51, -0.09, 0.31, -0.505)
n=1000000
doppie=0
triple=0

t0=time.time()
x0=sizex*np.random.random(n)
y0=sizey*np.random.random(n)
costh=(np.random.random(n))**(1/(1+expo))
phi=2*np.pi*np.random.random(n)
raggi=Ray(x0, y0, costh, phi)

hit1=cri1.hit(raggi)
hit2=cri3.hit(raggi)
doppie=np.sum(hit1)
triple=np.sum(hit1 & hit2)
t1=time.time()

print(f"efficienza è {(doppie/n):.2f}, {(doppie/triple):.2f}, ci ho messo {(t1-t0):.2f}")


