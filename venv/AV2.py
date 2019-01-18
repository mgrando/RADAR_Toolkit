# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 08:13:01 2015

@author: mbg
"""
#from parserbox import cstfssparser
from numba import jit
import numpy as np
from astropy import constants as const
from astropy import units as ureg
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import interpolate
import gc

class proposal_array():

    def __init__(self):
        self.freq = 500e6*ureg.Hz
        self.lambda0 = (const.c/self.freq).decompose()        
        self.Dcell = 0
        self.EIRP_passive = 0
        self.EIRP_active = 0
        self.Nele = 0
        self.D0 = 0 
        self.dx = 0
        self.dy = 0
    
    def Dcell_calc(self,lambda0,dx,dy,theta):
        self.Dcell = (np.pi/(lambda0**2))*(dx*dy)*np.cos(theta.ito(ureg.radians))

    def D0_calc(self):
        self.D0 = self.Nele*self.Dcell
    
    def EIRP_passive_calc(self,loss_e,Pin,Gamma):
        self.EIRP_passive = self.Nele*Pin*self.Dcell*(1-(np.abs(Gamma))**2)
        
    def EIRP_active_calc(self,loss_e,Pmodule,Gamma):
        self.EIRP_active = (self.Nele**2)*Pmodule*self.Dcell*(1-(np.abs(Gamma)**2))
        
        
        
class planar_array():
    
    def __init__(self):
        pass
    
    def to_meters(self, length):
        
        try:
            self.length = length.to(ureg.meter)
        except:
            self.length = length*ureg.meter
            
        return self.length
        
    def to_radians(self, angle_p):
        
        try:
            self.angle_p = angle_p.to(ureg.radian)
        except:
            self.angle_p = (angle_p*ureg.degree).to(ureg.radian)
            
        return self.angle_p

    def grid_array(self, dx, dy, rows, columns):
        
        self.dx = self.to_meters(dx)
        self.dy = self.to_meters(dy)            
            
        self.rows = int(rows)
        self.columns = int(columns)
        
        self.element_num = self.rows*self.columns
        self.x = np.zeros((self.rows,self.columns))*ureg.meter
        self.y = np.zeros((self.rows,self.columns))*ureg.meter
        

        for i in range(0,self.rows):
            for ii in range(0,self.columns):
                pass
                self.x[i,ii]=self.dx*ii
                self.y[i,ii]=self.dy*i
                
        maxy = self.y.max()
        miny = self.y.min()
        
        maxx = self.x.max()
        minx = self.x.min()
        
        delta_x = (maxx-minx)/2.0
        delta_y = (maxy-miny)/2.0

        self.x = self.x-delta_x
        self.y = self.y-delta_y

    def wg_slot_array(self,dx,dy,wg_offset,rows,columns):
        
        self.grid_array(dx, dy, rows, columns)
        self.wg_offset = self.to_meters(wg_offset)
        
        self.y[:,::2]=self.y[:,::2]+self.wg_offset
        self.y[:,1::2]=self.y[:,1::2]-self.wg_offset

        self.x = self.to_meters(self.x)
        self.y = self.to_meters(self.y)
        

    def hex_array(self,spacing,rows,columns,angle):
        
        self.dx = self.to_meters(spacing)     
        self.angle = self.to_radians(angle)
            
        self.dy = self.dx*np.sin(self.angle.to(ureg.radian))
            
        self.rows = rows
        self.columns = columns
        
        self.element_num = self.rows*self.columns
        self.x = np.zeros((self.rows,self.columns))*ureg.meter
        self.y = np.zeros((self.rows,self.columns))*ureg.meter
        
        for i in range(0,self.rows):
            for ii in range(0,self.columns):
                if i%2 == 0:
                    self.x[i,ii] = self.dx*ii
                    self.y[i,ii] = self.dy*i
                else:
                    self.x[i,ii] = self.dx*ii+(self.dx/2.0)
                    self.y[i,ii] = self.dy*i
                    
        maxy = self.y.max()
        miny = self.y.min()
        
        maxx = self.x.max()
        minx = self.x.min()
        
        delta_x = (maxx-minx)/2.0
        delta_y = (maxy-miny)/2.0

        self.x = self.x-delta_x
        self.y = self.y-delta_y
        
class array_params():
    
    def __init__(self,f,r_theta = 0.75,r_phi = 0.75):
        
        self.frequency = f*(1/ureg.second)
        self.lam = const.c/self.frequency
        self.k = np.pi*2*ureg.radian/(self.lam)     
        self.elements = 0
        
        self.steer_theta = 0 * ureg.degree
        self.steer_phi = 0 * ureg.degree
        
        self.sphere_flag = 1
        
        self.set_res_theta(r_theta)
        self.set_res_phi(r_phi)
        
    def to_radians(self, angle_p):
        
        try:
            self.angle_p = angle_p.to(ureg.radian)
        except:
            self.angle_p = (angle_p*ureg.degree).to(ureg.radian)
            
        return self.angle_p
        
    def set_res_theta(self,steps):
        
        self.res_theta = self.to_radians(steps)
        self.populate_thetas()        
        
    def set_res_phi(self,steps):
        
        self.res_phi = self.to_radians(steps)
        self.populate_phis()
        
    def populate_thetas(self):
        self.thetas = np.hstack([np.arange(0,np.pi,self.res_theta.to(ureg.radian).value)*ureg.radian,np.pi*ureg.radian])
    
    def populate_phis(self):
        self.phis = np.hstack([np.arange(0,2*np.pi,self.res_theta.to(ureg.radian).value)*ureg.radian,2*np.pi*ureg.radian]) 
        
    def set_position(self,positions):
        self.position_matrix = positions
    
    
    def set_alpha(self):
        self.alpha_x = (self.position_matrix[:,1]*np.sin(self.steer_theta.to(ureg.radian))*np.cos(self.steer_phi.to(ureg.radian)))
        self.alpha_y = (self.position_matrix[:,2]*np.sin(self.steer_theta.to(ureg.radian))*np.sin(self.steer_phi.to(ureg.radian)))
        self.alpha_z = (self.position_matrix[:,3]*np.sin(self.steer_theta.to(ureg.radian)))
        self.alpha = -self.k.value*(self.alpha_x.value+self.alpha_y.value+self.alpha_z.value) #technically this is k dot r
    
    
    def element_report(self):
        self.elem_rpt = pd.DataFrame(data=self.position_matrix,columns=['i','x','y','z'])
        self.set_alpha()
        self.elements,trash = np.shape(self.position_matrix)
        self.elem_rpt['alpha']= np.angle(np.exp(1j*self.alpha))
        
        ###This temporarily fills in the unit normals with zeros until this piece of the code works
        self.elem_rpt['nx'] = np.zeros((1,self.elements))[0]
        self.elem_rpt['ny'] = np.zeros((1,self.elements))[0]
        self.elem_rpt['nz'] = np.zeros((1,self.elements))[0]
        
        
    def element_pattern(self):
        
        try:
            self.TH, self.PH = np.meshgrid(self.thetas,self.phis)
            self.TH = self.TH.value
            self.PH = self.PH.value
            
        except:
            print('ERROR: You must first define theta and phi resolutions using "set_res_theta" and "set_res_phi')
            
        self.elements, trash = np.shape(self.position_matrix)
        self.THn,self.PHn = np.shape(self.TH)
        
        self.E = 1 #np.ones((self.elements,self.THn,self.PHn))*ureg.volt/ureg.meter

    @jit
    def AF_calc(self):
        
        self.elements, trash = np.shape(self.position_matrix)
        self.AF = np.zeros_like(self.TH)
        self.eta_temp = np.zeros_like(self.TH)   
        self.taper_total = (sum(self.position_matrix[:,0].value))
        
        for i in range(0,self.elements):
            self.eta_temp = 1j*self.k.value*(self.position_matrix[i,1].value*np.sin(self.TH)*np.cos(self.PH)+self.position_matrix[i,2].value*np.sin(self.TH)*np.sin(self.PH)+self.position_matrix[i,3].value*np.cos(self.TH))
            # self.AF = (self.AF+self.E[i,:,:]*(self.position_matrix[i,0].value/self.taper_total)*np.exp(1j*self.alpha[i])*np.exp(self.eta_temp))
            self.AF = (self.AF + 1*ureg.volt/ureg.meter * (self.position_matrix[i, 0].value / self.taper_total) * np.exp(1j * self.alpha[i]) * np.exp(self.eta_temp))

        self.F = self.AF
        
    def array_report(self):
        
        self.array_rpt = pd.DataFrame(columns=['Theta','Phi','Array_Patt','E_Theta','E_Phi'])
    
    @jit    
    def directivity(self):
        self.D_numerator = 4*np.pi*np.abs(self.F)**2       
        self.D_denominator = np.trapz(np.trapz((np.abs(self.F)**2)*np.sin(self.TH),axis=0,x = self.phis.value[np.newaxis].T[:,0]),x=self.thetas.value)
        self.D = self.D_numerator.value/self.D_denominator.value
        self.D_dB = 10*np.log10(self.D)



        
        

if __name__ == '__main__':
    x  = planar_array() 
    m = 5.37*ureg.mm
    offset = 0.28*ureg.mm
    
    x.wg_slot_array(m,m,offset,1,155)
    #x.hex_array(m,16,16,60)

    xs = x.x.ravel()[np.newaxis].T
    ys = x.y.ravel()[np.newaxis].T

    print(np.shape(xs))

    xm = []
    ym = []

    n = 50

    for i in range(0,n):
        xm.append(xs.value)
        ym.append(ys.value+(i*0.83/n)-(0.83/2))

    print(np.shape(xm))
    xs = np.ndarray.flatten(np.array(xm)) * ureg.meter
    ys = np.ndarray.flatten(np.array(ym)) * ureg.meter

    xs = xs[np.newaxis].T
    ys = ys[np.newaxis].T

    print(np.shape(xs))

    zs = np.ones((1,len(xs)))[np.newaxis].T[:,:,0]*ureg.meter

    print(np.shape(zs))
    i = np.ones((1,len(xs)))[np.newaxis].T[:,:,0]*ureg.volt
    position_matrix = np.concatenate((i,xs,ys,zs),1)
    g = array_params(35e9,r_theta = 0.5,r_phi = 0.5)
    g.steer_theta = 0*(np.pi/180)*ureg.radian
    g.steer_phi = 0*(np.pi/180)*ureg.radian
    g.set_position(position_matrix)
    g.element_report()
    g.element_pattern()
    print(g.elem_rpt)
    g.AF_calc()
    g.directivity()
#
# #
    s = proposal_array()
    p = g.AF.value
#
    plt.figure()

    print(x.x.unit, x.y.unit)
    plt.plot(xs,ys,'ko')
    plt.axis('equal')
#
    u = np.sin(g.TH)*np.cos(g.PH)
    v = np.sin(g.TH)*np.sin(g.PH)

    u_p = np.linspace(-1,1,500)
    v_p = np.linspace(-1,1,500)

    U,V = np.meshgrid(u_p,v_p)


    gc.collect()



    G_New = interpolate.griddata( (u.ravel(),v.ravel()),g.D_dB.ravel(), (U.ravel(), V.ravel()), method='linear')



    print(np.shape(U))

    indx,indy = np.shape(U)

    G_New = G_New.reshape((indx,indy))

    G_New[(np.sqrt(U ** 2 + V ** 2)) > 1] = np.nan

    plt.figure()
    plt.pcolor(U,V, G_New,vmin = g.D_dB.max()-50,vmax=g.D_dB.max(),cmap = cm.jet)
    plt.colorbar()
    plt.xlim([-.25,.25])
    plt.ylim([-.25,.25])
    plt.axis('equal')
    # #
    plt.show()