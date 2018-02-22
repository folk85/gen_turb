# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from scipy.fftpack import *

def plot_spectr(uin,vin,win):

  alpha = 1.339e0
  L = 1.0e-1
  sigma = 1.0e+1

  # x,y,z = np.genfromtxt('tests/spectr.dat',unpack=True)
  # x,y,z = np.genfromtxt('../hita/spectrum.dat',unpack=True)
  # x1,y1,z1 = np.genfromtxt('../hita/spectrum_32.dat',unpack=True)
  
  uvel,vvel,wvel = np.genfromtxt('./store.dat',unpack=True)
  nk = int(round(np.size(uvel)**(1./3.)))
  nel = nk
  ufft = fftn(uvel.reshape(nk,nk,nk))
  vfft = fftn(vvel.reshape(nk,nk,nk))
  wfft = fftn(wvel.reshape(nk,nk,nk))
  muu = ufft*np.conj(ufft) / nel**6
  mvv = vfft*np.conj(vfft) / nel**6
  mww = wfft*np.conj(wfft) / nel**6

  # calc std
  umean = np.array([np.mean(uvel),np.mean(vvel),np.mean(wvel)])
  std_i = np.array([np.std(uvel),np.std(vvel),np.std(wvel)])
  sigma = np.sqrt(np.sum(std_i[:]**2))
  print(std_i[0],np.sqrt(np.mean((uvel[:]-umean[0])**2)), sigma)
  dx = 10.
  k = np.arange(-nk//2,nk//2)*dx
  k = np.roll(k,nk//2)
  spectrum = np.zeros(nk)
  count = np.zeros(nk)
  # ?np.meshgrid(k,k,k)
  X,Y,Z = np.meshgrid(k,k,k)
  r = np.sqrt(X**2+Y**2+Z**2) #*dx
  # print(np.shape(r),r.min(),r.max(),k.max(),r[:,0,0])
  for i,ki in enumerate(k[:nk//2]):
      t = np.where((r<=ki+dx/2)&(r>ki-dx/2))
      spectrum[i] = np.sum(muu[t].real) + np.sum(mvv[t].real) + np.sum(mww[t].real)
      count[i] = np.size(t[0])  
      spectrum[i] *= 2.*np.pi*k[i]**2/dx**3/(count[i]+1.0e-30)

  font = {'family': 'Droid Sans',
          'weight': 'normal',
          'size': 12}
  m.rc('axes',linewidth=2)
  m.rc('font',**font)
  m.rc('lines',markeredgewidth=1.0)
  f,ax = plt.subplots()
  xf = np.linspace(np.log(k[1]/2),np.log(k[nk//2-1]*2.),100)
  xf = np.exp(xf)
  ax.loglog(xf,Ek(xf,alpha,L,sigma),c='g',lw=2)
  ax.loglog(k[:nk//2],spectrum[:nk//2],'bx-',lw=0.5,ms=8)
  # ax.loglog(x,y,'bx')
  # ax.loglog(x1,y1,'ro')
  ax.set_xlabel(u'$k, 1/м$',size='large')
  ax.set_ylabel(u'$E(k), м^3/с^2$',size='large')
  plt.grid()
  plt.tight_layout()
  plt.show()
  del(f)
  del(ax)
  plt.clf()

  Rij_x=(ufft*np.conj(ufft)) # compute velo. correlation tensor
  Rij_y=(vfft*np.conj(vfft))
  Rij_z=(wfft*np.conj(wfft))

  R1=ifftn(Rij_x)/np.std((uvel))**2/nel**3;
  R2=ifftn(Rij_y)/np.std((vvel))**2/nel**3;
  R3=ifftn(Rij_z)/np.std((wvel))**2/nel**3;
  
  NFFT=np.size(ufft,1)
  R11 = (R3[0,0,:]+R2[0,:,0]+R1[:,0,0])/3.
  # R11 = R11[:np.size(ufft)//2+1]
  R1_22 = (R1[0,:,0]+R3[0,:,0])/2.0e0
  R2_22 = (R2[:,0,0]+R3[:,0,0])/2.0e0
  R3_22 = (R1[0,0,:]+R2[0,0,:])/2.0e0

  R22 = (R1_22+R2_22+R3_22)/3.0e0
  # R22 = R22(1:size(u_fft)/2+1);
  Lx = 2.0*np.pi*1.0e-1
  r = np.linspace(0,Lx,NFFT)/(Lx/2);

  l11 = np.trapz(np.real(R11[:NFFT//2+1]),dx=r[1]-r[0])
  l22 = np.trapz(np.real(R22[:NFFT//2+1]),dx=r[1]-r[0])
  print("Integral Length Scale Longitudal: %g"%(l11))
  print("Integral Length Scale Tangent: %g"%(l22))

  f,ax = plt.subplots(1)
  ax.plot(r[:NFFT//2+1],R11[:NFFT//2+1],marker='>',mfc='w',lw=2,label=u'$R_{11}$')
  ax.plot(r[:NFFT//2+1],R22[:NFFT//2+1],marker='s',markerfacecolor='w',lw=2,label=u'$R_{22}$')
  ax.plot(r[:NFFT//2],np.exp(-r[:NFFT//2]/l11))
  ax.plot(r[:NFFT//2],1.e0+(1.0e0-R22[NFFT//2])*(np.exp(-r[:NFFT//2]/(l22-R22[NFFT//2]))-1.0e0))
  plt.legend()
  plt.tight_layout()
  ax.set_xlabel(u'$r$')
  ax.set_ylabel(u'$R_{11}, R_{22}$')
  plt.grid()
  plt.show()
  return [k[:nk//2],spectrum[:nk//2],r[:NFFT//2+1],R11[:NFFT//2+1],R22[:NFFT//2+1]]

def Ek(k,alpha=1.339,L=0.01,sigma=10.):
  tmp = (alpha * L * k) **2
  tmp =  sigma*sigma*L * tmp * tmp * 5.5e+1/ (27.0 * np.pi * (1.0 + tmp)**(1.7e+1/6.0e0))
  return tmp

if __name__ == '__main__':
  uvel,vvel,wvel = np.genfromtxt('./store.dat',unpack=True)
  [kx,e,rx,r11,r22] = plot_spectr(uvel,vvel,wvel)


