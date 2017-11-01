# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt

def plot_spectr():

  alpha = 1.339e0
  L = 7.0e-2
  sigma = 1.0e+1

  x,y,z = np.genfromtxt('tests/spectr.dat',unpack=True)


  font = {'family': 'Droid Sans',
          'weight': 'normal',
          'size': 12}
  m.rc('axes',linewidth=2)
  m.rc('font',**font)
  m.rc('lines',markeredgewidth=1.0)
  f,ax = plt.subplots()
  xf = np.linspace(x.min()/2,x.max(),100)
  ax.loglog(xf,Ek(xf,alpha,L,sigma),c='g',lw=2)
  ax.loglog(x,y,'bx')
  ax.set_xlabel(u'$k, 1/м$',size='large')
  ax.set_ylabel(u'$E(k), м^3/с^2$',size='large')
  plt.tight_layout()
  plt.show()
  return

def Ek(k,alpha=1.339,L=0.01,sigma=10.):
  tmp = (alpha * L * k) **2
  tmp =  sigma*sigma*L * tmp * tmp * 5.5e+1/ (27.0 * np.pi * (1.0 + tmp)**(1.7e+1/6.0e0))
  return tmp

if __name__ == '__main__':
  plot_spectr()


