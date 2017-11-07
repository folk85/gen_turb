# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from scipy.fftpack import *

def plot_rng(arr):
  '''Plot gauss distribution values'''
  # get std
  dstd = np.std(arr) 
  # get mean value
  dmean = np.mean(arr)

  font = {'family': 'Droid Sans',
          'weight': 'normal',
          'size': 12}
  m.rc('axes',linewidth=2)
  m.rc('font',**font)
  m.rc('lines',markeredgewidth=1.0)
  
  f,ax = plt.subplots()
  p,bins,v = ax.hist(arr,bins=20,normed=True)
  #ax.hist(arr,bins=40)
  drng = rng(bins,dmean,dstd)
  plt.twinx()
  ax.plot(bins,drng,lw=3)

  ax.set_xlabel(u'Value',size='large')
  ax.set_ylabel(u'Distribution',size='large')
  plt.grid()
  plt.tight_layout()
  plt.show()
  plt.clf()

def rng(x,mean,std):
  return np.exp(-(x-mean)**2*5.0e-1/std**2) / np.sqrt(2.0e0*np.pi*std**2) 

if __name__ == '__main__':
  ar = np.genfromtxt('./tests/gauss.dat',usecols=[0])
  plot_rng(ar)