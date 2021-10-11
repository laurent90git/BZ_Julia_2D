#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 18:08:33 2021

Ce script charge les données d'une simulation Fortran de BZ1D 2 ou 3 équations

@author: lfrancoi
"""
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def data_loader(datafile, bPrint=True):
  """ loads the solution history into Numpy arrays """
 
  if bPrint:
    print('----------- DATA LOADER -----------')

  rawdata = pd.read_csv(datafile, decimal='.', encoding='utf-8', dtype=np.float64, engine='c',
                    low_memory=True, skiprows=0, skipinitialspace=True, delimiter=' ', header=None).values
  
  nlines = rawdata.shape[0]
  nx_ny = np.where( rawdata[:,0]>rawdata[0,0] )[0][0]
  nx = np.where( rawdata[:,2]>rawdata[0,2] )[0][0] # x is enumerated first
  ny = nx_ny/nx
  assert int(nx)==nx
  assert int(ny)==ny
  
  nx, ny = int(nx), int(ny)
  
  nt = nlines / nx_ny
  assert int(nt)==nt
  nt = int(nt)
  
  neq = rawdata.shape[1]-3
  
  reshdata = rawdata.reshape((nt,ny,nx,neq+3))
  data_dict = {"t": reshdata[:,0,0,0],
          "xx": reshdata[0,:,:,1],
          "yy": reshdata[0,:,:,2],
          "nx": nx,
          "ny": ny,
          }
  
  # reshdata = rawdata.reshape((neq+3, nt, ny, nx), order="F")
  # data_dict = {"t": reshdata[0,:,0,0],
  #     "xx": reshdata[1,0,:,:],
  #     "yy": reshdata[2,0,:,:],
  #     "nx": nx,
  #     "ny": ny,
  #     }
  #data = rawdata.reshape((neq+3, nx,ny,nt))

  varnames = ('a','b','c')
  for i in range(neq):
      data_dict[varnames[i]] = reshdata[:,:,:,3+i]
      #data_dict[varnames[i]] = reshdata[3+i,:,:,:]

  if bPrint:
    print('\tstats: {} vars, {} mesh points ({}x{}), {} time steps)'.format(neq, nx*ny, nx, ny, nt))
  return data_dict

  

if __name__=='__main__':

  #%% Chargement de la solution 2D au cours du temps
  datafile = r'C:\Users\Laurent\Desktop\save_test_singlestep.txt'
  data = data_loader(datafile, bPrint=True)

  #%% Tracé 2D de l'évolution des champs
  bExport = False # export de PNG
  
  a_min, a_max = np.min(data['a']), np.max(data['a'])
  b_min, b_max = np.min(data['c']), np.max(data['b'])
  c_min, c_max = np.min(data['c']), np.max(data['c'])
  # a_min, a_max = 0, 1000
  # b_min, b_max = 0, 1
  # c_min,  c_max = 0, 0.2
  nlevels=100 # nombre de niveaux pour le contour
  
  if bExport:
      try:
        os.mkdir('export')
      except:
        pass

  npts = 1 # nb de pts équidistants en temps
  for iii in [int((len(data['t'])-1)/npts)*k for k in range(npts+1)]: #range(len(data['t'])):
    fig,ax = plt.subplots(3,1,sharex=True,sharey=True, figsize=np.array((5,10)))
    # for ix, var in ['a', 'b', 'c']
    ix=0
    levels = np.linspace(a_min, a_max, nlevels)
    cs = ax[ix].contourf(data['xx'], data['yy'], data['a'][iii,:,:], antialiased=True, cmap='jet', levels=levels)
    fig.colorbar(cs, ax=ax[ix], shrink=0.9, label='a', extend='both')

    ix+=1
    levels = np.linspace(b_min, b_max, nlevels)
    cs = ax[ix].contourf(data['xx'], data['yy'], data['b'][iii,:,:], antialiased=True, cmap='jet', levels=levels)
    fig.colorbar(cs, ax=ax[ix], shrink=0.9, label='b', extend='both')

    ix+=1
    levels = np.linspace(c_min, c_max, nlevels)
    cs = ax[ix].contourf(data['xx'], data['yy'], data['c'][iii,:,:], antialiased=True, cmap='jet', levels=levels)
    fig.colorbar(cs, ax=ax[ix], shrink=0.9, label='c', extend='both')

    for a in ax:
      a.set_xlabel('x')
      a.set_ylabel('y')
      a.axis("square")
      a.grid()
    fig.suptitle('t={:.3e}'.format(data['t'][iii]))
    if bExport:
       fig.savefig('export/{}.png'.format(iii), dpi=200)
    plt.show()
    
  #%% Create animation
  print('saving animation...')
  # os.system('convert -delay 15 -loop 0 export/*.png animation.gif')
  print('done')


