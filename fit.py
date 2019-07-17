
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plot
from scipy.optimize import curve_fit
import time


min_weight = -40
polynomial_degree = 2

def jacknife(vector):
  nsectors, nmeas = vector.shape
  nblocks = min(40, nmeas)

  nmeas = int(nmeas/nblocks)*nblocks
  vector_cut = vector[:,0:nmeas]

  vector_cut = vector_cut.reshape((nsectors, nblocks, int(nmeas/nblocks)))
  vector_cut = vector_cut.mean(2)

  idx = np.arange(nblocks)
  jkvector = np.array([ np.mean(vector_cut[:,idx!=i],1) for i in range(nblocks) ])
  return jkvector.T

def jackknife_mean(vector):
  return vector.mean(1)

def jackknife_sigma(vector):
  n = vector.shape[1]
  mean = jackknife_mean(vector)
  diff = vector.T-mean
  sigmasq = (n-1)/(n+0.0) * np.sum( diff**2, 0 )
  return np.sqrt(sigmasq)


def fit_function( x, *p ):
  r = np.polyval(p, x)
  return r


def plot_window(wl_f, center, width):
  mean = jackknife_mean(wl_f.T)
  mean = jackknife_sigma(wl_f.T)
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])

  window = np.logical_and(x > center-width-1, x < center+width+1)
  x = np.linspace( -width, width, 2*width+1)
  mean = mean[ window ]
  sigma = sigma[ window ]
  plot.errorbar( x, mean, sigma, fmt='o' , capsize=4 )


def fit_window( wl_f, center, width ):
  sigma = np.std(wl_f, axis=0)/np.sqrt((wl_f.shape[0]-1))
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])

  window = np.logical_and(x > center-width-50, x < center+width+50)
  window = np.logical_and(window, x >= 0)

  x = x[window]-center
  wl_f = wl_f[ :, window ]
  # Using a gaussian weight
  d = np.minimum( np.abs(x)/width, 10 ) 
  w = np.exp(d**2/2)
  sigma = sigma[ window ] * w

  par = []
  for m in range(wl_f.shape[0]):
    parameters, conv = curve_fit(fit_function, x, wl_f[m], [0]*(polynomial_degree+1), sigma=sigma)
    par.append(parameters)
  return par


def evaluate_fit(x, par):
  fit_evaluated = []
  for m in range(len(par)):
    fit_evaluated.append(fit_function(x, *par[m]))
  return np.array(fit_evaluated)


def average_fit(x, par):
  fit = evaluate_fit(x, par)
  mean = jackknife_mean( fit.T )
  return mean


def plot_fit(par, center, width):
  x = np.linspace( -width, width, 101)
  y = average_fit(x, par)
  plot.plot( x, y )


def read_data( datafilename ):
  weights = []
  with open(datafilename) as datafile:
    for line in datafile:
      if "SECTOR" in line:
        try:
          tag, sector, w, wf = line.split(' ')
          sector = int(sector)
          wf = float(wf)
          if sector >= len(weights):
            weights.append([wf])
          else:
            weights[sector].append(wf)
        except:
          pass

  nmeas = min( [len(m) for m in weights] )
  weights = [m[0:nmeas] for m in weights]

  wl_w = np.array(weights).astype(np.float)
  wl_w = jacknife( wl_w )

  wl_f = np.log(wl_w)
  wl_f[wl_f < min_weight-10] = min_weight-10
  weight_sums = np.sum(np.exp(wl_f), axis=0)
  energy_correction = np.log(weight_sums)
  wl_f = ( wl_f - energy_correction )

  return wl_f.T


def window_smooth( x, wl_f, width, eval = 0 ):
  if width > 0: 
    par = fit_window(wl_f, x, width)
    return evaluate_fit(0, par)
  else:
    return wl_f[:,x]


def plot_window_fit( datafilename, center, width ):
  wl_f = read_data( datafilename )
  plot_window(wl_f, center, width)

  par = fit_window(wl_f, center, width)
  plot_fit(par, center, width)
  
  plot.show()


def plot_smoothing( datafilename, width, max ):
  wl_f = read_data( datafilename )
  wl_f = wl_f[:,0:max]

  mean = jackknife_mean(wl_f.T)
  sigma = jackknife_sigma(wl_f.T)
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])
  plot.errorbar( x, mean, sigma, fmt='o' , capsize=4 )

  max_i = 0
  for i in range(x.shape[0]):
    if mean[i] > min_weight:
      max_i = i+1
      
  print("Max sector used", max_i)
  wl_f = wl_f[:,:max_i]
  x = x[:max_i]

  wl_f_fit = []
  intermediate_points = int(100/x.shape[0]+1)
  print(intermediate_points)
  x = np.linspace(x.min(), x.max(), x.shape[0]*intermediate_points)
  for i in range(x.shape[0]):
    point = x[i]-x.min()
    value = np.array([window_smooth(point, wl_f, width)])
    value = jackknife_mean(value)
    wl_f_fit.append(value)
  wl_f_fit = np.array(wl_f_fit)

  plot.plot( x, wl_f_fit )

  plot.xlabel('Negative loops')
  plot.ylabel('F')
  plot.ylim(int(1.1*min_weight), int(0.5*wl_f_fit.max()))
  plot.xlim(-1, int(1.1*x.max()))

  plot.show()


def average_sign( datafilename, width, max, print_weights = False ):
  wl_f = read_data( datafilename )
  wl_f = wl_f[:,:max]

  mean = jackknife_mean(wl_f.T)
  sigma = jackknife_sigma(wl_f.T)
  x = np.linspace(0, mean.shape[0]-1, mean.shape[0])
  window = mean > min_weight
  print(window.sum(), "sectors above minimal weight")

  mean = mean[window]
  sigma = sigma[window]
  wl_f = wl_f[:,window]
  x = x[window]

  weights = []
  free_energies = []
  for i in range(x.shape[0]):
    free_energy = window_smooth( i, wl_f, width )
    free_energies.append(free_energy)
  diffs = np.abs(jackknife_mean(np.array(free_energies)) - mean)/sigma
  print("mean diff:", diffs.mean())
  print("max diff:", diffs.max())

  weights = np.exp(np.array(free_energies).T)
  sign = -(x%2-0.5)*2
  weights = weights*sign
  
  sign = np.sum(weights, axis=1)
  mean = jackknife_mean(np.array([sign]))[0]
  sigma = jackknife_sigma(np.array([sign]))[0]
  if print_weights :
    print(jackknife_mean(weights.T))
  return [mean, sigma]
  


if __name__ == "__main__":
  if len(sys.argv) > 3 :
    datafilename = sys.argv[1]
    width = float(sys.argv[2])
    max_sector = int(sys.argv[3])
    min_weight = int(sys.argv[4])
    if len(sys.argv) > 5:
      do_plot = (sys.argv[5] == 'plot')
      print_weights = (sys.argv[5] == 'weights')
    else :
      do_plot = False
      print_weights = False
  else:
    print("usage: fit.py filename window_width max_sector ")

  if do_plot:
    plot_smoothing(datafilename, width, max_sector)
  else:
    mean, sigma = average_sign(datafilename, width, max_sector, print_weights)
    print(mean, sigma)

  
