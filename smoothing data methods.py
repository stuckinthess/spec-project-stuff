#!/usr/bin/env python
# coding: utf-8

# In[1]:


import astropy


# In[2]:


from astropy.io import fits


# In[3]:


import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.pyplot import axes


# In[4]:


import numpy as np


# In[5]:


#import pylab as pl
import statsmodels.api as sm


# In[6]:


fitspec = fits.open('spec-1373-53063-0583.fits')
fitspec.info()


# In[7]:


data =  fitspec[1].data


# In[8]:


sky = data['sky']
sky = -sky/100


# In[9]:


x = data['loglam']
y = data['flux'] - data['model'] - sky
x = 10**x
figure(figsize=(28, 12), dpi=90)
plt.plot(x, y, color = 'green')


# In[10]:


#бегающее среднее
def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')


# In[11]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x, y)
y_av = movingaverage(y, 10)
plt.plot(x, y_av, "r")
plt.show()


# In[12]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x, y-y_av, "r")


# In[13]:


#lowess
figure(figsize=(28, 12), dpi=90)
plt.plot(x,y)
ylowess = data['flux'] - data['model']
xlowess = data['loglam']
xlowess = 10**xlowess
lowess = sm.nonparametric.lowess(ylowess, xlowess, it=0, frac=0.0097)
plt.plot(lowess[:, 0], lowess[:, 1])


# In[14]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x, y-lowess[:, 1], 'r')


# In[15]:


#savitsky-golay filter
from scipy.signal import savgol_filter
yhat = savgol_filter(y, 51, 10) # window size 51, polynomial order 3
figure(figsize=(28, 12), dpi=90)
plt.plot(x, y, color = 'indigo')
plt.plot(x,y)
plt.plot(x, yhat, "r")
plt.show()


# In[16]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x, y-yhat, "r")
plt.show()


# In[ ]:




