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
import time
import math


# In[5]:


fitspec = fits.open('spec-2483-53852-0232.fits')
fitspec.info()


# In[6]:


cols = fitspec[1].columns
cols.info()
cols.names


# In[7]:


data =  fitspec[1].data


# In[8]:


x = data['loglam']
y = data['flux']
x = 10**x
figure(figsize=(28, 12), dpi=90)
plt.plot(x, y, color = 'indigo')


# In[9]:


sky = data['sky']
sky = -sky/100


# In[11]:


x = data['loglam']
x = 10**x
y = data['flux'] - data['model']
y = y-sky
figure(figsize=(28, 12), dpi=90)
plt.plot(x, y, color = 'indigo')
#plt.plot(x, y-sky, color = 'yellow')


# In[12]:


x_l = data['loglam']
x_l = 10**x_l
y_l = data['flux'] - data['model'] -sky


# In[13]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x_l, y_l, color = 'indigo')


# In[14]:


#следующая функция и класс были написаны другим человеком


# In[15]:


def tricubic(x0):
    y0 = np.zeros_like(x0)
    idx = (x0 >= -1) & (x0 <= 1)
    y0[idx] = np.power(1.0 - np.power(np.abs(x0[idx]), 3), 3)
    return y0


# In[16]:


class Loess(object):

    @staticmethod
    def normalize_array(array):
        min_val = np.min(array)
        max_val = np.max(array)
        return (array - min_val) / (max_val - min_val), min_val, max_val

    def __init__(self, xx, yy, degree=1):
        self.n_xx, self.min_xx, self.max_xx = self.normalize_array(xx)
        self.n_yy, self.min_yy, self.max_yy = self.normalize_array(yy)
        self.degree = degree

    @staticmethod
    def get_min_range(distances, window):
        min_idx = np.argmin(distances)
        n = len(distances)
        if min_idx == 0:
            return np.arange(0, window)
        if min_idx == n-1:
            return np.arange(n - window, n)

        min_range = [min_idx]
        while len(min_range) < window:
            i0 = min_range[0]
            i1 = min_range[-1]
            if i0 == 0:
                min_range.append(i1 + 1)
            elif i1 == n-1:
                min_range.insert(0, i0 - 1)
            elif distances[i0-1] < distances[i1+1]:
                min_range.insert(0, i0 - 1)
            else:
                min_range.append(i1 + 1)
        return np.array(min_range)

    @staticmethod
    def get_weights(distances, min_range):
        max_distance = np.max(distances[min_range])
        weights = tricubic(distances[min_range] / max_distance)
        return weights

    def normalize_x(self, value):
        return (value - self.min_xx) / (self.max_xx - self.min_xx)

    def denormalize_y(self, value):
        return value * (self.max_yy - self.min_yy) + self.min_yy

    def estimate(self, x, window, use_matrix=False, degree=1):
        n_x = self.normalize_x(x)
        distances = np.abs(self.n_xx - n_x)
        min_range = self.get_min_range(distances, window)
        weights = self.get_weights(distances, min_range)

        if use_matrix or degree > 1:
            wm = np.multiply(np.eye(window), weights)
            xm = np.ones((window, degree + 1))

            xp = np.array([[math.pow(n_x, p)] for p in range(degree + 1)])
            for i in range(1, degree + 1):
                xm[:, i] = np.power(self.n_xx[min_range], i)

            ym = self.n_yy[min_range]
            xmt_wm = np.transpose(xm) @ wm
            beta = np.linalg.pinv(xmt_wm @ xm) @ xmt_wm @ ym
            y = (beta @ xp)[0]
        else:
            xx = self.n_xx[min_range]
            yy = self.n_yy[min_range]
            sum_weight = np.sum(weights)
            sum_weight_x = np.dot(xx, weights)
            sum_weight_y = np.dot(yy, weights)
            sum_weight_x2 = np.dot(np.multiply(xx, xx), weights)
            sum_weight_xy = np.dot(np.multiply(xx, yy), weights)

            mean_x = sum_weight_x / sum_weight
            mean_y = sum_weight_y / sum_weight

            b = (sum_weight_xy - mean_x * mean_y * sum_weight) /                 (sum_weight_x2 - mean_x * mean_x * sum_weight)
            a = mean_y - b * mean_x
            y = a + b * n_x
        return self.denormalize_y(y)


# In[17]:


loess = Loess(x_l, y_l)

for i in range(len(x_l)):
    y_l[i] = loess.estimate(x_l[i], window=30) #длина участка (имеет ли смысл делать итерации?)


# In[18]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x, y, color = 'indigo')
plt.plot(x_l, y_l, color = 'yellow') #то что дальше вычитаем


# In[19]:


#сглаженная фунция (один раз)
figure(figsize=(28, 12), dpi=90)
plt.plot(x, y-y_l, color = 'green')


# In[20]:


data['ivar']


# In[21]:


ivars = data['ivar']
sigmas = (1/ivars) #без корня?
#sigmas_ = 3*sigmas три сигмы? с ними выглядит неправильно


# In[22]:


sigmas #похоже, что здесь уже взвешенные значения на три сигмы но это не точно


# In[23]:


figure(figsize=(28, 12), dpi=90)
plt.plot(x, y-y_l, color = 'green') #сглаженный спектр
plt.plot(x, sigmas) #сигмы
plt.plot(x, -sigmas)


# In[24]:


#для эмиссионных
ysmooth = y-y_l
lines_em = []
for k in range(len(ysmooth)):
    if ysmooth[k] > sigmas[k]:
        lines_em.append(x[k])


# In[26]:


lines_em #это не сами линни, а все точки, которые их образуют


# In[ ]:




