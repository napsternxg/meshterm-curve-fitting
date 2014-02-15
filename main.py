#!/usr/bin/env python

import pandas as pd
import pylab as pl
import numpy as np
from scipy import interpolate, optimize
import math

df = pd.read_csv("HIV.tsv",sep="\t",index_col="Year")
column_names = ["y"]
index_names = ["x"]
df.index.names = index_names

df.columns = column_names
df.y = df.y.apply(lambda y: math.log(y)+1)
min_year = 1970
max_year = 2014
x = range(min_year,max_year)
x1 = range(0,max_year-min_year)
y = [0 if i not in df.y else df.y[i] for i in x]

print x,x1,y

s = interpolate.PchipInterpolator(x1,y)
r = s(x1)

def func(x,a,b,c,d):
    return a+ (b/(1+np.exp(-(c+d*x))))

popt,pcov = optimize.curve_fit(func,x1,y)
# First derivative
def f_der1(x):
    b,c,d = popt[1:]
    e = 2+ np.exp(-(c+d*x)) + np.exp(c+d*x)
    return b*d/e

# Second derivative
def f_der2(x):
    b,c,d = popt[1:]
    t1 = b*d/np.power((2+ np.exp(-(c+d*x)) + np.exp(c+d*x)),2)
    t2 = (d*np.exp(-(c+d*x)) - d*np.exp(c+d*x))
    return t1*t2




y2 = [func(i,popt[0],popt[1],popt[2],popt[3]) for i in x1]
y3 = [f_der1(i) for i in x1]
y4 = [f_der2(i) for i in x1]

"""
yp = []
for i in range(3,8):
    z = np.polyfit(x1,y,i)
    f = np.poly1d(z)
    pl.plot(x,f(x1),"-",linewidth=2,label="{0}-D Poly Fit curve".format(i))
"""
pl.subplot(4,1,1)
pl.plot(x,y,"o",markersize=10,label="Data Points")
pl.plot(x,r,"-",label="Cubic Spline Interpolation")
pl.legend(loc="best")
pl.title("Curve fitting for HIV MeshTerm")
pl.xlabel("Year")
pl.ylabel("log(TotalCount)")

pl.subplot(4,1,2)
pl.plot(x,y,"o",markersize=10,label="Data Points")
pl.plot(x,y2,"--",linewidth=2,label="Regression Fit Curve")
pl.legend(loc="best")
pl.xlabel("Year")
pl.ylabel("log(TotalCount)")

pl.subplot(4,1,3)
pl.plot(x,y3,"-",linewidth=2,label="First Derivative")
pl.legend(loc="best")
pl.xlabel("Year")
pl.ylabel("f'(log(TotalCount))")

pl.subplot(4,1,4)
pl.plot(x,y4,"-",linewidth=2,label="Second Derivative")
pl.legend(loc="best")
pl.xlabel("Year")
pl.ylabel("f''(log(TotalCount))")

pl.show()

