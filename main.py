#!/usr/bin/env python

import pandas as pd
import matplotlib
matplotlib.use('Agg') 

import pylab as pl
import numpy as np
from scipy import interpolate, optimize
import math

"""
Year wise cumulative count of
pubmed ids
"""

def normalizeAndLogCount(k):
    global tyd
    normalizeF =  1.0*(tyd[tyd.year == k["x"]].TotalPMID.values[0])
    k["y"] = 1.0*k["y"]/normalizeF
    k["y"] = math.log(k["y"] + 1)
    return k

def prepareData(df):
    if normalized:
        print "Using Normalized routine"
        tyd = pd.read_csv("PMID_PER_YEAR.tsv", sep="\t")
        tyd = tyd[(tyd.year > 1901) & (tyd.year <= 2011)] 
        df = df.apply(normalizeAndLogCount, axis=1)
        df[["x"]] = df[["x"]].astype(int)
    else:
        df.y = df.y.apply(lambda y: math.log(y+1))

    print df

    min_year = df.x.min() - 5
    max_year = df.x.max() -1
    x = range(min_year,max_year)
    x1 = range(0,max_year-min_year)
    y = []
    for i in x:
        if i not in df.x.values:
            y.append(0)
        else:
            y.append(df[df.x == i].y.values[0])

    #y = [0 if i not in df.x else df[df.x == i].y for i in x]

    print x,x1,y
    return x,x1,y

def getData(filename, normalized):
    global tyd
    column_names = ["x","y"]
    df = pd.read_csv(filename,sep="\t", dtype = {'Year':int, 'TotalPMID':float})
    df.columns = column_names
    return prepareData(df)


def func(x, a, b, c, d):
    return a+ (b/(1.0+np.exp(-(c+d*x))))

# First derivative
def f_der1(x, popt):
    b,c,d = popt[1:]
    e = 2+ np.exp(-(c+d*x)) + np.exp(c+d*x)
    return b*d/e

# Second derivative
def f_der2(x, popt):
    b,c,d = popt[1:]
    t1 = b*d/np.power((2+ np.exp(-(c+d*x)) + np.exp(c+d*x)),2)
    t2 = (d*np.exp(-(c+d*x)) - d*np.exp(c+d*x))
    return t1*t2

def plot(x,y,x2,r,y2,y3,y4, saveFigure=False, filename="None"):
    print x,y,x2,r,y2,y3,y4
    title  = "Curve fitting for MeshTerm"
    if normalized:
        title += " [Normalized]"
        filename += "_Normalized"
    pl.subplot(4,1,1)
    pl.plot(x,y,"o",markersize=10,label="Data Points")
    pl.plot(x,r,"-",label="Cubic Spline Interpolation")
    pl.legend(loc="best")
    pl.title(title)
    pl.xlabel("Year")
    pl.ylabel("log(TotalCount)")

    pl.subplot(4,1,2)
    pl.plot(x,y,"o",markersize=10,label="Data Points")
    pl.plot(x2,y2,"--",linewidth=2,label="Regression Fit Curve")
    pl.legend(loc="best")
    pl.xlabel("Year")
    pl.ylabel("log(TotalCount)")

    pl.subplot(4,1,3)
    pl.plot(x2,y3,"-",linewidth=2,label="First Derivative")
    pl.legend(loc="best")
    pl.xlabel("Year")
    pl.ylabel("f'(log(TotalCount))")

    pl.subplot(4,1,4)
    pl.plot(x2,y4,"-",linewidth=2,label="Second Derivative")
    pl.legend(loc="best")
    pl.xlabel("Year")
    pl.ylabel("f''(log(TotalCount))")
    if saveFigure:
        filename = filename + ".png"
        pl.savefig(filename, bbox_inches = "tight")
        print "Figure saved as: %s" % filename
    else:
        pl.show()

def main(filename, normalized=False, saveFigure=False):
    x,x1,y = getData(filename, normalized)
    #Interpolate using spline
    s = interpolate.PchipInterpolator(x,y)
    r = s(x)
    
    #Fit curve to logistic regression and get coefficients
    popt,pcov = optimize.curve_fit(func,x1,y)
    print "Coefficients: %s" % popt
    x2 = np.arange(min(x1),max(x1)+1, 0.2)
    y2 = [func(i,popt[0],popt[1],popt[2],popt[3]) for i in x2]
    y3 = [f_der1(i,popt) for i in x2]
    y4 = [f_der2(i,popt) for i in x2]
    x2 = [min(x)+k for k in x2]
    filename = filename.split('.')[0]
    plot(x,y,x2,r,y2,y3,y4, saveFigure = saveFigure, filename = filename)





"""
yp = []
for i in range(3,8):
    z = np.polyfit(x1,y,i)
    f = np.poly1d(z)
    pl.plot(x,f(x1),"-",linewidth=2,label="{0}-D Poly Fit curve".format(i))
"""

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Create Plots for Mesh Term\
                                     Data")
    parser.add_argument("-f",help="Name of file to get data from. Should\
                        contain\
                        2 columns only seperated by tab. Year and Count")
    parser.add_argument("-n", type=bool, help="Normalize counts. Set to activate")
    parser.add_argument("-s", type=bool, help="Save figure and not show. Set to activate")
    args = parser.parse_args()
    filename = args.f
    normalized = args.n
    saveFigure  = args.s
    print "Normalized:", normalized
    if filename == None:
        print parser.print_help()
        exit(-1)
    if normalized == None:
        normalized == False
    if saveFigure == None:
        saveFigure == False
    main(filename, normalized, saveFigure)

