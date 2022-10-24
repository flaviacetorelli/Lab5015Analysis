import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import numpy as np
from sklearn import datasets, linear_model
from scipy.optimize import curve_fit
import pandas as pd

import matplotlib.pyplot as plt


#parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('--savePlots',          action='store_true',       default=False,      help='save plots')
parser.add_argument("-o", "--outdir",    action="store",      type=str,    help="output directory")
parser.add_argument("-i", "--inputs",      action="store",      type=int, nargs = '+' ,   help="txt files inputs from last step")


args = parser.parse_args()


outputPlot = args.outdir #"/eos/user/f/fcetorel/www/MTD/TBjune22/summaryOffset/run5335/"
runs = args.inputs#['offset4thr_5192.txt', 'offset4thr_5194.txt', 'offset4thr_5335.txt']

vovdict = {  
         5192 : 3.5,  
         5194 : 1.5, 
         5335 : 1.8, 

}

filenames = []
for run in runs:
    filenames.append('offset4thr_'+str(run)+'.txt')

dfs = []
for filename in filenames:
    dfs.append(pd.read_csv(filename,delimiter=','))

dftot = pd.concat(dfs)

# some cuts to eliminte some buggy points
dftot = dftot[dftot["mean"]>0]
dftot = dftot[dftot["mean"]<1]
dftot = dftot[dftot["emean"]<0.05]

bars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]



def func(x,b):
    return 0*x + b

formats = [".pdf", ".png"]

with open('offsetFit.txt', 'w') as f:
    for run in runs:
        vov = vovdict[run]
        for b in bars:
            df1 = dftot [dftot["bar"]==b]
            dfL = df1[df1["lab"]=="L"]
            dfR = df1[df1["lab"]=="R"]
    
            ax1 = dfL[dfL["vov"]==vov].plot.scatter(x="thr", y="mean", yerr = "emean", label='L', color = 'b', )   
            ax2 = dfR[dfR["vov"]==vov].plot.scatter(x="thr", y="mean", yerr = "emean",  label=' R', color = 'r', ax=ax1)
    
    
            # save mean in a txt file        
            #f.write("{b} 3.5 L {mean}  \n".format(b = b, mean = dfL[dfL["vov"]==3.5]["mean"].mean()))
            #f.write("{b} 1.5 L {mean}  \n".format(b = b, mean = dfL[dfL["vov"]==1.5]["mean"].mean()))
            #f.write("{b} 3.5 R {mean}  \n".format(b = b, mean = dfR[dfR["vov"]==3.5]["mean"].mean()))
            #f.write("{b} 1.5 R {mean}  \n".format(b = b, mean = dfR[dfR["vov"]==1.5]["mean"].mean()))
    
            
            
            #fit with a constant
            pars1, cov1 = curve_fit(func, xdata=dfL[dfL["vov"]==vov]["thr"].values, ydata=dfL[dfL["vov"]==vov]["mean"].values)
            pars2, cov2 = curve_fit(func, xdata=dfR[dfR["vov"]==vov]["thr"].values, ydata=dfR[dfR["vov"]==vov]["mean"].values)
   
    
            plt.plot(dfL[dfL["vov"]==vov]["thr"].values, func(dfL[dfL["vov"]==vov]["thr"].values,pars1[0]), '--', color ='b', label = "fit L ")
            plt.plot(dfR[dfR["vov"]==vov]["thr"].values, func(dfR[dfR["vov"]==vov]["thr"].values,pars2[0]), '--', color ='r', label ="fit R ")
   
            f.write("{b} {vov} L {mean}  \n".format(b = b, vov = str(vov) , mean = str(pars1).replace("[","").replace("]","")))
            f.write("{b} {vov} R {mean}  \n".format(b = b, vov = str(vov) ,  mean = str(pars2).replace("[","").replace("]","")))
   
            plt.ylim(0.,1.)
    
            plt.legend()
            plt.xlabel('ithr1')
            plt.ylabel('mean')
            plt.title("bar "+ str(b))
            plt.grid()
    
            if args.savePlots:
                if not os.path.exists(outputPlot+'/run'+str(run)):
                    os.makedirs(outputPlot+'/run'+str(run))   
 
                for form in formats:
                    plt.savefig(outputPlot+'/run'+str(run)+"/offset_bar_"+ str(b)+form)
            plt.clf()
            plt.close()
#if drawAllTogheter:
#    for b in bars:
#        df1 = dftot [dftot["bar"]==b]
#        dfL = df1[df1["lab"]=="L"]
#        dfR = df1[df1["lab"]=="R"]
#
#        ax1 = dfL[dfL["vov"]==3.5].plot.scatter(x="thr", y="mean", yerr = "emean", label='3.5 L')
#        ax3 = dfR[dfR["vov"]==3.5].plot.scatter(x="thr", y="mean", yerr = "emean", label='3.5 R', color = 'g', ax=ax1)
#        ax2 = dfL[dfL["vov"]==1.5].plot.scatter(x="thr", y="mean", yerr = "emean", label='1.5 L', color = 'r', ax=ax1)   
#        ax4 = dfR[dfR["vov"]==1.5].plot.scatter(x="thr", y="mean", yerr = "emean",  label='1.5 R', color = 'y', ax=ax1)
#        ax5 = dfL[dfL["vov"]==1.8].plot.scatter(x="thr", y="mean", yerr = "emean", label='1.8 L', color = 'c', )   
#        ax6 = dfR[dfR["vov"]==1.8].plot.scatter(x="thr", y="mean", yerr = "emean",  label='1.8 R', color = 'm', ax=ax5)
#
#
#        # save mean in a txt file        
#        #f.write("{b} 3.5 L {mean}  \n".format(b = b, mean = dfL[dfL["vov"]==3.5]["mean"].mean()))
#        #f.write("{b} 1.5 L {mean}  \n".format(b = b, mean = dfL[dfL["vov"]==1.5]["mean"].mean()))
#        #f.write("{b} 3.5 R {mean}  \n".format(b = b, mean = dfR[dfR["vov"]==3.5]["mean"].mean()))
#        #f.write("{b} 1.5 R {mean}  \n".format(b = b, mean = dfR[dfR["vov"]==1.5]["mean"].mean()))
#
#        
#        
#        #fit with a constant
#        pars1, cov1 = curve_fit(func, xdata=dfL[dfL["vov"]==3.5]["thr"].values, ydata=dfL[dfL["vov"]==3.5]["mean"].values)
#        pars2, cov2 = curve_fit(func, xdata=dfR[dfR["vov"]==3.5]["thr"].values, ydata=dfR[dfR["vov"]==3.5]["mean"].values)
#        pars3, cov3 = curve_fit(func, xdata=dfL[dfL["vov"]==1.5]["thr"].values, ydata=dfL[dfL["vov"]==1.5]["mean"].values)
#        pars4, cov4 = curve_fit(func, xdata=dfR[dfR["vov"]==1.5]["thr"].values, ydata=dfR[dfR["vov"]==1.5]["mean"].values)
#        pars5, cov5 = curve_fit(func, xdata=dfL[dfL["vov"]==1.8]["thr"].values, ydata=dfL[dfL["vov"]==1.8]["mean"].values)
#        pars6, cov6 = curve_fit(func, xdata=dfR[dfR["vov"]==1.8]["thr"].values, ydata=dfR[dfR["vov"]==1.8]["mean"].values)
#
#
#        plt.plot(dfL[dfL["vov"]==3.5]["thr"].values, func(dfL[dfL["vov"]==3.5]["thr"].values,pars1[0]), '--', color ='b', label = "fit L 3.5")
#        plt.plot(dfR[dfR["vov"]==3.5]["thr"].values, func(dfR[dfR["vov"]==3.5]["thr"].values,pars2[0]), '--', color ='g', label ="fit R 3.5")
#        plt.plot(dfL[dfL["vov"]==1.5]["thr"].values, func(dfL[dfL["vov"]==1.5]["thr"].values,pars3[0]), '--', color ='r', label ="fit L 1.5")
#        plt.plot(dfR[dfR["vov"]==1.5]["thr"].values, func(dfR[dfR["vov"]==1.5]["thr"].values,pars4[0]), '--', color ='y', label ="fit R 1.5")
#        plt.plot(dfL[dfL["vov"]==1.8]["thr"].values, func(dfL[dfL["vov"]==1.8]["thr"].values,pars5[0]), '--', color ='c', label ="fit L 1.8")
#        plt.plot(dfR[dfR["vov"]==1.8]["thr"].values, func(dfR[dfR["vov"]==1.8]["thr"].values,pars6[0]), '--', color ='m', label ="fit R 1.8")
#
#        f.write("{b} 3.5 L {mean}  \n".format(b = b, mean = str(pars1).replace("[","").replace("]","")))
#        f.write("{b} 3.5 R {mean}  \n".format(b = b, mean = str(pars2).replace("[","").replace("]","")))
#        f.write("{b} 1.5 L {mean}  \n".format(b = b, mean = str(pars3).replace("[","").replace("]","")))
#        f.write("{b} 1.5 R {mean}  \n".format(b = b, mean = str(pars4).replace("[","").replace("]","")))
#        f.write("{b} 1.8 L {mean}  \n".format(b = b, mean = str(pars5).replace("[","").replace("]","")))
#        f.write("{b} 1.8 R {mean}  \n".format(b = b, mean = str(pars6).replace("[","").replace("]","")))
#
#        plt.ylim(0.,1.)
#
#        plt.legend()
#        plt.xlabel('ithr1')
#        plt.ylabel('mean')
#        plt.title("bar "+ str(b))
#        plt.grid()
#
#        if opttions.savePlots: 
#            for form in formats:
#                plt.savefig(outputPlot+"offset_bar_all_"+ str(b)+form)
#    























































