
### This script works for analyzing Alessandros rat data in a more comples way.
### Here we analyze clusters of channels, and get one number for each channel, ### from each rat for each measure in each condition, averaged over trials.
### We analyze 6 random channels (or other number) from the set of all channels ### The script also filters the data in different bands if needed.

#%% IMPORTS, INITIALIZING, ETC

## Starting stuff
#%load_ext autoreload
#%reload_ext autoreload
#%autoreload 2

import read_write as rw
import func_ale_data as fad
import func_various as fv
import func_plots as fp

import numpy as np
import pyconscious as pc

# import functions as func
import matplotlib as mpl
import time
import pandas as pd
import scipy.stats as sct
import matplotlib.pyplot as plt
import autoreload
from importlib import reload
import scipy.io <http://scipy.io>  as sio import numpy.random as ra import scipy.signal as ss import copy import scipy.interpolate as sci import collections import os import pickle import mne import scipy import math from copy import deepcopy as cp

import networkx as netx
import bct

filepath = "Z:\\" #r"Z:\results\Nilsen_2021_rat_spont_analysis"
filepath_raw = filepath + r"Data\Rat_anesthesia"#r"Z:\results\Nilsen_2021_rat_spont_analysis"
filepath = filepath + r"results\Nilsen_2021_rat_spont_analysis"#r"Z:\results\Nilsen_2021_rat_spont_analysis"

condition = [r"W&K", r"W&P", r"W&S"]
wake = 3*["Wake"]
anes = ["Keta","Prop","Sevo"]

#General properties
trials=80
times=[300,2300]#[300,2300] # normal is 300, 2300
iters=10
tau=1
chE = 6  # nr chans
frqbands = [0.1,45.] # added a bit more
measures = ["LZ","LZ"]#["pow"]#["wPLI"]# "alt_star","alt_geo","alt_SI","alt_MI","alt_M1"] #,"alt_geo","alt_SI","alt_MI","alt_M1"]#"LZ","ACE","SCE"] #,"alt_star","alt_geo","alt_SI","alt_MI","alt_M1","phi","ents"]
chanlist = [1,2,5,6,9,10,13,14]
hz = 500
chans = 14

# Reloads
reload(rw)
reload(fad)
reload(fv)
reload(fp)
reload(plt)




#%% LOAD DATA
# Make nestable dict
dataA = fv.nested_dict()

#### Load data
for co in range(len(condition)): #iterate through folders
    
    files = sorted(os.listdir(os.path.join(filepath_raw, condition[co])))
    its = -1
    rat = "j"
    
    for fy in files:

        y = copy.deepcopy(fy)
          
        if not fy[:3] == rat:
            rat = fy[:3]
            its = its + 1
            
        con = anes[co] if anes[co] in fy else "Wake"
        
        d = rw.readMatFile(os.path.join(filepath_raw, condition[co]), fy, wait = False) # Load
          
        d = fad.badchan(d) # Remove chans marked bad and rename properly
        hz = d["fs"]
        
        dataA[co][con][its] = copy.deepcopy(d)

rw.save(dataA, os.path.join(filepath, "raw.pkl")) dataA = rw.load(os.path.join(filepath, "raw.pkl"))

#%% PREPROCESSING        
#### Select trials and channels
for co in dataA.keys(): #iterate through conditions
    
    for con in dataA[co].keys(): #iterate through states
    
        for its in dataA[co][con].keys(): #iterate through rats
            
            d = dataA[co][con][its]
            
            a = np.shape(d["data"])[2]
            
            tr = a if a < trials else trials  #find max number of trials
            tr = np.random.choice(list(range(a)), tr, replace=False)               # Selecting random trials
            channels = sorted(np.random.choice(d["gCh"], chans, replace=False))               # Selecting 15 random channels
            chs = [d["gCh"].index(ch) for ch in channels]                          # Selecting from chanlist
            dataA[co][con][its]["channels"] = channels
            dataA[co][con][its]["t_data"] = copy.deepcopy(d["data"][:, chs][times[0]:times[1], :, tr]) # Store it again (trimmed data)
            print("Rat {} - Con {}".format(its,anes[co]))
            print(channels)
rw.save(dataA, os.path.join(filepath, "trimmed.pkl")) dataA = rw.load(os.path.join(filepath, "trimmed.pkl"))

#### FILTER DATA
#  Create filter params and plot the response b, a = fv.filter_params(4, hz, frqbands[0], frqbands[1], "bandpass", "butter", analog = False)

fp.freqrespplot(b, a, frqbands, "Digital Butterworth filter") fp.savePlot(filepath, "freq_resp{}-{}.png".format(int(frqbands[0]), int(frqbands[1])), close = True)

# Filter the data itself
for co in dataA.keys(): #iterate through conditions
    
    for con in dataA[co].keys(): #iterate through states
    
        for its in dataA[co][con].keys(): #iterate through rats
            
            d = dataA[co][con][its]
            k = []
                       
            for z in range(np.shape(d["t_data"])[2]):
                k.append([fv.filter_it(b, a, d["t_data"][:, zz, z], 150, "even") for zz in range(len(d["t_data"][0, :, 0]))])
            k
            dataA[co][con][its]["f_data"] = np.transpose(k) # Store it again (filtered data)
            dataA[co][con][its].pop("data") rw.save(dataA, os.path.join(filepath, "filtered.pkl"))
dataA = rw.load(os.path.join(filepath, "filtered.pkl"))     







#%% ANALYSIS 

results = fv.nested_dict()

#### Generate a simulated dataset
dataB = []
import pdc_dtf as pd
reload(pd)
t=10
p = 3
for i in range(5000):
    try:
        # for t,p in zip([10,10,10,10],[1,2,3,4]):
        print("Running number: " + str(i))
        temp = {"a" : 0}#np.random.rand()} #dropping connections %
        temp["c"] = np.random.uniform(0.0,1.0)
        temp["weight"] = [np.random.uniform(0.001,10.1), np.random.uniform(0.001,10.1)]
        #temp["weight"] = [temp["weight"], temp["weight"]/np.random.randint(1,20)] #normal(poisson) dist for uni(bi)directional
        d, temp["constr"], temp["cons"], temp["modules"], temp["mod"], temp["mat"], temp["noise"], temp["weight"] = simmod(10000, 14, temp["a"], temp["c"], temp["weight"], direct= "bi", t_thick = t, delay = p)
        temp["gaussp"]=np.mean(sct.normaltest(d,axis=1)[1])
        temp["LZs"]=np.mean([pc.LZc([e]) for e in d])
        temp["ACE"]=pc.ACE(d)
        temp["adf"]=pc.stationarity(d)
        temp["corr"]=np.mean(np.abs(np.corrcoef(d)))
        for meas in ["alt_star","alt_SI", "alt_MI", "alt_M1", "alt_geo"]: #
            # while True:
                # try:
            rw.saveMatFile(r"Y:\Coding\MATLAB_Packages\phi_toolbox", r"data.mat", {"ls":d, "m":meas, "frBand":[1.0, 40.0] ,"sRate":500})
    
            time.sleep(1.0)
                    # break
                # except:
                #     print("Issue writing, trying again...")
                #     time.sleep(1.0)
            r = rw.readMatFile(r"Y:\Coding\MATLAB_Packages\phi_toolbox", r"result.mat", wait = True)
            temp[meas] = copy.deepcopy(r["res"][0][0])
            os.remove(r"Y:\Coding\MATLAB_Packages\phi_toolbox\result.mat")
    
        temp_dtf = []
        D, P, f = pd.run_all(d)
        print(t,p)
        # plt.title("{}-{}".format(t,p))
        temp_dtf.append(copy.deepcopy(D))
        x = np.arange(0,40)#len(temp_dtf[0]))
        temp["dtf"] = np.trapz(D[x,:,:],x,axis=0)/len(x)
        
        data=1-temp["dtf"]
        data[np.isinf(data)] = 0
        data[np.isnan(data)] = 0
        #data[data==m] = 0
    
        temp["dtf_int"] = net_integration(data)
    
        # m = dataB[i]["int"]*1 if dataB[i]["int"]>m else m
        # data = np.abs(np.log2(1-dataB[i]["mat"]))
        # np.fill_diagonal(data, 0)
        #data = np.abs(np.log2(1-dataB[i]["mat"]/m))
        data = temp["dtf"]*1
        data[np.isinf(data)] = 0
        data[np.isnan(data)] = 0
    
        temp["dtf_mods"], temp["dtf_mod"] = net_mod_mods(data)
        dataB.append(copy.deepcopy(temp))
        print(i)
        #print(temp["LZs"],temp["alt_geo"],temp["adf"],temp["modules"],temp["cons"],temp["constr"])
        # print("LZ - geo: ", sct.spearmanr([x["LZs"] for x in dataB],[x["alt_geo"] for x in dataB]))
        # print("geo - modules: ", sct.spearmanr([x["alt_geo"] for x in dataB],[x["modules"] for x in dataB]))
        # print("geo - con*str: ", sct.spearmanr([x["alt_geo"] for x in dataB],[x["constr"]*x["cons"] for x in dataB]))
        # print("lz - con*str: ", sct.spearmanr([x["gaussp"] for x in dataB],[x["LZs"] for x in dataB]))
        rw.save(dataB, os.path.join(filepath, "simulated3.pkl"))
    except:
        print("Something went wrong, skipping to next")
        os.remove(r"Y:\Coding\MATLAB_Packages\phi_toolbox\data.mat")
        try:
            os.remove(r"Y:\Coding\MATLAB_Packages\phi_toolbox\result.mat")
        except:
            print("No result file found")
        break

  

#### PCI ST numbers
pcist1 = {"Keta":{"Keta":[62.98,77.51,48.10,70.77,58.67,17.34,45.15], # 0-600ms
                  "Wake":[107.23,101.53,60.68,67.65,60.73,42.71,66.42]},
          "Prop":{"Prop":[19.34,63.50,36.48,34.77,60.48,32.66,43.13,21.53,50.17],
                  "Wake":[65.21,93.11,80.29,64.24,82.25,71.99,61.98,43.39,64.84]},
          "Sevo":{"Sevo":[15.87,72.91,51.61,48.60,31.78,63.82,36.12,21.71,24.33],
                  "Wake":[62.39,76.01,81.40,98.18,77.18,65.19,53.99,46.47,68.05]}}
pcist2 = {"Keta":{"Keta":[32.47,32.71,21.41,38.68,20.74,4.96,9.12], # 80-600ms
                  "Wake":[53.08,66.14,30.28,54.81,30.11,24.39,36.04]},
          "Prop":{"Prop":[1.74,10.71,5.89,5.61,10.60,1.36,5.46,4.85,13.43],
                  "Wake":[35.18,54.91,50.96,42.40,53.14,27.21,50.55,31.37,35.42]},
          "Sevo":{"Sevo":[4.34,15.73,8.61,14.60,2.21,18.02,10.60,3.05,3.64],
                  "Wake":[31.61,31.50,51.70,55.28,34.66,33.32,29.70,7.86,36.33]}}

for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            results[con][state][rat]["PCIst1"] = pcist1[anes[con]][state][rat]
            results[con][state][rat]["PCIst2"] = pcist2[anes[con]][state][rat]


#### Spectral exponent, for each channel, PSD combined over epoch (OLD does PSD for each epoch) frwq= [20, 40] c = 12 for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            x = exponent(d, frwq, dataA[con][state][rat]["fs"], method = "linear", ea = "mean", ca = "mean", rat = rat, con = state, plotit=c)
#            rw.saveMatFile(r"Y:\Coding\MATLAB_Packages\SpecExp_Colombo", r"data.mat", {"ls":d, "frBand":[1.0, 20.0] ,"sRate":dataA[con][state][rat]["fs"]})
#            r = rw.readMatFile(r"Y:\Coding\MATLAB_Packages\SpecExp_Colombo", r"result.mat", wait = True)
            results[con][state][rat]["SpecExp_20_40"] = copy.deepcopy(x)
#            results[con][state][rat]["SpecExp2_1_20"] = copy.deepcopy(r["result"][0][0])
#            os.remove(r"Y:\Coding\MATLAB_Packages\SpecExp_Colombo\result.mat")
            title = "SpecExp-rat{}-ch{}-frq{}_{}.png".format(rat,c,frwq[0],frwq[1])
            fp.savePlot(filepath, title, close = True)
            
#### LZ single channel, mean, self normed for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            
            temp = []
            
            for c in range(np.shape(d)[1]):
                temp.append([pc.LZc(np.transpose(d[:, c, :, np.newaxis], (0, 2, 1)), threshold = "mean", shuffles = 3, ea="mean")])
            
            results[con][state][rat]["LZs"] = copy.deepcopy(temp) #### Stationarity test for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            temp = pc.stationarity(d)
            results[con][state][rat]["adf"] = copy.deepcopy(temp)
            print(temp)


#### LZ all channels, mean, self normed
for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            
            temp = pc.LZc(d, threshold = "mean", shuffles = 3, ea="mean")
            
            results[con][state][rat]["LZc"] = copy.deepcopy(temp)

#### ACE all channels, mean, self normed for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            
            temp = pc.ACE(d, threshold = "mean", shuffles = 3, ea="mean")
            
            results[con][state][rat]["ACE"] = copy.deepcopy(temp)
            
#### mean correlation
for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            c = []
            for i in d:
                c.append(np.mean(np.abs(np.corrcoef(i))))
            print(np.mean(c))
            results[con][state][rat]["meanC"] = np.mean(c)
            
#### SCE all channels, mean, self normed for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            
            temp = pc.SCE(d, threshold = 0.8, shuffles = 3, ea="mean")
            
            results[con][state][rat]["SCE"] = copy.deepcopy(temp)
  
            
#### PHI MEASURES
import mne
for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["f_data"])
            for meas in ["alt_geo"]: #["alt_star","alt_SI", "alt_MI", "alt_M1"]: #
                rw.saveMatFile(r"Y:\Coding\MATLAB_Packages\phi_toolbox", r"data.mat", {"ls":d, "m":meas, "frBand":[1.0, 40.0] ,"sRate":dataA[con][state][rat]["fs"]})
                r = rw.readMatFile(r"Y:\Coding\MATLAB_Packages\phi_toolbox", r"result.mat", wait = True)
                results[con][state][rat]["PHI_"+meas+"_high"] = copy.deepcopy(r["res"][0])
                os.remove(r"Y:\Coding\MATLAB_Packages\phi_toolbox\result.mat")
        rw.save(results, os.path.join(filepath, "results/results3.pkl")) 
        rw.save([con, state, rat], os.path.join(filepath, "results/temp.pkl")) 
        
#### PSI measure (and other cons)
for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            d = np.transpose(dataA[con][state][rat]["t_data"])
            for meas in [[tuple(np.arange(0,35,1)),tuple(np.arange(10,45,1))]]:#[[0.5,4]]:
                # r = mne.connectivity.spectral_connectivity(d,sfreq=hz,method=meas,cwt_freqs=np.arange(1,45,1) ,faverage = True, mode="cwt_morlet", fmin=1,fmax=45)
                # results[con][state][rat]["con_"+meas] = copy.deepcopy(r[0])
                # results[con][state][rat]["con_mean_"+meas] = np.nanmean(r[0])
                # results[con][state][rat]["con_std_"+meas] = np.nanstd(r[0])
                x = mne.connectivity.phase_slope_index(d, sfreq=hz,  mode="multitaper", fmin = meas[0], fmax=meas[1], n_jobs=4)[0]
                xjk = [mne.connectivity.phase_slope_index(d[np.arange(len(d))!=i,:,:], sfreq=hz,  mode="multitaper",fmin = meas[0], fmax=meas[1], n_jobs=4)[0] for i in range(len(d))]
                #x = np.squeeze(np.mean(mne.connectivity.phase_slope_index(d, sfreq=hz,  mode="multitaper", fmin = meas[0], fmax=meas[1], n_jobs=4)[0],axis=2)) #, fmin=1, fmax=45
                #xjk = [np.squeeze(np.mean(mne.connectivity.phase_slope_index(d[np.arange(len(d))!=i,:,:], sfreq=hz,  mode="multitaper",fmin = meas[0], fmax=meas[1], n_jobs=4)[0],axis=2)) for i in range(len(d)-40)]
                xjk = np.sqrt(len(d))*np.std(xjk,axis=0)
                xjk = x/xjk
                xjk[np.isnan(xjk)] = 0
                xjk = np.trapz(xjk,np.arange(0,len(meas[0]),1),axis=2)/len(meas[0])
                xjk = 0.5 - sct.norm.sf(check_empty_and_transform(xjk, 0))
                results[con][state][rat]["con_psi"] = copy.deepcopy(xjk)
                # results[con][state][rat]["con_psi_mt"] = copy.deepcopy(x)

#### Calc absolute power in bands
for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():
            for f in [[0.1, 5.0], [5.0, 8.0],[8.0,12.0],[12.0,25.0],[25.0,45.0]]:
                results[con][state][rat]["pow_"+str(f)] = bandpow(dataA[con][state][rat]["t_data"], hz, f[0],f[1])

#### Calc mean and relative power over all bands (can then be used for adjustment or control) for con in results.keys():
    for state in list(results[con].keys())[:2]:
        for rat in results[con][state].keys():
            results[con][state][rat]["pow_avg"] = np.mean([results[con][state][rat]["pow_delta"],
                                                           results[con][state][rat]["pow_theta"],
                                                           results[con][state][rat]["pow_alpha"],
                                                           results[con][state][rat]["pow_beta"],
                                                           results[con][state][rat]["pow_gamma"]])
            results[con][state][rat]["pow_delta_rel"] = results[con][state][rat]["pow_delta"] / results[con][state][rat]["pow_avg"]
            results[con][state][rat]["pow_theta_rel"] = results[con][state][rat]["pow_theta"] / results[con][state][rat]["pow_avg"]
            results[con][state][rat]["pow_alpha_rel"] = results[con][state][rat]["pow_alpha"] / results[con][state][rat]["pow_avg"]
            results[con][state][rat]["pow_beta_rel"] = results[con][state][rat]["pow_beta"] / results[con][state][rat]["pow_avg"]
            results[con][state][rat]["pow_gamma_rel"] = results[con][state][rat]["pow_gamma"] / results[con][state][rat]["pow_avg"]

#### Calc fb flow ratio - negative = back to front for con in results.keys():
    for state in list(results[con].keys())[:2]:
        for rat in results[con][state].keys():
            for m in ["dtf", "pdc", "psi"]:
                results[con][state][rat]["con_"+m+"_fb_ratio"] = np.log2(results[con][state][rat]["con_"+m+"_fb_flow"] / results[con][state][rat]["con_"+m+"_bf_flow"])


#### DTF and PDC test
import pdc_dtf as pd
reload(pd)

for con in dataA.keys():
    for state in dataA[con].keys():
        
        for rat in dataA[con][state].keys():
            data = np.transpose(dataA[con][state][rat]["t_data"], (2,1,0))
            data = data.astype('float32')
            temp_dtf = []
            temp_pdc = []
            for i in range(np.shape(data)[0]):
                # mu = np.mean(data[i], axis=1)
                # X = (data[i] - mu[:, None])/np.std(data[i],axis=1)[:, None]
                # p, bic = pd.compute_order(X + np.random.randn(np.shape(X)[0],np.shape(X)[1])/10, p_max=20)
                # temp_dtf.append(copy.deepcopy(bic))
                
                # #plt.figure()
                # plt.plot(np.arange(1,22,1), bic)#np.nanmean(bic,axis=0))
                
                D, P, f = pd.run_all(data[i])
                temp_dtf.append(copy.deepcopy(D))
                x = np.arange(0,40)#len(temp_dtf[0]))
                
                temp_pdc.append(calc_graph_stuff(np.trapz(D[x,:,:],x,axis=0)/len(x),results[con][state][rat]["channels"]))
                
                #temp_pdc.append(copy.deepcopy(P))
            temp_pdc = np.median(temp_pdc,axis=0)
            results = place_graph_results(results, con, state, rat, temp_pdc, "con_dtf_t_")
            x = np.arange(0,40)#len(temp_dtf[0]))
            results[con][state][rat]["con_dtf"] = np.trapz(np.nanmean(temp_dtf, axis=0)[x,:,:],x,axis=0)/len(x)
            #results[con][state][rat]["con_pdc"] = np.trapz(np.nanmean(temp_pdc, axis=0)[x,:,:],x,axis=0)/len(x)
            # np.fill_diagonal(results[con][state][rat]["con_dtf"],0)
            # np.fill_diagonal(results[con][state][rat]["con_pdc"],0)
        
del D, P, temp_dtf, temp_pdc, data, f, axes, m, N, i, j

#### Multivar transfer entropy test with idtxl import idtxl import pyopencl as cl from idtxl.multivariate_te import MultivariateTE from idtxl.data import Data as idData import idtxl.estimators_opencl

data = idData()
data.generate_mute_data(n_samples=800, n_replications=5)

network_analysis = MultivariateTE()
settings = {'cmi_estimator': "JidtGaussianCMI",#'JidtGaussianCMI',
            'gpuid':"10DE 128B 118C196E",
            'max_lag_sources': 25,
            'min_lag_sources': 1,
            'n_perm_max_stat':200,
            'n_perm_min_stat': 200,
            'n_perm_omnibus': 200,
            'n_perm_max_seq': 200,}
for con in dataA.keys():
    for state in dataA[con].keys():
        for rat in dataA[con][state].keys():    
            results[con][state][rat]["gaussp"] = np.mean(sct.normaltest(dataA[con][state][rat]["t_data"], axis=0)[1])
            print(np.mean(sct.normaltest(dataA[con][state][rat]["t_data"], axis=0)[1]))
            d = idData(dataA[con][state][rat]["t_data"][:,:,0:5], dim_order="spr")
            results[con][state][rat]["cmi_net"] = network_analysis.analyse_network(settings=settings, data=d) del d, r import pyopencl as cl
reload(cl)


cl.release()

platform = cl.get_platforms()
gpus = platform[0].get_devices(device_type=cl.device_type.GPU)
ctx = cl.Context(devices=gpus)



platform = cl.get_platforms()[0]
device = platform.get_devices()[0]
platform.get_devices(device_type=cl.device_type.GPU) != [] context = cl.Context([device]) # del context, platform, device # b) Initialise analysis object and define settings network_analysis = MultivariateTE() settings = {'cmi_estimator': 'OpenCLKraskovCMI', #JidtGaussianCMI OpenCLKraskovCMI JidtKraskovCMI
            'gpuid':0,
            'max_lag_sources': 5,
            'min_lag_sources': 1,
            "debug": True}

r = network_analysis.analyse_network(settings=settings, data=data)



#### Test graph theory connectivity stuff    


def calc_graph_stuff(mat, channels):
    data = 0 + mat
    data[np.isinf(data)] = 0
    data[np.isnan(data)] = 0
    bf_flow = total_information_flow(data,channels, direction="back_front")
    fb_flow = total_information_flow(data, channels, direction="front_back")
                    
    temp = net_basic_mes(data)
    weight = np.mean(temp[0])
    cons = np.mean(temp[1])
    cons_str = np.mean(temp[2])
    temp = net_mod_mods(data)
    mods = np.mean(temp[0])
    mod = np.mean(temp[1])
    
    data = 1-data
    integ = net_integration(data)
    
    return [bf_flow, fb_flow, weight, cons, cons_str, mods, mod, integ]

def place_graph_results(results, con, state, rat, temp, m):
    results[con][state][rat][m+"_bf_flow"] = temp[0]
    results[con][state][rat][m+"_fb_flow"] = temp[1]
            
    results[con][state][rat][m+"_weight"] = temp[2]
    results[con][state][rat][m+"_cons"] = temp[3]
    results[con][state][rat][m+"_cons*str"] = temp[4]
    
    results[con][state][rat][m+"_mods"] = temp[5]
    results[con][state][rat][m+"_mod"] = temp[6]

    results[con][state][rat][m+"_int"] = temp[7]
    
    return results
    
meas = ["con_dtf"] # d for delta, p for pci band, g for gamma for con in results.keys():
    for state in results[con].keys():
        for rat in results[con][state].keys():
            # f,ax=plt.subplots(1,3)
            # [ax[i].imshow(np.abs(np.log2(results[con][state][rat][meas[i]][:,:]))) for i in range(len(meas))]
            for m in meas:
                temp = calc_graph_stuff(results[con][state][rat][m], results[con][state][rat]["channels"])
                results = place_graph_results(results, con, state, rat, temp, m)
