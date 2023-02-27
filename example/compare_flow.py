import os
import sys
import numpy as np
import scipy.stats as st
import scipy.optimize as op
from matplotlib import pyplot as plt
import matplotlib as mat
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from scipy import interpolate

def fit(x,b,k): 
    return k*x+b

def fit2(x, b, k, a):
    # return k*x+b+a*10**(-x)
    return (k*x+b)/(a*x+1)

# inverse of R0=1/(klnT+b)+aTlnT
def fitR0(x, b, k, a):
    return (k*x + b) / (1 + (k*x + b)*a*x*10**(-x))

# inverse of chi0=1/(klnT+b)+c+aTlnT
def fitchi0(x, b, k, a, c):
    return (k*x + b) / (1 + (k*x+b)*(a*x*10**(-x) + c))

# generic for fit functions
# assume lamus go to 0 at Tc
def new_fit_flow(lnbetas, lamus, fitfunc=fit, init=0, fin=None, pguess=np.zeros(2)):
    print("fit function=", fitfunc)
    popt, pcov = op.curve_fit(fitfunc, lnbetas[init:fin], lamus[init:fin], p0=pguess)
    print(popt)
    log10tc = -op.fsolve(fitfunc, -popt[0]/popt[1], args=tuple(popt))
    print("log10tc=", log10tc)
    return log10tc, popt

# assume fname has betas, chis, invR0s, rss, ls
def new_extract_flow(fname, ischi=False, fitfunc=fit, init=0, fin=None, pguess=np.zeros(2), lnwd = 0.0):
    print("fit function=", fitfunc)
    fnamesp = fname.split("_")
    rs = float(fnamesp[2][2:])
    l = int(fnamesp[3][1:])
    betas, chis, invR0s, rss, ls = np.loadtxt(fname, unpack = True)
    lamus = invR0s
    lnbetas = np.log10(betas)
    lamus1 = [lamus[i] for i in range(len(lamus)) if (abs(lamus[i])>1e-4) and lnbetas[i]>lnwd]
    clamus1 = [chis[i] for i in range(len(lamus)) if (abs(lamus[i])>1e-4) and lnbetas[i]>lnwd]
    lnbetas1 = [lnbetas[i] for i in range(len(lamus)) if (abs(lamus[i])>1e-4) and lnbetas[i]>lnwd]
    if ischi:
        log10tc, popt = new_fit_flow(lnbetas1, clamus1, fitfunc=fitfunc, init=init, fin=fin, pguess=pguess)
    else:
        log10tc, popt = new_fit_flow(lnbetas1, lamus1, fitfunc=fitfunc, init=init, fin=fin, pguess=pguess)
    return log10tc, popt, np.array(lnbetas1), np.array(lamus1), np.array(clamus1), lnbetas, invR0s, chis

lnwd = np.log10(10)

fnamerpa = "gap3D_rpachi_rs2.0_l0_v1.txt"
fitfunc, nfitpara = fitR0, 3
log10tcrpa, poptrpa, lnbetasrpa, lamusrpa, clamusrpa, lnbetascrpa, invR0scrpa, chiscrpa = new_extract_flow(fnamerpa, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamusrpa = fitfunc(np.array(lnbetasrpa), *poptrpa)

fnameph0 = "gap3D_phchi_rs2.0_l0_v0.txt"
fitfunc, nfitpara = fitR0, 3
log10tcph0, poptph0, lnbetasph0, lamusph0, clamusph0, lnbetascph0, invR0scph0, chiscph0 = new_extract_flow(fnameph0, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamusph0 = fitfunc(np.array(lnbetasph0), *poptph0)

fnameph1 = "gap3D_phchi_rs2.0_l0_v1.txt"
fitfunc, nfitpara = fitR0, 3
log10tcph1, poptph1, lnbetasph1, lamusph1, clamusph1, lnbetascph1, invR0scph1, chiscph1 = new_extract_flow(fnameph1, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamusph1 = fitfunc(np.array(lnbetasph1), *poptph1)

fnameph2 = "gap3D_phchi_rs2.0_l0_v2.txt"
fitfunc, nfitpara = fitR0, 3
log10tcph2, poptph2, lnbetasph2, lamusph2, clamusph2, lnbetascph2, invR0scph2, chiscph2 = new_extract_flow(fnameph2, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamusph2 = fitfunc(np.array(lnbetasph2), *poptph2)

fname0 = "gap3D_phrpachi_rs2.0_l0_v0.txt"
fitfunc, nfitpara = fitR0, 3
log10tc0, popt0, lnbetas0, lamus0, clamus0, lnbetasc0, invR0sc0, chisc0 = new_extract_flow(fname0, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus0 = fitfunc(np.array(lnbetas0), *popt0)

fname1 = "gap3D_phrpachi_rs2.0_l0_v1.txt"
fitfunc, nfitpara = fitR0, 3
log10tc1, popt1, lnbetas1, lamus1, clamus1, lnbetasc1, invR0sc1, chisc1 = new_extract_flow(fname1, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus1 = fitfunc(np.array(lnbetas1), *popt1)

fname2 = "gap3D_phrpachi_rs2.0_l0_v2.txt"
fitfunc, nfitpara = fitR0, 3
log10tc2, popt2, lnbetas2, lamus2, clamus2, lnbetasc2, invR0sc2, chisc2 = new_extract_flow(fname2, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus2 = fitfunc(np.array(lnbetas2), *popt2)

fname3 = "gap3D_phrpachi_rs2.0_l0_v3.txt"
fitfunc, nfitpara = fitR0, 3
log10tc3, popt3, lnbetas3, lamus3, clamus3, lnbetasc3, invR0sc3, chisc3 = new_extract_flow(fname3, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus3 = fitfunc(np.array(lnbetas3), *popt3)

fname4 = "gap3D_phrpachi_rs2.0_l0_v4.txt"
fitfunc, nfitpara = fitR0, 3
log10tc4, popt4, lnbetas4, lamus4, clamus4, lnbetasc4, invR0sc4, chisc4 = new_extract_flow(fname4, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus4 = fitfunc(np.array(lnbetas4), *popt4)

fname5 = "gap3D_phrpachi_rs2.0_l0_v5.txt"
fitfunc, nfitpara = fitR0, 3
log10tc5, popt5, lnbetas5, lamus5, clamus5, lnbetasc5, invR0sc5, chisc5 = new_extract_flow(fname5, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus5 = fitfunc(np.array(lnbetas5), *popt5)

fname6 = "gap3D_phrpachi_rs2.0_l0_v6.txt"
fitfunc, nfitpara = fitR0, 3
log10tc6, popt6, lnbetas6, lamus6, clamus6, lnbetasc6, invR0sc6, chisc6 = new_extract_flow(fname6, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus6 = fitfunc(np.array(lnbetas6), *popt6)

fname7 = "gap3D_phrpachi_rs2.0_l0_v7.txt"
fitfunc, nfitpara = fitR0, 3
log10tc7, popt7, lnbetas7, lamus7, clamus7, lnbetasc7, invR0sc7, chisc7 = new_extract_flow(fname7, ischi=False, fitfunc=fitfunc, pguess=np.zeros(nfitpara), lnwd=lnwd)
fitlamus7 = fitfunc(np.array(lnbetas7), *popt7)

# plot setting
zoomin = 2
plt.switch_backend('TkAgg')
plt.style.use(['science'])
mat.rcParams.update({'font.size': 34/zoomin})
mat.rcParams["font.family"] = "Times New Roman"
size, sizein = 38/zoomin, 22/zoomin
Colormap = ["blue", "orange", "green", "red", "purple", "brown",
            "pink", "gray", "olive", "cyan",  "black"]
pt_L = ['o', '^', '8', 's', 'v', 'p', '<', 'h', '>', 'H', 'D']

fig, ax2 = plt.subplots()
fig.set_size_inches(13.5/zoomin, 10/zoomin)
plt.xscale("log")
plt.xlabel(r"$T/T_F$")
plt.ylabel(r"$1/R_0$")

plt.plot(1/10**lnbetascrpa, -invR0scrpa, "k+-", fillstyle="none", ms=2, label=r"rpa only")

plt.plot(1/10**lnbetasc0, -invR0sc0, "r+-", fillstyle="none", ms=2, label=r"$g=0.4,\omega_D={10}^{-8}$")
plt.plot(1/10**lnbetasc1, -invR0sc1, "ro", fillstyle="none", ms=2, label=r"$g=0.4,\omega_D=0.05$")
plt.plot(1/10**lnbetasc2, -invR0sc2, "rx", fillstyle="none", ms=2, label=r"$g=0.4,\omega_D=0.005$")
plt.plot(1/10**lnbetasc3, -invR0sc3, "rv", fillstyle="none", ms=2, label=r"$g=0.4,\omega_D=0.0005$")

plt.plot(1/10**lnbetasc7, -invR0sc7, "g+-", fillstyle="none", ms=2, label=r"$g=0.2,\omega_D={10}^{-8}$")
plt.plot(1/10**lnbetasc4, -invR0sc4, "go", fillstyle="none", ms=2, label=r"$g=0.2,\omega_D=0.005$")

plt.plot(1/10**lnbetasc6, -invR0sc6, "b+-", fillstyle="none", ms=2, label=r"$g=0.8,\omega_D={10}^{-8}$")
plt.plot(1/10**lnbetasc5, -invR0sc5, "bo", fillstyle="none", ms=2, label=r"$g=0.8,\omega_D=0.005$")

plt.plot(1/10**lnbetascph0, -invR0scph0, "c+-", fillstyle="none", ms=2, label=r"ph only,$g=0.4,\omega_D={10}^{-8}$")
plt.plot(1/10**lnbetascph1, -invR0scph1, "co", fillstyle="none", ms=2, label=r"ph only,$g=0.4,\omega_D=0.05$")
plt.plot(1/10**lnbetascph2, -invR0scph2, "cx", fillstyle="none", ms=2, label=r"ph only,$g=0.4,\omega_D=0.005$")

plt.plot(1/10**lnbetasc1, lnbetasc1*0.0, "k--")
# plt.plot(1/10**lnbetasc1, -invR0sc1, "ro", fillstyle="none", ms=2, label=r"$RPA$")
# plt.plot(1/10**lnbetas1, -fitlamus1, "r-")

# plt.plot(1/10**lnbetasc2, -invR0sc2, "bo", fillstyle="none", ms=2, label=r"$RPA+PH$")
# plt.plot(1/10**lnbetas2, -fitlamus2, "b-")

# plt.plot(1/10**lnbetasc1, invR0sc2-invR0sc1, "go", fillstyle="none", ms=2, label=r"$RPAPH-RPA$")
# plt.plot(1/10**lnbetasc1, invR0sc3, "co", fillstyle="none", ms=2, label=r"$PH$")
# plt.plot(1/10**lnbetasc1, -invR0sc1-invR0sc3+invR0sc2, "mo", fillstyle="none", ms=2, label=r"$RPAPH-RPA-PH$")

# plt.plot(1/10**lnbetas3, lamus3, "bo:", fillstyle="none", ms=2, label=r"$1/\chi_0(T)$")
#plt.plot(1/10**lnbetas2, 1-lamus2, "gv--", fillstyle="none", ms=2, label=r"$1-\lambda(T)$")
# plt.plot([10 ** log10tc,], [0.0,], "+", ms=8)
plt.legend(fontsize=6)
plt.savefig("compareflow.pdf")

