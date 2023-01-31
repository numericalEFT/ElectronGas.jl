#!/home/pchou/miniconda3/bin/python
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from scipy import interpolate

plt.switch_backend('TkAgg')
plt.style.use(['science'])
mat.rcParams.update({'font.size': 32})
mat.rcParams["font.family"] = "Times New Roman"
size, sizein = 37, 22

Colormap = ["blue", "orange", "green", "red",  "purple", "brown",
            "pink", "gray", "olive", "cyan",  "black"]
pt_L = ['o', '^', '8', 's', 'v', 'p', '<', 'h', '>', 'H', 'D']
# fig, axes = plt.subplots(
# 1, 2, sharex='col', gridspec_kw={'wspace': 0})
# fig.set_size_inches(22, 10)
# fig, axes = plt.subplots(1, 1, sharex='col', gridspec_kw={'wspace': 0})
# fig.set_size_inches(13, 10)
# axes = axes.reshape(-1)
# Rss = [0.2, 0.4, 0.8]
fig, ax2 = plt.subplots()
fig.set_size_inches(13.5, 10)
# channels = range(6, 10+1)
ell = 0  # 0: Explicit; 1: Renorm
sigtype = 1  # 1: RPA;  2: G0W0

D = 3

msize, msizein = 12, 6


def readdata(rs, channel, type):
    if type == 0:
        fname = "rs{0}_sigtype{3}/gap{1}DExplicit_rs{0}_l{2}.txt".format(
            rs, D, channel, sigtype)
    elif type == 1:
        fname = "rs{0}_sigtype{3}/gap{1}D_rs{0}_l{2}.txt".format(
            rs, D, channel, sigtype)
    print('Loading '+fname)
    dat = np.loadtxt(fname)
    bet, rss, chan = dat[:, 0], dat[:, -1], dat[:, -2]
    lam = dat[:, 1]
    outp = (rss == rs)*(chan == channel)  # *(lam < 0.999)
    if True in outp:
        beta, data = dat[outp, 0], dat[outp, 1]
        return beta, data
    else:
        return [0], [0]


def readchiR(rs, channel, is_ph):
    if is_ph:
        fname = "gap{1}D_phchi_rs{0}_l{2}.txt".format(rs, D, channel)
    else:
        fname = "gap{1}D_chi_rs{0}_l{2}.txt".format(rs, D, channel)
    print('Loading '+fname)
    dat = np.loadtxt(fname)
    bet, rss, chan = dat[:, 0], dat[:, -1], dat[:, -2]
    chi, lam = dat[:, 1], dat[:, 2]

    return bet, chi, lam


def fitfunc(x, a, b):
    return a*x+b


def fitfunc0(x, a, tc):
    return a*np.log(x/tc)


def fitfunc1(x, a, b, tc):
    return a*np.log(x/tc) / (1 + b*np.log(x/tc))
    # return a*np.log(x/tc) + b*(x-tc)


############### figures ############
rs, ell = 2.0, 0
# xgrid, ygrid = readdata(rs, ell, 0)
# func = interpolate.PchipInterpolator(np.log(1/xgrid[::-1]), ygrid[::-1])
# ax2.plot(1/xgrid, ygrid-1, 's', color="blue",
#          ms=msize, label=r"$\lambda_{\rm max}(T)-1$")
# grid = np.linspace(-14.5, 0.0, num=200)
# ax2.plot(np.exp(grid), func(grid)-1, color="blue", lw=2, ls='--')

# K, data = readdata(rs, ell, 1)
K, chi0_inv, R0_inv = readchiR(rs, ell, True)

x, y = 1/K, R0_inv
x1, y1 = 1/K, -chi0_inv

err = y*1e-8
popt, pcov = curve_fit(
    fitfunc0, x[K > 100], y[K > 100], p0=[-0.32, 5e-5])
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
# Tc = np.exp(-popt[1]/popt[0])
Tc = popt[1]
if Tc > 1:
    Tc = 0
print(Tc)
print(Tc*np.exp(-1/popt[0]))

popt1, pcov1 = curve_fit(
    fitfunc1, x1[K >= 100], y1[K >= 100], p0=[-0.328, 0.04, 5e-5])
perr1 = np.sqrt(np.diag(pcov1))
ax2.plot(x1, y1, '^', ms=msize, color="green", label=r"$1/\chi_0(T)$")

print(popt1)
print(perr1)
grid = np.linspace(-10.5, 0.0, num=500)
ylam1 = fitfunc1(np.exp(grid), popt1[0], popt1[1], popt1[2])
ax2.plot(np.exp(grid), ylam1, color="green", lw=2)


ax2.plot(x, y, 'o', ms=10, color="red", label=r"$1/R_0(T)$")
t = np.array([1e-5, 1])
# ylam = np.log(t)*popt[0] + popt[1]
# ylam0 = fitfunc0(t, -0.225, 1e-5)
# ax2.plot(t, ylam0, color="green", lw=2, ls='--')
ylam = fitfunc0(t, popt[0], popt[1])
ax2.plot(t, ylam, color="red", lw=2)

# without phonon
K, chi0_inv, R0_inv = readchiR(rs, ell, False)

x, y = 1/K, R0_inv
x1, y1 = 1/K, -chi0_inv

err = y*1e-8
# popt, pcov = curve_fit(
# fitfunc0, x[K > 100], y[K > 100], p0=[-0.1, 1e-7])
# perr = np.sqrt(np.diag(pcov))
# print(popt)
# print(perr)
# Tc = popt[1]
# if Tc > 1:
#     Tc = 0
# print(Tc)
# print(Tc*np.exp(-1/popt[0]))

# popt1, pcov1 = curve_fit(
#     fitfunc1, x1[K >= 100], y1[K >= 100], p0=[-0.1, 0.01, 1e-7])
# perr1 = np.sqrt(np.diag(pcov1))
ax2.plot(x1, y1, '^', ms=msize, color="blue", label=r"$1/\chi_0(T)$")

# print(popt1)
# print(perr1)
# grid = np.linspace(-10.5, 0.0, num=500)
# ylam1 = fitfunc1(np.exp(grid), popt1[0], popt1[1], popt1[2])
# ax2.plot(np.exp(grid), ylam1, color="blue", lw=2)


ax2.plot(x, y, 'o', ms=10, color="black", label=r"$1/R_0(T)$")
t = np.array([1e-5, 1])
# ylam = np.log(t)*popt[0] + popt[1]
# ylam0 = fitfunc0(t, -0.225, 1e-5)
# ax2.plot(t, ylam0, color="green", lw=2, ls='--')
# ylam = fitfunc0(t, popt[0], popt[1])
# ax2.plot(t, ylam, color="black", lw=2)


ax2.axhline(y=0, ls=':', lw=2, color='black')
ax2.set_xlim(2e-5, 1)
ax2.set_ylim(-2.55, 0.2)
ax2.set_xscale("log")
ax2.set_xlabel(r"$T/T_{\textsc f}$", size=size)

# ax2.text(0.75, 0.76, r"$\lambda_{\rm max}(T)-1$",
#          transform=ax2.transAxes, size=35)
ax2.text(0.32, 0.28, r"$1/R_0(T)$",
         transform=ax2.transAxes, size=35)
ax2.text(0.65, 0.52, r"$1/\chi_0(T)$",
         transform=ax2.transAxes, size=35)

ax2.tick_params(which='major', length=8, width=1.0)
ax2.tick_params(which='minor', length=4, width=1.0)
for tick in ax2.xaxis.get_major_ticks():
    tick.set_pad(8)
ax2.set_title(r"$r_s=2, \ell=0$, 3D RPA-$G_0$ with el-ph")

####################### INSET ###############
# ax2ins = inset_axes(ax2, width="42%", height="50%", bbox_transform=ax2.transAxes,
#                     loc=3, bbox_to_anchor=(0.0, 0.025, 0.9, 0.94))
# labels = ['a', 'b', 'c', 'd']

# rss = [0.5, 0.6, 0.7, 0.8]
# # ells = [0, 1, 2, 3]
# ell = 3
# gs = [-0.028, -0.043, -0.058, -0.075]
# Tcs = [3e-22, 1.45e-15, 3.86e-12, 4.66e-10]
# for (i, rs) in enumerate(rss):
#     K, data = readdata(rs, ell, 1)
#     x, y = 1/K, data-1
#     # popt, pcov = curve_fit(fitfunc, np.log(x[K >= 800]), y[K >= 800])
#     # Tc = np.exp((-popt[1])/popt[0])
#     popt, pcov = curve_fit(
#         fitfunc0, x[K > 200], y[K > 200], p0=[gs[i], Tcs[i]])
#     print(popt)
#     print(pcov)
#     perr = np.sqrt(np.diag(pcov))
#     ax2ins.plot(x, y, 'o', color=Colormap[i],
#                 # ax2.plot(x, y-1, pt_L[i], color=Colormap[i],
#                 ms=msizein, label=r"${0}$".format(rs))
#     Tc = popt[1]
#     print(Tc)
#     print(Tc*np.exp(-1/popt[0]))
#     t = np.array([Tc, 1])
#     # ylam = np.log(t)*popt[0] + popt[1]
#     ylam = fitfunc0(t, popt[0], popt[1])
#     ax2ins.plot(t, ylam, color=Colormap[i], lw=1.5)
# ax2ins.axhline(y=0, ls=':', lw=1.5, color='black')

# ax2ins.set_xlim(1.1e-22, 2)
# ax2ins.set_ylim(-1.42, 0.08)
# # ax2.set_xlim(1e-22, 0.01)
# # ax2.set_ylim(-1.25, 0.06)
# ax2ins.set_xscale("log")

# ax2ins.legend(loc=(0.07, 0.08), title="$r_s$", fontsize=sizein,
#               title_fontsize=24, handletextpad=0.1)
# ax2ins.set_xlabel(r"$T/T_{\textsc f}$", size=24)
# ax2ins.set_ylabel(r"$1/R_0(\ell=3; T)$", size=24)

# ax2ins.tick_params(labelsize=sizein)
# ax2ins.yaxis.set_ticks_position('right')
# ax2ins.xaxis.set_label_coords(0.5, 0.1)
# ax2ins.yaxis.set_label_coords(0.1, 0.5)

fig.tight_layout()
if sigtype == 1:
    # title = "(c1) $r_s={0}$, RPA".format(rs)
    # ax2.text(0.3, 0.9, title, size=size, transform=ax2.transAxes)
    plt.savefig("./{0}dRPA_phchi.pdf".format(D, rs))
elif sigtype == 2:
    # title = "(c2) $r_s={0}$, G0W0".format(rs)
    # ax2.text(0.3, 0.9, title, size=size, transform=ax2.transAxes)
    plt.savefig("./lam{0}dvsT_G0W0_fig.pdf".format(D, rs))

# plt.show()
exit()
