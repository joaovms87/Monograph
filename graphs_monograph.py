import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
plt.rcParams['errorbar.capsize'] = 3
# plt.rcParams['mathtext.fontset'] = 'cm'
# plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['grid.color'] = 'k'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.5
plt.rcParams['font.size'] = 12
plt.rcParams['legend.fontsize'] = 'large'
plt.rcParams['figure.titlesize'] = 'medium'

plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


beta = 0.1
DELTA = 2
for B0 in [-1, -2]:
    Bf = B0 + DELTA
    # reading the previously calculated values for W
    Wpoints = np.loadtxt(f'W; B0={B0}; beta={beta}', complex)
    WADpoints = np.loadtxt(f'W_ad; B0={B0}; beta={beta}', complex)
    n = len(Wpoints)

    Dpoints = np.loadtxt(f'D; B0={B0}; beta={beta}', complex)
    D2points = np.loadtxt(f'D2; B0={B0}; beta={beta}', complex)

    # list of the process durations that will be investigated
    taupoints = np.logspace(-2, 4, n)

    br = 100
    for i, tau in enumerate(taupoints):
        if tau > br:
            break

    # fit for log scale (W)
    fit = np.polyfit(np.log(taupoints[i:]), np.log(np.real(Wpoints[i:]-WADpoints[i:])), deg=1, cov=True)
    print(f'Wex... For B0={B0}:')
    print(ufloat(fit[0][0], np.sqrt(fit[1][0, 0])))
    W_fit = taupoints**fit[0][0] * np.exp(fit[0][1])

    # fit for log scale (D2)
    fit = np.polyfit(np.log(taupoints[i:]), np.log(np.real(D2points[i:])), deg=1, cov=True)
    print(f'D2... For B0={B0}:')
    print(ufloat(fit[0][0], np.sqrt(fit[1][0, 0])))
    D2_fit = taupoints**fit[0][0] * np.exp(fit[0][1])


    # plot W in linear scale
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.plot(taupoints, np.real(Wpoints), 'r.', label=r'$W_\tau(\tau)$')
    ax.plot(taupoints, np.real(WADpoints), 'k--', label=r'$W_{ad}(\tau)$')
    ax.legend(fancybox=False, edgecolor='black')
    ax.set_xlabel(r'$\tau J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1}$')
    ax.set_title(r'$\tau$-dependence of the total work')
    ax.set_xlim(0, 100)
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\LZ\B0={B0}/W(tau),B0={B0},beta={beta}.png', bbox_inches='tight')

    # plot Wex in log scale
    label = r'Fit for $\tau J/\hbar > $'+f'{br}'
    Wex = np.real(Wpoints-WADpoints)
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.loglog(taupoints, Wex, 'r.', label=r'$W_{\tau,ex}(\tau)$')
    ax.loglog(taupoints, W_fit, 'k--', label=label)
    ax.legend(fancybox=False, edgecolor='black')
    ax.set_xlabel(r'$\tau J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1}$')
    ax.set_title(r'$\tau$-dependence of the excess work')
    ax.set_ylim(min(Wex), 1.5*max(Wex))
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\LZ\B0={B0}/log(Wex)_log(tau),B0={B0},beta={beta}.png', bbox_inches='tight')

    # plot D in log scale with Wex
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.loglog(taupoints, Wex, 'r.', label=r'$W_{\tau,ex}(\tau)$')
    ax.loglog(taupoints, np.real(Dpoints), 'k.', label=r'$D[\rho_\tau (\tau) || \rho_{eq} (\tau)]$')
    ax.legend(fancybox=False, edgecolor='black', loc='lower left')
    ax.set_xlabel(r'$\tau J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1} \ $ or $\ \ D$')
    ax.set_title(r'Comparison between $W_{ex}$ and $D$')
    ax.set_xlim(1e-2, 1e4)
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\LZ\B0={B0}/Wex vs D,B0={B0},beta={beta}.png', bbox_inches='tight')

    # plot D2 in log scale together with Wex
    label_W = r'$W_{ex}$ fit for $\tau \cdot J/\hbar > $'+f'{br}'
    label_D = r'$D_{\tau,ad}$ fit for $\tau \cdot J/\hbar > $'+f'{br}'
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.loglog(taupoints, Wex, 'r.', label=r'$W_{\tau,ex}(\tau)$')
    ax.loglog(taupoints, W_fit, 'k-', label=label_W)
    ax.loglog(taupoints, np.real(D2points), 'b.', label=r'$D_{\tau,ad}(\tau)$')
    ax.loglog(taupoints, D2_fit, 'k--', label=label_D)
    ax.legend(fancybox=False, edgecolor='black', loc='lower left', fontsize=11)
    ax.set_xlabel(r'$\tau J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1} \ $ or $\ \ D$')
    ax.set_title(r'Comparison between $W_{\tau,ex}(\tau)$ and $D_{\tau,ad}(\tau)$')
    ax.set_xlim(1e-2, 1e4)
    ax.set_ylim(top=1)
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\LZ\B0={B0}/Wex vs D(rho,rho_ad),B0={B0},beta={beta}.png', bbox_inches='tight')


