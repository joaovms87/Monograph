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


betapoints = [0.1, 10]
h0 = 0.1
DELTA = 2

# list of the process durations that will be investigated
taupoints = np.logspace(-2, 4, 200)

br = 100
for i, tau in enumerate(taupoints):
    if tau > br:
        break


for j, beta in enumerate(betapoints):
    j += 1
    hf = h0 + DELTA
    # reading the previously calculated values for W
    Wpoints = np.loadtxt(f'new W beta={beta} h0={h0}', complex)
    WADpoints = np.loadtxt(f'new W_ad beta={beta} h0={h0}', complex)
    n = len(Wpoints)

    Wex = np.real(Wpoints - WADpoints)

    Dpoints = np.loadtxt(f'new D beta={beta} h0={h0}', complex)
    D2points = np.loadtxt(f'new D2 beta={beta} h0={h0}', complex)


    fit1, cov1 = np.polyfit(np.log(taupoints[i:]), np.log(np.real(Wpoints[i:]-WADpoints[i:])), deg=1, cov=True)
    print(f'For h0={h0}, Wex fit...')
    print(ufloat(fit1[0], np.sqrt(cov1[0, 0])))
    W_fit1 = taupoints**fit1[0] * np.exp(fit1[1])

    # plot Wex in log scale
    fig, ax = plt.subplots()
    label = r'Fit for $\tau J/\hbar > $' + f'{br}'
    ax.loglog(taupoints, Wex, 'k.', label=r'$W_\tau(\tau)-W_{ad}(\tau)$')
    ax.loglog(taupoints, W_fit1, 'b--', label=label)
    ax.legend(fancybox=False, edgecolor='black')
    ax.set_xlabel(r'$\tau\cdot J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1}$')
    ax.set_title(r'$\tau$-dependence of the excess work')
    ax.set_xlim(1e-2, 1e4)
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\spins3\case{j}/log(W)_log(tau).png', bbox_inches='tight')


    # plot W in linear scale
    fig, ax = plt.subplots()
    ax.plot(taupoints[:i], np.real(Wpoints[:i]), 'r.', label=r'$W_\tau(\tau)$')
    ax.plot(taupoints[:i], np.real(WADpoints[:i]), 'k--', label=r'$W_{ad}(\tau)$')
    ax.legend(fancybox=False, edgecolor='black')
    ax.set_xlabel(r'$\tau\cdot J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1}$')
    ax.set_title(r'Dependência com $\tau$ do trabalho total')
    ax.set_box_aspect(1)
    ax.set_xlim(1e-2, 100)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\spins3\case{j}/W_linear.png', bbox_inches='tight')

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
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\spins3\case{j}/D vs Wex.png', bbox_inches='tight')

    # fit for D2
    fit1, cov1 = np.polyfit(np.log(taupoints[i:]), np.log(np.real(D2points[i:])), deg=1, cov=True)
    print(f'For h0={h0}, D2 fit...')
    print(ufloat(fit1[0], np.sqrt(cov1[0, 0])))
    D2_fit1 = taupoints**fit1[0] * np.exp(fit1[1])


    # plot D2 in log scale
    fig, ax = plt.subplots()
    label = r'Fit for $\tau \cdot J/\hbar >$'+f'{br}'
    ax.loglog(taupoints, np.real(D2points), 'k.', label=r'$D_{\tau,ad}(\tau)$')
    ax.loglog(taupoints, D2_fit1, 'b--', label=label)
    ax.legend(fancybox=False, edgecolor='black')
    ax.set_xlabel(r'$\tau J/\hbar$')
    ax.set_ylabel(r'$D_{\tau,ad}(\tau)$')
    ax.set_title(r'$\tau$-dependence of $D_{\tau,ad}(\tau)$')
    ax.set_xlim(1e-2, 1e4)
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\spins3\case{j}/log(D2)_log(tau).png', bbox_inches='tight')

    # plot D2 in log scale together with Wex
    fig, ax = plt.subplots()
    ax.grid(True)
    label1 = r'$W_{ex}$ fit for $\tau J/\hbar > $'+f'{br}'
    label2 = r'$D_{\tau, ad}(\tau)$ fit for $\tau J/\hbar > $'+f'{br}'
    ax.loglog(taupoints, Wex, 'r.', label=r'$W_{\tau,ex}(\tau)$')
    ax.loglog(taupoints, np.real(D2points), 'b.', label=r'$D_{\tau,ad}(\tau)$')
    ax.loglog(taupoints, W_fit1, 'k-', label=label1)
    ax.loglog(taupoints, D2_fit1, 'k--', label=label2)
    ax.legend(fancybox=False, edgecolor='black', loc='lower left', fontsize=10)
    ax.set_xlabel(r'$\tau J/\hbar$')
    ax.set_ylabel('$W \cdot J^{-1} \ $ or $\ \ D$')
    ax.set_title(r'Comparison between $W_{\tau,ex}(\tau)$ and $D_{\tau,ad}(\tau)$')
    ax.set_xlim(1e-2, 1e4)
    ax.set_box_aspect(1)
    plt.savefig(f'D:\GRÁFICOS IC\MONOGRAFIA\spins3\case{j}/Wex vs D2.png', bbox_inches='tight')

