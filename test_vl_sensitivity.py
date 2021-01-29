'''
Exploring the effect of viral load profiles and their relationship to test sensitivity, and infectiousness
'''

# Imports and settings
import os
import sciris as sc
import covasim as cv
import covasim.utils as cvu
import pylab as pl
import seaborn as sns
import numpy as np
import scipy.stats as stats

do_plot   = 1
do_save   = 1

#%% Analyses

def test_vl_testing(fig_path=None):
    sc.heading('Test viral loads with different types of test')
    sc.heading('Setting up...')
    sc.tic()

    if fig_path is None:
        fig_path = f'results/test_with_vl.png'

    pars = dict(
        pop_size=50e3,
        pop_infected=100,
        start_day='2020-02-01',
        end_day='2020-03-30',
    )
    sim = cv.Sim(pars)
    sim.initialize()

    # Create interventions
    # Set parameters
    n_people = sim['pop_size']
    n_tests = 0.02 * n_people
    delay = 1
    start_day = 20

    tn3 = cv.test_num(daily_tests=n_tests, lod={3.0:0.9, 5.0:0.1}, symp_test=1.0, start_day=start_day, test_delay=delay)
    tn5 = cv.test_num(daily_tests=n_tests, lod={3.0:0.6, 5.0:0.4}, symp_test=1.0, start_day=start_day, test_delay=delay)
    tp3 = cv.test_prob(symp_prob=0.1, asymp_prob=0.01, lod=3.0, start_day=start_day, test_delay=delay)
    tp5 = cv.test_prob(symp_prob=0.1, asymp_prob=0.01, lod=5.0, start_day=start_day, test_delay=delay)

    # Define the simulations
    sim3n = cv.Sim(interventions=tn3, label='90% PCR tests, test_num')
    sim5n = cv.Sim(interventions=tn5, label='60% PCR tests, test_num')
    sim3p = cv.Sim(interventions=tp3, label='PCR tests, test_prob')
    sim5p = cv.Sim(interventions=tp5, label='Rapid tests, test_prob')

    # Run and plot results
    to_plot = ['new_infections', 'cum_infections', 'new_inf_neg', 'cum_inf_neg', 'new_noninf_pos', 'cum_noninf_pos', 'new_diagnoses', 'cum_diagnoses']

    msim = cv.MultiSim([sim3n, sim5n])#, sim3p, sim5p])
    msim.run(keep_people=True)
    msim.plot(to_plot=to_plot, do_save=True, fig_path=fig_path, do_show=False)

    return msim

#%% Run as a script
if __name__ == '__main__':

    sc.tic()
    msim = test_vl_testing()
    sc.toc()


print('Done.')



###
# Relevant things to plot:
#   1. Number of infectious negatives
#   2. Number of noninfectious positives
#   3. %presymptomatic/asymptomatic infections



### OLD
# def sim_vl(x, t_0, t_peak, t_fade, v_peak):
#
#     npts = len(x)
#
#     slope1 = (v_peak - 3) / (t_peak - t_0)
# #    eqn1 = 3 + slope1*(x - t_0)
#     slope2 = (6 - v_peak) / (t_fade - t_peak)
# #    eqn2 = v_peak + slope2*(x - t_peak)
#
#     inds1 = sc.findinds(x<t_0)
#     inds2 = x[(x>=t_0) * (x<t_peak)]
#     inds3 = x[(x>=t_peak) * (x<28)]
#     inds4 = x[(x>=28)]
#
#     y = np.zeros(npts)
#     y[inds1] = np.nan
#     y[inds2] = 3 + slope1*(x[inds2] - t_0)
#     y[inds3] = v_peak + slope2*(x[inds3] - t_peak)
#     y[inds4] = np.nan
#     y[y<3] = np.nan
#
#     return y
#
# def msim_vls(n, fig_path='results/plot_vl_larremore2.png'):
#
#     x = np.arange(0,30)
#     n_pts = len(x)
#     y = {}
#     for inf_type in ['symp', 'asymp']: y[inf_type] = np.zeros((n, n_pts))
#
#     for p in range(n):
#         # Define parameters
#         t_0 = np.random.uniform(2.5, 3.5)
#         t_peak = 0.5 + t_0 + np.minimum(np.random.gamma(1.5), 3)
#         v_peak = np.random.uniform(7, 11)
#         t_fade = {}
#         t_fade['symp'] = t_peak + np.random.uniform(4, 9)
#         t_fade['asymp'] = t_peak + np.random.uniform(0, 3) + np.random.uniform(4, 9)
#
#         for inf_type in ['symp','asymp']:
#             y[inf_type][p, :] = sim_vl(x, t_0, t_peak, t_fade[inf_type], v_peak)
#
#     pl.figure(figsize=(24, 20))
#     font_size = 36
#     font_family = 'Libertinus Sans'
#     pl.rcParams['font.size'] = font_size
#     pl.rcParams['font.family'] = font_family
#     xl, xr, yb, ym, yt = 0.05, 0.05, 0.05, 0.05, 0.05
#     dx, dy = 1 - xl - xr, (1 - yb - ym - yt) / 2
#     ax = {}
#
#     # Symptomatic trajectories
#     ax[0] = pl.axes([xl, yb, dx, dy])
#     ax[0].plot(x, y['symp'].T, color='grey', linewidth=2)
#     ax[0].set_xlim(0, 28)
#
#     # Asymptomatic trajectories
#     ax[1] = pl.axes([xl, yb + ym + dy, dx, dy])
#     ax[1].plot(x, y['asymp'].T, color='grey', linewidth=2)
#     ax[1].set_xlim(0, 28)
#
#     if do_save:
#         cv.savefig(fig_path, dpi=100)
#
#
# def sim_vls(n, n_days, n_pts=None, fig_path='results/plot_vl_larremore.png'):
#     '''
#     n: number of people
#     n_days: number of days
#     n_pts: number of time points, typically more than the number of days, e.g. 10 points/day
#     '''
#
#     if n_pts is None: n_pts = n_days*10 # Evaluate 10 time points per day
#     pts_per_day = n_pts/n_days
#
#     # Process randomly generated numbers to indices
#     def process_to_inds(arr):
#         return np.round(arr*pts_per_day).astype(int)
#
#     # Define parameters
#     raw_t_0 = np.random.uniform(2.5,3.5,size=n)
#     raw_t_peak = 0.5 + raw_t_0 + np.minimum(np.random.gamma(1.5,size=n),3)
#     v_peak = np.random.uniform(7,11,size=n)
#     raw_t_f_asymp = raw_t_peak + np.random.uniform(4,9,size=n)
#     raw_t_f_symp  = raw_t_peak + np.random.uniform(0,3,size=n) + np.random.uniform(4,9,size=n)
#
#     # Round them for processing
#     t_0 = process_to_inds(raw_t_0)
#     t_peak = process_to_inds(raw_t_peak)
#     t_f_asymp = process_to_inds(raw_t_f_asymp)
#     t_f_symp = process_to_inds(raw_t_f_symp)
#
#     # Generate the hinge function - any way to do this without looping over people??
#     x = np.arange(n_pts + 50)
#     y = {}
#     y['asymp'] = np.zeros((n, n_pts+50))*np.nan
#     y['symp']  = np.zeros((n, n_pts+50))*np.nan
#     polys = {'asymp': {}, 'symp': {}}
#     for inf_type in ['asymp', 'symp']:
#         polys[inf_type]['xvals'] = []
#         polys[inf_type]['yvals'] = []
#         polys[inf_type]['pars']  = []
#
#     for p in range(n):
#         # Asymptomatic infections
#         y['asymp'][p, t_0[p]:t_peak[p]] = np.linspace(3, v_peak[p], t_peak[p] - t_0[p])
#         y['asymp'][p, t_peak[p]:2*t_f_asymp[p]-t_peak[p]] = np.linspace(v_peak[p], 3, (t_f_asymp[p] - t_peak[p])*2)
#         # Symptomatic infections
#         y['symp'][p, t_0[p]:t_peak[p]] = np.linspace(3, v_peak[p], t_peak[p] - t_0[p])
#         y['symp'][p, t_peak[p]:2*t_f_symp[p]-t_peak[p]] = np.linspace(v_peak[p], 3, (t_f_symp[p] - t_peak[p])*2)
#
#         for inf_type in ['asymp', 'symp']:
#             inds = ~np.isnan(y[inf_type][p])
#             polys[inf_type]['xvals'].append(x[inds])
#             thesepars = np.polyfit(x[inds], y[inf_type][p][inds], 3)
#             polys[inf_type]['pars'].append(thesepars)
#             polyfn = np.poly1d(thesepars)
#             polys[inf_type]['yvals'].append(polyfn(x[inds]))
#
#     if do_plot:
#         pl.figure(figsize=(24, 20))
#         font_size = 36
#         font_family = 'Libertinus Sans'
#         pl.rcParams['font.size'] = font_size
#         pl.rcParams['font.family'] = font_family
#         xl, xr, yb, ym, yt = 0.05, 0.05, 0.05, 0.05, 0.05
#         dx, dy = 1-xl-xr, (1-yb-ym-yt)/2
#         X = np.array([x]*n)
#         ax = {}
#
#         # Symptomatic trajectories
#         ax[0] = pl.axes([xl, yb, dx, dy])
#         ax[0].scatter(X/10, y['symp'], color='grey', linewidth=2)
#         for p in range(n):
#             ax[0].plot(polys['symp']['xvals'][p] / 10, polys['symp']['yvals'][p], color='red', linewidth=2)
#         ax[0].set_xlim(0, 28)
#
#         # Asymptomatic trajectories
#         ax[1] = pl.axes([xl, yb+(ym+dy), dx, dy])
#         ax[1].scatter(X/10, y['asymp'], color='grey', linewidth=2)
#         for p in range(n):
#             ax[1].plot(polys['asymp']['xvals'][p] / 10, polys['asymp']['yvals'][p], color='red', linewidth=2)
#         ax[1].set_xlim(0, 28)
#
#         if do_save:
#             cv.savefig(fig_path, dpi=100)
#
#     return y
#
#
#
# def test_vl_sensitivity(method, beta=None, do_plot=True, do_save=True, fig_path=None, verbose=0.1):
#     sc.heading('Test viral loads over time')
#     sc.heading('Setting up...')
#     sc.tic()
#
#     if fig_path is None:
#         fig_path = f'results/plot_vl_{method}.png'
#
#     if beta is not None:
#         pars = {'beta': beta}
#         sim = cv.Sim(pars) # create sim object
#     else:
#         sim = cv.Sim()
#
#     # Create interventions
#     # Set parameters
#     n_people = sim['pop_size']
#     n_tests = 0.1 * n_people
#     delay = 1
#     start_day = 40
#
#     tn = cv.test_num(daily_tests=n_tests, symp_test=1.0, start_day=start_day, test_delay=delay)
#     tp = cv.test_prob(symp_prob=0.1, asymp_prob=0.1, start_day=start_day, test_delay=delay)
#
#     # Finally, update the parameters
#     sim.update_pars(interventions=[tn])
#
#     sim.run(verbose=verbose, keep_people=True)
#
#     vl = sc.promotetoarray(sim.people.viral_log)
#     vl[vl == 0.] = np.nan
#
#     if method=='new2':
#         # Get useful indices
#         exp_inds  = cvu.defined(sim.people.date_exposed)
#         symp_inds = cvu.defined(sim.people.date_symptomatic)
#         asymp_inds = np.setdiff1d(exp_inds, symp_inds)
#
#         vls2plot = {'symp': [],'asymp': []}
#         shedding_dur = {'symp': [],'asymp': []}
#         for new_ind, orig_ind in enumerate(symp_inds):
#             thisvl = vl[:, orig_ind]
#             inds2plot = np.arange(sim.people.date_exposed[orig_ind].astype(int), cvu.defined(thisvl)[-1]+2)
#             vls2plot['symp'].append(vl[inds2plot, orig_ind])
#             shedding_dur['symp'].append(len(cvu.defined(thisvl)))
#         for p in asymp_inds:
#             thisvl = vl[:, p]
#             inds2plot = np.arange(sim.people.date_exposed[p].astype(int), cvu.defined(thisvl)[-1]+2)
#             vls2plot['asymp'].append(vl[inds2plot, p])
#             shedding_dur['asymp'].append(len(cvu.defined(thisvl)))
#
#         ss = np.array(shedding_dur['symp'])
#         sa = np.array(shedding_dur['asymp'])
#
#     if do_plot:
#         # Make a heat map
#         pl.figure(figsize=(24, 10))
#         font_size = 36
#         font_family = 'Libertinus Sans'
#         pl.rcParams['font.size'] = font_size
#         pl.rcParams['font.family'] = font_family
#         ax = pl.axes([0.05, 0.11, 0.95, 0.85])
#         sns.heatmap(vl[0:30,0:60].T, cmap=sns.cm.rocket_r, vmin=0, cbar_kws={'label': 'Log10 viral load (copies/mL)'})
#         ax.set_yticklabels([])
#         ax.set_xlabel('Simulation days')
#         ax.set_ylabel('People')
#         if do_save:
#             cv.savefig(fig_path, dpi=100)
#
#         if method=='new2':
#             # Make individual plots
#             pl.figure(figsize=(24, 20))
#             xl, xr, yb, ym, yt = 0.07, 0.03, 0.07, 0.07, 0.05
#             dx, dy = 1-xl-xr, (1-yb-ym-yt)/2
#             ax = {}
#
#             # Symptomatic trajectories
#             ax[0] = pl.axes([xl, yb, dx, dy])
#             ax[0].set_title('Symptomatic infections')
#             for count_ind, person_ind in enumerate(symp_inds):
#                 x = np.arange(len(vls2plot['symp'][count_ind]))
#                 ax[0].plot(x, vls2plot['symp'][count_ind], color='grey', linewidth=1)
#             ax[0].set_xlim(0, 20)
#             ax[0].set_ylim(0, 11)
#             ax[0].set_xlabel('Days post exposure')
#             ax[0].set_ylabel('Log10 viral load (copies/mL)')
#             ax[0].text(14, 8, f'Median (95%): {np.median(ss)}', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
#
#             # Asymptomatic trajectories
#             ax[1] = pl.axes([xl, yb+(ym+dy), dx, dy])
#             ax[1].set_title('Asymptomatic infections')
#             for count_ind, person_ind in enumerate(asymp_inds):
#                 x = np.arange(len(vls2plot['asymp'][count_ind]))
#                 ax[1].plot(x, vls2plot['asymp'][count_ind], color='grey', linewidth=1)
#             ax[1].set_xlim(0, 20)
#             ax[1].set_ylim(0, 11)
#             ax[1].set_xlabel('Days post exposure')
#             ax[1].set_ylabel('Log10 viral load (copies/mL)')
#             ax[1].text(14, 8, f'Median (95%): {np.median(sa)}', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
#             if do_save:
#                 cv.savefig('results/sample_vls', dpi=100)
#
#
#     return sim
