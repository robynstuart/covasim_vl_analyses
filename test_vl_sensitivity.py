'''
Exploring the effect of viral load profiles and their relationship to test sensitivity, and infectiousness
'''

# Imports and settings
import os
import sciris as sc
import covasim as cv
import pylab as pl
import seaborn as sns

do_plot   = 1
do_save   = 1

#%% Analyses
def test_vl_sensitivity(do_plot=True, do_save=True, fig_path='results/plot_vl.png', verbose=0.1):
    sc.heading('Test viral loads over time')

    sc.heading('Setting up...')

    sc.tic()

    sim = cv.Sim(pop_size=100) # create sim object
    sim.run(verbose=verbose, keep_people=True)

    vl = sc.promotetoarray(sim.people.viral_log)

    if do_plot:
        pl.figure(figsize=(24, 20))
        font_size = 20
        font_family = 'Libertinus Sans'
        pl.rcParams['font.size'] = font_size
        pl.rcParams['font.family'] = font_family
        ax = pl.axes([0.05, 0.05, 0.9, 0.9])
        ax = sns.heatmap(vl, cmap=sns.cm.rocket_r, vmin=0)
        if do_save:
            cv.savefig(fig_path, dpi=100)

    return sim



#%% Run as a script
if __name__ == '__main__':

    sc.tic()
    sim = test_vl_sensitivity(do_plot=do_plot, do_save=do_save)
    sc.toc()


print('Done.')
