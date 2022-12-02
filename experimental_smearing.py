import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy import signal
import b2plot as bp
plt.style.use("belle2")

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

smearing_dist = stats.crystalball(loc=0, scale=0.03,beta=0.7,m=3)
smearing_dist2 = stats.norm(loc=0, scale=0.025,)
std = 0.01
normal_dist = stats.cauchy(loc=0, scale=std)

## "True (experimental) smearing
smearing_dist12 = stats.crystalball(loc=-0.01, scale=0.035,beta=0.6,m=3.2)
smearing_dist22 = stats.norm(loc=-0.01, scale=0.03,)


delta = 1e-4
big_grid = np.arange(-0.5,0.5,delta)

pmf1 = smearing_dist.pdf(big_grid)*delta
pmf12 = smearing_dist12.pdf(big_grid)*delta
print("Sum of uniform pmf: "+str(sum(pmf1)))

pmf2 = normal_dist.pdf(big_grid)*delta
print("Sum of normal pmf: "+str(sum(pmf2)))

pmf3 = smearing_dist2.pdf(big_grid)*delta
pmf32 = smearing_dist22.pdf(big_grid)*delta

conv_pmf = signal.fftconvolve(pmf1,pmf2,'same')
conv_pmf2 = signal.fftconvolve(pmf12,pmf2,'same')
conv_pmf = signal.fftconvolve(conv_pmf,pmf3,'same')
conv_pmf2 = signal.fftconvolve(conv_pmf2,pmf32,'same')
print("Sum of convoluted pmf: "+str(sum(conv_pmf)))

pdf1 = pmf1/delta
pdf2 = pmf2/delta
pdf3 = pmf3/delta
conv_pdf = conv_pmf/delta
conv_pdf2 = conv_pmf2/delta
conv_pdf2[np.abs(conv_pdf2)<1e-8]=0
print("Integration of convoluted pdf: " + str(np.trapz(conv_pdf, big_grid)))


plt.plot(big_grid,pdf2, label="True distribution", color=CB_color_cycle[0])
plt.plot(big_grid,pdf1, label="Smearing effect 1",ls="dashed", color=CB_color_cycle[2])
plt.plot(big_grid,pdf3, label='Smearing effect 2',ls="dotted", color=CB_color_cycle[4])
plt.plot(big_grid,conv_pdf, label="Experimental distribution", color=CB_color_cycle[6])
plt.xlabel("Observable")
plt.ylabel("Probability density")
plt.xlim(-0.5,0.5)
plt.ylim(-0,)
plt.legend(loc='upper left')
plt.savefig("experimental_smearing.pdf", bbox_inches="tight")

random_draw1 = np.random.choice(big_grid, size=100000, p=pdf2/np.sum(pdf2))
random_draw2 = np.random.choice(big_grid, size=100000, p=conv_pdf/np.sum(conv_pdf))
random_draw25 = np.random.choice(big_grid, size=100, p=(conv_pdf2)/np.sum(conv_pdf2))
random_draw3 = np.random.choice(big_grid, size=10000, p=(conv_pdf2)/np.sum(conv_pdf2))

bins=np.linspace(-0.1, 0.1, 100)
plt.close()

# Make 2 axes objects with height ratio 10/3, a space of 0.1 between them, that share the x axis.
# fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [10, 3], 'hspace': 0.1}, sharex=True)



# n_matrix, _, _ = np.histogram2d(x=random_draw1, y=random_draw2, bins=bins)

# x = np.linalg.norm(n_matrix, ord=1, axis=1)
# mig_matrix = n_matrix/np.expand_dims(x, axis=1)
# plt.close()
# plt.imshow(mig_matrix)
# plt.show()
#
#
# binned1,_,_ = bp.hist(random_draw1, bins=bins, ax=ax[0], scale=len(random_draw25)/len(random_draw1))
# binned2,_,_ = bp.hist(random_draw2, bins=bins, ax=ax[0], scale=len(random_draw25)/len(random_draw2))
# binned3,_,_ = bp.hist(random_draw3, bins=bins, ax=ax[0], scale=len(random_draw25)/len(random_draw3))
# binned25,_ = np.histogram(random_draw25, bins=bins) 
#
# unfolded = mig_matrix.dot(binned25)
# unfolded_err = mig_matrix.dot(np.sqrt(binned25))
#
# _ = ax[0].errorbar(bp.bc(bins), unfolded, yerr=unfolded_err, fmt="k.")
# _ = ax[1].errorbar(bp.bc(bins), unfolded/binned1, yerr=np.sqrt(0), fmt="k.")
#
#
# # Afterwards you can style your axes the way you like, here I just do some basic formatting:
# ax[0].set_xlim(-0.1,0.1)
# ax[1].set_ylim(0.9,1.1)
# ax[0].legend(loc='best')
# ax[0].set_ylabel('Candidates')
# ax[1].set_xlabel('Variable')
# ax[1].set_ylabel('Ratio')
#
#
# plt.show()
