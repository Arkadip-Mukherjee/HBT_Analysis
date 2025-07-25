import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit


kTmin = 0.15
kTmax = 0.55
kT_num_bins = 4
kT_bin_width = (kTmax - kTmin) / kT_num_bins

hbar_c = 0.197326980;

Rout_store, Rside_store, Rlong_store, Ro_Rs_store = [], [], [], []
Routerr_store, Rsideerr_store, Rlongerr_store = [], [], []

chi2_ndf = []
kT_val = []

for ii in range(kT_num_bins):
	kT_val.append(kTmin + (ii + 0.5)*kT_bin_width)


def correlation_function(q, lambda_, Ro, Rs, Rl):
	qo, qs, ql = q
	return 1 + lambda_ * np.exp((-(qo**2 * Ro**2) - (qs**2 * Rs**2) - (ql**2 * Rl**2)) / (hbar_c**2))


def load_data(input_filename):
	data = np.loadtxt(input_filename)
	qo = data[:, 0]
	qs = data[:, 1]
	ql = data[:, 2]
	Cq = data[:, 3]
	Cq_err = data[:, 4]
	return qo, qs, ql, Cq, Cq_err


def fit_correlation_function(input_filename):	
	qo, qs, ql, Cq, Cq_err = load_data(input_filename)
	initial_guess = [1.0, 4.0, 4.0, 5.0]  # [lambda, Rout, Rside, Rlong]
	popt, pcov = curve_fit(correlation_function, (qo, qs, ql), Cq, sigma = Cq_err, p0 = initial_guess, absolute_sigma = True)
	return popt, pcov


def calculate_errors(pcov):
	perr = np.sqrt(np.diag(pcov))
	return perr


def calculate_chi2_ndf(qo, qs, ql, Cq, Cq_err, popt):
	Cq_fit = correlation_function((qo, qs, ql), *popt)
	chi2 = np.sum(((Cq - Cq_fit) / Cq_err)**2)
	ndf = len(Cq) - len(popt)
	return chi2, ndf, chi2/ndf


def plot_results(input_filename, popt):

	qo, qs, ql, Cq, Cq_err = load_data(input_filename)

	Cq_fit = correlation_function((qo, qs, ql), *popt)

	fig, ax = plt.subplots(1, 3, figsize = (15, 5))

	ax[0].scatter(qo, Cq, label = 'Data', color = 'blue', alpha = 0.5)
	ax[0].scatter(qo, Cq_fit, label = 'Fit', color = 'red', alpha = 0.5)
	ax[0].set_xlabel('q1')
	ax[0].set_ylabel('C(q)')
	ax[0].legend()
	ax[0].set_title('q1 projection')

	ax[1].scatter(qs, Cq, label = 'Data', color = 'blue', alpha = 0.5)
	ax[1].scatter(qs, Cq_fit, label = 'Fit', color = 'red', alpha = 0.5)
	ax[1].set_xlabel('q2')
	ax[1].legend()
	ax[1].set_title('q2 projection')

	ax[2].scatter(ql, Cq, label = 'Data', color = 'blue', alpha = 0.5)
	ax[2].scatter(ql, Cq_fit, label = 'Fit', color = 'red', alpha = 0.5)
	ax[2].set_xlabel('q3')
	ax[2].legend()
	ax[2].set_title('q3 projection')

	plt.savefig('Cq3D_Plot.png', format = 'png', bbox_inches = 'tight', dpi = 800)


if __name__ == "__main__":
	
	for ii in range(kT_num_bins):
		input_filename = f"Cq3D_kT_{ii}.data"  
		
		qo, qs, ql, Cq, Cq_err = load_data(input_filename)

		popt, pcov = fit_correlation_function(input_filename)

		perr = calculate_errors(pcov)

		lambda_opt, r1_opt, r2_opt, r3_opt = popt
		lambda_err, r1_err, r2_err, r3_err = perr

		plot_results(input_filename, popt)
		
		Rout_store.append(r1_opt)
		Rside_store.append(r2_opt)
		Rlong_store.append(r3_opt-1)
		Ro_Rs_store.append((r1_opt/r2_opt))
		
		Routerr_store.append(r1_err)
		Rsideerr_store.append(r2_err)
		Rlongerr_store.append(r3_err)

		chisq, ndf, chisq_ndf = calculate_chi2_ndf(qo, qs, ql, Cq, Cq_err, popt)
		
		chi2_ndf.append(chisq_ndf)
		
		print(f"\nFitting results for kTbin-{ii}:")
		print(f"λ = {lambda_opt:.6f} ± {lambda_err:.6f}")
		print(f"R_out  = {r1_opt:.6f} ± {r1_err:.6f}")
		print(f"R_side = {r2_opt:.6f} ± {r2_err:.6f}")
		print(f"R_long = {r3_opt:.6f} ± {r3_err:.6f}")
		print(f"χ² = {chisq:.2f}")
		print(f"ndf = {ndf:.2f}")
		print(f"χ²/ndf = {chisq_ndf:.2f}")
		
	fig, ax = plt.subplots(2, 3, figsize = (15, 10))
	ax = ax.flatten()
	
	
	ax[0].minorticks_on()
	
	ax[0].tick_params(axis = 'both', which = 'major', labelsize = 10, direction = 'inout',
		              left = True, right = True, top = True, width = 1.5, length = 4.0)
	ax[0].tick_params(axis = 'both', which = 'minor', labelsize = 10, direction = 'inout',
		              left = True, right = True, top = True, width = 1.5, length = 2.0)
	
	ax[0].plot(kT_val, Rout_store, marker = 'o', linestyle = '-', linewidth = 2, color = 'red', label = "Python")	
	ax[0].errorbar(kT_val, Rout_store, yerr= Routerr_store, fmt='.', markersize = 1,
	               ecolor='red', elinewidth = 0.8, capsize=2.5, capthick=0.7)
	
	ax[0].set_xlabel(r"$k_{T}$ [GeV/c]", fontsize = 20)
	ax[0].set_ylabel(r"$R_{out}$ [fm]", fontsize = 20)
	ax[0].set_xlim(0.15, 0.55)
	ax[0].set_ylim(3.55, 6.55)
	ax[0].legend(loc = "best", fontsize = 15)
	
	
	ax[1].minorticks_on()

	ax[1].tick_params(axis = 'both', which = 'major', labelsize = 10, direction = 'inout',
		         left = True, right = True, top = True, width = 1.5, length = 4.0)
	ax[1].tick_params(axis = 'both', which = 'minor', labelsize = 10, direction = 'inout',
		         left = True, right = True, top = True, width = 1.5, length = 2.0)
	
	ax[1].plot(kT_val, Rside_store, marker = 'o', linestyle = '-', linewidth = 2, color = 'blue', label = "Python")	
	ax[1].errorbar(kT_val, Rside_store, yerr= Rsideerr_store, fmt='.', markersize = 1,
	               ecolor='blue', elinewidth = 0.8, capsize=2.5, capthick=0.7)

	ax[1].set_xlabel(r"$k_{T}$ [GeV/c]", fontsize = 20)
	ax[1].set_ylabel(r"$R_{side}$ [fm]", fontsize = 20)
	ax[1].set_xlim(0.15, 0.55)
	ax[1].set_ylim(3.55, 6.55)
	ax[1].legend(loc = "best", fontsize = 15)
	
	ax[2].minorticks_on()
	
	ax[2].tick_params(axis = 'both', which = 'major', labelsize = 10, direction = 'inout',
		         left = True, right = True, top = True, width = 1.5, length = 4.0)
	ax[2].tick_params(axis = 'both', which = 'minor', labelsize = 10, direction = 'inout',
		         left = True, right = True, top = True, width = 1.5, length = 2.0)
	
	ax[2].plot(kT_val, Rlong_store, marker = 'o', linestyle = '-', linewidth = 2, color = 'green', label = "Python")
	ax[2].errorbar(kT_val, Rlong_store, yerr= Rlongerr_store, fmt='.', markersize = 1,
	               ecolor='green', elinewidth = 0.8, capsize=2.5, capthick=0.7)

	ax[2].set_xlabel(r"$k_{T}$ [GeV/c]", fontsize = 20)
	ax[2].set_ylabel(r"$R_{long}$ [fm]", fontsize = 20)
	ax[2].set_xlim(0.15, 0.55)
	ax[2].set_ylim(3.55, 8.55)
	ax[2].legend(loc = "best", fontsize = 15)
	
	
	ax[3].minorticks_on()
	
	ax[3].tick_params(axis = 'both', which = 'major', labelsize = 10, direction = 'inout',
		         left = True, right = True, top = True, width = 1.5, length = 4.0)
	ax[3].tick_params(axis = 'both', which = 'minor', labelsize = 10, direction = 'inout',
		         left = True, right = True, top = True, width = 1.5, length = 2.0)
	
	ax[3].plot(kT_val, Ro_Rs_store, marker = 'o', linestyle = '-', linewidth = 2, color = 'purple', label = "Python")

	ax[3].axhline(y = 1, color = 'k', linestyle = '--', alpha = 0.5)
	ax[3].set_xlabel(r"$k_{T}$ [GeV/c]", fontsize = 20)
	ax[3].set_ylabel(r"$R_{out}$ / $R_{side}$ [fm]", fontsize = 20)
	ax[3].set_xlim(0.15, 0.55)
	ax[3].set_ylim(0.8, 1.2)
	ax[3].legend(loc = "best", fontsize = 15)


  ax[4].minorticks_on()
	
	ax[4].tick_params(axis = 'both', which = 'major', labelsize = 10, direction = 'inout',
		              left = True, right = True, top = True, width = 1.5, length = 4.0)
	ax[4].tick_params(axis = 'both', which = 'minor', labelsize = 10, direction = 'inout',
		              left = True, right = True, top = True, width = 1.5, length = 2.0)
	
	ax[4].plot(kT_val, chi2_ndf, marker = 'o', linestyle = '-', linewidth = 2, color = 'purple', label = "Python")		
	
	ax[4].set_xlabel(r"$k_{T}$ [GeV/c]", fontsize = 20)
	ax[4].set_ylabel(r"$\chi^{2}$/ndf", fontsize = 20)
	ax[4].set_xlim(0.15, 0.55)
	ax[4].set_ylim(0.0, 2.0)
	ax[4].legend(loc = "best", fontsize = 15)
	

  fig.delaxes(ax[5])


	plt.savefig('HBT_Plot_.png', format = 'png', bbox_inches = 'tight', dpi = 800)
	print("\nPlots generated succesfully!\n")
