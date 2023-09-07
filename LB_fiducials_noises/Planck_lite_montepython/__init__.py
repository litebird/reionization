import numpy as np
from montepython.likelihood_class import Likelihood
import os

class Planck_lite_montepython(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        self.need_cosmo_arguments(
            data, {'lensing': 'yes', 'output': 'tCl lCl pCl'})

        try:
            self.planck_bin_file
        except:
            raise io_mp.LikelihoodError("For reading planck bins from file, you must provide planck_bin_file")

        try:
            self.cov_file
        except:
            raise io_mp.LikelihoodError("For reading covariance from file, you must provide cov_file")

        try:
            self.ell_cut
        except:
            raise io_mp.LikelihoodError("you must provide ell_cut")

        
        try:
            self.T_cmb
        except:
            raise io_mp.LikelihoodError("you must provide T_cmb")

        if os.path.exists(os.path.join(self.data_directory, self.planck_bin_file)):

        # Read Planck lite bin centers
            av = np.loadtxt(os.path.join(
                    self.data_directory, self.planck_bin_file))
            ltt = av[:214]    #last bin of planck dropped, we stay below ell = 2500
            lte = av[214:413]
            lee = av[413:]

            # Deduce Planck lite ell bins from centers, and keep only those above ell_cut
            ixs_ells = []
            # TT bins
            ell_bins_tt = []
            start_ell = 30
            for i, mid_ell in enumerate(ltt):
                end_ell = 2 * int(mid_ell) - start_ell
                if start_ell >= self.ell_cut:
                    ell_bins_tt.append(list(range(start_ell, end_ell + 1)))
                    ixs_ells.append(i)
                start_ell = end_ell + 1
            # TE bins
            ell_bins_te = []
            start_ell = 30
            for i, mid_ell in enumerate(lte):
                end_ell = 2 * int(mid_ell) - start_ell
                if start_ell >= self.ell_cut:
                    ell_bins_te.append(list(range(start_ell, end_ell + 1)))
                    ixs_ells.append(len(ltt) + i)
                start_ell = end_ell + 1
            # EE bins
            ell_bins_ee = []
            start_ell = 30
            for i, mid_ell in enumerate(lee):
                end_ell = 2 * int(mid_ell) - start_ell
                if start_ell >= self.ell_cut:
                    ell_bins_ee.append(list(range(start_ell, end_ell + 1)))
                    ixs_ells.append(len(ltt) + len(lte) + i)
                start_ell = end_ell + 1

            # Derive weight matrices for binning Cells
            n_bins_tt, n_bins_te, n_bins_ee = len(ell_bins_tt), len(ell_bins_te), len(ell_bins_ee)
            n_bins_tot = n_bins_tt + n_bins_te + n_bins_ee
            self.weights_tt = np.zeros((2476, n_bins_tt))  #np.zeros((2509, n_bins_tt)) we dropped the bin [2476,2508]
            self.weights_te = np.zeros((2476, n_bins_te))
            self.weights_ee = np.zeros((2476, n_bins_ee))
            for i, ell_bins in enumerate(ell_bins_tt):
                self.weights_tt[ell_bins, i] = 1. / len(ell_bins)
            for i, ell_bins in enumerate(ell_bins_te):
                self.weights_te[ell_bins, i] = 1. / len(ell_bins)
            for i, ell_bins in enumerate(ell_bins_ee):
                self.weights_ee[ell_bins, i] = 1. / len(ell_bins)


        if os.path.exists(os.path.join(self.data_directory, self.cov_file)):

            # Read, cut, and invert Planck lite covariance matrix
            cov = np.loadtxt(os.path.join(
                        self.data_directory, self.cov_file))
            self.icov = np.linalg.inv(cov[ixs_ells, :][:, ixs_ells])


        # impose that the cosmological code computes Cl's up to maximum l
        # needed by the window function
        self.need_cosmo_arguments(data, {'l_max_scalars': self.l_max})

        ###########################################################################
        # implementation of default settings for flags describing the likelihood: #
        ###########################################################################

        try:
            self.delensing
        except:
            self.delensing = False
        # - do not include lensing extraction by default:
        try:
            self.LensingExtraction
        except:
            self.LensingExtraction = False
        # - neglect TD correlation by default:
        try:
            self.neglect_TD
        except:
            self.neglect_TD = True
        # - use lthe lensed TT, TE, EE by default:
        try:
            self.unlensed_clTTTEEE
        except:
            self.unlensed_clTTTEEE = False
        # - do not exclude TTEE by default:

        ###############################################################
        # Read data for TT, EE, TE, [eventually BB or phi-phi, phi-T] #
        ###############################################################

        #default
        numCls = 3


        # deal with pp, pT (p = CMB lensing potential):
        if self.LensingExtraction:
            self.index_pp = numCls
            numCls += 1


        # deal with fiducial model:
        # If the file exists, initialize the fiducial values
        self.Cl_fid = np.zeros((numCls, 2476), 'float64')
        self.fid_values_exist = False
        if os.path.exists(os.path.join(
                self.data_directory, self.fiducial_file)):
            self.fid_values_exist = True
            fid_file = open(os.path.join(
                self.data_directory, self.fiducial_file), 'r')
            line = fid_file.readline()
            while line.find('#') != -1:
                line = fid_file.readline()
            while (line.find('\n') != -1 and len(line) == 1):
                line = fid_file.readline()
            for l in range(self.l_min, self.l_max+1):
                ll = int(line.split()[0])
                self.Cl_fid[0, ll] = float(line.split()[1])   #tt
                self.Cl_fid[1, ll] = float(line.split()[2])   #ee
                self.Cl_fid[2, ll] = float(line.split()[3])   #te
                # read DD, TD (D = deflection field):
                if self.LensingExtraction:
                    try:
                        self.Cl_fid[self.index_pp, ll] = float(line.split()[self.index_pp+1])
                        if not self.ExcludeTTTEEE:
                            self.Cl_fid[self.index_tp, ll] = float(line.split()[self.index_tp+1])
                    except:
                        raise io_mp.LikelihoodError(
                            "The fiducial model does not have enough columns.")

                line = fid_file.readline()


            l = np.arange(2476)
            self.ell_factor = l*(l+1)/2./np.pi

            self.Cl_fid[0, 2:] /= self.ell_factor[2:]
            self.Cl_fid[1, 2:] /= self.ell_factor[2:]
            self.Cl_fid[2, 2:] /= self.ell_factor[2:]

            # Bin mock data vector removing the ell_factor
            self.binned_fid_cl = np.hstack((
                self.Cl_fid[0] @ self.weights_tt,
                self.Cl_fid[2] @ self.weights_te,
                self.Cl_fid[1] @ self.weights_ee,
            )) * 1e12 * self.T_cmb**2.



        # Else the file will be created in the loglkl() function.

        # Explicitly display the flags to be sure that likelihood does what you expect:
        print("Initialised likelihood_mock_cmb with following options:")
        if self.unlensed_clTTTEEE:
            print("  unlensed_clTTTEEE is True")
        else:
            print("  unlensed_clTTTEEE is False")
        if self.delensing:
            print("  delensing is True")
        else:
            print("  delensing is False")
        if self.LensingExtraction:
            print("  LensingExtraction is True")
        else:
            print("  LensingExtraction is False")
        if self.neglect_TD:
            print("  neglect_TD is True")
        else:
            print("  neglect_TD is False")
        print("")

        # end of initialisation
        return

    def loglkl(self, cosmo, data):

        # get Cl's from the cosmological code (returned in muK**2 units)

        # if we want unlensed Cl's
        if self.unlensed_clTTTEEE:
            cl = self.get_unlensed_cl(cosmo)
            # exception: for non-delensed B modes we need the lensed BB spectrum
            # (this case is usually not useful/relevant)
            # if self.Bmodes and (not self.delensing):
            #         cl_lensed = self.get_cl(cosmo)
            #         for l in range(self.l_min,self.l_max+1):
            #             cl['bb'][l]=cl_lensed['bb'][l]

        # if we want lensed Cl's
        else:
            cl = self.get_cl(cosmo)
            # exception: for delensed B modes we need the unlensed spectrum
            # if self.Bmodes and self.delensing:
            #     cl_unlensed = self.get_unlensed_cl(cosmo)
            #     for l in range(self.l_min,self.l_max+1):
            #             cl['bb'][l]=cl_unlensed['bb'][l]

        # get likelihood
        lkl = self.compute_lkl(cl, cosmo, data)

        return lkl

    def compute_lkl(self, cl, cosmo, data):


        # Theoretical Cells coming from CLASS or CAMB
        # Note: expected to be dimensioneless, no l(l+1) factor, ells=[0, 2508]

       
        cl_theo_tt = np.zeros(2476)
        cl_theo_te = np.zeros(2476)
        cl_theo_ee = np.zeros(2476)


        cl_theo_tt[:self.l_max+1]  = cl['tt']
        cl_theo_te[:self.l_max+1]  = cl['te']
        cl_theo_ee[:self.l_max+1]  = cl['ee']

        # Compute likelihood
        binned_th_cl = np.hstack((
            cl_theo_tt @ self.weights_tt,
            cl_theo_te @ self.weights_te,
            cl_theo_ee @ self.weights_ee,
        )) # cl by class produced as dimensional spectra, without ell factor 


        diff = binned_th_cl - self.binned_fid_cl

        lnlike = -0.5 * diff @ self.icov @ diff

        return lnlike

