#Likelihoods for LB reionisation studies
#
# Nov 2025   - M. Tristram -
import os
from typing import Optional

import astropy.io.fits as fits
import numpy as np
from scipy.linalg import sqrtm

from cobaya.likelihoods.base_classes.InstallableLikelihood import InstallableLikelihood
from cobaya.likelihoods.base_classes.DataSetLikelihood import DataSetLikelihood
from cobaya.likelihoods.base_classes.cmblikes import CMBlikes
from cobaya.likelihoods.base_classes.bao import BAO
from cobaya.log import LoggedError

Tcmb = 2.7255

class _Gaussian(InstallableLikelihood):
    """
    Gaussian likelihood
    based on Planck lite covariance matrix

    Units: Cl, muK^2
    """

    data_folder: Optional[str] = "litebird"
    lmin: Optional[int] = 2
    lmax: Optional[int] = 30

    install_options = {}

    type = "CMB"

    def initialize(self):
        # Data folder
        data_file_path = os.path.normpath(
            getattr(self, "path", None) or os.path.join(self.packages_path, "data")
        )
        self.data_folder = os.path.join(data_file_path, self.data_folder)
        
        # Likelihood modes
        self.cl_used = self.cl_used.split(" ")  # TT,EE,BB,TE,EB
        self.log.info( f"Spectra used: {self.cl_used}")

        # data ranges
        self.data_vector = []
        self.cl_range = {spec:{'lmin':0,'lmax':0} for spec in ['TT','EE','TE']}
        for spec in self.cl_used:
            lmin,lmax = self.data_ranges[spec].split()
            self.cl_range[spec] = {'lmin':int(lmin),'lmax':int(lmax),'nl':int(lmax)-int(lmin)+1}
            self.log.debug(f"l-range {spec}: {lmin}-{lmax}")

        # Compute bin window
        av = np.loadtxt(os.path.join(self.data_folder, self.bin_ref))
        ltt = av[:215]
        lte = av[215:414]
        lee = av[414:]
        ix_tt,self.window_tt = self.BinWindow( ltt, self.cl_range['TT'])
        ix_te,self.window_te = self.BinWindow( lte, self.cl_range['TE'])
        ix_ee,self.window_ee = self.BinWindow( lee, self.cl_range['EE'])

        # Fiducial spectrum Dl (l,TT,EE,TE,BB,phiphi,Tphi,Ephi) from CLASS
        self.log.debug(f"Reading fiducial model ({self.cl_fiducial})")
        clfid = np.loadtxt(os.path.join(self.data_folder, self.cl_fiducial)).T
        self.clfid = dict( (s,np.zeros(int(max(clfid[0]))+1)) for s in ['TT','EE','TE','BB','pp','Tp','Ep'])
        for i,s in enumerate(['TT','EE','TE','BB','pp','Tp','Ep']):
            self.clfid[s][np.array(clfid[0],int)] = clfid[i+1] * (1e6*Tcmb)**2 / (clfid[0]*(clfid[0]+1)/2/np.pi) #Cl, muK2

        #Bin fiducial
        self.binned_data = np.hstack((
            self.clfid['TT'][:self.cl_range['TT']['lmax']+1] @ self.window_tt,
            self.clfid['TE'][:self.cl_range['TE']['lmax']+1] @ self.window_te,
            self.clfid['EE'][:self.cl_range['EE']['lmax']+1] @ self.window_ee,
            ))

        # Covariance
        self.log.debug(f"Reading covariance ({self.cl_cov})")
        clcov = np.loadtxt( os.path.join(self.data_folder, self.cl_cov))
        ixs_ells = np.hstack((ix_tt,len(ltt)+ix_te,len(ltt)+len(lte)+ix_ee))
        self.icov = np.linalg.inv( clcov[ixs_ells,:][:,ixs_ells])
        
        self.log.info("Initialized!")


    def BinWindow( self, ells, lra):
        ixs_ells = []
        ell_bins = []

        start_ell = 30
        for i, mid_ell in enumerate(ells):
            end_ell = 2 * int(mid_ell) - start_ell
            if start_ell >= lra['lmin'] and end_ell <= lra['lmax']:
                ell_bins.append(list(range(start_ell, end_ell + 1)))
                ixs_ells.append(i)
            start_ell = end_ell + 1

        # Derive weight matrices for binning Cells
        n_bins = len(ell_bins)
        weights = np.zeros((lra['lmax']+1, n_bins))
        for i, ell_bins in enumerate(ell_bins):
            weights[ell_bins, i] = 1. / len(ell_bins)

        return np.array(ixs_ells,int), weights
        

    def chi_squared(self, clth):
#        clth = {spec.lower():cl for spec,cl in self.clfid.items()}

        binned_clth = np.hstack((
            clth['tt'][:self.cl_range['TT']['lmax']+1] @ self.window_tt,
            clth['te'][:self.cl_range['TE']['lmax']+1] @ self.window_te,
            clth['ee'][:self.cl_range['EE']['lmax']+1] @ self.window_ee,
            ))

        chi2 = self._fast_chi_squared(self.icov, self.binned_data - binned_clth)

        return chi2

    def logp(self, **data_params):
        Cls = self.provider.get_Cl(ell_factor=False,units='muK2')
        return -0.5 * self.chi_squared(Cls)

    def get_requirements(self):
        # State requisites to the theory code
        return {"Cl": {spec:clra['lmax'] for spec,clra in self.cl_range.items()}}





class _Wishart(InstallableLikelihood):
    """
    Wishart 'Exact' likelihood

    Units: Cl, muK^2
    """

    data_folder: Optional[str] = "litebird"
    fsky: Optional[float] = 0.8
    lmin: Optional[int] = 2
    lmax: Optional[int] = 126

    install_options = {}

    type = "CMB"

    def initialize(self):
        # Data folder
        data_file_path = os.path.normpath(
            getattr(self, "path", None) or os.path.join(self.packages_path, "data")
        )
        self.data_folder = os.path.join(data_file_path, self.data_folder)
        
        # Likelihood approx
        self.cl_used = self.cl_used.split()  # TT,EE,BB,TE,EB
        self.log.info( f"Spectra used: {self.cl_used}")
        self.log.info( f"l-range: {self.lmin}-{self.lmax}")
        self.log.info( f"fsky: {self.fsky}")
        
        # Noise Cl (ell | NlTT | NlEE)
        self.log.debug(f"Reading noise: {self.cl_noise}")
        clnoise = np.loadtxt(os.path.join(self.data_folder, self.cl_noise)).T  #Cl, muK2
        if max( clnoise[0]) < self.lmax:
            raise LoggedError( self.log, f"lmax from noise file too low ({int(max(clnoise[0]))})")
        self.clnoi = dict( (s,np.zeros(int(max(clnoise[0]))+1)) for s in ['TT','EE','BB','TE','TB','EB'])
        for i,s in enumerate(['TT','EE']):
            self.clnoi[s][np.array(clnoise[0],int)] = clnoise[i+1]
        self.clnoi['BB'] = clnoise[2]
        
        # Fiducial spectrum Dl (l,TT,EE,TE,BB,phiphi,Tphi,Ephi) from CLASS
        self.log.debug(f"Reading model: {self.cl_fiducial}")
        clfid = np.loadtxt(os.path.join(self.data_folder, self.cl_fiducial)).T
        if max( clfid[0]) < self.lmax:
            raise LoggedError( self.log, f"lmax from fiducial file too low ({int(max(clfid[0]))})")
        self.clfid = dict( (s,np.zeros(int(max(clfid[0]))+1)) for s in ['TT','EE','BB','TE','TB','EB'])
        for i,s in enumerate(['TT','EE','TE','BB']):
            self.clfid[s][np.array(clfid[0],int)] = clfid[i+1] * (1e6*Tcmb)**2 / (clfid[0] * (clfid[0]+1) / 2 / np.pi) #Cl, muK2
        
        self.data = []
        self.noise= []
        for ell in range(self.lmax+1):
            self.data.append( self.X_to_M([self.clfid[spec][ell] for spec in self.cl_used]))
            self.noise.append( self.X_to_M([self.clnoi[spec][ell] for spec in self.cl_used]))
        
        self.log.info("Initialized!")

    def X_to_M(self, X):
        tags = ['T','E','B']

        M = np.zeros( (3,3) )
        for spec,cl in zip(self.cl_used,X):
            i,j = [tags.index(s) for s in spec]
            M[i,j] = M[j,i] = cl

        icols = [0,1,2]
        if 'TT' not in self.cl_used: icols.remove(0)
        if 'EE' not in self.cl_used: icols.remove(1)
        if 'BB' not in self.cl_used: icols.remove(2)
        M = M[icols,:][:,icols]
        
        return M

    def exact_chi_sq(self, C, Chat, L):
        if C.shape[0] == 1:
            return (
                (2 * L + 1)
                * self.fsky
                * (Chat[0, 0] / C[0, 0] - np.log(Chat[0, 0] / C[0, 0]) - 1)
            )
        else:
            M = np.linalg.inv(C).dot(Chat)
            return (
                (2 * L + 1)
                * self.fsky
                * (np.trace(M) - np.linalg.slogdet(M)[1] - len(C))
            )

    def log_likelihood(self, clth):
        """
        Compute Wishart Likelihood for TT,EE,BB
        """
#        clth = {spec.lower():cl for spec,cl in self.clfid.items()}

        chisq = 0.
        for l in range(self.lmin,self.lmax+1):
            C = self.X_to_M( [clth[spec.lower()][l] for spec in self.cl_used])
#            print( l, self.exact_chi_sq(C+self.noise[l], self.data[l]+self.noise[l], l))
            chisq += self.exact_chi_sq(C+self.noise[l], self.data[l]+self.noise[l], l)
        
        return -0.5 * chisq


    def logp(self, **data_params):
        Dls = self.provider.get_Cl(ell_factor=False,units='muK2')
        return self.log_likelihood(Dls)

    def get_requirements(self):
        # State requisites to the theory code
        return {"Cl": {spec:self.lmax for spec in self.cl_used}}




class _HM(CMBlikes):
    """
    Hamimeche and Lewis likelihood

    Units: Cl, muK^2
    """

    data_folder: Optional[str] = "litebird"
    lmin: Optional[int] = 2
    lmax: Optional[int] = 126
    fsky: Optional[float] = 0.8

    install_options = {}

    type = "CMB"

    def initialize(self):
        # Data folder
        data_file_path = os.path.normpath(
            getattr(self, "path", None) or os.path.join(self.packages_path, "data")
        )
        self.data_folder = os.path.join(data_file_path, self.data_folder)
        
        # Likelihood approx
        self.cl_used = self.cl_used.split(" ")  # TT,EE,BB,TE,EB
        self.log.info( f"Spectra used: {self.cl_used}")
        self.log.info( f"l-range: {self.lmin}-{self.lmax}")
        
        # Fiducial spectrum Dl (l,TT,EE,TE,BB,phiphi,Tphi,Ephi) from CLASS
        self.log.debug(f"Reading fiducial model ({self.cl_fiducial})")
        clfid = np.loadtxt(os.path.join(self.data_folder, self.cl_fiducial)).T
        self.clfid = dict( (s,np.zeros(int(max(clfid[0]))+1)) for s in ['TT','EE','BB','TE','TB','EB'])
        for i,s in enumerate(['TT','EE','TE','BB']):
            self.clfid[s][np.array(clfid[0],int)] = clfid[i+1] * (1e6*Tcmb)**2 / (clfid[0] * (clfid[0]+1) / 2 / np.pi)  #muK2
        
        self.data_matrix = []
        self.sqrt_fid_matrix = []
        for ell in range(self.lmax+1):
            data  = self.X_to_M([self.clfid[spec][ell] for spec in self.cl_used])
            self.data_matrix.append( data)
            self.sqrt_fid_matrix.append( sqrtm(data) )
        
        # Covariance
        clcov = fits.getdata( os.path.join(self.data_folder, self.cl_cov))  #EE,BB,EB, start at l=2
        ndim = len(clcov)//3
        nls = self.lmax-self.lmin+1
        ixs = []
        if 'EE' in self.cl_used: ixs.append( np.arange(self.lmin,self.lmax+1)-2 )
        if 'BB' in self.cl_used: ixs.append( np.arange(self.lmin,self.lmax+1)-2 )
        if 'EB' in self.cl_used: ixs.append( np.arange(self.lmin,self.lmax+1)-2 )
        ixs = np.array(ixs).flatten()
        self.icov = np.linalg.inv( clcov[ixs,:][:,ixs])
        
        self.log.info("Initialized!")


    def X_to_M(self, X):
        tags = ['T','E','B']

        M = np.zeros( (3,3) )
        for spec,cl in zip(self.cl_used,X):
            i,j = [tags.index(s) for s in spec]
            M[i,j] = M[j,i] = cl

        icols = [0,1,2]
        if 'TT' not in self.cl_used: icols.remove(0)
        if 'EE' not in self.cl_used: icols.remove(1)
        if 'BB' not in self.cl_used: icols.remove(2)
        M = M[icols,:][:,icols]
        
        return M.reshape( (len(icols),len(icols)) )

    def M_to_X(self, M):
        tags = ['T','E','B']
        icols = [0,1,2]
        if 'TT' not in self.cl_used: icols.remove(0)
        if 'EE' not in self.cl_used: icols.remove(1)
        if 'BB' not in self.cl_used: icols.remove(2)

        X = np.zeros( len(icols)*(len(icols)+1)//2)
        for i in range(len(icols)):
            X[i] = M[i,i]
        if len(icols) > 1:
            for k,(i,j) in enumerate(np.triu_indices(len(icols),1)):
                X[len(icols)+k] = M[i,j]
          
        return X

    def log_likelihood(self, clth):
        """
        Compute Hamimeche&Lewis Likelihood for TT,EE,BB
        """
#        clth = {spec.lower():cl for spec,cl in self.clfid.items()}

        X = []
        for l in range(self.lmin,self.lmax+1):
            model = [clth[spec.lower()][l] for spec in self.cl_used]
            C = self.X_to_M( model)
            try:
                self.transform(
                    C, self.data_matrix[l], self.sqrt_fid_matrix[l]
                    )
            except np.linalg.LinAlgError:
                self.log.debug("Likelihood computation failed.")
                return -np.inf

            X.append( self.M_to_X(C))
        
#        print( np.array(X).flatten())
        return -0.5 * self._fast_chi_squared(self.icov, np.array(X).flatten())

    def logp(self, **data_params):
        Dls = self.provider.get_Cl(ell_factor=False,units='muK2')
        return self.log_likelihood(Dls)

    def get_requirements(self):
        # State requisites to the theory code
        return {"Cl": {spec:self.lmax for spec in self.cl_used}}




class _AnalyticWishart(InstallableLikelihood):
    """
    Analytical Wishart Distribution
    from a fit on Planck lowT

    Units: Cl, K^2
    """

    data_folder: Optional[str] = "litebird"
    lmin: Optional[int] = 2
    lmax: Optional[int] = 29

    install_options = {}

    type = "CMB"

    def initialize(self):
        # Data folder
        data_file_path = os.path.normpath(
            getattr(self, "path", None) or os.path.join(self.packages_path, "data")
        )
        self.data_folder = os.path.join(data_file_path, self.data_folder)

        # Read Planck PR3 low-ell TT posterior fit parameters
        self.log.debug("Reading Analytical Distribution")
        ell, alpha_ref, beta_ref = np.loadtxt(os.path.join(self.data_folder, self.Wishart_distrib), unpack=True)

        # Read mock data
        # Fiducial spectrum Dl (l,TT,EE,TE,BB,phiphi,Tphi,Ephi) from CLASS
        self.log.debug("Reading model")
        clfid = np.loadtxt(os.path.join(self.data_folder, self.cl_fiducial)).T
        self.fid_cl = np.zeros( int(max(clfid[0]))+1)
        self.fid_cl[np.array(clfid[0],int)] = clfid[1] / (clfid[0] * (clfid[0]+1) / 2 / np.pi)  * Tcmb**2

        # Compute adjusted fit parameters
        self.alpha = alpha_ref.copy()
        self.beta = (self.alpha + 1) * self.fid_cl[self.lmin:self.lmax+1]
    
    def logp(self, **data_params):
        Cls = self.provider.get_Cl(ell_factor=False,units='K2')
        cltt = Cls['tt'][self.lmin:self.lmax+1]
#        cltt = self.fid_cl[self.lmin:self.lmax+1]

        lnl = (-self.alpha-1)*np.log(cltt*(self.alpha + 1)/self.beta) - self.beta/cltt + self.alpha + 1

        return lnl.sum()

    def get_requirements(self):
        # State requisites to the theory code
        return {"Cl": {'tt':self.lmax}}




class _Lensing(InstallableLikelihood):
    """
    Wishart 'Exact' likelihood for Cl_phiphi

    C_l^dd (deflection) = l(l+1) C_l^phi-phi
    C_l^gg/kk (shear/convergence) = 1/4 (l(l+1))^2 C_l^phi-phi

    Units: Cl, muK^2
    """

    data_folder: Optional[str] = "litebird"
    fsky: Optional[float] = 0.4
    lmin: Optional[int] = 2
    lmax: Optional[int] = 5000

    install_options = {}

    type = "CMB"

    def initialize(self):
        # Data folder
        data_file_path = os.path.normpath(
            getattr(self, "path", None) or os.path.join(self.packages_path, "data")
        )
        self.data_folder = os.path.join(data_file_path, self.data_folder)
        
        # Likelihood approx
        self.log.info( f"l-range: {self.lmin}-{self.lmax}")
        self.log.info( f"fsky: {self.fsky}")
        
        # Noise Cl (ell | Nlkk)
        self.log.debug("Reading noise")
        cl = np.loadtxt(os.path.join(self.data_folder, self.cl_noise)).T  #Cl, muK2
        self.clnoi = np.zeros(int(max(cl[0]))+1)
        self.clnoi[np.array(cl[0],int)] = cl[1] / ( (cl[0]*(cl[0]+1))**2/4 )
        
        # Fiducial spectrum Dl (l,TT,EE,TE,BB,phiphi,Tphi,Ephi) from CLASS
        self.log.debug("Reading model")
        clfid = np.loadtxt(os.path.join(self.data_folder, self.cl_fiducial)).T
        self.clfid = np.zeros(int(max(clfid[0]))+1)
        self.clfid[np.array(clfid[0],int)] = clfid[5] / (clfid[0]*(clfid[0]+1)/2/np.pi) #Dl_phiphi -> Cl_phiphi
        
        self.log.info("Initialized!")

    def exact_chi_sq(self, Cl, Clhat, L):
            return (
                (2 * L + 1)
                * self.fsky
                * (Clhat / Cl - np.log(Clhat / Cl) - 1)
            )

    def log_likelihood(self, clth):
        """
        Compute Wishart Likelihood for TT,EE,BB
        """
#        clth = {spec.lower():cl for spec,cl in self.clfid.items()}

        chisq = 0.
        for l in range(self.lmin,self.lmax+1):
            clpp = clth['pp'][l]
            chisq += self.exact_chi_sq(clpp+self.clnoi[l], self.clfid[l]+self.clnoi[l], l)
        
        return -0.5 * chisq


    def logp(self, **data_params):
        Dls = self.provider.get_Cl(ell_factor=False)
        return self.log_likelihood(Dls)

    def get_requirements(self):
        # State requisites to the theory code
        lmax = 10000 #self.lmax
        return {"Cl": {'tt':lmax,'pp':lmax}}



class _BAO(BAO):
    """
    Gaussian likelihood for BAO data
    """

    data_folder: Optional[str] = "litebird"

    install_options = {}

    data_folder: Optional[str] = "litebird"
    rs_fid: Optional[float] = None

    type = "BAO"

    def initialize(self):
        # Data folder
        data_file_path = os.path.normpath(
            getattr(self, "path", None) or os.path.join(self.packages_path, "data")
        )
        self.data_folder = os.path.join(data_file_path, self.data_folder)

        # Fiducial BAO data (Hz*rs, DA/rs)
        self.log.debug(f"Reading fiducial BAO data ({self.bao_data})")
        z,data1,data2,icov11,icov22,icov12 = np.loadtxt(os.path.join(self.data_folder, self.bao_data), usecols=[1,2,3,4,5,6],unpack=True) #dataset,z,data1,data2,icov11,icov22,icov12,8
        self.redshift = z
        self.data = [ np.array([d1,d2]) for d1,d2 in zip(data1,data2)]
        self.icov = [ np.array([i11,i12,i12,i22]).reshape(2,2) for i11,i22,i12 in zip(icov11,icov22,icov12)]

        self.observables = ["Hz_rs","DA_over_rs"]
        self.nobs = len(self.redshift)
        # Rescaling by a fiducial value of the sound horizon

        # Rescaling by a fiducial value of the sound horizon
        if self.rs_fid is not None:
            self.rs_rescale = 1 / self.rs_fid
        else:
            self.rs_rescale = 1
        
        self.log.info("Initialized!")

    def logp(self, **data_params):
        chi2 = 0.
        theo_values = {obs:self.theory_fun( self.redshift, obs) for obs in self.observables}

        for i in range(self.nobs):
            theo = np.array([theo_values[obs][i] for obs in self.observables])
            chi2 += self._fast_chi_squared(self.icov[i], theo - self.data[i])

        return -0.5 * chi2

    def get_requirements(self):
        # State requisites to the theory code
        req = {
            "Hz_rs": {"Hubble": {"z": self.redshift}, "rdrag": None},
            "DA_over_rs": {
                "angular_diameter_distance": {"z": self.redshift},
                "rdrag": None,
                }
            }
        requirements = {}
        for obs in self.observables:
            for k,v in req[obs].items(): requirements[k] = v
        return requirements






class PLK(_Gaussian):
    """
    Planck High-L Likelihood for TT+EE+TE
    Gaussian based on planck_lite covariance
    """

class PLK_lowlT(_AnalyticWishart):
    """
    Planck High-L Likelihood for TT
    Analytical Wishart based on lowT Planck fit
    """

class PLK_lowlE(_HM):
    """
    Planck High-L Likelihood for EE
    Hamimeche&Lewis based on lollipop covariance
    """

class LB(_Wishart):
    """
    LiteBIRD likelihood for TT+EE+BB+TE
    """

class SO(_Wishart):
    """
    Simons-Observatory likelihood for TT+EE+BB+TE
    """

class S4(_Wishart):
    """
    CMB-S4 likelihood for TT+EE+BB+TE
    """

class SOlens(_Lensing):
    """
    Simons-Observatory likelihood for TT+EE+BB+TE
    """

class S4lens(_Lensing):
    """
    CMB-S4 likelihood for TT+EE+BB+TE
    """

class bao(_BAO):
    """
    forecast BAO for DESI+EUCLID
    ref: ??
    """
