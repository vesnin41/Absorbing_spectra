import numpy as np
import scipy.constants as consts
import pandas as pd
from numpy.polynomial import polynomial as P
from scipy.interpolate import interp1d

class MyBulkMaterial:
    '''
    Class for Bulk film
    [T, rho, M, sigma_0, wavelength, N_m]
    '''

    # Constants | Константы
    c = consts.c # Speed of light | Скорость света [м/с]
    Na = consts.Avogadro # Avogadro constant | Постоянная Авагадро [1/моль]
    e = consts.e # The electron charge | Заряд электрона
    eps_0 = consts.epsilon_0 # Vacuum permittivity | Диэлектрическая проницаемость

    def __init__(self, 
                 rho: 'Density array [kg/m^3]', 
                 M: 'Molar mass [kg/mol]', 
                 sigma_0: "The material's electrical conductivity", 
                 wavelength: 'Wavelength array [m]', 
                 N_m: "Elemental N/m * 10^-59",
                 pathfile_n: "n [file csv]",
                 pathfile_k: "k [file csv]",
                 ): 
        self.rho = rho
        self.M = M
        self.sigma_0 = sigma_0
        self.N_m = N_m
        self.wavelength = wavelength
        self.pathfile_n = pathfile_n
        self.pathfile_k = pathfile_k

        self.N = self.free_e_density()
        self.m = self.eff_mass_e()
        self.omega = self.angular_freq()
        self.omega_p = self.plasma_freq()
        self.omega_dless = self.omega_omega_p()

        self.n_1 = self.n_1_LD()
        self.k_1 = self.k_1_LD()
        self.A_bulk = self.A_bulk_LD()

    
    def free_e_density(self):
        '''
        The free electron density
        |
        Плотность свободных электронов
        '''
        
        return (self.rho * self.Na) / self.M

    def eff_mass_e(self):


        '''
        The optical effective mass of the electron
        |
        Оптически эффективная масса электрона
        '''

        return self.free_e_density() / self.N_m

    def angular_freq(self):
        '''
        [omega]
        The angular frequency of the incident light
        | 
        Угловая частота падующего света


        Params:
        wavelenght,

        '''
    

        return 2 * np.pi * (self.c / self.wavelength)


    def plasma_freq(self):
        '''
        [omega_p]
        The plasma frequency of the material
        |
        Плазменная частота материала
        '''


        return ((self.N * self.e**2) / (self.m * self.eps_0))**0.5


    def e_relax_time(self):
        '''
        [tau] 
        The electron relaxation time
        |
        Время релаксации электрона
        '''

        #N = N(self)
        #m = m(self, N)

        return (self.m * self.sigma_0) / (self.N * self.e**2)

    def omega_omega_p(self):
        '''
        ω / ωp
        '''

        return self.omega / self.omega_p[0]


    #------------From Rakic------------------------------------------------
    def n_1_LD(self):
        '''
        The index of refraction
        '''
        # Import from a file
        data = pd.read_csv(self.pathfile_n)

        # Interpolation
        wavelength_p = np.array(data['wl']) * (10**-6)
        n_1 = np.array(data['n'])
        n_1 = interp1d(wavelength_p, n_1, kind='cubic')

        return n_1(self.wavelength)

    def k_1_LD(self):
        '''
        
        '''
        # Import from a file
        data = pd.read_csv(self.pathfile_k)

        # Interpolation
        wavelength_p = np.array(data['wl']) * (10**-6)
        k_1 = np.array(data['k'])
        k_1 = interp1d(wavelength_p, k_1, kind='cubic')

        return k_1(self.wavelength)
    
    def A_bulk_LD(self):
        '''
        Absorptivity for the bulk material from Palik
        '''
        return 1 - ((1 - self.n_1)**2 + self.k_1**2) / ((1 + self.n_1)**2 + self.k_1**2)


class MyThinMaterial(MyBulkMaterial):
    '''
    Class for Thin film
    [h, pathsubstrate]
    '''

    # Constants | Константы
    #c = consts.c # Speed of light | Скорость света [м/с]
    #Na = consts.Avogadro # Avogadro constant | Постоянная Авагадро [1/моль]
    #e = consts.e # The electron charge | Заряд электрона
    #eps_0 = consts.epsilon_0 # Vacuum permittivity | Диэлектрическая проницаемость

    def __init__(self,
                 rho: 'Density array [kg/m^3]', 
                 M: 'Molar mass [kg/mol]', 
                 sigma_0: "The material's electrical conductivity", 
                 wavelength: 'Wavelength array [m]', 
                 N_m: "Elemental N/m * 10^-59",
                 pathfile_n: "n [file csv]",
                 pathfile_k: "k [file csv]",
                 h: 'Film thickness',
                 pathsubstrate: 'Substrate material',
                 averaging_n2_k2: 'True or False'
                 ):
        
        super().__init__(
                rho, 
                M, 
                sigma_0, 
                wavelength,  
                N_m,  
                pathfile_n,
                pathfile_k,
                )
        self.h = h
        self.pathsubstrate = pathsubstrate

        self.h_dless = self.h_dless_()
        self.k_2 = self.k_2_()
        self.n_2 = self.n_2_()
        if averaging_n2_k2 == True:
            self.k_2 = self.k_2_average()
            self.n_2 = self.n_2_average()

        self.mu = self.mu_LD()
        self.beta = self.beta_LD()


        self.A1p = self.A_1_plus_LD()
        self.A1m = self.A_1_minus_LD()
        
        
        self.A2p = self.A_2_plus_LD()
        self.A2m = self.A_2_minus_LD()

        
        self.A3p = self.A_3_plus_LD()
        self.A3m = self.A_3_minus_LD()

        
        self.A4p = self.A_4_plus_LD()
        self.A4m = self.A_4_minus_LD()

        self.A_thin = self.A_thinfilm_LD()

    def h_dless_(self):
        '''
        
        '''
        return self.h * self.omega_p[0] / (2 * np.pi * self.c)

    def k_2_(self):
        # Import data
        dataSiO2 = pd.read_csv(self.pathsubstrate)

        wavelength_2p = np.array(dataSiO2['Wavelength'])
        k_2p = np.array(dataSiO2['k'])

        # Interpolation
        k_2p = interp1d(wavelength_2p, k_2p, kind='cubic')

        return  k_2p(self.wavelength)

    def n_2_(self):
        # Import data
        dataSiO2 = pd.read_csv(self.pathsubstrate)

        wavelength_2p = np.array(dataSiO2['Wavelength'])
        n_2p = np.array(dataSiO2['n'])

        # Interpolation
        n_2p = interp1d(wavelength_2p, n_2p, kind='cubic')

        return  n_2p(self.wavelength)

    def k_2_average(self):
        return np.ones(len(self.k_2))*(self.k_2.sum() / np.size(self.k_2))

    def n_2_average(self):
        return np.ones(len(self.n_2)) * (self.n_2.sum() / np.size(self.n_2))

    def mu_LD(self):
        '''
        
        '''
        mu = 4*np.pi*self.k_1*self.h_dless * self.omega_dless
        return mu

    def beta_LD(self):
        '''
        
        '''
        beta = (4*np.pi*self.n_1*self.h_dless) * self.omega_dless
        return beta


    #Coefficient A1:

    def A_1_plus_LD(self):
        '''
        Coefficient A1+
        |
        Коэффициент A1+
        '''
        A_1p = ((1 + self.n_1)**2 + self.k_1**2) * ((self.n_1 + self.n_2)**2 + (self.k_1 + self.k_2)**2)
        return  A_1p

    def A_1_minus_LD(self):
        '''
        Coefficient A1-
        |
        Коэффициент A1-
        '''
        A_1o = ((1 - self.n_1)**2 + self.k_1**2) * ((self.n_1 + self.n_2)**2 + (self.k_1 + self.k_2)**2)      
        return A_1o
    

    #Coefficient A2:

    def A_2_plus_LD(self):
        '''
        Coefficient A2+
        |
        Коэффициент A2+
        '''
        A_2p = ((1 + self.n_1)**2 + self.k_1**2) * ((self.n_1 - self.n_2)**2 + (self.k_1 - self.k_2)**2)       
        return A_2p

    def A_2_minus_LD(self):
        '''
        Coefficient A2-
        |
        Коэффициент A2-       
        '''
        A_2o = ((1 - self.n_1)**2 + self.k_1**2) * ((self.n_1 - self.n_2)**2 + (self.k_1 - self.k_2)**2)
        return A_2o


    #Coefficient A3:

    def A_3_plus_LD(self):
        '''
        Coefficient A3+
        |
        Коэффициент A3+  
        '''
        A_3p = 2*((1 - self.n_1**2 - self.k_1**2) * (self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2) + 4*self.k_1*(self.n_1*self.k_2 - self.n_2*self.k_1))        
        return A_3p

    def A_3_minus_LD(self):
        '''
        Coefficient A3-
        |
        Коэффициент A3-
        '''
        A_3o = 2*((1 - self.n_1**2 - self.k_1**2) * (self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2) - 4*self.k_1*(self.n_1*self.k_2 - self.n_2*self.k_1))
        return A_3o


    # Coefficient A4:

    def A_4_plus_LD(self):
        '''
        Coefficient A4+
        |
        Коэффициент A4+
        '''
        A_4p = 4*((1 - self.n_1**2 - self.k_1**2) * (self.n_1*self.k_2 - self.n_2*self.k_1) + self.k_1*(self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2))
        return A_4p

    def A_4_minus_LD(self):
        '''
        Coefficient A4- 
        |
        Коэффициент A4-
        '''
        A_4o = 4*((1 - self.n_1**2 - self.k_1**2) * (self.n_1*self.k_2 - self.n_2*self.k_1) - self.k_1*(self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2))
    
        return A_4o

    def A_thinfilm_LD(self):
        R1 = self.A1m*np.exp(self.mu) + self.A2p*np.exp(-self.mu) + self.A3p*np.cos(self.beta) + self.A4m*np.sin(self.beta)
        R2 = self.A1p*np.exp(self.mu) + self.A2m*np.exp(-self.mu) + self.A3m*np.cos(self.beta) + self.A4p*np.sin(self.beta)
        R = R1 / R2
        return 1 - R


def main():
    Titan_10mkm_P = MyThinMaterial(
        rho = np.array([4500, 4490, 4470, 4460, 4450, 4430, 4420, 4400]),
        M = 0.047867,
        sigma_0 = np.array([0.0207e8, 0.0158e8, 0.01227e8, 0.01007e8, 0.00861e8, 0.00762e8, 0.00699e8, 0.00657e8]),
        wavelength =  np.linspace(2e-6, 30e-6, 29),
        N_m = 1.254e59,
        pathfile_n = 'Data/Ti_Rakic-LD_n.csv',
        pathfile_k = 'Data/Ti_Rakic-LD_k.csv',
        h = 10e-9, 
        pathsubstrate = 'Data/SiO2.csv',
        averaging_n2_k2 = False,
        )
    print(Titan_10mkm_P.A_thin)

if __name__ == '__main__':
    main()