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
                 T: 'Temperature array [K]', 
                 rho: 'Density array [kg/m^3]', 
                 M: 'Molar mass [kg/mol]', 
                 sigma_0: "The material's electrical conductivity", 
                 wavelength: 'Wavelength array [m]', 
                 N_m: "Elemental N/m * 10^-59",
                 ): 
        
        self.T = T
        self.rho = rho
        self.M = M
        self.sigma_0 = sigma_0
        self.wavelength = wavelength
        self.N_m = N_m

        self.N = self.free_e_density()
        self.m = self.eff_mass_e()
        self.omega = self.angular_freq()
        self.omega_p = self.plasma_freq()
        self.tau = self.e_relax_time()
        self.omega_dless = self.omega_omega_p()
        self.tau_dless = self.tau_omega_p()
        self.Q = self.Q_()
        self.n_1 = self.n_()
        self.k_1 = self.k_()

        self.A_bulk = self.A_bulk_()

    #------------Block for Bulk-------------------------------------
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

    def tau_omega_p(self):
        '''
        τ*ω_p
        '''

        return self.tau * self.omega_p[0]


    def Q_(self):
        '''
        
        '''
        Q = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
        for i in range(len(self.omega_dless)):
            Q[i,:] = 1 / (self.omega_dless[i]**2 + self.tau_dless**-2)
        return Q

    def n_(self):
        '''
        
        '''

        n = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
        for i in range(len(self.omega_dless)):
            n[i,:] = (1/np.sqrt(2)) * (((1 - self.Q[i,:])**2 + (self.Q[i,:] / (self.omega[i] * self.tau))**2)**0.5 - self.Q[i,:] + 1)**0.5
        return n

    def k_(self):
        '''
        
        '''

        k = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
        for i in range(len(self.omega_dless)):
            k[i,:] = (1/np.sqrt(2)) * (((1 - self.Q[i,:])**2 + (self.Q[i,:] / (self.omega[i] * self.tau))**2)**0.5 + self.Q[i,:] - 1)**0.5 
        return k

    def A_bulk_(self):
        '''
        Absorptivity
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
                 T: 'Temperature array [K]', 
                 rho: 'Density array [kg/m^3]', 
                 M: 'Molar mass [kg/mol]', 
                 sigma_0: "The material's electrical conductivity", 
                 wavelength: 'Wavelength array [m]', 
                 N_m: "Elemental N/m * 10^-59",
                 h: 'Film thickness',
                 pathsubstrate: 'Substrate material',
                 averaging_n2_k2: 'True or False'
                 ):
        
        super().__init__(T, 
                 rho, 
                 M, 
                 sigma_0, 
                 wavelength, 
                 N_m,)
        self.h = h
        self.pathsubstrate = pathsubstrate

        self.h_dless = self.h_dless_()
        self.k_2 = self.k_2_()
        self.n_2 = self.n_2_()
        if averaging_n2_k2 == True:
            self.k_2 = self.k_2_average()
            self.n_2 = self.n_2_average()
        self.mu = self.mu_()
        self.beta = self.beta_()

        self.A1p = self.A_1_plus()
        self.A1m = self.A_1_minus()
        
        self.A2p = self.A_2_plus()
        self.A2m = self.A_2_minus()

        self.A3p = self.A_3_plus()
        self.A3m = self.A_3_minus()

        self.A4p = self.A_4_plus()
        self.A4m = self.A_4_minus()

        self.A_thin = self.A_thinfilm()

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

    def mu_(self):
        '''
        
        '''
        mu = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
        for i in range(len(self.omega_dless)):
            mu[i,:] = 4*np.pi*self.k_1[i,:]*self.h_dless * self.omega_dless[i]
        return mu

    def beta_(self):
        '''
        
        '''
        beta = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
        for i in range(len(self.omega_dless)):
            beta[i,:] = (4*np.pi*self.n_1[i,:]*self.h_dless) * self.omega_dless[i]
        return beta

    #Coefficient A1:

    def A_1_plus(self):
        '''
        Coefficient A1+
        |
        Коэффициент A1+
        '''
        A_1p = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            #for j in range(len(T)):
            A_1p[i,:] = ((1 + self.n_1[i,:])**2 + self.k_1[i,:]**2) * ((self.n_1[i,:] + self.n_2[i])**2 +  (self.k_1[i,:] + self.k_2[i])**2)
        return  A_1p

    def A_1_minus(self):
        '''
        Coefficient A1-
        |
        Коэффициент A1-
        '''
        A_1o = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_1o[i,j] = ((1 - self.n_1[i,j])**2 + self.k_1[i,j]**2) * ((self.n_1[i,j] + self.n_2[i])**2 + (self.k_1[i,j] + self.k_2[i])**2)
        
        return A_1o

    #Coefficient A2:

    def A_2_plus(self):
        '''
        Coefficient A2+
        |
        Коэффициент A2+
        '''
        A_2p = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_2p[i,j] = ((1 + self.n_1[i,j])**2 + self.k_1[i,j]**2) * ((self.n_1[i,j] - self.n_2[i])**2 + (self.k_1[i,j] - self.k_2[i])**2)
                
        return A_2p


    def A_2_minus(self):
        '''
        Coefficient A2-
        |
        Коэффициент A2-
        '''
        A_2o = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_2o[i,j] = ((1 - self.n_1[i,j])**2 + self.k_1[i,j]**2) * ((self.n_1[i,j] - self.n_2[i])**2 +  (self.k_1[i,j] - self.k_2[i])**2)
        return A_2o


    #Coefficient A3:

    def A_3_plus(self):
        '''
        Coefficient A3+
        |
        Коэффициент A3+
        '''
        A_3p = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_3p[i,j] = 2*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2) + 4*self.k_1[i,j]*(self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]))
                
        return A_3p

    def A_3_minus(self):
        '''
        Coefficient A3-
        |
        Коэффициент A3-
        '''
        A_3o = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_3o[i,j] = 2*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2) - 4*self.k_1[i,j]*(self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]))
        return A_3o



    # Coefficient A4:

    def A_4_plus(self):
        '''
        Coefficient A4+
        |
        Коэффициент A4+
        '''
        A_4p = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_4p[i,j] = 4*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]) + self.k_1[i,j]*(self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2))
                     
        return A_4p


    def A_4_minus(self):
        '''
        Coefficient A4-
        |
        Коэффициент A4-
        '''
        A_4o = np.zeros((len(self.omega_dless), len(self.tau_dless)), dtype=float)
    
        for i in range(len(self.omega_dless)):
            for j in range(len(self.tau_dless)):
                A_4o[i,j] = 4*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]) - self.k_1[i,j]*(self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2))
        
        return A_4o


    def A_thinfilm(self):
        '''
        Absorptivity of the thin film
        '''
        R1 = self.A1m*np.exp(self.mu) + self.A2p*np.exp(-self.mu) + self.A3p*np.cos(self.beta) + self.A4m*np.sin(self.beta)
        R2 = self.A1p*np.exp(self.mu) + self.A2m*np.exp(-self.mu) + self.A3m*np.cos(self.beta) + self.A4p*np.sin(self.beta)
        R = R1 / R2
        return 1 - R

def main():
    Titan_10mkm = MyThinMaterial(
        T = np.array([300, 400, 500, 600, 700, 800, 900, 1000]),
        rho = np.array([4500, 4490, 4470, 4460, 4450, 4430, 4420, 4400]),
        M = 0.047867,
        sigma_0 = np.array([0.0207e8, 0.0158e8, 0.01227e8, 0.01007e8, 0.00861e8, 0.00762e8, 0.00699e8, 0.00657e8]),
        wavelength =  np.linspace(2e-6, 30e-6, 29),
        N_m = 1.254e59,
        h = 10e-9, 
        pathsubstrate = 'Data/SiO2.csv',
        averaging_n2_k2 = False,
        )
    print(Titan_10mkm.n_1)

if __name__ == '__main__':
    main()