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

    #------------Block for Bulk-------------------------------------
    @property
    def N_e(self):
        '''
        The free electron density
        |
        Плотность свободных электронов
        '''
        
        return (self.rho * self.Na) / self.M

    @property
    def m_e(self):


        '''
        The optical effective mass of the electron
        |
        Оптически эффективная масса электрона
        '''

        return self.N_e / self.N_m

    @property
    def omega(self):
        '''
        [omega]
        The angular frequency of the incident light
        | 
        Угловая частота падующего света


        Params:
        wavelenght,

        '''
    

        return 2 * np.pi * (self.c / self.wavelength)

    @property
    def omega_p(self):
        '''
        [omega_p]
        The plasma frequency of the material
        |
        Плазменная частота материала
        '''


        return ((self.N_e * self.e**2) / (self.m_e * self.eps_0))**0.5

    @property
    def tau(self):
        '''
        [tau] 
        The electron relaxation time
        |
        Время релаксации электрона
        '''

        #N = N(self)
        #m = m(self, N)

        return (self.m_e * self.sigma_0) / (self.N_e * self.e**2)

    @property
    def Q(self):
        '''
        
        '''
        Q = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
        for i in range(len(self.wavelength)):
            Q[i,:] = self.omega_p**2 / (self.omega[i]**2 + self.tau**-2)
        return Q

    @property
    def n_1(self):
        '''
        The index of refraction
        '''

        n = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
        for i in range(len(self.wavelength)):
            n[i,:] = (1/np.sqrt(2)) * (((1 - self.Q[i,:])**2 + (self.Q[i,:] / (self.omega[i] * self.tau))**2)**0.5 - self.Q[i,:] + 1)**0.5
        return n

    @property
    def k_1(self):
        '''
        The extinction coefficient
        '''

        k = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
        for i in range(len(self.wavelength)):
            k[i,:] = (1/np.sqrt(2)) * (((1 - self.Q[i,:])**2 + (self.Q[i,:] / (self.omega[i] * self.tau))**2)**0.5 + self.Q[i,:] - 1)**0.5 
        return k

    @property
    def A_bulk(self):
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
                 ):
        
        super().__init__(T, 
                 rho, 
                 M, 
                 sigma_0, 
                 wavelength, 
                 N_m,
                 )
        self.h = h
        self.pathsubstrate = pathsubstrate


    @property
    def k_2(self):
        # Import data
        dataSiO2 = pd.read_csv(self.pathsubstrate)

        wavelength_2p = np.array(dataSiO2['Wavelength'])
        k_2p = np.array(dataSiO2['k'])

        # Interpolation
        k_2p = interp1d(wavelength_2p, k_2p, kind='cubic')

        return  k_2p(self.wavelength)

    @property
    def n_2(self):
        # Import data
        dataSiO2 = pd.read_csv(self.pathsubstrate)

        wavelength_2p = np.array(dataSiO2['Wavelength'])
        n_2p = np.array(dataSiO2['n'])

        # Interpolation
        n_2p = interp1d(wavelength_2p, n_2p, kind='cubic')

        return  n_2p(self.wavelength)

    @property
    def mu(self):
        '''
        
        '''
        mu = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
        for i in range(len(self.wavelength)):
            mu[i,:] = (4*np.pi*self.k_1[i,:]*self.h) / self.wavelength[i]
        return mu

    @property
    def beta(self):
        '''
        
        '''
        beta = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
        for i in range(len(self.wavelength)):
            beta[i,:] = (4*np.pi*self.n_1[i,:]*self.h) / self.wavelength[i]
        return beta



    #Coefficient A1:

    @property
    def A_1p(self):
        '''
        Coefficient A1+
        |
        Коэффициент A1+
        '''
        A_1p = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            #for j in range(len(T)):
            A_1p[i,:] = ((1 + self.n_1[i,:])**2 + self.k_1[i,:]**2) * ((self.n_1[i,:] + self.n_2[i])**2 +  (self.k_1[i,:] + self.k_2[i])**2)
        return  A_1p

    @property
    def A_1m(self):
        '''
        Coefficient A1-
        |
        Коэффициент A1-
        '''
        A_1o = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_1o[i,j] = ((1 - self.n_1[i,j])**2 + self.k_1[i,j]**2) * ((self.n_1[i,j] + self.n_2[i])**2 + (self.k_1[i,j] + self.k_2[i])**2)
        
        return A_1o

    #Coefficient A2:

    @property
    def A_2p(self):
        '''
        Coefficient A2+
        |
        Коэффициент A2+
        '''
        A_2p = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_2p[i,j] = ((1 + self.n_1[i,j])**2 + self.k_1[i,j]**2) * ((self.n_1[i,j] - self.n_2[i])**2 + (self.k_1[i,j] - self.k_2[i])**2)
                
        return A_2p

    @property
    def A_2m(self):
        '''
        Coefficient A2-
        |
        Коэффициент A2-
        '''
        A_2o = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_2o[i,j] = ((1 - self.n_1[i,j])**2 + self.k_1[i,j]**2) * ((self.n_1[i,j] - self.n_2[i])**2 +  (self.k_1[i,j] - self.k_2[i])**2)
        return A_2o

    #Coefficient A3:

    @property
    def A_3p(self):
        '''
        Coefficient A3+
        |
        Коэффициент A3+
        '''
        A_3p = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_3p[i,j] = 2*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2) + 4*self.k_1[i,j]*(self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]))
                
        return A_3p

    @property
    def A_3m(self):
        '''
        Coefficient A3-
        |
        Коэффициент A3-
        '''
        A_3o = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_3o[i,j] = 2*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2) - 4*self.k_1[i,j]*(self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]))
        return A_3o

    # Coefficient A4:

    @property
    def A_4p(self):
        '''
        Coefficient A4+
        |
        Коэффициент A4+
        '''
        A_4p = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_4p[i,j] = 4*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]) + self.k_1[i,j]*(self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2))
                     
        return A_4p

    @property
    def A_4m(self):
        '''
        Coefficient A4-
        |
        Коэффициент A4-
        '''
        A_4o = np.zeros((len(self.wavelength), len(self.T)), dtype=float)
    
        for i in range(len(self.wavelength)):
            for j in range(len(self.T)):
                A_4o[i,j] = 4*((1 - self.n_1[i,j]**2 - self.k_1[i,j]**2) * (self.n_1[i,j]*self.k_2[i] - self.n_2[i]*self.k_1[i,j]) - self.k_1[i,j]*(self.n_1[i,j]**2 + self.k_1[i,j]**2 - self.n_2[i]**2 - self.k_2[i]**2))
        
        return A_4o


    @property
    def R_thinfilm(self):
        '''
        Absorptivity of the thin film
        '''
        R1 = self.A_1m*np.exp(self.mu) + self.A_2p*np.exp(-self.mu) + self.A_3p*np.cos(self.beta) + self.A_4m*np.sin(self.beta)
        R2 = self.A_1p*np.exp(self.mu) + self.A_2m*np.exp(-self.mu) + self.A_3m*np.cos(self.beta) + self.A_4p*np.sin(self.beta)
        R = R1 / R2
        return R

    @property
    def A_thinfilm(self):
        '''
        Absorptivity of the thin film
        '''
        return 1 - self.R_thinfilm

    


def main():
    Titan_Bulk = MyBulkMaterial(
        T = np.array([300, 400, 500, 600, 700, 800, 900, 1000]),
        rho = np.array([4500, 4490, 4470, 4460, 4450, 4430, 4420, 4400]),
        M = 0.047867,
        sigma_0 = np.array([0.0207e8, 0.0158e8, 0.01227e8, 0.01007e8, 0.00861e8, 0.00762e8, 0.00699e8, 0.00657e8]),
        wavelength =  np.linspace(2e-6, 30e-6, 29),
        N_m = 1.254e59,
    )
    print('Hello')
    print(Titan_Bulk.omega)
    print('***************************************************')
    
    print(Titan_Bulk.omega_p)
    print('***************************************************')
    print(Titan_Bulk.tau)
    print('***************************************************')
    print(Titan_Bulk.N_e)
    print('***************************************************')
    print(Titan_Bulk.m_e)
    print("Q = {}".format(Titan_Bulk.Q))
    print('***************************************************')

    print('Omega_p = {}'.format(6.582e-16 * Titan_Bulk.omega_p))
    
    Titan_10mkm = MyThinMaterial(
        T = np.array([300, 400, 500, 600, 700, 800, 900, 1000]),
        rho = np.array([4500, 4490, 4470, 4460, 4450, 4430, 4420, 4400]),
        M = 0.047867,
        sigma_0 = np.array([0.0207e8, 0.0158e8, 0.01227e8, 0.01007e8, 0.00861e8, 0.00762e8, 0.00699e8, 0.00657e8]),
        wavelength =  np.linspace(2e-6, 30e-6, 29),
        N_m = 1.254e59,
        h = 10e-9, 
        pathsubstrate = 'Data/SiO2.csv')
    print(Titan_10mkm.A_thinfilm)



if __name__ == '__main__':
    main()