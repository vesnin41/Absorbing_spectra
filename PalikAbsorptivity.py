import numpy as np
import scipy.constants as consts
import pandas as pd
from numpy.polynomial import polynomial as P
from scipy.interpolate import interp1d

class MyBulkMaterial:
    '''
    Class for Bulk film
    '''

    # Constants | Константы
    c = consts.c # Speed of light | Скорость света [м/с]
    Na = consts.Avogadro # Avogadro constant | Постоянная Авагадро [1/моль]
    e = consts.e # The electron charge | Заряд электрона
    eps_0 = consts.epsilon_0 # Vacuum permittivity | Диэлектрическая проницаемость

    def __init__(self, 
                 wavelength: 'Wavelength array [m]', 
                 pathfile_n_k: "n, k from Palik [file csv]",
                 ): 
        
        self.wavelength = wavelength
        self.pathfile_n_k = pathfile_n_k


        self.n_1 = self.n_1_Palik()
        self.k_1 = self.k_1_Palik()
        self.A_bulk = self.A_bulk_Palik()

    #------------Block for Bulk-------------------------------------
    
    #------------From Palik------------------------------------------------
    def n_1_Palik(self):
        '''
        The index of refraction
        '''
        # Import from a file
        data = pd.read_csv(self.pathfile_n_k)

        # Interpolation
        wavelength_p = np.array(data['wl'])
        n_1 = np.array(data['n'])
        n_1 = interp1d(wavelength_p, n_1, kind='cubic')

        return n_1(self.wavelength)

    def k_1_Palik(self):
        '''
        
        '''
        # Import from a file
        data = pd.read_csv(self.pathfile_n_k)

        # Interpolation
        wavelength_p = np.array(data['wl'])
        k_1 = np.array(data['k'])
        k_1 = interp1d(wavelength_p, k_1, kind='cubic')

        return k_1(self.wavelength)
    
    def A_bulk_Palik(self):
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
                 wavelength: 'Wavelength array [m]', 
                 pathfile_n_k: "n, k from Palik [file csv]",
                 h: 'Film thickness',
                 pathsubstrate: 'Substrate material',
                 ):
        
        super().__init__( 
                 wavelength,
                 pathfile_n_k)

        self.h = h
        self.pathsubstrate = pathsubstrate

        self.k_2 = self.k_2_()
        self.n_2 = self.n_2_()
        self.mu_p = self.mu_Palik()
        self.beta_p = self.beta_Palik()

        self.A1p_p = self.A_1_plus_Palik()
        self.A1m_p = self.A_1_minus_Palik()
        
        self.A2p_p = self.A_2_plus_Palik()
        self.A2m_p = self.A_2_minus_Palik()

        self.A3p_p = self.A_3_plus_Palik()
        self.A3m_p = self.A_3_minus_Palik()


        self.A4p_p = self.A_4_plus_Palik()
        self.A4m_p = self.A_4_minus_Palik()

        self.A_thin = self.A_thinfilm_Palik()

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


    def mu_Palik(self):
        '''
        
        '''
        mu = (4*np.pi*self.k_1*self.h) / self.wavelength
        return mu


    def beta_Palik(self):
        '''
        
        '''
        beta = (4*np.pi*self.n_1*self.h) / self.wavelength
        return beta


    #Coefficient A1:

    def A_1_plus_Palik(self):
        '''
        Coefficient A1+ for Palik
        |
        Коэффициент A1+ c n и k из Палика
        '''
        A_1p = ((1 + self.n_1)**2 + self.k_1**2) * ((self.n_1 + self.n_2)**2 + (self.k_1 + self.k_2)**2)
        return  A_1p

    def A_1_minus_Palik(self):
        '''
        Coefficient A1- for Palik
        |
        Коэффициент A1- c n и k из Палика
        '''
        A_1o = ((1 - self.n_1)**2 + self.k_1**2) * ((self.n_1 + self.n_2)**2 + (self.k_1 + self.k_2)**2)      
        return A_1o
    

    #Coefficient A2:

    def A_2_plus_Palik(self):
        '''
        Coefficient A2+ for Palik
        |
        Коэффициент A2+ c n и k из Палика
        '''
        A_2p = ((1 + self.n_1)**2 + self.k_1**2) * ((self.n_1 - self.n_2)**2 + (self.k_1 - self.k_2)**2)       
        return A_2p

    def A_2_minus_Palik(self):
        '''
        Coefficient A2- for Palik
        |
        Коэффициент A2- c n и k из Палика        
        '''
        A_2o = ((1 - self.n_1)**2 + self.k_1**2) * ((self.n_1 - self.n_2)**2 + (self.k_1 - self.k_2)**2)
        return A_2o


    #Coefficient A3:

    def A_3_plus_Palik(self):
        '''
        Coefficient A3+ for Palik
        |
        Коэффициент A3+ c n и k из Палика     
        '''
        A_3p = 2*((1 - self.n_1**2 - self.k_1**2) * (self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2) + 4*self.k_1*(self.n_1*self.k_2 - self.n_2*self.k_1))        
        return A_3p

    def A_3_minus_Palik(self):
        '''
        Coefficient A3- for Palik
        |
        Коэффициент A3- c n и k из Палика
        '''
        A_3o = 2*((1 - self.n_1**2 - self.k_1**2) * (self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2) - 4*self.k_1*(self.n_1*self.k_2 - self.n_2*self.k_1))
        return A_3o


    # Coefficient A4:

    def A_4_plus_Palik(self):
        '''
        Coefficient A4+ for Palik
        |
        Коэффициент A4+ c n и k из Палика
        '''
        A_4p = 4*((1 - self.n_1**2 - self.k_1**2) * (self.n_1*self.k_2 - self.n_2*self.k_1) + self.k_1*(self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2))
        return A_4p

    def A_4_minus_Palik(self):
        '''
        Coefficient A4- for Palik
        |
        Коэффициент A4- c n и k из Палика
        '''
        A_4o = 4*((1 - self.n_1**2 - self.k_1**2) * (self.n_1*self.k_2 - self.n_2*self.k_1) - self.k_1*(self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2))
    
        return A_4o


    def A_thinfilm_Palik(self):
        R1 = self.A1m_p*np.exp(self.mu_p) + self.A2p_p*np.exp(-self.mu_p) + self.A3p_p*np.cos(self.beta_p) + self.A4m_p*np.sin(self.beta_p)
        R2 = self.A1p_p*np.exp(self.mu_p) + self.A2m_p*np.exp(-self.mu_p) + self.A3m_p*np.cos(self.beta_p) + self.A4p_p*np.sin(self.beta_p)
        R = R1 / R2
        return 1 - R