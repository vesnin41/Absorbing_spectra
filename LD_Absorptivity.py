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
                 pathfile_n: "n from Rakic [file csv]",
                 pathfile_k: "k from Rakic [file csv]",
                 ): 
        
        self.wavelength = wavelength
        self.pathfile_n = pathfile_n
        self.pathfile_k = pathfile_k


        self.n_1 = self.n_1_LD()
        self.k_1 = self.k_1_LD()
        self.A_bulk = self.A_bulk_LD()

    #------------Block for Bulk-------------------------------------
    
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
                 wavelength: 'Wavelength array [m]', 
                 pathfile_n: "n from Rakic [file csv]",
                 pathfile_k: "k from Rakic [file csv]",
                 h: 'Film thickness',
                 pathsubstrate: 'Substrate material',
                 ):
        
        super().__init__( 
                 wavelength,
                 pathfile_n,
                 pathfile_k,
                 )

        self.h = h
        self.pathsubstrate = pathsubstrate

        self.k_2 = self.k_2_()
        self.n_2 = self.n_2_()
        self.mu = self.mu_()
        self.beta = self.beta_()

        self.A1p = self.A_1_plus_()
        self.A1m = self.A_1_minus_()
        
        self.A2p = self.A_2_plus_()
        self.A2m = self.A_2_minus_()

        self.A3p = self.A_3_plus_()
        self.A3m = self.A_3_minus_()


        self.A4p = self.A_4_plus_()
        self.A4m = self.A_4_minus_()

        self.A_thin = self.A_thinfilm_()

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


    def mu_(self):
        '''
        
        '''
        mu = (4*np.pi*self.k_1*self.h) / self.wavelength
        return mu


    def beta_(self):
        '''
        
        '''
        beta = (4*np.pi*self.n_1*self.h) / self.wavelength
        return beta


    #Coefficient A1:

    def A_1_plus_(self):
        '''
        Coefficient A1+ for LD
        |
        Коэффициент A1+ 
        '''
        A_1p = ((1 + self.n_1)**2 + self.k_1**2) * ((self.n_1 + self.n_2)**2 + (self.k_1 + self.k_2)**2)
        return  A_1p

    def A_1_minus_(self):
        '''
        Coefficient A1-
        |
        Коэффициент A1- 
        '''
        A_1o = ((1 - self.n_1)**2 + self.k_1**2) * ((self.n_1 + self.n_2)**2 + (self.k_1 + self.k_2)**2)      
        return A_1o
    

    #Coefficient A2:

    def A_2_plus_(self):
        '''
        Coefficient A2+
        |
        Коэффициент A2+
        '''
        A_2p = ((1 + self.n_1)**2 + self.k_1**2) * ((self.n_1 - self.n_2)**2 + (self.k_1 - self.k_2)**2)       
        return A_2p

    def A_2_minus_(self):
        '''
        Coefficient A2-
        |
        Коэффициент A2-       
        '''
        A_2o = ((1 - self.n_1)**2 + self.k_1**2) * ((self.n_1 - self.n_2)**2 + (self.k_1 - self.k_2)**2)
        return A_2o


    #Coefficient A3:

    def A_3_plus_(self):
        '''
        Coefficient A3+
        |
        Коэффициент A3+     
        '''
        A_3p = 2*((1 - self.n_1**2 - self.k_1**2) * (self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2) + 4*self.k_1*(self.n_1*self.k_2 - self.n_2*self.k_1))        
        return A_3p

    def A_3_minus_(self):
        '''
        Coefficient A3-
        |
        Коэффициент A3-
        '''
        A_3o = 2*((1 - self.n_1**2 - self.k_1**2) * (self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2) - 4*self.k_1*(self.n_1*self.k_2 - self.n_2*self.k_1))
        return A_3o


    # Coefficient A4:

    def A_4_plus_(self):
        '''
        Coefficient A4+
        |
        Коэффициент A4+
        '''
        A_4p = 4*((1 - self.n_1**2 - self.k_1**2) * (self.n_1*self.k_2 - self.n_2*self.k_1) + self.k_1*(self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2))
        return A_4p

    def A_4_minus_(self):
        '''
        Coefficient A4-
        |
        Коэффициент A4-
        '''
        A_4o = 4*((1 - self.n_1**2 - self.k_1**2) * (self.n_1*self.k_2 - self.n_2*self.k_1) - self.k_1*(self.n_1**2 + self.k_1**2 - self.n_2**2 - self.k_2**2))
    
        return A_4o


    def A_thinfilm_(self):
        R1 = self.A1m*np.exp(self.mu) + self.A2p*np.exp(-self.mu) + self.A3p*np.cos(self.beta) + self.A4m*np.sin(self.beta)
        R2 = self.A1p*np.exp(self.mu) + self.A2m*np.exp(-self.mu) + self.A3m*np.cos(self.beta) + self.A4p*np.sin(self.beta)
        R = R1 / R2
        return 1 - R


def main():
    Titan_Bulk_LD = MyBulkMaterial(
        wavelength =  np.linspace(2e-6, 30e-6, 29),
        pathfile_n = 'Data/Ti_Rakic-LD_n.csv',
        pathfile_k = 'Data/Ti_Rakic-LD_k.csv',
)
    print(Titan_Bulk_LD.n_1)
if __name__ == '__main__':
    main()