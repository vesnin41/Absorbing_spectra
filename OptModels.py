import numpy as np
import scipy.constants as consts
import pandas as pd
from numpy.polynomial import polynomial as P
from scipy.interpolate import interp1d
import materials

class DrudeModel():

    # Constants| Константы
    c = consts.c # Speed of light | Скорость света [м/с]
    Na = consts.Avogadro # Avogadro constant | Постоянная Авагадро [1/моль]
    e = consts.e # The electron charge | Заряд электрона
    eps_0 = consts.epsilon_0 # Vacuum permittivity | Диэлектрическая проницаемость

    def __init__(self,
                 wl: 'Wavelength/s | Длины волн',
                 T: 'Temperature | Температура',
                 metal_props: '''M - Molar mass [kg/mol], 
                                 rho - Density array [kg/m^3], 
                                 sigma_0 - Electrical conductivity,
                                 N_m - N / m''',
                from_data=None,
                ):

        self.wl = wl
        self.T = T
        self.M = metal_props['M']
        self.N_m = metal_props['N_m']
        self.rho = self._interp_props(metal_props['T'], metal_props['rho'], self.T)
        self.sigma_0 = self._interp_props(metal_props['T'], metal_props['sigma_0'], self.T)

        self.N = self._get_free_e_density()
        self.m = self._get_eff_mass_e()
        self.omega = self._get_angular_freq()
        self.omega_p = self._get_plasma_freq()
        self.tau = self._get_relax_time()

        if from_data == None:
            self.eps_1 = self._get_real_dielectric_drude()
            self.eps_2 = self._get_imag_dielectric_drude()
            self.n_1 = self._get_n_1()
            self.k_1 = self._get_k_1()
        else:
            self.path = from_data
            self.n_1 = self._get_n_1_fromfile()[:, np.newaxis]
            self.k_1 = self._get_k_1_fromfile()[:, np.newaxis]


        self.R = self._get_R()
        self.A = self._get_A()







    def _interp_props(self, arg_0, func, arg_1):
        """
        Performs density interpolation
        |
        Выполняет интерполяцию плотности 
        """
        return interp1d(arg_0, func, kind='cubic')(arg_1)

    def _get_free_e_density(self):
        '''
        The free electron density
        |
        Плотность свободных электронов
        '''

        return (self.rho * self.Na) / self.M

    def _get_eff_mass_e(self):
        '''
        The optical effective mass of the electron
        |
        Оптически эффективная масса электрона
        '''
        return self.N / self.N_m

    def _get_plasma_freq(self):
        '''
        [omega_p]
        The plasma frequency of the material
        |
        Плазменная частота материала
        '''

        return ((self.N * self.e**2) / (self.m * self.eps_0))**0.5

    def _get_angular_freq(self):
        '''
        [omega]
        The angular frequency of the incident light
        | 
        Угловая частота падующего света
        '''

        return 2 * np.pi * (self.c / self.wl)

    def _get_relax_time(self):
        '''
        [tau] 
        The electron relaxation time
        |
        Время релаксации электрона
        '''

        return (self.m * self.sigma_0) / (self.N * self.e**2)

    def _get_real_dielectric_drude(self):
        '''
        Real part of the complex dielectric function
        |
        Вещественная часть комплексной диэлектрической проницаемости
        '''

        return np.array([1 - self.omega_p**2 / (self.omega[i]**2 + self.tau**-2) for i in range(len(self.wl))])

    def _get_imag_dielectric_drude(self):
        '''
        Imaginary part of the complex dielectric function
        |
        Мнимая часть комплексной диэлектрической проницаемости
        '''

        return np.array([(self.omega_p**2 * self.tau**-1) / (self.omega[i] * (self.omega[i] + self.tau**-2)) for i in range(len(self.wl))])

    def _get_n_1(self):
        '''
        The index of refraction
        |
        Показатель преломления
        '''

        return (0.5 * (self.eps_1 + (self.eps_1**2 + self.eps_2**2)**0.5))**0.5

    def _get_k_1(self):
        '''
        The extinction coefficient
        |
        Коэффициент экстинкции
        '''

        return (0.5 * (-1*self.eps_1 + (self.eps_1**2 + self.eps_2**2)**0.5))**0.5

    def _get_n_1_fromfile(self):
        '''
        The index of refraction
        '''
        # Import from a file
        data = pd.read_csv(self.path)

        return self._interp_props(data['wl'], data['n'], self.wl)

    def _get_k_1_fromfile(self):
        '''
        The index of refraction
        '''
        # Import from a file
        data = pd.read_csv(self.path)

        return self._interp_props(data['wl'], data['k'], self.wl)

    def _get_R(self):
        '''
        Reflectivity
        |
        Коэффициент отражения
        '''

        return ((self.n_1 - 1)**2 + self.k_1**2) / ((self.n_1 + 1)**2 + self.k_1**2)

    def _get_A(self):
        '''
        Absorptivity
        |
        Коэффициент пропускания
        '''

        return 1 - self.R

    def result_table(self):
         
        table = pd.DataFrame()
        for i in range(len(self.wl)):
            index = pd.MultiIndex.from_product([self.wl[i:i+1], self.T])
            df = pd.DataFrame({'Conductivity': self.sigma_0,
                                'Plasma frequency': self.omega_p,
                                'Relaxation time': self.tau,
                                'n': self.n_1[i,:],
                                'k': self.k_1[i,:],
                                'Reflectivity': self.R[i,:],
                                'Absorptivity': self.A[i,:]}, 
                                index=index, )
            table = table.append(df)


        index = pd.MultiIndex.from_tuples(table.index)
        table = table.reindex(index)
        table.index.names = ['Wavelength, m', 'Temperature, K']
        table = table[['Conductivity', 'Plasma frequency',
                    'Relaxation time','n', 'k', 
                    'Reflectivity', 'Absorptivity']]

        return table

class LorentzDrudeModel(DrudeModel):

    def __init__(self,
                 wl: 'Wavelength/s | Длины волн',
                 T: 'Temperature | Температура',
                 metal_props: '''M - Molar mass [kg/mol], 
                                 rho - Density array [kg/m^3], 
                                 sigma_0 - Electrical conductivity,
                                 N_m - N / m''',
                omega_0,
                from_data=None,
                ):
        super().__init__(wl, T, metal_props, from_data)
        self.wl = wl
        self.T = T
        self.M = metal_props['M']
        self.N_m = metal_props['N_m']
        self.omega_0 = omega_0

        self.eps_1 = self._get_real_dielectric()
        self.eps_2 = self._get_imag_dielectric()
        

        self.n_1 = self._get_n_1()
        self.k_1 = self._get_k_1()
        
    def _get_real_dielectric_lorentz(self):
        '''
        Real part of the complex dielectric function
        |
        Вещественная часть комплексной диэлектрической проницаемости
        '''

        return np.array([1+((self.e**2 * self.N) / (self.m * self.eps_0)) * ((self.omega_0**2 - self.omega[i]**2) / ((self.omega_0**2 - self.omega[i]**2)**2 + self.omega[i]**2 * self.tau**-2)) for i in range(len(self.wl))])

    def _get_imag_dielectric_lorentz(self):
        '''
        Imaginary part of the complex dielectric function
        |
        Мнимая часть комплексной диэлектрической проницаемости
        '''

        return np.array([((self.e**2 * self.N) / (self.m * self.eps_0)) * ((self.omega[i] * self.tau**-1)  / ((self.omega_0**2 - self.omega[i]**2)**2 + self.omega[i]**2 * self.tau**-2)) for i in range(len(self.wl))])

    def _get_real_dielectric(self):


        return self.eps_1 + self._get_real_dielectric_lorentz()

    def _get_imag_dielectric(self):

        return self.eps_2 + self._get_imag_dielectric_lorentz()

    def _get_n_1(self):
        '''
        The index of refraction
        |
        Показатель преломления
        '''
        return (0.5 * (self.eps_1 + (self.eps_1**2 + self.eps_2**2)**0.5))**0.5

    def _get_k_1(self):
        '''
        The extinction coefficient
        |
        Коэффициент экстинкции
        '''

        return (0.5 * (-1*self.eps_1 + (self.eps_1**2 + self.eps_2**2)**0.5))**0.5

    

def main():
    titan_bulk = DrudeModel(wl=np.arange(2e-6, 31e-6, 1e-6),
                            T=np.arange(300, 1100, 100),
                            metal_props=materials.Titanium_prop,
                            )

    omega_0=2 * np.pi * 1.878778621627e14
    titan_bulk_norm = DrudeModel(wl=np.arange(2e-6, 31e-6, 1e-6),
                                T=np.array([300]),
                                metal_props=materials.Titanium_prop,
                                )
    print(titan_bulk.omega)

    titan_bulk_Palik = DrudeModel(wl=np.arange(2e-6, 31e-6, 1e-6),
                                T=np.array([300]),
                                metal_props=materials.Titanium_prop,
                                from_data='Data/Ti_Rakic_LD.csv')
    
    ti_LD = LorentzDrudeModel(wl=np.arange(2e-6, 31e-6, 1e-6),
                              T=np.arange(300, 1100, 100),
                              omega_0=omega_0,
                              metal_props=materials.Titanium_prop,)

    print(ti_LD.eps_1)
    print(ti_LD.n_1)

    # print(titan_bulk.R)
    # table = titan_bulk.result_table()
    # print(table.loc[(2e-6, 300)])


    # table1 = titan_bulk_norm.result_table()
    # print(table1)

    # print(titan_bulk_Palik.n_1)

    # table2 = titan_bulk_Palik.result_table()
    # print(table2)



    

if __name__ == '__main__':
    main()
    
    


    

