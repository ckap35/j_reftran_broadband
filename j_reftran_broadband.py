import numpy as np
import gc
from numba import jit
from numba import njit, float32, int32,complex64
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

#@jit
def sind(x):
    y = np.sin(np.radians(x))
    I = x/180.
    if type(I) is not np.ndarray:   #Returns zero for elements where 'X/180' is an integer
        if(I == np.trunc(I) and np.isfinite(I)):
            return 0
    return y


'''
https://ww2.mathworks.cn/matlabcentral/fileexchange/50923-jreftran-a-layered-thin-film-transmission-and-reflection-coefficient-calculator?requestedDomain=zh
    l = free space wavelength, nm
    d = layer thickness vector, nm
    n = layer complex refractive index vector
    t0= angle of incidence
    polarization should be 0 for TE (s-polarized), otherwise TM (p-polarized)
Example: Finding the coefficients for a 200nm gold layer surrounded by air, using the Johnson and Christy data
    input:
        [r,t,R,T,A]=jreftran_rt(500,[NaN,200,NaN],[1,0.9707 + 1.8562i,1],0,0)
    output:
        r = -0.4622 - 0.5066i               reflection coefficient
        t = -0.0047 + 0.0097i               transmission coefficient
        R = 0.4702                          reflectivity (fraction of intensity that is reflected)
        T = 1.1593e-04                      transmissivity (fraction of intensity that is transmitted)
        A = 0.5296                          absorptivity (fraction of intensity that is absorbed in the film stack)
'''

def jreftran_rt(wavelength, d, n, t0, polarization):

    M = np.zeros((2, 2, d.shape[0]), dtype=complex)
    M_t = np.identity(2, dtype=complex)
    Z0 = 376.730313     #impedance of free space, Ohms
    Y = n / Z0
    g = 1j * 2 * np.pi * n / wavelength        #propagation constant in terms of free space wavelength and refractive index
    t = (n[0] / n * sind(t0))
    t2 = t*t        #python All arithmetic operates elementwise,Array multiplication is not matrix multiplication!!!
    ct = np.sqrt(1 -t2);        # ct=sqrt(1-(n(1)./n*sin(t0)).^2); %cosine theta
    if polarization == 0:       # tilted admittance
        eta = Y * ct;       # tilted admittance, TE case
    else:
        eta = Y / ct;      # tilted admittance, TM case
    delta = 1j * g * d * ct
    ld = d.shape[0]
    #M = np.zeros((2, 2, ld),dtype=complex)
    for j in range(ld):
        a = delta[j]
        M[0, 0, j] = np.cos(a);
        M[0, 1, j] = 1j / eta[j] * np.sin(a);
        M[1, 0, j] = 1j * eta[j] * np.sin(a);
        M[1, 1, j] = np.cos(a)
        # ("M(:,:,{})={}\n\n".format(j,M[:,:,j]))
    #M_t = np.identity(2,dtype=complex)        #total charateristic matrix
    for j in range(1,ld - 1):
        M_t = np.matmul(M_t,M[:,:, j])

    e_1, e_2 = eta[0], eta[-1]
    De = M_t[0, 0] + M_t[0, 1] * eta[-1]
    Nu = M_t[1, 0] + M_t[1, 1] * eta[-1]
    if False:  # Add by Dr.zhu
        Y_tot = (M_t(2, 1) + M_t(2, 2) * eta(len(d))) / (M_t(1, 1) + M_t(1, 2) * eta(len(d)))
        eta_one = eta[0]
        Re = Y_tot.real
        Im = Y_tot.imag
        fx = 2 * Im * eta[0]

    e_de_nu = e_1 * De + Nu
    r = (e_1 * De - Nu) / e_de_nu;
    t = 2 * e_1 / e_de_nu;

    R = abs(r) * abs(r);
    T = (e_2.real / e_1) * abs(t) * abs(t)
    T = T.real
    a = De * np.conj(Nu) - e_2
    A = (4 * e_1 * a.real) / (abs(e_de_nu)*abs(e_de_nu));
    A = A.real

    return r, t, R, T, A


#plot wavelength dependence
max_pairs=10
all_T=np.array([])
all_R=np.array([])
all_A=np.array([])
all_z=np.array([])
d = np.array([np.NaN])     #%d = layer thickness vector, nm
n = np.array([1])


for pairs in range(0,max_pairs,1):
    n=np.append(n,np.array([1,0.9707/1+1.8562j/1, 1.0, 0.9707+1.8562j/1]))
    n_end=np.append(n,np.array([1]))
    d=np.append(d,np.array([60*40,10,60*40,10]))
    d_end=np.append(d,np.array([np.NaN]))



wavelengths=np.logspace(2,4,num=1000)
for wavelength in wavelengths:
    r, t, R, T, A = jreftran_rt(wavelength, d_end, n_end, 0, 0)
    all_R=np.append(all_R,R)
    all_T=np.append(all_T,T)
    all_A=np.append(all_A,A)

plt.plot(wavelengths,all_R,"-b", label="R")
plt.plot(wavelengths,all_T,"-k", label="T")
plt.plot(wavelengths,all_A,"-r", label="A")
plt.gca().legend()
plt.ylabel('Reflectance')
plt.xlabel('wavelength')
plt.show()
