# Compute reflection and transmission coefficients for arbitrary lossy dielectric filmstack using transfer matrix formalism
Compute reflection and transmission coefficients for arbitrary lossy dielectric filmstack using transfer matrix formalism

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
