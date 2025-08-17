Program Description: 

Computes the normalized intensity SDF for a plane wave that traverses a phase screen

I(mu) = int_{-inf}^{inf} exp[-gamma(eta,mu)] exp(-i eta mu) d eta

The normalized screen SDF is specified as a piecewise power law

           {   muo^(-p1),          if  0<= mu <= muo,
P(mu) = U1 {    mu^(-p1),          if muo < mu <= mub,
           { mub^(p2-p1) mu^(-p2), if mub < mu <= mui,

where U1 = Cp rhof^(p1-1) and 

    Cp   = phase spectral strength
    rhof = sqrt(z/wavk) is the Fresnel scale
    z    = distance past the screen
    wavk = signal wavenumber
    k    = transverse wavenumber
    mu   = rhof*k  normalized transverse wavenumber
    muo  = rhof*ko normalized outer scale wavenumber
    mub  = rhof*kb normalized break scale wavenumber
    mui  = rhof*ki normalized inner scale wavenumber

We assume that 0 < muo < mub < mui and also that muo << muf << mui, where
muf = 2*pi is the normalized wavenumber corresponding to the Fresnel scale.

The structure interaction function gamma(eta,mu) is given by

gamma(eta,mu) = 
   16 U1 int_{  0,muo} muo^(-p1)   sin^2(chi eta/2) sin^2(chi mu/2) d chi/(2 pi)
 + 16 U1 int_{muo,mub}             sin^2(chi eta/2) sin^2(chi mu/2) d chi/(2 pi)
 + 16 U1 int_{muo,mui} mub^(p2-p1) sin^2(chi eta/2) sin^2(chi mu/2) d chi/(2 pi)

Universal scattering strength U equals U1 if mub>=1 and U1 mub^(p2-p1) otherwise

This program fully supports the limiting cases muo->0 and/or mui->infinity

--------------------------------------------------------------------------------

Program Usage: ispectrum U p1 p2 mub muo mui [mu] | ([mu_min] [mu_max] [mu_num])

Calling options: 
 1: ispectrum U p1 p2 mub muo mui
 2: ispectrum U p1 p2 mub muo mui mu
 3: ispectrum U p1 p2 mub muo mui mu_min mu_max mu_num

Notes:
 * Setting muo and/or mui to zero omits them from the model
 * Option1 uses adaptive quadrature to compute S4
 * Option2 computes I(mu) at the specified value of mu
 * Option3 computes I(mu) at log spaced mu values between mu_min and mu_max
 * The intensity SDF I(mu) is written to ispectrum.dat
 * Parameters and moments are written to ispectrum.log

Written by Charles Carrano, Boston College (charles.carrano@bc.edu

For additional details, see

   Carrano, C. S., and C. L. Rino (2016), A theory of scintillation for two-component
      power law irregularity spectra: Overview and numerical results, Radio Sci., 
      51, 789–813, doi:10.1002/2015RS005903.

