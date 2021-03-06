function IRI_IO
%
%   
%   -TIME: Times to compute IRI model either in MATLAB serial date number
%   format or a string that can be converted into MATLAB serial date number
%   format using DATENUM with no format specified (see documentation of
%   DATENUM for more information). Whether the times are local or UTC are
%   determined by the input UTC. Valid range is from year 1958 to year 2016
%   currently (optional, default is January 1, 2000 at 01:30).
%   
%   -LATITUDE: Latitude in degrees to compute IRI model. Whether this is
%   geodetic, geocentric, or geomagnetic latitude is determined by the
%   input COORD. Valid range is -90 degrees to 90 degrees (optional,
%   default is 50 degrees).
%   
%   -LONGITUDE: Longitude in degrees to compute IRI model. Whether this is
%   geodetic, geocentric, or geomagnetic longitude is determined by the
%   input COORD. Valid range is 0 degrees to 360 degrees (optional, default
%   is 40 degrees).
%   
%   -HEIGHT: For geodetic or geomagnetic coordiates, the height in km above
%   the Earth's surface. For geocentric coordiates, the radius in km from
%   the center of the Earth. Valid range for altitude is 60 km to 2000 km,
%   although the recommended upper limit is 1500 km (optional, default when
%   all other inputs are scalars is to sweep from 100:50:2000 km and when
%   any other input is an array is 100 km).
%   
%   -UTC: Set to true (or 'UTC', 'U') if the times in TIME are in
%   Coordinated Universal Time (UTC) or false (or 'Local', 'LT', 'L') if
%   the times in TIME are in local time (optional, default is true).
%   
%   -COORD: String specifying the coordinate system to use. Can be
%   geodetic ('geodetic', 'geod', or 'gd'), geomagnetic ('geomagnetic',
%   'geom', or 'gm'), or geocentric ('geocentric', 'geom', or 'gm')
%   (optional, default is geodetic).
%   
%   -CURLDIR: Directory where curl.exe can be found (optional, only
%   necessary for Windows computers, and default for those is the same
%   directory that this function is located).
%   
%   -Rz12: 13-month running mean of sunspot number index. This index is
%   typically calculated automatically from the time you specify, but you
%   can enter your own index here and not rely on the internal index file.
%   Valid range is 0 to 400.
%   
%   -IG12: Index based on foF2 measurements from a dozen ionosondes
%   correlated with the CCIR foF2 maps [1].  This index is typically
%   calculated automatically from the time you specify, but you can enter
%   your own index here and not rely on the internal index file. Valid
%   range is -50 to 400.
% 
%   -F10_7_DAILY: F10.7 radio flux, a daily index based on full disc solar
%   flux measurements at 2800 MHz (10.7 cm wavelength) first made at the
%   Algonquin Radio Observatory, near Ottawa, Canada (1947 until May 31,
%   1991) and then from the Dominion Radio Astrophysical Observatory, near
%   Penticton, British Columbia at local noon (1700 UT at Ottawa and 2000
%   UT at Penticton). The flux values are expressed in solar flux units 
%   (1 s.f.u. = 10-22 W*m-2*Hz-1 ). Indices are tabulated in two forms: the
%   "observed flux" (S), and the "adjusted flux" (Sa). The former are the
%   actual measured values, and are affected by the changing distance
%   between the Earth and Sun throughout the year, whereas the latter are
%   scaled to a standard distance of 1 AU. The "observed flux" values are
%   used in IRI. This index is typically calculated for you, but you can
%   enter your own and not rely on the internal index file. Valid range is
%   0 to 400.
%   
%   -F10_7_81DAY: The 81-day running mean of the daily F10.7 index
%   described above. 81 days are about 3 solar rotations and statistical
%   studies have found good correlation between this parameter or
%   PF10.7=(F10.7_daily + F10.7_81day)/2 and ionospheric parameters. This
%   index is typically calculated for you, but you can enter your own and
%   not rely on the internal index file. Valid range is 0 to 400.
%   
%   -TEC_HMAX: Electron content upper boundary in km. The electron content
%   is the vertical electron content (vTEC) calculated using numerical
%   integration from a lower to an upper boundary. You must enter a value
%   for the upper boundary if you want to obtain electron content values.
%   The value should not be higher than 2000 km because this is the upper
%   limit of the IRI validity range for electron density profiles. The
%   lower boundary is set to the starting height of the IRI validity range
%   (60 km during the day and 80 km during the night).
% 
%   -Ne_TOP: Ne Topside. There are three options here:
%     1. IRI-2001: The IRI-2001 model that was based primarily on Alouette
%     1 topside sounder data with some AE-C, AEROS, and DE-2 in situ data
%     [2]. Selected by inputting one of 3, 'IRI01', 'IRI2001', or '2001'.
%     2. IRI01-corr: A correction of the 2001 model with the help of
%     Aloutte 2, ISIS 1 and 2 topside sounder data [3]. Selected by
%     inputting one of 2, 'IRI01-corr', 'IRI2001-corr', 'corr', or 'c'.
%     3. NeQuick: The model developed by Radicella and his group at ICTP
%     using Intercosmos 19 topside sounder data in addition to the ISIS 1
%     and 2 data [4]. Selected by inputing one of 1, 'NeQuick', 'Quick',
%     'Qu', 'Q', or 'Ne'.
%   The default option is NeQuick.
% 
%   -FPEAK: F peak model. There are two options here:
%     1. CCIR: The CCIR model for the F peak plasma frequency foF2 was
%     developed by Jones and Gallet [5] using data from the worldwide
%     network of ionosondes. It is the model recommended by the Comite
%     Consultatif International des Radiocommunications (CCIR) of the
%     International Telecommunications Union [6]. Selected by inputting one
%     of 2, 'CCIR', or 'C'.
%     2. URSI: The URSI model was developed by Rush et al. [7] using a
%     physical model to obtain screen points over regions not covered by
%     ionosondes instead of the extrapolation along magnetic field lines
%     employed by Jones and Gallet [5]. Selected by inputting one of 1,
%     'URSI', or 'U'.
%   The default option is URSI. The CCIR model is recommended over the
%   continents and the URSI model over the oceans. NOTE: Changes here will
%   also affect the height of the F peak hmF2 because the hmF2 model
%   depends on the ratio foF2/foE where foE is the plasma frequency at the
%   E peak.
%   
%   -F2STORM: foF2 Storm model on (true or 'on') or off (false or 'off')
%   (optional, default is on). The F peak storm model describes the average
%   storm behavior in terms of the ratio foF2_storm/foF2_quiet based on the
%   ap history over the preceding 33 hours [8]. A large volume of ionosonde
%   data for storms during the 1980-1990 time period were used to describe
%   the most coherent and repeatable features of the ionospheric storm
%   response.
% 
%   -BOTTOM: Bottomside thickness. The bottomside thickness is the height
%   difference between the F peak height hmF2 and the height where the
%   electron density profile has dropped down to half the F peak value
%   (NmF2/2). There are three options here:
%     1. Bilitza-2000: This (table-)option is based on incoherent scatter
%     data [9]. Selected by inputting one of 1, 'B0 Table', 'B0Table',
%     'B0', or 'B'.
%     2. Gulyaeva: This option is based on ionosonde data mostly from mid-
%     latitudes [10]. Selected by inputting one of 2, 'Gulyaeva', 'Guly',
%     or 'G'.
%     3. ABT-2009: Altadill et al. [11] used a large volume of ionosonde
%     data to develop a much improved representation of latitudinal and
%     solar cycle variation of B0 and B1. B1 is a parameter describing the
%     bottomside profile shape. Selected by inputting one of 3, 'ABT-2009',
%     or 'ABT'.
%   The default option is ABT-2009.
% 
%   -F1_PROB: This parameter describes the occurrence probability of an F1
%   layer. There are three options here:
%     1. IRI-95: This option uses the ITU-recommended Ducharme et al. model
%     [12] that applies a simple cutoff solar zenith angle, so probability
%     is either 0 or 1. Selected by inputting one of 3, 'IRI-95', 'IRI95',
%     'IRI', or '95'.
%     2. Scotto-1997 no L: This option uses the model developed by Scotto
%     et al. [13] using only ionograms with clear F1 layer presence and
%     excluding the more uncertain L condition cases. Ionograms often
%     exhibit a Fl ledge rather than a fully developed cusp, primarily
%     during the time period just before the Fl layer disappears. These
%     cases are described as L condition according to the URSI standard
%     nomenclature. Selected by inputting one of 1, 'Scotto-1997 no L', 
%     'no L', or 'nL'.
%     3. Scotto-1997 with L: Here L-condition cases were included. Selected
%     by inputting one of 2, 'Scotto-1997 with L', 'with L', or 'wL'.
%   The default option is Scotto-1997 no L.
%   
%   -AURORAL_BOUNDARY: Turn auroral boundary on (true or 'on') or off
%   (false or 'off') (optional, default is off). Leaving this off will
%   produce faster results.
%   
%   -foE_STORM: Turn the E peak auroral storm model on (true or 'on') or
%   off (false or 'off') (optional, default is on). The model was developed
%   by Mertens et al. [14] and Fernandez et al. [15] and describes the
%   average storm behavior in terms of the ratio foE_storm/foE_quite based
%   on the ap index history.
% 
%   -Ne_Dreg: D-region electron density. There are two options here:
%     1. FT-2001: Model based on Friedrich's compilation of reliable rocket
%     data [16]. Selected by inputting one of 2, 'FPT-2000', 'FPT', or
%     '2000'.
%     2. IRI-95: Model based on a much smaller selection of typical rocket
%     profiles [17]. Selected by inputting one of 1, 'IRI-95', 'IRI95',
%     'IRI', or '95'.
%   The default option is IRI-95.
% 
%   -Te_TOP: Topside electron temperature. There are two options here:
%     1. IRI-95 is using the global models at fixed altitudes developed by
%     Brace and Theis [18] based on their ISIS and AE-C electron
%     temperature measurements. The model is described in [19] and also in
%     the IRI-90 report that is available as PDF document from the
%     references section of the IRI homepage. Selected by inputting one of
%     2, 'IRI-95', 'IRI95', 'IRI', or '95'.
%     2. TBT-2012 model is newer and uses a large volume of satellite in
%     situ measurements and includes variations with solar activity [20].
%     Selected by inputting one of 1, 'TBT-2012', 'TBT', or '2012'.
%   The default option is TBT-2012.
% 
%   -IONCOMP: Ion composition. There are two options here:
%     1. DS95/DY85: Uses the Danilov and Smirnova model [21] based on their
%     compilation of rocket data in the region below the F-peak and the
%     Danilov and Yaichnikov model [22] based on their compilation of
%     Russian high-altitude rocket measurements. Selected by inputting one
%     of 2, 'DS95/DY85', or 'DS95'.
%     2. RBV10/TTS03: Based on an adjustment of the ion composition from
%     the FLIP-model photochemistry to the IRI electron density profile in
%     region below the F peak [23]. In the region above the F-peak the
%     model of Triskova et al. [24] based on satellite ion mass
%     spectrometer measurements. Selected by inputting one of 1,
%     'RBV10/TTS03', or 'RBV10'.
%   The default option is RBV10/TTS03.
%   
%   -NmF2_foF2: A measured value to update the IRI profile to actual
%   conditions. Only valid when varying HEIGHT linearly. If this number is
%   between 10^3 and 10^8, it represents the F2 peak density (NmF2) in
%   cm^-3. If instead it is between 2 and 14, it represents the F2 plasma
%   frequency foF2 in MHz.
%   
%   -hmF2_M3000F2: A measured value to update the IRI profile to actual
%   conditions. Only valid when varying HEIGHT linearly. If this number is
%   between 100 and 1000, it represents the F2 peak height (hmF2) in km. If
%   instead it is between 1.5 and 4, it represents the propagation factor
%   M(3000)F2.
%   
%   -NmE_foE: A measured value to update the IRI profile to actual
%   conditions. Only valid when varying HEIGHT linearly. If this number is
%   between 20 and 10^8, it represents the E peak density (NmE) in cm^-3.
%   If instead it is between 0.1 and 14, it represents the E plasma
%   frequency in MHz.
% 
%   -hmE: Measured E peak height in km to update the IRI profile to actual
%   conditions. Valid range is between 70 km and 200 km.
% 
% Outputs:
%   
%   -OUT: Array with the same number of rows as elements in the position
%   and time inputs and the following data in each column:
%     1. Electron density (Ne) in m^-3.
%     2. Ratio of Ne to the F2 peak density (NmF2).
%     3. Neutral temperature (Tn) in K.
%     4. Ion temperature (Ti) in K.
%     5. Electron temperature (Te) in K.
%     6. Atomic oxygen ions (O+) percentage.
%     7. Atomic hydrogen ions (H+) percentage.
%     8. Atomic helium ions (He+) percentage.
%     9. Molecular oxygen ions (02+) percentage.
%     10. Nitric oxide ions (NO+) percentage.
%     11. Cluster ions percentage.
%     12. Atomic nitrogen ions (N+) percentage.
%     13. Total electron content (TEC) in 10^16 m^-2.
%     14. TEC top percentage.
%     15. Height of the F2 peak (hmF2) in km.
%     16. Height of the F1 peak (hmF1) in km.
%     17. Height of the E peak (hmE) in km.
%     18. Height of the D peak (hmD) in km.
%     19. Density of the F2 peak (NmF2) in m^-3.
%     20. Density of the F1 peak (NmF1) in m^-3.
%     21. Density of the E peak (NmE) in m^-3.
%     22. Density of the D peak (NmD) in m^-3.
%     23. Propagation factor M(3000)F2.
%     24. Bottomside thickness (B0) in km.
%     25. Bottomside shape (B1).
%     26. E-valley width in km.
%     27. E-valley depth (Nmin/NmE) in km.
%     28. F2 plasma frequency (foF2) in MHz.
%     29. F1 plasma frequency (foF1) in MHz.
%     30. E plasma frequency (foE) in MHz.
%     31. D plasma frequency (foD) in MHz.
%     32. Equatorial vertical ion drift in m/s.
%     33. Ratio of foF2 storm to foF2 quiet.
%     34. F1 probability.
%     35. CGM latitude of auroral oval boundary.
%     36. Spread F probability.
%     37. Ratio of foE storm to foE quiet.
%     38. 12-month running mean of sunspot number Rz12 used by the model.
%     39. Ionospheric index IG12 used by the model.
%     40. Daily solar radio flux F107D used by the model.
%     41. 81 day solar radio flux F107_81D used by the model.
%     42. 3 hour ap index used by the model.
%     43. Daily ap index used by the model.
%     44. 3 hour Kp index used by the model.
%   A value of -1 for any output indicates that the parameter is not
%   available for the specified range. TEC = -1 means you have not entered
%   an upper boundary height for TEC_HMAX.'


