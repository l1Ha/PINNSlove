MODULE Propag_Module
   USE Numeric_Kinds_Module
   USE Regins_Module
   USE DBS_Module
   USE VIVS_Module
   SAVE
   NAMELIST/propag/ SwitchToVIVS, nsteps
   ! nsteps: Number of steps per asymptotic wavelength used in the Log-derivative method.
   ! SwitchToVIVS:  Distance to switch from the Log-derivative propagator to the VIVAS propagator.
END MODULE Propag_Module
