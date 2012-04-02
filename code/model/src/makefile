#FFLAGS= -fast -fpp
##FFLAGS= -fast -fpp -prof-gen
##FFLAGS= -fpp -prof-use -ipo

#FFLAGS= -fpp 
FFLAGS= -fpp -real-size 64 
#FFLAGS= -cpp
CFLAGS= -c
#CFLAGS= -O3 -c -g
FC=ifort





#LIB_DIR           = /usr/local/netcdf-4.1.3/lib/libnetcdff.a
#INC_DIR           = -I/usr/local/netcdf-4.1.3/include
#LNK_FLAGS         = -lnetcdff
LIB_DIR           = -L/Users/jgilet/Library/netcdf/lib
INC_DIR           = -I/Users/jgilet/Library/netcdf/include
LNK_FLAGS         = -lnetcdf
#LNK_FLAGS         = 
#TRACMASS_FLAGS    = -Dzgrid3d -Druntracmass -Dfull_wflux
#TRACMASS_FLAGS    = -Dzgrid3d  -Dfull_wflux


%.o : %.f90
	$(FC) $(FFLAGS) $(CFLAGS) $(ARG_FLAGS) $(TRACMASS_FLAGS) $(LNK_FLAGS) $< -o $@ 
	#$(FC) $(FFLAGS) $(CFLAGS) $(ARG_FLAGS) $(TRACMASS_FLAGS) $< -o $@


NHO  = mymodules.o mod_particles.o mod_grids.o relaxation.o main.o ana_uvw.o init_tr.o  ini_rho.o inith_thermalwind.o ini_tracer.o mod_io.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rpflattopog.o correctbc.o surfaceflux.o conadjust_sT.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o writeslice.o writeksurf.o writeyslice.o writeisopycnal.o writestrain.o density.o geostroph.o ran3.o  mprove.o dgtsv.o pcorrect.o  diag_vort.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o srcface_nopy.o tracersource_nut.o windstress.o calc_dvdy.o streamfunction.o n2budget_topbot.o calcn2.o pvcalc.o analytic_eval.o sigma2z.o msg_error.o stprofile_X.o  mod_grids.o utils.o
#NHO  = mymodules.o mod_particles.o main.o mod_grids.o ana_uvw.o init_tr.o ini_rho.o inith_thermalwind.o ini_tracer.o sub_io.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_rho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rgrad_song_bndry0dsdy.o correctbc.o surfaceflux.o conadjust_sT.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o  density.o geostroph.o ran3.o  mprove.o dgtsv.o pcorrect.o  diag_vort.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o srcface_nopy.o tracersource_nut.o windstress.o calc_dvdy.o streamfunction.o n2budget_topbot.o calcn2.o pvcalc.o analytic_eval.o sigma2z.o
#NHO  = mymodules.o mod_particles.o main.o mod_grids.o ana_uvw.o init_tr.o ini_rho.o ini_tracer.o sub_io.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rgrad_song_bndry0dsdy.o correctbc.o surfaceflux.o conadjust_sT.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o  density.o geostroph.o ran3.o  mprove.o dgtsv.o pcorrect.o  inith_thermalwind.o diag_vort.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o srcface_nopy.o tracersource_nut.o windstress.o calc_dvdy.o streamfunction.o n2budget_topbot.o calcn2.o pvcalc.o analytic_eval.o 
nh: $(NHO)
	$(FC) $(FFLAGS)  -o nh $(NHO) $(LIB_DIR) $(LNK_FLAGS) 

$(NHO) : size.h ini_param.f90



# NHTRMASSO  = mymodules.o tracmass_modules.o tracmass_interface.o tracmass.o seed_particles.o time_subs.o main.o init_tr.o outcdf.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_rho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rgrad_song_bndry0dsdy.o correctbc.o surfaceflux.o conadjust_rho.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o writeslice.o density.o geostroph.o ran3.o tracerinit_dots.o mprove.o dgtsv.o pcorrect.o stprofile_lightintrusion.o topog.o inith_fixdh.o vort.o writeksurf.o writeyslice.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o srcface_nopy.o tracersource_isopyc.o windstress.o calc_dvdy.o streamfunction.o n2budget_topbot.o calcn2.o pvcalc.o writen2budget.o outpv.o writeisopycnal.o analytic_eval.o tracerrelease_dots.o

# NHCNO  =  maincn.o initcn_tr.o outcdf.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_rho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rgrad_song_bndry0dsdy.o correctbc.o surfaceflux.o conadjust_tr.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o writeframe.o writeslice.o density.o geostroph.o ran3.o tracerinit.o mprove.o dgtsv.o pcorrect.o topog.o inith_fixdh.o vort.o writeksurf.o writeyslice.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o readcdfcn_tr.o tracersource.o srcface_nopy.o windstress.o calc_dvdy.o streamfunction.o n2budget.o calcn2.o pvcalc.o writen2budget.o

# VORO  =  vort_writecdf.o initcn_tr.o seval.o sigma_toplayer.o staticsigma.o findzall_toplayer_moves.o density.o readcdfcn_vor.o evalrho_tr.o smooth.o ran3.o meanh.o potdens.o

# VORXTO  =  vorticity_xt.o initcn_tr.o seval.o sigma.o findzall.o density.o readcdfcn_vor.o evalrho_tr.o smooth.o ran3.o meanh.o potdens.o

# VELXTO  =  vel_xt.o initcn_tr.o seval.o sigma.o findzall.o density.o readcdfcn_vor.o evalrho_tr.o smooth.o ran3.o meanh.o potdens.o


nhcn: $(NHCNO)
	$(FC) $(FFLAGS)  -o nhcn $(NHCNO) -L/usr/local/lib -lnetcdf

writevor: $(VORO)
	gfortran $(FFLAGS)  -o writevor $(VORO) -L/usr/local/lib -lnetcdf 

vorxt: $(VORXTO)
	gfortran $(FFLAGS)  -o vorxt $(VORXTO) -L/usr/local/lib -lnetcdf

velxt: $(VELXTO)
	gfortran $(FFLAGS)  -o velxt $(VELXTO) -L/usr/local/lib -lnetcdf

clean :
	rm *.o *.mod
