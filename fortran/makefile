# $Id$
# Makefile for rrtmg - LINUX platform

SHELL  = /bin/sh
KGSRC = dat
################################
# Set variables
################################

PLATFORM = linux_gfortran
VERSION = v5.00

##################################
# Options for the FORTRAN compiler
##################################

FC = gfortran
   FCOPTS = -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -g -fbacktrace -ffree-line-length-none -std=legacy
# Compile absorption coefficient data without optimization
   FCOPTK = -O0 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -ffree-line-length-none -std=legacy

#################
# Executable name
#################

# OUTPUT = rrtmg_lw_$(VERSION)_$(PLATFORM)
OUTPUT = cfram_rrtmg

###########################################
# Source, module and object names and paths
###########################################

MPATH_LW = ./RRTMG_LW-master/modules
CPATH_LW = ./RRTMG_LW-master/src
CPATH_UTL_LW = ./RRTMG_LW-master/aer_rt_utils
MPATH_SW = ./RRTMG_SW-master/modules
CPATH_SW = ./RRTMG_SW-master/src
CPATH_UTL_SW = ./RRTMG_SW-master/aer_rt_utils
CDOTO = $(PLATFORM)_$(FC)_dbl.obj

MAIN    = cfram_rrtmg
MATH    = math
DRDT    = drdt
DRIVER  = rad_driver
DRIVERL = rad_driver_lw
FOUT    = writeout
MPRK_LW = parkind
MPAR_LW = parrrtm
MCLD_LW = rrlw_cld
MCON_LW = rrlw_con
MK01_LW = rrlw_kg01
MK02_LW = rrlw_kg02
MK03_LW = rrlw_kg03
MK04_LW = rrlw_kg04
MK05_LW = rrlw_kg05
MK06_LW = rrlw_kg06
MK07_LW = rrlw_kg07
MK08_LW = rrlw_kg08
MK09_LW = rrlw_kg09
MK10_LW = rrlw_kg10
MK11_LW = rrlw_kg11
MK12_LW = rrlw_kg12
MK13_LW = rrlw_kg13
MK14_LW = rrlw_kg14
MK15_LW = rrlw_kg15
MK16_LW = rrlw_kg16
MREF_LW = rrlw_ref
MTBL_LW = rrlw_tbl
MVSN_LW = rrlw_vsn
MWVN_LW = rrlw_wvn
CRTM_LW = rrtmg_lw_cfram
CCLD_LW = rrtmg_lw_cldprop
CCMC_LW = rrtmg_lw_cldprmc
CINI_LW = rrtmg_lw_init
CKGS_LW = rrtmg_lw_k_g
CRTR_LW = rrtmg_lw_rtrn
CRMR_LW = rrtmg_lw_rtrnmr
CRMC_LW = rrtmg_lw_rtrnmc
CSET_LW = rrtmg_lw_setcoef
CTAU_LW = rrtmg_lw_taumol
CMRN_LW = mcica_random_numbers
CMSC_LW = mcica_subcol_gen_lw.1col
CATM_LW = rrtatm
CEXT_LW = extra
CUTL_LW = util_gfortran

MPAR_SW = parrrsw
MAER_SW = rrsw_aer
MCLD_SW = rrsw_cld
MCON_SW = rrsw_con
MK16_SW = rrsw_kg16
MK17_SW = rrsw_kg17
MK18_SW = rrsw_kg18
MK19_SW = rrsw_kg19
MK20_SW = rrsw_kg20
MK21_SW = rrsw_kg21
MK22_SW = rrsw_kg22
MK23_SW = rrsw_kg23
MK24_SW = rrsw_kg24
MK25_SW = rrsw_kg25
MK26_SW = rrsw_kg26
MK27_SW = rrsw_kg27
MK28_SW = rrsw_kg28
MK29_SW = rrsw_kg29
MREF_SW = rrsw_ref
MTBL_SW = rrsw_tbl
MVSN_SW = rrsw_vsn
MWVN_SW = rrsw_wvn
CRTM_SW = rrtmg_sw_cfram
CCLD_SW = rrtmg_sw_cldprop
CCMC_SW = rrtmg_sw_cldprmc
CINI_SW = rrtmg_sw_init
CKGS_SW = rrtmg_sw_k_g
CRFT_SW = rrtmg_sw_reftra
CSET_SW = rrtmg_sw_setcoef
CSVT_SW = rrtmg_sw_spcvrt
CSMC_SW = rrtmg_sw_spcvmc
CTAU_SW = rrtmg_sw_taumol
CVRT_SW = rrtmg_sw_vrtqdr
CMSC_SW = mcica_subcol_gen_sw.1col

OBJFIL = $(MAIN).o \
		 $(MATH).o \
		 $(DRDT).o   \
		 $(DRIVER).o \
		 $(DRIVERL).o\
		 $(FOUT).o   \
         $(CCLD_LW).o \
         $(CCMC_LW).o \
         $(CKGS_LW).o \
         $(CRTR_LW).o \
         $(CRMR_LW).o \
         $(CRMC_LW).o \
         $(CSET_LW).o \
         $(CTAU_LW).o \
         $(CMRN_LW).o \
         $(CMSC_LW).o \
         $(CINI_LW).o \
         $(CRTM_LW).o \
         $(CATM_LW).o \
         $(CEXT_LW).o \
         $(CUTL_LW).o \
		 $(CCLD_SW).o \
		 $(CCMC_SW).o \
		 $(CKGS_SW).o \
		 $(CRFT_SW).o \
		 $(CSET_SW).o \
		 $(CTAU_SW).o \
		 $(CVRT_SW).o \
		 $(CSVT_SW).o \
		 $(CSMC_SW).o \
		 $(CMSC_SW).o \
		 $(CINI_SW).o \
		 $(CRTM_SW).o 

OBJPTH = $(CDOTO)/$(CMRN_LW).o \
         $(CDOTO)/$(CMSC_LW).o \
         $(CDOTO)/$(CCLD_LW).o \
         $(CDOTO)/$(CCMC_LW).o \
         $(CDOTO)/$(CKGS_LW).o \
         $(CDOTO)/$(CRTR_LW).o \
         $(CDOTO)/$(CRMR_LW).o \
         $(CDOTO)/$(CRMC_LW).o \
         $(CDOTO)/$(CSET_LW).o \
         $(CDOTO)/$(CTAU_LW).o \
         $(CDOTO)/$(CINI_LW).o \
         $(CDOTO)/$(CRTM_LW).o \
         $(CDOTO)/$(CATM_LW).o \
         $(CDOTO)/$(CEXT_LW).o \
         $(CDOTO)/$(CUTL_LW).o \
		 $(CDOTO)/$(CMSC_SW).o \
         $(CDOTO)/$(CCLD_SW).o \
		 $(CDOTO)/$(CCMC_SW).o \
         $(CDOTO)/$(CKGS_SW).o \
         $(CDOTO)/$(CRFT_SW).o \
         $(CDOTO)/$(CSET_SW).o \
         $(CDOTO)/$(CTAU_SW).o \
         $(CDOTO)/$(CVRT_SW).o \
         $(CDOTO)/$(CSVT_SW).o \
         $(CDOTO)/$(CSMC_SW).o \
         $(CDOTO)/$(CINI_SW).o \
         $(CDOTO)/$(CRTM_SW).o \
		 $(CDOTO)/$(FOUT).o    \
		 $(CDOTO)/$(DRDT).o   \
		 $(CDOTO)/$(MATH).o\
		 $(CDOTO)/$(DRIVER).o \
		 $(CDOTO)/$(DRIVERL).o\
		 $(CDOTO)/$(MAIN).o


OBMPTH = $(CDOTO)/$(MPRK_LW).o \
         $(CDOTO)/$(MPAR_LW).o \
         $(CDOTO)/$(MCLD_LW).o \
         $(CDOTO)/$(MCON_LW).o \
         $(CDOTO)/$(MK01_LW).o \
         $(CDOTO)/$(MK02_LW).o \
         $(CDOTO)/$(MK03_LW).o \
         $(CDOTO)/$(MK04_LW).o \
         $(CDOTO)/$(MK05_LW).o \
         $(CDOTO)/$(MK06_LW).o \
         $(CDOTO)/$(MK07_LW).o \
         $(CDOTO)/$(MK08_LW).o \
         $(CDOTO)/$(MK09_LW).o \
         $(CDOTO)/$(MK10_LW).o \
         $(CDOTO)/$(MK11_LW).o \
         $(CDOTO)/$(MK12_LW).o \
         $(CDOTO)/$(MK13_LW).o \
         $(CDOTO)/$(MK14_LW).o \
         $(CDOTO)/$(MK15_LW).o \
         $(CDOTO)/$(MK16_LW).o \
         $(CDOTO)/$(MREF_LW).o \
         $(CDOTO)/$(MTBL_LW).o \
         $(CDOTO)/$(MVSN_LW).o \
         $(CDOTO)/$(MWVN_LW).o \
	     $(CDOTO)/$(MPAR_SW).o \
		 $(CDOTO)/$(MAER_SW).o \
         $(CDOTO)/$(MCLD_SW).o \
         $(CDOTO)/$(MCON_SW).o \
		 $(CDOTO)/$(MK16_SW).o \
         $(CDOTO)/$(MK17_SW).o \
         $(CDOTO)/$(MK18_SW).o \
         $(CDOTO)/$(MK19_SW).o \
         $(CDOTO)/$(MK20_SW).o \
         $(CDOTO)/$(MK21_SW).o \
         $(CDOTO)/$(MK22_SW).o \
         $(CDOTO)/$(MK23_SW).o \
         $(CDOTO)/$(MK24_SW).o \
         $(CDOTO)/$(MK25_SW).o \
         $(CDOTO)/$(MK26_SW).o \
		 $(CDOTO)/$(MK27_SW).o \
         $(CDOTO)/$(MK28_SW).o \
         $(CDOTO)/$(MK29_SW).o \
         $(CDOTO)/$(MREF_SW).o \
         $(CDOTO)/$(MTBL_SW).o \
         $(CDOTO)/$(MVSN_SW).o \
         $(CDOTO)/$(MWVN_SW).o 

#################
# Load line
#################

#$(OUTPUT) : \
#	intro makedir $(OBMPTH) $(OBJPTH) conclude
#	$(FC) $(FCOPTS) -I$(CDOTO) -o $(OUTPUT) $(OBMPTH) $(OBJPTH)

$(OUTPUT) : \
	makedir $(OBMPTH) $(OBJPTH)
	$(FC) $(FCOPTS) -I$(CDOTO) -o $(OUTPUT) $(OBMPTH) $(OBJPTH) /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3


#################
# Compile modules
#################

$(CDOTO)/$(MPRK_LW).o : $(MPATH_LW)/$(MPRK_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MPRK_LW).f90
	\mv $(MPRK_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MPAR_LW).o : $(MPATH_LW)/$(MPAR_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MPAR_LW).f90
	\mv $(MPAR_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MCLD_LW).o : $(MPATH_LW)/$(MCLD_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MCLD_LW).f90
	\mv $(MCLD_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MCON_LW).o : $(MPATH_LW)/$(MCON_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MCON_LW).f90
	\mv $(MCON_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK01_LW).o : $(MPATH_LW)/$(MK01_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK01_LW).f90
	\mv $(MK01_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK02_LW).o : $(MPATH_LW)/$(MK02_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK02_LW).f90
	\mv $(MK02_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK03_LW).o : $(MPATH_LW)/$(MK03_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK03_LW).f90
	\mv $(MK03_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK04_LW).o : $(MPATH_LW)/$(MK04_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK04_LW).f90
	\mv $(MK04_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK05_LW).o : $(MPATH_LW)/$(MK05_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK05_LW).f90
	\mv $(MK05_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK06_LW).o : $(MPATH_LW)/$(MK06_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK06_LW).f90
	\mv $(MK06_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK07_LW).o : $(MPATH_LW)/$(MK07_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK07_LW).f90
	\mv $(MK07_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK08_LW).o : $(MPATH_LW)/$(MK08_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK08_LW).f90
	\mv $(MK08_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK09_LW).o : $(MPATH_LW)/$(MK09_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK09_LW).f90
	\mv $(MK09_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK10_LW).o : $(MPATH_LW)/$(MK10_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK10_LW).f90
	\mv $(MK10_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK11_LW).o : $(MPATH_LW)/$(MK11_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK11_LW).f90
	\mv $(MK11_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK12_LW).o : $(MPATH_LW)/$(MK12_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK12_LW).f90
	\mv $(MK12_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK13_LW).o : $(MPATH_LW)/$(MK13_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK13_LW).f90
	\mv $(MK13_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK14_LW).o : $(MPATH_LW)/$(MK14_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK14_LW).f90
	\mv $(MK14_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK15_LW).o : $(MPATH_LW)/$(MK15_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK15_LW).f90
	\mv $(MK15_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MK16_LW).o : $(MPATH_LW)/$(MK16_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MK16_LW).f90
	\mv $(MK16_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MREF_LW).o : $(MPATH_LW)/$(MREF_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MREF_LW).f90
	\mv $(MREF_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MTBL_LW).o : $(MPATH_LW)/$(MTBL_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MTBL_LW).f90
	\mv $(MTBL_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MVSN_LW).o : $(MPATH_LW)/$(MVSN_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MVSN_LW).f90
	\mv $(MVSN_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MWVN_LW).o : $(MPATH_LW)/$(MWVN_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_LW)/$(MWVN_LW).f90
	\mv $(MWVN_LW).o *.mod $(CDOTO)

$(CDOTO)/$(MPAR_SW).o : $(MPATH_SW)/$(MPAR_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MPAR_SW).f90
	\mv $(MPAR_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MAER_SW).o : $(MPATH_SW)/$(MAER_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MAER_SW).f90
	\mv $(MAER_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MCLD_SW).o : $(MPATH_SW)/$(MCLD_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MCLD_SW).f90
	\mv $(MCLD_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MCON_SW).o : $(MPATH_SW)/$(MCON_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MCON_SW).f90
	\mv $(MCON_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK16_SW).o : $(MPATH_SW)/$(MK16_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK16_SW).f90
	\mv $(MK16_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK17_SW).o : $(MPATH_SW)/$(MK17_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK17_SW).f90
	\mv $(MK17_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK18_SW).o : $(MPATH_SW)/$(MK18_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK18_SW).f90
	\mv $(MK18_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK19_SW).o : $(MPATH_SW)/$(MK19_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK19_SW).f90
	\mv $(MK19_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK20_SW).o : $(MPATH_SW)/$(MK20_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK20_SW).f90
	\mv $(MK20_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK21_SW).o : $(MPATH_SW)/$(MK21_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK21_SW).f90
	\mv $(MK21_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK22_SW).o : $(MPATH_SW)/$(MK22_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK22_SW).f90
	\mv $(MK22_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK23_SW).o : $(MPATH_SW)/$(MK23_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK23_SW).f90
	\mv $(MK23_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MK24_SW).o : $(MPATH_SW)/$(MK24_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK24_SW).f90
	\mv $(MK24_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MK25_SW).o : $(MPATH_SW)/$(MK25_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK25_SW).f90
	\mv $(MK25_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MK26_SW).o : $(MPATH_SW)/$(MK26_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK26_SW).f90
	\mv $(MK26_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MK27_SW).o : $(MPATH_SW)/$(MK27_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK27_SW).f90
	\mv $(MK27_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MK28_SW).o : $(MPATH_SW)/$(MK28_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK28_SW).f90
	\mv $(MK28_SW).o *.mod $(CDOTO)
	 
$(CDOTO)/$(MK29_SW).o : $(MPATH_SW)/$(MK29_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MK29_SW).f90
	\mv $(MK29_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MREF_SW).o : $(MPATH_SW)/$(MREF_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MREF_SW).f90
	\mv $(MREF_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(MTBL_SW).o : $(MPATH_SW)/$(MTBL_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MTBL_SW).f90
	\mv $(MTBL_SW).o *.mod $(CDOTO)

$(CDOTO)/$(MVSN_SW).o : $(MPATH_SW)/$(MVSN_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MVSN_SW).f90
	\mv $(MVSN_SW).o *.mod $(CDOTO)
	 
$(CDOTO)/$(MWVN_SW).o : $(MPATH_SW)/$(MWVN_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MPATH_SW)/$(MWVN_SW).f90
	\mv $(MWVN_SW).o *.mod $(CDOTO)

#################
# Compile source
#################

$(CDOTO)/$(CMRN_LW).o : $(CPATH_LW)/$(CMRN_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CMRN_LW).f90
	\mv $(CMRN_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CMSC_LW).o : $(CPATH_LW)/$(CMSC_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CMSC_LW).f90
	\mv $(CMSC_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CCLD_LW).o : $(CPATH_LW)/$(CCLD_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CCLD_LW).f90
	\mv $(CCLD_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CCMC_LW).o : $(CPATH_LW)/$(CCMC_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CCMC_LW).f90
	\mv $(CCMC_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CKGS_LW).o : $(CPATH_LW)/$(CKGS_LW).f90
	$(FC) $(FCOPTK) -I$(CDOTO) -c $(CPATH_LW)/$(CKGS_LW).f90
	\mv $(CKGS_LW).o $(CDOTO)

$(CDOTO)/$(CRTR_LW).o : $(CPATH_LW)/$(CRTR_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CRTR_LW).f90
	\mv $(CRTR_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CRMR_LW).o : $(CPATH_LW)/$(CRMR_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CRMR_LW).f90
	\mv $(CRMR_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CRMC_LW).o : $(CPATH_LW)/$(CRMC_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CRMC_LW).f90
	\mv $(CRMC_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CSET_LW).o : $(CPATH_LW)/$(CSET_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CSET_LW).f90
	\mv $(CSET_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CTAU_LW).o : $(CPATH_LW)/$(CTAU_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CTAU_LW).f90
	\mv $(CTAU_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CINI_LW).o : $(CPATH_LW)/$(CINI_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CINI_LW).f90
	\mv $(CINI_LW).o *.mod $(CDOTO)

$(CDOTO)/$(CRTM_LW).o : $(CPATH_LW)/$(CRTM_LW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_LW)/$(CRTM_LW).f90
	\mv $(CRTM_LW).o $(CDOTO)

$(CDOTO)/$(CATM_LW).o : $(CPATH_LW)/$(CATM_LW).f
	$(FC) $(FCOPTS) -c $(CPATH_LW)/$(CATM_LW).f
	\mv $(CATM_LW).o $(CDOTO)

$(CDOTO)/$(CEXT_LW).o : $(CPATH_LW)/$(CEXT_LW).f
	$(FC) $(FCOPTS) -c $(CPATH_LW)/$(CEXT_LW).f
	\mv $(CEXT_LW).o $(CDOTO)

$(CDOTO)/$(CUTL_LW).o : $(CPATH_UTL_LW)/$(CUTL_LW).f90
	$(FC) $(FCOPTS) -c $(CPATH_UTL_LW)/$(CUTL_LW).f90
	\mv $(CUTL_LW).o $(CDOTO)

$(CDOTO)/$(CMSC_SW).o : $(CPATH_SW)/$(CMSC_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CMSC_SW).f90
	\mv $(CMSC_SW).o *.mod $(CDOTO)

$(CDOTO)/$(CCLD_SW).o : $(CPATH_SW)/$(CCLD_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CCLD_SW).f90
	\mv $(CCLD_SW).o *.mod $(CDOTO)
 
$(CDOTO)/$(CCMC_SW).o : $(CPATH_SW)/$(CCMC_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CCMC_SW).f90
	\mv $(CCMC_SW).o *.mod $(CDOTO)

$(CDOTO)/$(CKGS_SW).o : $(CPATH_SW)/$(CKGS_SW).f90
	$(FC) $(FCOPTK) -I$(CDOTO) -c $(CPATH_SW)/$(CKGS_SW).f90
	\mv $(CKGS_SW).o $(CDOTO)
	 
$(CDOTO)/$(CRFT_SW).o : $(CPATH_SW)/$(CRFT_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CRFT_SW).f90
	\mv $(CRFT_SW).o *.mod $(CDOTO)
	 
$(CDOTO)/$(CSET_SW).o : $(CPATH_SW)/$(CSET_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CSET_SW).f90
	\mv $(CSET_SW).o *.mod $(CDOTO)
	
$(CDOTO)/$(CTAU_SW).o : $(CPATH_SW)/$(CTAU_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CTAU_SW).f90
	\mv $(CTAU_SW).o *.mod $(CDOTO)
	 
$(CDOTO)/$(CVRT_SW).o : $(CPATH_SW)/$(CVRT_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CVRT_SW).f90
	\mv $(CVRT_SW).o *.mod $(CDOTO)

$(CDOTO)/$(CSVT_SW).o : $(CPATH_SW)/$(CSVT_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CSVT_SW).f90
	\mv $(CSVT_SW).o *.mod $(CDOTO)
	
$(CDOTO)/$(CSMC_SW).o : $(CPATH_SW)/$(CSMC_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CSMC_SW).f90
	\mv $(CSMC_SW).o *.mod $(CDOTO)
	 
$(CDOTO)/$(CINI_SW).o : $(CPATH_SW)/$(CINI_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CINI_SW).f90
	\mv $(CINI_SW).o *.mod $(CDOTO)
	 
$(CDOTO)/$(CRTM_SW).o : $(CPATH_SW)/$(CRTM_SW).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(CPATH_SW)/$(CRTM_SW).f90
	\mv $(CRTM_SW).o $(CDOTO)

$(CDOTO)/$(MAIN).o : $(MAIN).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MAIN).f90
	\mv $(MAIN).o $(CDOTO)

$(CDOTO)/$(DRIVER).o : $(DRIVER).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(DRIVER).f90
	\mv $(DRIVER).o  $(CDOTO)

$(CDOTO)/$(DRIVERL).o : $(DRIVERL).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(DRIVERL).f90
	\mv $(DRIVERL).o  $(CDOTO)

$(CDOTO)/$(FOUT).o : $(FOUT).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(FOUT).f90
	\mv $(FOUT).o *.mod $(CDOTO)

$(CDOTO)/$(DRDT).o : $(DRDT).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(DRDT).f90
	\mv $(DRDT).o $(CDOTO)

$(CDOTO)/$(MATH).o : $(MATH).f90
	$(FC) $(FCOPTS) -I$(CDOTO) -c $(MATH).f90
	\mv $(MATH).o *.mod $(CDOTO)


###########################################
# Create object file directory if necessary
###########################################

makedir :
	if [ ! -d $(CDOTO) ] ; \
	then \
	   mkdir $(CDOTO) ; \
	fi

###################
# Intro information
###################

#.SILENT:
intro :
	echo
	echo '-----------------'
	echo '  $(OUTPUT) makefile'
	echo '-----------------'
	echo
	echo 'This makefile was designed for the $(PLATFORM) platform.'
	echo
	echo 'It uses the $(FC) compiler, with the following options:'
	echo '      $(FCOPTS)'
	echo 'except for the absorption coefficient data file which is'
	echo 'compiled without optimization (-O0) to improve compile time.'
	echo
	echo 'The object files used are as follows:'
	echo
	echo '      $(OBJFIL)'
	echo
	echo
	echo '$(OUTPUT) make is in progress ...'
	echo

###############################
# Concluding information
###############################

conclude :

	echo
	echo '================='
	echo '  makefile done  '
	echo '================='
	echo

print :
	echo '-----stop 1---------------'

