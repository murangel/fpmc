#-------------------------------------------------------------------------------------
# * 29/01/2009, Oldrich Kepka 
# * Makefile for the inclusive or exclusive production with gaps
# 
# * module				 : module that accept parameters from the Format-Free card
#							   syntax to run: module < datacard.ffread
#							   histograms with HBOOK - cernlib assumed installed
# 
# * module_reco:			module with pxcone algorithm interface
#							   syntax to run: module < datacard.ffread
#							   histograms with HBOOK - cernlib assumed installed
# 
#-------------------------------------------------------------------------------------
# old specific examples
#-------------------------------------------------------------------------------------
#
# * example_excDPE    :standard example with exclusive bbbar production, in the BL model
#							  no external libraries required
#
# * example_incDPE    : an example with inclusive diphoton production, in the BPR model
#							  no external libraries required
#
# * example_hbook  	 : a standard example with exclusive bbbar production, in the BL model
#							  histograms with HBOOK - cernlib assumed installed
#
# * example_sd  		 :  dijets in single diffraction
#							  histograms with HBOOK - cernlib assumed installed
#
# * example_ddis		 :  Diffractive DIS
#							  histograms with HBOOK - cernlib assumed installed
#
# * example_qedww		 : pp --> p WW p via QED
#							  histograms with HBOOK - cernlib assumed installed
#
# * example_excChi    : an example of exclusive Chi_c production. This needs a special phase 
#							  space routine, prepared by Murilo Rangel
#							  histograms with HBOOK - cernlib assumed installed
#
# * example_incChi	 :  inclusive Chi_c production.
#							  histograms with HBOOK - cernlib assumed installed
#
#-------------------------------------------------------------------------------------


#########################################################################
# IMPORTANT USER SETUP
#########################################################################

CERNLIB=`cernlib mathlib pawlib packlib` -L$(PWD) -lgfortran
#CERNLIB=`cernlib mathlib pawlib packlib pdflib804` -L$(PWD) -lgfortran -L$(LHAPDFLIB) -lLHAPDF
#CERNLIB=`cernlib` -L$(PWD) -lgfortran
GSLLIB=`gsl-config --cflags --libs`

LHAPDF_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/lhapdf/5.9.1-cms3
LHAPDFLIB=-L$(LHAPDF_BASE)/lib -lLHAPDF
LHAPDF_INCLUDE=-I$(LHAPDF_BASE)/include

#########################################################################

#-------------------------------------------------------------------------------------
# Utilities
#-------------------------------------------------------------------------------------

all: allModules 
#allModules: module module_reco fpmc
allModules: fpmc
examples:oldExamples
allApp: Fpmc Herwig Externals examples allModules 
	
clean:clean_sqme clean_excl_aaaa
	@find ./ -name "*~" -exec rm -v {} \;
	@find ./ -name ".*.swp" -exec rm -v {} \;
	rm -f Objects/* module* fort.* *.hbook last.kumac *.ntp example_* *.mod lib/*.so lib/*.a
	
# FLAGS
# -------
  
# debuging flag to turn off anoying messages
#SPEC_FL=-Wno-globals -Wno-implicit

# g77 - setup: 
# flags very important for simulation interfaces!
#F_FLAGS = -g -O1 -Wno-all -fno-f2c -finit-local-zero -fno-automatic -Iinc
#F_COMP = g77 $(F_FLAGS) $(SPEC_FL)

# gforan - setup
F_FLAGS = -g -O1 -fno-automatic -Iinc
F_COMP = gfortran $(F_FLAGS) $(SPEC_FL)

# other 
CC=g++

# Directories
# -----------

OBJDIR  = Objects
DPEDIR	= Fpmc
HERDIR	= Herwig
EXTDIR	= External

# FPMC additions to Herwig
# -------------------------

Fpmc: $(OBJDIR)/fpmc.o


$(OBJDIR)/fpmc.o: $(DPEDIR)/fpmc.f 
	$(F_COMP) -c $^ -o $@

$(OBJDIR)/fpmc_chi.o: $(DPEDIR)/fpmc_chi.f
	$(F_COMP) -c $< -o $@


# HERWIG
# ------

Herwig: $(OBJDIR)/herwig6500.o 

$(OBJDIR)/herwig6500.o: $(HERDIR)/herwig6500.f
	$(F_COMP) -c $< -o $@


# External contributions
# ----------------------

ext=tools.f lininter.f
ext_obj=$(ext:%.f=%.o)
ext_obj_dest=$(ext_obj:%=$(OBJDIR)/%)
$(ext_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/%.f
	$(F_COMP) -c $< -o $@

# ----- pxcone algorithm and reconstruction 
reco=reco_ok/reco_ok.f reco_ok/pxcone_mod_new.f

reco_obj=$(reco:reco_ok/%.f=%.o)
reco_obj_dest=$(reco_obj:%=$(OBJDIR)/%)
$(reco_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/reco_ok/%.f
	$(F_COMP) -c $< -o $@

# ---- kmr luminosity
ext_kmr=kmrlumi/KMR.f kmrlumi/kmrLumi_ExHume_Tev.f kmrlumi/kmrLumi_ExHume_LHC.f kmrlumi/kmrLumi_Lonnblad_Tev.f kmrlumi/kmrLumi_Lonnblad_LHC.f

ext_kmr_obj=$(ext_kmr:kmrlumi/%.f=%.o)
ext_kmr_obj_dest=$(ext_kmr_obj:%=$(OBJDIR)/%)
$(ext_kmr_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/kmrlumi/%.f
	$(F_COMP) -c $< -o $@ 

# ----- KMR2 - direct implementation of KMR (based on ExHuME)
ext_KMR2=kmrlumi2/kmr2.f kmrlumi2/mrst2002.f 

ext_KMR2_obj=$(ext_KMR2:kmrlumi2/%.f=%.o)
ext_KMR2_obj_dest=$(ext_KMR2_obj:%=$(OBJDIR)/%)
$(ext_KMR2_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/kmrlumi2/%.f
	$(F_COMP) -c $< -o $@ 


################# CHIDe MODEL ##################

# -----  common files
ext_CHIDeCommon = \
CHIDe/Common/CHIDedcadredo.f  \
CHIDe/Common/CHIDedgd2008.f  \
CHIDe/Common/CHIDeFunctions.f \
CHIDe/Common/CHIDedgd2007.f  \
CHIDe/Common/CHIDedgdforward.f

ext_CHIDeCommon_obj=$(ext_CHIDeCommon:CHIDe/Common/%.f=%.o)
ext_CHIDeCommon_obj_dest=$(ext_CHIDeCommon_obj:%=$(OBJDIR)/%)
$(ext_CHIDeCommon_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/CHIDe/Common/%.f
	$(F_COMP) -c $< -o $@ 

# ----- Higgs
ext_CHIDeHiggs= \
CHIDe/Higgs/CHIDeHiggs.f \
CHIDe/Higgs/CHIDeHiggsInit.f

ext_CHIDeHiggs_obj=$(ext_CHIDeHiggs:CHIDe/Higgs/%.f=%.o)
ext_CHIDeHiggs_obj_dest=$(ext_CHIDeHiggs_obj:%=$(OBJDIR)/%)
$(ext_CHIDeHiggs_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/CHIDe/Higgs/%.f
	$(F_COMP) -c $< -o $@ 

# ----- GG 
ext_CHIDeGG=\
CHIDe/GG/CHIDeGG.f \
CHIDe/GG/CHIDeGGAmplitudes.f \
CHIDe/GG/CHIDeGGFunctions.f \
CHIDe/GG/CHIDeGGInit.f \
CHIDe/GG/CHIDeGGDurhamlike.f 

ext_CHIDeGG_obj=$(ext_CHIDeGG:CHIDe/GG/%.f=%.o)
ext_CHIDeGG_obj_dest=$(ext_CHIDeGG_obj:%=$(OBJDIR)/%)
$(ext_CHIDeGG_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/CHIDe/GG/%.f
	$(F_COMP) -c $< -o $@ 

# ----- Diphoton
ext_CHIDeDiphoton=\
CHIDe/Diphoton/CHIDeDiphoton.f \
CHIDe/Diphoton/CHIDeDiphotonFunctions.f \
CHIDe/Diphoton/CHIDeDiphotonAmplitudes.f \
CHIDe/Diphoton/CHIDeDiphotonInit.f \

ext_CHIDeDiphoton_obj=$(ext_CHIDeDiphoton:CHIDe/Diphoton/%.f=%.o)
ext_CHIDeDiphoton_obj_dest=$(ext_CHIDeDiphoton_obj:%=$(OBJDIR)/%)
$(ext_CHIDeDiphoton_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/CHIDe/Diphoton/%.f
	$(F_COMP) -c $< -o $@ 

################################################


# ---- softc 
ext_softc=softc/getsoftc.f softc/soft.2TeV.f softc/soft.14TeV.f softc/soft.2TeV.effopa.f softc/soft.14TeV.effopa.f

ext_softc_obj=$(ext_softc:softc/%.f=%.o)
ext_softc_obj_dest=$(ext_softc_obj:%=$(OBJDIR)/%)
$(ext_softc_obj_dest): $(OBJDIR)/%.o: $(EXTDIR)/softc/%.f
	$(F_COMP) -c $< -o $@ 

#---- Parton densities
ext_pdf = pdf/i_2006_fita.f pdf/i_2006_fitb.f pdf/qcd_2006.f pdf/h1qcd.f pdf/i_nlotables.f
ext_pdf_obj=$(ext_pdf:pdf/%.f=%.o)
ext_pdf_dest=$(ext_pdf_obj:%=$(OBJDIR)/%)
$(ext_pdf_dest): $(OBJDIR)/%.o: $(EXTDIR)/pdf/%.f
	$(F_COMP) -c $< -o $@

# interface to comphep routines 
ext_comphep_dest=$(OBJDIR)/comphep_wraper.o $(OBJDIR)/sqme_aaww.a $(OBJDIR)/sqme_aazz.a $(OBJDIR)/sqme_aaaa.a

$(OBJDIR)/sqme_aaww.a:
	$(MAKE) -C $(EXTDIR)/comphep_interface/sqme_aaww	
	cp -f $(EXTDIR)/comphep_interface/sqme_aaww/sqme_aaww.a $(OBJDIR)

$(OBJDIR)/sqme_aazz.a:
	$(MAKE) -C $(EXTDIR)/comphep_interface/sqme_aazz	
	cp -f $(EXTDIR)/comphep_interface/sqme_aazz/sqme_aazz.a $(OBJDIR)

$(OBJDIR)/sqme_aaaa.a:
	$(MAKE) -C $(EXTDIR)/comphep_interface/sqme_aaaa	
	cp -f $(EXTDIR)/comphep_interface/sqme_aaaa/sqme_aaaa.a $(OBJDIR)

$(OBJDIR)/comphep_wraper.o:$(EXTDIR)/comphep_interface/comphep_wraper.cpp
	$(CC) -c -o $@ $<

clean_sqme:
	$(MAKE) -C $(EXTDIR)/comphep_interface/sqme_aaww clean
	$(MAKE) -C $(EXTDIR)/comphep_interface/sqme_aazz clean
	$(MAKE) -C $(EXTDIR)/comphep_interface/sqme_aaaa clean


# interface to excl aa->aa routines
ext_excl_aaaa_dest=$(OBJDIR)/excl_aaaa_wraper.o $(OBJDIR)/sm_sqme_aaaa.a $(OBJDIR)/bsmf_sqme_aaaa.a $(OBJDIR)/bsmv_sqme_aaaa.a $(OBJDIR)/eft_sqme_aaaa.a

$(OBJDIR)/sm_sqme_aaaa.a:
	$(MAKE) -C $(EXTDIR)/excl_aaaa/sm_sqme_aaaa	
	cp -f $(EXTDIR)/excl_aaaa/sm_sqme_aaaa/sm_sqme_aaaa.a $(OBJDIR)

$(OBJDIR)/bsmf_sqme_aaaa.a:
	$(MAKE) -C $(EXTDIR)/excl_aaaa/bsmf_sqme_aaaa	
	cp -f $(EXTDIR)/excl_aaaa/bsmf_sqme_aaaa/bsmf_sqme_aaaa.a $(OBJDIR)

$(OBJDIR)/bsmv_sqme_aaaa.a:
	$(MAKE) -C $(EXTDIR)/excl_aaaa/bsmv_sqme_aaaa	
	cp -f $(EXTDIR)/excl_aaaa/bsmv_sqme_aaaa/bsmv_sqme_aaaa.a $(OBJDIR)

$(OBJDIR)/eft_sqme_aaaa.a:
	$(MAKE) -C $(EXTDIR)/excl_aaaa/eft_sqme_aaaa	
	cp -f $(EXTDIR)/excl_aaaa/eft_sqme_aaaa/eft_sqme_aaaa.a $(OBJDIR)

$(OBJDIR)/excl_aaaa_wraper.o:$(EXTDIR)/excl_aaaa/excl_aaaa_wraper.cpp
	$(CC) -c -o $@ $<

clean_excl_aaaa:
	$(MAKE) -C $(EXTDIR)/excl_aaaa/sm_sqme_aaaa clean
	$(MAKE) -C $(EXTDIR)/excl_aaaa/bsmf_sqme_aaaa clean
	$(MAKE) -C $(EXTDIR)/excl_aaaa/bsmv_sqme_aaaa clean
	$(MAKE) -C $(EXTDIR)/excl_aaaa/eft_sqme_aaaa clean




# ----- user objects
# ----------------------
# ntuple
$(OBJDIR)/ntuple.o:External/ntuple.f 
	$(F_COMP) -c $< -o $@

# reading datacards
$(OBJDIR)/ffcard.o:Examples/ffcard.f Examples/ffcard.inc
	$(F_COMP) -c $< -o $@

# LHE functions
$(OBJDIR)/fpmc_lhe.o:Examples/fpmc_lhe.f
	$(F_COMP) -c $< -o $@

# ----------------
# Objects variables
# ----------------
OBJSTAND = $(OBJDIR)/herwig6500.o  $(OBJDIR)/fpmc.o $(OBJDIR)/ffcard.o $(OBJDIR)/fpmc_lhe.o
OBJEXT   = $(ext_obj_dest) $(ext_pdf_dest) $(ext_comphep_dest) $(ext_excl_aaaa_dest) $(ext_kmr_obj_dest)  $(ext_softc_obj_dest) \
	$(ext_CHIDeCommon_obj_dest) $(ext_CHIDeHiggs_obj_dest) $(ext_KMR2_obj_dest) \
	$(ext_CHIDeGG_obj_dest) $(ext_CHIDeDiphoton_obj_dest) 
OBJUSR   = $(OBJDIR)/ntuple.o
LIBS     = $(CERNLIB) $(LHAPDFLIB) $(GSLLIB) $(LIB_OMEGA)
OBJRECO  = $(reco_obj_dest)


# ----------------
# Example programs
# ----------------

# old examples
pgm=example_excDPE example_incDPE example_hbook example_sd example_ddis \
	 example_qedww example_excChi example_incChi
$(pgm): % : Examples/%.f $(OBJSTAND) $(OBJEXT)
	$(F_COMP) $< -o $@ $(OBJSTAND) $(OBJEXT) $(LIBS)  -lstdc++
oldExamples:$(pgm)

### modules to be used by the end user###
fpmc: Examples/fpmc_main.f $(OBJEXT) $(OBJSTAND)  $(OBJUSR) $(OBJRECO) Examples/ffcard.inc
	$(F_COMP) -o $@ $< $(OBJSTAND) $(OBJEXT) $(OBJUSR)  $(OBJRECO) $(LIBS) -lstdc++  	
module_reco: Examples/module_reco.f $(OBJEXT) $(OBJSTAND)  $(OBJUSR) $(OBJRECO) Examples/ffcard.inc
	$(F_COMP) -o $@ $< $(OBJSTAND) $(OBJEXT) $(OBJUSR)  $(OBJRECO) $(LIBS) -lstdc++  	
module: Examples/module.f $(OBJEXT) $(OBJSTAND)  $(OBJUSR) \
  	Examples/ffcard.inc
	$(F_COMP) -o $@ $< $(OBJSTAND) $(OBJEXT) $(OBJUSR)  $(LIBS) -lstdc++ 

#---- Wrapper/HepMC
#----
HEPMC_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07-cms4
HEPMCLIB=-L$(HEPMC_BASE)/lib -lHepMCfio -lHepMC
HEPMC_INCLUDE=-I$(HEPMC_BASE)/include

BOOST_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.51.0-cms2
BOOSTLIB=-L$(BOOST_BASE)/lib -lboost_thread -lboost_signals -lboost_date_time
BOOST_INCLUDE=-I$(BOOST_BASE)/include

LIBDIR  = lib
CFLAGS  = -g -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
LDFLAGS = -g -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC

wrapper_f=HepMCWrapper/f77out.f HepMCWrapper/hwaend_dummy.f
wrapper_f_obj=$(wrapper_f:HepMCWrapper/%.f=%.o)
wrapper_f_obj_dest=$(wrapper_f_obj:%=$(OBJDIR)/%)

$(wrapper_f_obj_dest): $(OBJDIR)/%.o: HepMCWrapper/%.f
	$(F_COMP) -fPIC -c $< -o $@

wrapper=HepMCWrapper/fostream.cc HepMCWrapper/Fpmc.cc
wrapper_obj=$(wrapper:HepMCWrapper/%.cc=%.o)
wrapper_obj_dest=$(wrapper_obj:%=$(OBJDIR)/%)

$(wrapper_obj_dest): $(OBJDIR)/%.o: HepMCWrapper/%.cc
	$(CC) $(CFLAGS) $(HEPMC_INCLUDE) $(LHAPDF_INCLUDE) -c $< -o $@

$(LIBDIR)/FPMCHepMCWrapper.so:$(wrapper_f_obj_dest) $(wrapper_obj_dest)
	mkdir -p $(LIBDIR); \
	$(CC) $(LDFLAGS) -shared $(wrapper_f_obj_dest) $(wrapper_obj_dest) $(HEPMCLIB) -o $@

$(LIBDIR)/FPMCHepMCWrapper.a:$(wrapper_f_obj_dest) $(wrapper_obj_dest)
	mkdir -p $(LIBDIR); \
	ar -r $@ $(wrapper_f_obj_dest) $(wrapper_obj_dest)

$(OBJDIR)/fpmc-hepmc.o: HepMCWrapper/main.cc
	$(CC) $(CFLAGS) $(BOOST_INCLUDE) $(HEPMC_INCLUDE) $(LHAPDF_INCLUDE) -c $< -o $@

fpmc-hepmc: \
$(OBJDIR)/herwig6500.o \
$(OBJDIR)/fpmc.o \
$(OBJDIR)/ffcard.o \
$(OBJEXT) \
$(wrapper_f_obj_dest) $(wrapper_obj_dest) \
$(OBJDIR)/fpmc-hepmc.o
	$(CC) $(LDFLAGS) $(OBJDIR)/herwig6500.o $(OBJDIR)/fpmc.o $(OBJDIR)/ffcard.o $(OBJEXT) \
	$(wrapper_f_obj_dest) $(wrapper_obj_dest) $(OBJDIR)/fpmc-hepmc.o \
	$(CERNLIB) $(LHAPDFLIB) $(GSLLIB) $(LIB_OMEGA) $(HEPMCLIB) -o $@
#----
