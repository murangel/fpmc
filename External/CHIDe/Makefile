objdir=Objects

Common_code=Common/CHIDedcadredo.f  Common/CHIDedgd2007.f  Common/CHIDedgd2008.f\
	    Common/CHIDedgdforward.f  Common/CHIDeFunctions.f

GG_code=GG/CHIDeGG.f            GG/CHIDeGGAmplitudes.f      \
	GG/CHIDeGGFunctions.f   GG/CHIDeGGInit.f  \
      	GG/CHIDeGGDurhamlike.f       

Diphoton_code=Diphoton/CHIDeDiphoton.f            Diphoton/CHIDeDiphotonAmplitudes.f      \
	      GG/CHIDeGGFunctions.f Diphoton/CHIDeDiphotonFunctions.f  Diphoton/CHIDeDiphotonInit.f

Higgs_code=Higgs/CHIDeHiggs.f  Higgs/CHIDeHiggsInit.f

Standalone_code=Standalone/Standalone.f  Standalone/StandaloneGG.f  Standalone/StandaloneHiggs.f\
		Standalone/vegas11.f Standalone/VegasGG.f Standalone/VegasHiggs.f Standalone/Vegas.f\
		Standalone/StandaloneDiphoton.f Standalone/VegasDiphoton.f Standalone/StandaloneDiphoton2.f


Common_obj=$(Common_code:%.f=$(objdir)/%.o)
GG_obj=$(GG_code:%.f=$(objdir)/%.o)
Diphoton_obj=$(Diphoton_code:%.f=$(objdir)/%.o)
Higgs_obj=$(Higgs_code:%.f=$(objdir)/%.o)
Standalone_obj=$(Standalone_code:%.f=$(objdir)/%.o)

obj=$(Common_obj) $(GG_obj) $(Diphoton_obj) $(Higgs_obj) $(Standalone_obj) 
code=$(Common_code) $(GG_code) $(Diphoton_code) $(Higgs_code) $(Standalone_code) 

#FF = g77 -ggdb
FF = gfortran -ggdb

all: CHIDeGG CHIDeHiggs CHIDeDiphoton CHIDe VegasGG VegasDiphoton CHIDeDiphoton2

$(obj): $(objdir)/%.o : %.f
	$(FF) $< -c -o $@

VegasHiggs: $(Common_obj) $(Higgs_obj) Objects/Standalone/VegasHiggs.o Objects/Standalone/vegas11.o
	$(FF) $(Common_obj) $(Higgs_obj) Objects/Standalone/VegasHiggs.o Objects/Standalone/vegas11.o -o VegasHiggs
	
VegasDiphoton: $(Diphoton_obj) $(Common_obj) Objects/Standalone/vegas11.o Objects/Standalone/VegasDiphoton.o 
	$(FF) $(Diphoton_obj) $(Common_obj) Objects/Standalone/vegas11.o  Objects/Standalone/VegasDiphoton.o -o VegasDiphoton

VegasGG: $(GG_obj) $(Common_obj) Objects/Standalone/vegas11.o Objects/Standalone/VegasGG.o 
	$(FF) $(GG_obj) $(Common_obj) Objects/Standalone/vegas11.o  Objects/Standalone/VegasGG.o -o VegasGG
	
CHIDeDiphoton: $(Common_obj) $(Diphoton_obj) Objects/Standalone/StandaloneDiphoton.o
	$(FF) $(Common_obj) $(Diphoton_obj) Objects/Standalone/StandaloneDiphoton.o -o CHIDeDiphoton

CHIDeGG: $(Common_obj) $(GG_obj) Objects/Standalone/StandaloneGG.o
	$(FF) $(Common_obj) $(GG_obj) Objects/Standalone/StandaloneGG.o -o CHIDeGG
	
CHIDeHiggs: $(Common_obj) $(Higgs_obj) Objects/Standalone/StandaloneHiggs.o
	$(FF) $(Common_obj) $(Higgs_obj) Objects/Standalone/StandaloneHiggs.o -o CHIDeHiggs
	
CHIDe: $(Common_obj) $(Higgs_obj) $(GG_obj)  Objects/Standalone/Standalone.o
	$(FF) $(Common_obj) $(Higgs_obj) $(GG_obj)  Objects/Standalone/Standalone.o -o CHIDe

Vegas: $(Common_obj) $(Higgs_obj) $(GG_obj) Objects/Standalone/vegas11.o Objects/Standalone/Vegas.o  
	$(FF) $(Common_obj) $(Higgs_obj) $(GG_obj) Objects/Standalone/vegas11.o Objects/Standalone/Vegas.o -o Vegas

clean:
	rm -f *~ */*~ $(objdir)/*/*
	rm -f CHIDe CHIDeHiggs CHIDeGG CHIDeDiphoton Vegas VegasGG VegasHiggs VegasDiphoton
	rm -f fort.*
