/*
 *  $Date: $
 *  $Revision: $
 *  
 */

#include "Fpmc.h"
#include "fostream.h"

#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cstring>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/PdfInfo.h>
#include <HepMC/HerwigWrapper.h>
#include <HepMC/HEPEVT_Wrapper.h>
#include <HepMC/IO_HERWIG.h>

#include "herwig.h"
#include "fpmc.h"

using namespace std;
using namespace fpmc;

Fpmc::Fpmc(double comEnergy, long int seed, vector<string> const& params):
  event_(0), 
  comEnergy_(comEnergy),
  seed_(seed), 
  params_(params) {
  
  herwigVerbosity_ = 1; 
  hepMCVerbosity_ = true;
  maxEventsToPrint_ = 2;
  hadronize_ = true;
  debug_ = false; 
}

Fpmc::~Fpmc(){
}

void Fpmc::write(std::ostream& out){
  hepMCEvt_->write(out);
}

void Fpmc::begin() {
  event_ = 0;

  // Write datacard  
  int iunit = 5;	
  fostream config(iunit,"datacard.txt");	
  /* ===============================    
  Example:    
  config << "TYPEPR      'INC'";
  config << "TYPINT      'QED'";
  config << "IPROC       16010";
  config << "NFLUX       15";
  config << "NRN1        33781";
  config << "NRN2        11776";
  config << "YJMAX       6.";
  config << "YJMIN      -6.";
  config << "PTMIN       10.";
  config << "IFIT        10";
  config << "ISOFTM      1";
  =============================== */

  std::vector<std::string> invalidParams;

  // Loop over all parameters and stop in case of error
  for( vector<string>::const_iterator itPar = params_.begin(); itPar != params_.end(); ++itPar ) {
     // Check for invalid parameters
     for(std::vector<std::string>::const_iterator itInvPar = invalidParams.begin(); itInvPar != invalidParams.end(); ++itInvPar){
        if( 0 == itPar->compare(0,itInvPar->size(),*itInvPar) ){
           stringstream oss;
           oss << "FpmcError: The following parameter is not accepted in this mode: " << endl
               << *itInvPar << endl;      
           throw runtime_error( oss.str() );
        }
     }
     // Pass string to datacard
     config << itPar->c_str();
  }

  // Use random seeds from datacard
  /* ===============================    
  Using CLHEP engines: 
  fRandomEngine = new CLHEP::HepJamesRandom(seed_);
  fRandomGenerator = new CLHEP::RandFlat(fRandomEngine);
 
  long seed0 = fRandomGenerator->fireInt(1L,179L);
  long seed1 = fRandomGenerator->fireInt(1L,179L);
  =============================== */  

  config.rewind();

  cout << "[FPMC Wrapper] Read datacard" << endl;

  int read = 1;
  fpmc_var_ini(&read); 

  if(debug_) cout << "[FPMC Wrapper] UTYPEPR = " << string(cc1.UTYPEPR) << endl; 
  if(debug_) cout << "[FPMC Wrapper] UTYPINT = " << string(cc2.UTYPINT) << endl; 
  if(debug_) cout << "[FPMC Wrapper] UTMASS  = " << myffread1.UTMASS << endl;

  cout << "[FPMC Wrapper] Initializing HERWIG/FPMC" << endl;

  // Call hwudat to set up HERWIG block data
  //hwudat();

  //PART1=UPART1
  //PART2=UPART2
  //
  for(int i=0;i<8;++i){
     hwbmch.PART1[i] = ' ';
     hwbmch.PART2[i] = ' ';
  } 
  memcpy(hwbmch.PART1, cc3.UPART1, 4);
  memcpy(hwbmch.PART2, cc4.UPART2, 4);

  hwproc.PBEAM1 = comEnergy_/2.;
  hwproc.PBEAM2 = comEnergy_/2.;

  //TYPEPR = UTYPEPR
  //TYPINT = UTYPINT
  //
  for(int i=0;i<3;++i){
     prtype.TYPEPR[i] = ' ';
     prtype.TYPINT[i] = ' ';
  }
  memcpy(prtype.TYPEPR, cc1.UTYPEPR, 3);
  memcpy(prtype.TYPINT, cc2.UTYPINT, 3);

  // 
  //IPROC = UIPROC
  //
  hwproc.IPROC = myffread3.UIPROC;

  //ANSWER=UHADR 
  //
  hadronize_ = strcmp(cc0.UHADR, "Y") == 0;  
  cout << "[FPMC Wrapper] Run hadronization/showering: " << hadronize_ << endl; 
  
  if(debug_) cout << "[FPMC Wrapper] PART1  = '" << hwbmch.PART1 << "'" << endl; 
  if(debug_) cout << "[FPMC Wrapper] PART2  = '" << hwbmch.PART2 << "'" << endl; 
  if(debug_) cout << "[FPMC Wrapper] TYPEPR = " << prtype.TYPEPR << endl; 
  if(debug_) cout << "[FPMC Wrapper] TYPINT = " << prtype.TYPINT << endl; 
  if(debug_) cout << "[FPMC Wrapper] IPROC  = " << hwproc.IPROC << endl; 

  //CALL HWIGIN
  //
  hwigin();

  // Read random seeds from datacard  
  //NRN(1) = UNRN1 ! set again later in the code
  //NRN(2) = UNRN2 ! set again later in the code
  //
  hwevnt.NRN[0] = myffread3.UNRN1;
  hwevnt.NRN[1] = myffread3.UNRN2;

  //EFFMIN = 1d-6
  //
  hwpram.EFFMIN = 1E-06;

  hwevnt.MAXER = 100000000; // O(inf)
  hwpram.LWSUD = 0;         // don't write Sudakov form factors
  hwdspn.LWDEC = 0;         // don't write three/four body decays
                            // (no fort.77 and fort.88 ...)
  // Init LHAPDF glue
  std::memset(hwprch.AUTPDF, ' ', sizeof(hwprch.AUTPDF));
  for(unsigned int i = 0; i < 2; i++) {
     hwpram.MODPDF[i] = -111;
     std::memcpy(hwprch.AUTPDF[i], "HWLHAPDF", 8);
  }

  hwevnt.MAXPR = maxEventsToPrint_;
  hwpram.IPRINT = herwigVerbosity_;

  hwpram.MODPDF[0] = cc5.UMODPDF1;
  hwpram.MODPDF[1] = cc5.UMODPDF2;

  //    PARAMETER(NMXRES=500)
  //    COMMON/HWPROP/RLTIM(0:NMXRES),RMASS(0:NMXRES),RSPIN(0:NMXRES),
  //   & ICHRG(0:NMXRES),IDPDG(0:NMXRES),IFLAV(0:NMXRES),NRES,
  //   & VTOCDK(0:NMXRES),VTORDK(0:NMXRES),
  //   & QORQQB(0:NMXRES),QBORQQ(0:NMXRES) */
  //RMASS(201) = UHMASS
  //RMASS(6)   = UTMASS ! top mass
  //RMASS(198) = UWMASS ! W mass
  //RMASS(406) = UMST1 ! stop1 mass
  //RMASS(405) = UMSB1 ! stop2 mass
  //
  hwprop.RMASS[201] = myffread1.UHMASS;
  hwprop.RMASS[6]   = myffread1.UTMASS;
  hwprop.RMASS[198] = myffread1.UWMASS;
  hwprop.RMASS[406] = myffread1.UMST1;
  hwprop.RMASS[405] = myffread1.UMSB1;

  //    COMMON/HWHARD/ASFIXD,CLQ(7,6),COSS,COSTH,CTMAX,DISF(13,2),EMLST,
  //   & EMMAX,EMMIN,EMPOW,EMSCA,EPOLN(3),GCOEF(7),GPOLN,OMEGA0,PHOMAS,
  //   & PPOLN(3),PTMAX,PTMIN,PTPOW,Q2MAX,Q2MIN,Q2POW,Q2WWMN,Q2WWMX,QLIM,
  //   & SINS,THMAX,Y4JT,TMNISR,TQWT,XX(2),XLMIN,XXMIN,YBMAX,YBMIN,YJMAX,
  //   & YJMIN,YWWMAX,YWWMIN,WHMIN,ZJMAX,ZMXISR,IAPHIG,IBRN(2),IBSH,
  //   & ICO(10),IDCMF,IDN(10),IFLMAX,IFLMIN,IHPRO,IPRO,MAPQ(6),MAXFL,
  //   & BGSHAT,COLISR,FSTEVT,FSTWGT,GENEV,HVFCEN,TPOL,DURHAM   */
  //Q2WWMN=UQ2WWMN
  //Q2WWMX=UQ2WWMX
  //
  hwhard.Q2WWMN = myffread2.UQ2WWMN;
  hwhard.Q2WWMX = myffread2.UQ2WWMX;

  //YWWMIN=UYWWMIN
  //YWWMAX=UYWWMAX
  //
  hwhard.YWWMIN = myffread2.UYWWMIN;
  hwhard.YWWMAX = myffread2.UYWWMAX;

  //YJMAX = UYJMAX
  //YJMIN = UYJMIN
  //
  hwhard.YJMAX = myffread1.UYJMAX;
  hwhard.YJMIN = myffread1.UYJMIN;
  
  //PTMIN=UPTMIN
  //PTMAX=UPTMAX
  //
  hwhard.PTMIN = myffread1.UPTMIN;
  hwhard.PTMAX = myffread1.UPTMAX;

  //EMMIN=UEMMIN
  //
  hwhard.EMMIN = myffread1.UEMMIN;

  //NFLUX = UNFLUX
  //
  xsect.NFLUX = myffread3.UNFLUX;

  //IFITPDF = UIFIT
  //
  pdfs.IFITPDF = myffread3.UIFIT;

  //AAANOM = UAAANOM
  //D_KAPPA = UDKAPPA
  //LAMBDA = UDLAMBDA
  //A0W = UA0W
  //ACW = UACW
  //A0Z = UA0Z
  //ACZ = UACZ
  //A1A = UA1A
  //A2A = UA2A
  //ANOMCUTOFF = UANOMCUTOFF
  //
  //AAEXOTIC = UAAEXOTIC
  //AAM = UAAM
  //AAQ = UAAQ
  //AAN = UAAN
  //
  aaanomal.AAANOM = myffread3.UAAANOM;
  aaanomal.D_KAPPA = myffread1.UDKAPPA;
  aaanomal.LAMBDA = myffread2.UDLAMBDA;
  aaanomal.A0W = myffread1.UA0W;
  aaanomal.ACW = myffread1.UACW;
  aaanomal.A0Z = myffread1.UA0Z;
  aaanomal.ACZ = myffread1.UACZ;
  aaanomal.A1A = myffread1.UA1A;
  aaanomal.A2A = myffread1.UA2A;
  aaanomal.ANOMCUTOFF = myffread2.UANOMCUTOFF;

  aaexotical.AAEXOTIC = myffread3.UAAEXOTIC;
  aaexotical.AAM = myffread1.UAAM;
  aaexotical.AAQ = myffread1.UAAQ;
  aaexotical.AAN = myffread1.UAAN;

  //CHIDeIGLU = UCHIDeIGLU
  //CHIDeX   =  UCHIDeX
  //CHIDeXP  =  UCHIDeXP
  //CHIDeS2  =  UCHIDeS2
  //CHIDeS   =  UECMS*UECMS
  //XI1MIN = UXI1MIN
  //XI1MAX = UXI1MAX
  //XI2MIN = UXI2MIN
  //XI2MAX = UXI2MAX
  //CHIDeGapMin = UCHIDeGapMin
  //CHIDeGapMax = UCHIDeGapMax
  //CHIDePATH = UCHIDePATH
  //
  chidefpmc.CHIDeIGLU = myffread3.UCHIDeIGLU;
  chidefpmc.CHIDeX = myffread1.UCHIDeX;
  chidefpmc.CHIDeXP = myffread1.UCHIDeXP;
  chidefpmc.CHIDeS2 = myffread1.UCHIDeS2;
  chidefpmc.CHIDeS = comEnergy_*comEnergy_;
  
  //KMR2Q2CUT=UKMR2Q2CUT
  //KMR2SURV=UKMR2SURV
  //KMR2SCALE=UKMR2SCALE
  //KMR2DELTA=UKMR2DELTA
  //
  kmr2fpmc.KMR2DELTA = myffread3.UKMR2DELTA;
  kmr2fpmc.KMR2Q2CUT = myffread1.UKMR2Q2CUT;
  kmr2fpmc.KMR2SURV = myffread1.UKMR2SURV;
  kmr2fpmc.KMR2SCALE = myffread1.UKMR2SCALE;

  //ISOFTM = UISOFTM
  //  
  xsect.ISOFTM = myffread3.UISOFTM;

  //ZION = UZION
  //AION = UAION
  //RBMIN = UBMIN
  ion.ZION = myffread3.UZION;
  ion.AION = myffread3.UAION;
  ion.RBMIN = myffread3.UBMIN;

  //c---Initialize model/pdf dependant parameters
  //CALL HWMODINI
  //
  hwmodini();
  
  //c---Compute parameter dependent constants
  //CALL HWUINC
  //
  hwuinc();
  
  //c---Check POMWIG Settings + Initialisations for consistency
  //CALL HWCHEK
  //
  hwchek();

  //c---Call HWUSTA to make any particle stable
  //CALL HWUSTA('PI0     ')      
  //
  int iopt = 1;
  int iwig = 0;
  char nwig[9] = "        ";

  int ipdg = 111;
  hwuidt(&iopt, &ipdg, &iwig, nwig);
  if(ipdg) hwusta(nwig, 1);
 
  //c---Initialize elementary process
  //CALL HWEINI
  //
  hweini();

  //c---Initialize event record fixing : this will replace the beam 
  //c   electrons by protons, radiated photons by pomerons/reggeons etc
  //CALL HWFXER(.TRUE.,IPROC)
  //
  int init = 1;
  hwfxer(&init); 

}

void Fpmc::end() {
 
}

bool Fpmc::run() {

   //c---Loop over events
   //      DO 100 N=1,MAXEV
   //c...Initialize event
   //         CALL HWUINE
   //c...Generate hard subprocesses
   //         CALL HWEPRO
   //c...Include showering and hadronization
   //         IF (ANSWER.EQ.'Y') THEN
   //            CALL HWBGEN
   //            CALL HWDHOB
   //            CALL HWCFOR
   //            CALL HWCDEC
   //            CALL HWDHAD
   //            CALL HWDHVY
   //            CALL HWMEVT
   //         END IF
   //c...Finish event
   //         CALL HWUFNE
   //c...Fix event record (i.e. restore correct intermediate states); print result
   //         CALL HWFXER(.FALSE.,IPROC)
   //         IF(N.LE.MAXPR) THEN
   //           PRINT*, ' '
   //           PRINT*, ' '
   //           PRINT*, ' '
   //           PRINT*, ' '
   //           PRINT*, 'AFTER EVENT RECORD FIXING:'
   //           CALL HWUEPR
   //         ENDIF
   // 100  CONTINUE

   // Call herwig routines to create HEPEVT
   //
   hwuine();	// initialize event

   hwepro();

   if(hadronize_){
      hwbgen();	// parton cascades

      hwdhob();	// heavy quark decays
      hwcfor();	// cluster formation
      hwcdec();	// cluster decays

      hwdhad();	// unstable particle decays
      hwdhvy();	// heavy flavour decays
      hwmevt();	// soft underlying event		
   }

   hwufne();	// finalize event

   int init = 0;
   hwfxer(&init);

   if (hwevnt.IERROR) {
      return false;
   }

   hepMCEvt_ = conv_.read_next_event();
   ++event_;

   hepMCEvt_->set_event_number(event_ - 1);
   hepMCEvt_->set_signal_process_id( hwproc.IPROC );
   hepMCEvt_->weights().push_back( hwevnt.EVWGT );

   hepMCEvt_->set_event_scale( -1. );

   HepMC::PdfInfo pdfInfo;
   pdfInfo.set_x1( hwhard.XX[0] );
   pdfInfo.set_x2( hwhard.XX[1] );
   pdfInfo.set_scalePDF( hwhard.EMSCA );
   hepMCEvt_->set_pdf_info( pdfInfo );

   /*HepMC::GenParticle* incomingParton = NULL;
   HepMC::GenParticle* targetParton = NULL;
   // find incoming parton (first entry with IST=121)
   for(HepMC::GenEvent::particle_const_iterator it = event()->particles_begin(); 
      (it != event()->particles_end() && incomingParton==NULL); it++)
      if((*it)->status()==121) incomingParton = (*it);
  
   // find target parton (first entry with IST=122)
   for(HepMC::GenEvent::particle_const_iterator it = event()->particles_begin(); 
      (it != event()->particles_end() && targetParton==NULL); it++)
      if((*it)->status()==122) targetParton = (*it);*/


   //******** Verbosity ********
   if(event_ <= maxEventsToPrint_ && hepMCVerbosity_) {
      // Prints HepMC event
      if(hepMCVerbosity_) {
	 stringstream oss;
	 oss << "\n----------------------" << endl	
	     << "Event process id = " << hepMCEvt_->signal_process_id() << endl; 
	 cout << oss.str();
	 hepMCEvt_->print();
      }
   }
   
   return true;
}
