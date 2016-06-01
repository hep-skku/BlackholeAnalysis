C----------------------------------------------------------------------
C----------------------------------------------------------------------
C              Initilisation and Event Generation Routines
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, UPINIT.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Chris Harris and Peter Richardson
C----------------------------------------------------------------------
      SUBROUTINE UPINIT
C----------------------------------------------------------------------
C     Les Houches initialisation routine
C----------------------------------------------------------------------
#if !defined(HERWIG)
      IMPLICIT NONE
#else
      INCLUDE 'HERWIG65.INC'
#endif
#if defined(CHARYBDIS)
      INTEGER NRN
      COMMON/CHRANDOM/NRN(2)
#endif
      INCLUDE 'charybdis2.inc'
C--Les Houches run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &                IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &                XMAXUP(MAXPUP),LPRUP(MAXPUP)
C--event common block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)

C--local variables
      INTEGER NCHSRCH,I
      DOUBLE PRECISION CHWSUM,CHWSQ,HWRSET,DUMMY
      EXTERNAL HWRSET
C--parameters for the initial maximum weight search
      DATA NCHSRCH /100000/
C--Default paramter values are set in BLOCK DATA CHDATA,
C  so that they can be reset in main program.
C--You can still change them here instead if you want.
C--Beam particles and energies (only proton=2212 and antiproton=-2212)
C--N.B.for HERWIG these are overidden by the values in the main program
C      IDBMUP(1) = 2212
C      IDBMUP(2) = 2212
C      EBMUP(1)  = 7000.0D0
C      EBMUP(2)  = 7000.0D0
C--Set MPLNCK and define what is meant by it:
C--MSSDEF=1 means M_GT, MSSDEF=2 means M_DL and MSSDEF=3 means M_D
C      MPLNCK=1000.0D0
C      MSSDEF=3
C--Set number of dimensions (number of EXTRA dimensions is n=TOTDIM-4)
C--TOTDIM can be in the range 6-11
C      TOTDIM=6
C--Use Giddings+Thomas momentum scale for calling pdfs?
C--(See page 12 of hep-ph/0106219)
C      GTSCA=.FALSE.
C--Set mass window for black holes produced
C      MINMSS=5000.0D0
C      MAXMSS=EBMUP(1)+EBMUP(2)
C--Set NBODY decay of BH remnant - NBODY can be in range 2-5.
C      NBODY=2
C--Turn time variation of BH Hawking temperature in decay on or off
C      TIMVAR=.TRUE.
C--Set which decay products are to be used:
C--MSSDEC=1 gives no heavy particles
C--MSSDEC=2 gives t, W and Z as well
C--MSSDEC=3 gives t, W, Z and Higgs as well  
C      MSSDEC=3
C--Turn grey-body factors on/off
C      GRYBDY=.TRUE.
C--Turn kinematic cut-off of decay on (rather than M=MPLANCK cut-off)
C      KINCUT=.FALSE.
C--Use Yoshino-Rychkov factors in cross-section
C      YRCSEC= .FALSE.
C--Max Hawking temperature
C      THWMAX =1000.0D0
C--Remnant decay model: boiling at THWMAX
C      RMBOIL = .FALSE.
C--Min BH mass, ends boiling
C      RMMINM = 1000.D0
C--New parameters for version 1.004
C--Include effects of BH spin
C      BHSPIN = .FALSE.
C--Vary BH spin axis during decay (if BHSPIN)
C      BHJVAR = .TRUE.
C--Include anisotropy of decay products (if BHSPIN)
C      BHANIS = .TRUE.
C--New parameters for version 2.001
C--/Total lepton number conservation, Electron lepton number conservation, Tau lepton number conservation, Mu lepton number conservation/
C     NLEPTONCSV(0)=.FALSE.
C     NLEPTONCSV(1)=.FALSE.
C     NLEPTONCSV(2)=.FALSE.
C     NLEPTONCSV(3)=.FALSE.
C--Let BH lose mass and angular momentum in production
C      MJLOST = .TRUE.
C--Use MINMSS for the minimum mass of the black holes generated (default .TRUE.). Otherwise use it as minimum partonic centre of mass energy (.FALSE.)
C      USEMINMSSBH=.TRUE.
C--Mean mass fraction lost: <MLOST> = FMLOST*(ECM-MBH) (power law)
C      FMLOST=0.99D0 
C--Bias towards constant horizon angular velocity lines in the allowed M,J plane for black hole masses
C      CVBIAS = .FALSE.
C--Stable remnant model
C      RMSTAB = .TRUE.
C--Choose the recoil model. Default is emission in rest frame of initial BH
C      RECOIL = 1
C--Use the estimate for <N> Nto terminate evaporation. Default is NBODYAVERAGE=.TRUE.
C      NBODYAVERAGE = .TRUE.
C--Variable NBODY decay. Default is NBODYVAR=.TRUE.
C      NBODYVAR = .TRUE.
C--Use phase space or hawking spectrum for NBODY decay. Default is NBODYPHASE=.FALSE.
C      NBODYPHASE = .FALSE.
C--Skip evaporation and go to remnant. Default is SKIP2REMNANT=.FALSE.
C      SKIP2REMNANT=.FALSE.
C--print out option (0=no printout, 1 =errors only, 2=errors+info)
C      IBHPRN = 1
C--Les houches event file name and histories file name (must be 8 characters long)
C      LHEFILENAME = 'lhouches'
C      HISFILENAME = 'histfile'
C--pdf's for the beams (use the generator default)
C--MUST BE THE SAME FOR BOTH BEAMS
C      PDFGUP(1) = -1
C      PDFGUP(2) = -1
C      PDFSUP(1) = -1
C      PDFSUP(2) = -1
C--LHAPDF set number
C      LHAPDFSET = 10000
C--NMAXEV the maximum number of bh events to generate
C      CHNMAXEV = 100
C--------------------------------------------------------------------
C     END OF USER VARIABLES DON'T TOUCH ANYTHING BELOW HERE
C--------------------------------------------------------------------

#if defined(CHARYBDIS)
C--SET SEED OF RANDOM NUMBER GENERATOR
      DUMMY=HWRSET(NRN)
#endif

C--what to do with the weights(here are generating weighted events)
      IDWTUP = 1
C--only one process
      NPRUP  = 1
C--communication code
      LPRUP(1) = 1
C--calculate the maximum weight
C--Initialise to zero
      CHWSUM = 0.0D0
      CHWSQ = 0.0D0
      XMAXUP(1) = 0.0D0
      DO I=1,NCHSRCH
C--Generate weigthed events, determine maximum and statistical error
         CALL CHEVNT(.FALSE.)
         CHWSUM = CHWSUM+XWGTUP
         CHWSQ  = CHWSQ+XWGTUP**2
         XMAXUP(1) = MAX(XMAXUP(1),XWGTUP)
      ENDDO
      CHWSUM = CHWSUM/DBLE(NCHSRCH)
      CHWSQ  = MAX(CHWSQ/DBLE(NCHSRCH)-CHWSUM**2,0.0D0)
      CHWSQ  = SQRT(CHWSQ/ DBLE(NCHSRCH))
C--cross section
      XSECUP(1) = CHWSUM
C--error on the cross section
      XERRUP(1) = CHWSQ
C--output initialisation information
C--header
      WRITE(*,10)
C--beam parameters
      WRITE(*,11) IDBMUP(1),EBMUP(1),IDBMUP(2),EBMUP(2)
C--basic parameters
#if !defined(GINGRICH)
      WRITE(*,12) MINMSS,MAXMSS,MPLNCK,THWMAX,RMMINM,MSSDEF,TOTDIM,
     &     RECOIL,YRCSEC,TIMVAR,GRYBDY,KINCUT,RMBOIL,RMSTAB,BHSPIN,
     &     NLEPTONCSV(0),NLEPTONCSV(1),NLEPTONCSV(2),NLEPTONCSV(3),
     &     MJLOST,USEMINMSSBH,CVBIAS,NBODYVAR,NBODYPHASE,NBODY,FMLOST,
     &     SKIP2REMNANT
#else
      WRITE(*,12) MINMSS,MAXMSS,MPLNCK,DGMS,DGGS,THWMAX,RMMINM,
     &            MSSDEF,TOTDIM,RECOIL,DGSB,YRCSEC,TIMVAR,GRYBDY,
     &            KINCUT,RMBOIL,RMSTAB,BHSPIN,NLEPTONCSV(0),
     &            NLEPTONCSV(1),NLEPTONCSV(2),NLEPTONCSV(3),MJLOST,
     &            USEMINMSSBH,CVBIAS,NBODYVAR,NBODYPHASE,FMLOST,
     &            SKIP2REMNANT
#endif
      IF (RMBOIL) THEN
         IF (KINCUT) THEN
            WRITE (*,8)
 8          FORMAT(/10X,'KINEMATIC CUT INCONSISTENT'/
     &              10X,'WITH BOILING REMNANT MODEL:'/
     &              10X,'RESETTING KINCUT = .FALSE.'/)
            KINCUT=.FALSE.
         ENDIF
      ELSE
         IF (RMMINM.LT.MPLNCK) THEN
#if !defined(GINGRICH)
            WRITE (*,9)
 9          FORMAT(/10X,'BOILING REMNANT MODEL NOT USED:'/ 
     &              10X,'RESETTING MIN REMNANT MASS = MPLNCK'/)
            RMMINM=MPLNCK
#else
            IF (DGSB) THEN 
              WRITE (*,91)
 91           FORMAT(/10X,'BOILING REMNANT MODEL NOT USED:'/ 
     &                10X,'RESETTING MIN REMNANT MASS = DGMS'/)
              RMMINM=DGMS
            ELSE
              WRITE (*,9)
 9            FORMAT(/10X,'BOILING REMNANT MODEL NOT USED:'/ 
     &                10X,'RESETTING MIN REMNANT MASS = MPLNCK'/)
              RMMINM=MPLNCK
            ENDIF
#endif
         ENDIF
      ENDIF
      IF ((.NOT.TIMVAR).AND.(BHJVAR)) THEN
         WRITE(*,21) 
         BHJVAR=.FALSE.
      END IF
      IF (BHSPIN) THEN
         WRITE (*,20) BHJVAR,BHANIS
      ENDIF
C--choice of outgoing particles
      IF(MSSDEC.EQ.1) THEN
         WRITE(*,13)
      ELSEIF(MSSDEC.EQ.2) THEN
         WRITE(*,14)
      ELSEIF(MSSDEC.EQ.3) THEN
         WRITE(*,15)
#if defined (GINGRICH)
      ELSEIF(MSSDEC.EQ.4) THEN
         WRITE(*,151)
         IF (BHSPIN) THEN
           WRITE(*,*) 'Brane tension model does not work for BHSPIN=T'
           STOP
         ENDIF
         IF (NBODY.NE.2) THEN
           WRITE(*,*) 'Brane tension model only works for NBODY=2'
           STOP
         ENDIF
         IF (NBODYVAR) THEN
           WRITE(*,*) 'Brane tension model does not work for NBODYVAR=T'
           STOP
         ENDIF
#endif
      ENDIF
      IF(MINMSS.LT.MPLNCK) WRITE(*,22)
#if defined (GINGRICH)
      IF (DGSB) DGTENSION =(DGMS/MPLNCK)**4/(DGGS*(2.0D0*ACOS(-1.0))**3)
      WRITE(*,152) DGTENSION
#endif
C--choice of scale
      IF(GTSCA) THEN
         WRITE(*,17)
      ELSE
         WRITE(*,18)
      ENDIF
C--information on the cross section
      WRITE(*,19) XSECUP(1),XERRUP(1),XMAXUP(1)
C--particle production probabilites
      WRITE(*,*)
      WRITE(*,*)'         Production probabilities according to '
      WRITE(*,*)'         integrated fluxes and a_star dependent'
      WRITE(*,*)

      RETURN
C--format statements for the output
 10   FORMAT(//10X,'       CHARYBDIS2-1.0.3   February 2011      ',//,
     &         10X,'Please reference: C.M.Harris, P.Richardson &',/,
     &         10X,'B.R.Webber,JHEP0308(2003)033 [hep-ph/0307305]',/,
     &         10X,'and J.A.Frost, J.R.Gaunt, M.O.P.Sampaio,',/,
     &         10X,'M. Casals, S.R.Dolan, M.A.Parker & B.R. Webber',/,
#if !defined(GINGRICH)
     &         10X,'ArXiv:0904.0979 [hep-th]'/)
#else
     &         10X,'ArXiv:0904.0979 [hep-th]'//,
     &         10X,'Graviton and recoil added by D.M.Gingrich'/,
     &         10X,'JHEP0711(2007)064; 0706.0623 [hep-ph]'/,
     &         10X,'String balls added by D.M.Gingrich & K.Martell'/,
     &         10X,'Phys.Rev.D 78(2007)115009; 0808.2512 [hep-ph]'/)
#endif
 11   FORMAT(/10X,'INPUT CONDITIONS FOR THIS RUN'//
     &        10X,'BEAM 1 (',I8,') ENG. =',F10.2/
     &        10X,'BEAM 2 (',I8,') ENG. =',F10.2/)
#if !defined(GINGRICH)
 12   FORMAT(
     &     10X,'MIN PARTONIC CM ENERGY=',F11.2/
     &     10X,'MAX PARTONIC CM ENERGY=',F11.2/
     &     10X,'PLANCK MASS           =',F11.2/
     &     10X,'MAX HAWKING TEMP      =',F11.2/
     &     10X,'MIN REMNANT MASS      =',F11.2/
     &     10X,'DEFN OF PLANCK MASS   =',I5/
     &     10X,'NO. OF DIMENSIONS     =',I5/
     &     10X,'RECOIL OPTION         =',I5/
     &     10X,'YOSHINO-RYCHKOV C-S   =',L5/
     &     10X,'TIME VARIATION        =',L5/
     &     10X,'GREY BODY EFFECTS     =',L5/
     &     10X,'KINEMATIC CUT         =',L5/
     &     10X,'BOILING REMN. MODEL   =',L5/
     &     10X,'STABLE  REMN. MODEL   =',L5/
     &     10X,'BH SPIN INCLUDED      =',L5/
     &     10X,'TOTAL NLEPTON CONSERV.=',L5/
     &     10X,'NELECTLEPTON CONSERV. =',L5/
     &     10X,'NTAULEPTON CONSERV.   =',L5/
     &     10X,'NMULEPTON CONSERV.    =',L5/
     &     10X,'MASS & ANG MOM LOST   =',L5/
     &     10X,'MINMSS->MBH THRESHOLD =',L5/
     &     10X,'CONST. V. BIAS        =',L5/
     &     10X,'NBODY VARIABLE        =',L5/
     &     10X,'NBODY PHASE           =',L5/
     &     10X,'NBODY                 =',I5/
     &     10X,'FMLOST                =',F11.2/
     &     10X,'SKIP TO REMNANT       =',L5)
#else
 12   FORMAT(
     &     10X,'MIN PARTONIC CM ENERGY=',F11.2/
     &     10X,'MAX PARTONIC CM ENERGY=',F11.2/
     &     10X,'PLANCK MASS           =',F11.2/
     &     10X,'STRING SCALE          =',F11.2/
     &     10X,'STRING COUPLING       =',F11.2/
     &     10X,'MAX HAWKING TEMP      =',F11.2/
     &     10X,'MIN REMNANT MASS      =',F11.2/
     &     10X,'DEFN OF PLANCK MASS   =',I5/
     &     10X,'NO. OF DIMENSIONS     =',I5/
     &     10X,'RECOIL OPTION         =',I5/
     &     10X,'STRING BALLS          =',L5/
     &     10X,'YOSHINO-RYCHKOV C-S   =',L5/
     &     10X,'TIME VARIATION        =',L5/
     &     10X,'GREY BODY EFFECTS     =',L5/
     &     10X,'KINEMATIC CUT         =',L5/
     &     10X,'BOILING REMN. MODEL   =',L5/
     &     10X,'STABLE  REMN. MODEL   =',L5/
     &     10X,'BH SPIN INCLUDED      =',L5/
     &     10X,'TOTAL NLEPTON CONSERV.=',L5/
     &     10X,'NELECTLEPTON CONSERV. =',L5/
     &     10X,'NTAULEPTON CONSERV.   =',L5/
     &     10X,'NMULEPTON CONSERV.    =',L5/
     &     10X,'MASS & ANG MOM LOST   =',L5/
     &     10X,'MINMSS->MBH THRESHOLD =',L5/
     &     10X,'CONST. V. BIAS        =',L5/
     &     10X,'NBODY VARIABLE        =',L5/
     &     10X,'NBODY PHASE           =',L5/
     &     10X,'FMLOST                =',F11.2/
     &     10X,'SKIP TO REMNANT       =',L5)
#endif
 13   FORMAT(
     &     10X,'ONLY LIGHT SM PARTICLES PRODUCED')
 14   FORMAT(
     &     10X,'ALL SM PARTICLES BUT HIGGS PRODUCED')
 15   FORMAT(
     &     10X,'ALL SM PARTICLES PRODUCED')
#if defined (GINGRICH)
 151  FORMAT(
     &     10X,'GRAVITONS INCLUDED')
 152  FORMAT(
     &     10X,'BRANE TENSION',E12.4)
#endif
 17   FORMAT(
     &     10X,'USING BLACK HOLE RADIUS AS SCALE')
 18   FORMAT(
     &     10X,'USING BLACK HOLE MASS AS SCALE')
 19   FORMAT(
     &     10X,'CROSS SECTION (PB) =',G12.4/
     &     10X,'ERROR IN C-S  (PB) =',G12.4/
     &     10X,'MAXIMUM WEIGHT     =',E12.4)
 20   FORMAT(
     &     10X,'RECOILING BH SPIN  =',L5/
     &     10X,'ANISOTROPIC DECAY  =',L5)
 21   FORMAT(
     &     10X,'TIMVAR OFF. RESETING BHJVAR TO .FALSE.')
 22   FORMAT(/
     &     10X,'WARNING: PARTONIC CM ENERGY SMALLER ',/,
     &     10X,'THAN PLANCK MASS. MODEL MAY NOT BE ',/,
     &     10X,'VALID IN THIS REGION. ',/,
     &     10X,'RUNS CONTINUES AT USER`S RISK',/)
      END

CDECK  ID>, UPEVNT.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Chris Harris and Peter Richardson
C----------------------------------------------------------------------
      SUBROUTINE UPEVNT
C----------------------------------------------------------------------
C     Les Houches event routine
C----------------------------------------------------------------------
      IMPLICIT NONE
      CALL CHEVNT(.TRUE.)
      END

CDECK  ID>, CHDATA.
*CMZ :-        -21/08/06  10.33.53  by  Peter Richardson
*-- Author :    
C-----------------------------------------------------------------------
      BLOCK DATA CHDATA
C-----------------------------------------------------------------------
C    BLOCK DATA TO SET INITIAL VALUES OF PARAMETERS
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
C--Les Houches run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &                IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &                XMAXUP(MAXPUP),LPRUP(MAXPUP)
C--event common block
C--Set MPLNCK and define what is meant by it:
C--MSSDEF=1 means M_GT, MSSDEF=2 means M_DL and MSSDEF=3 means M_D
      DATA MPLNCK /1000.0D0/
      DATA MSSDEF /3/
C--Set number of dimensions (number of EXTRA dimensions is n=TOTDIM-4)
C--TOTDIM can be in the range 6-11
      DATA TOTDIM/6/
C--Use Giddings+Thomas momentum scale for calling pdfs?
C--(See page 12 of hep-ph/0106219)
      DATA GTSCA /.FALSE./
C--Set mass window for black holes produced
      DATA  MINMSS /5000.0D0/
      DATA  MAXMSS /14000.0D0/
C--Set NBODY decay of BH remnant - NBODY can be in range 2-5.
      DATA NBODY /2/
C--Turn time variation of BH Hawking temperature in decay on or off
      DATA TIMVAR /.TRUE./
C--Set which decay products are to be used:
C--MSSDEC=1 gives no heavy particles
C--MSSDEC=2 gives t, W and Z as well
C--MSSDEC=3 gives t, W, Z and Higgs as well  
#if defined(GINGRICH)
C--MSSDEC=4 gives graviton in D-dimensions as well
#endif
      DATA MSSDEC /3/
C--Turn grey-body factors on/off
      DATA GRYBDY /.TRUE./
C--Turn kinematic cut-off of decay on (rather than M=MPLANCK cut-off)
      DATA KINCUT /.FALSE./
C--Use Yoshino-Rychkov factors in cross-section
      DATA YRCSEC /.TRUE./
C--Max Hawking temperature
      DATA THWMAX /1000.0D0/
C--Remnant decay model: boiling at THWMAX
      DATA RMBOIL /.FALSE./
C--Min BH mass, ends boiling
      DATA RMMINM /100.D0/
C--New parameters for version 1.004
C--Include effects of BH spin
      DATA BHSPIN /.TRUE./
C--Vary BH spin axis during decay (if BHSPIN)
      DATA BHJVAR /.TRUE./
C--Include anisotropy of decay products (if BHSPIN)
      DATA BHANIS /.TRUE./
C--New parameters for version 2.001
C--Switches to control lepton numbers conservation
C--/Total lepton number conservation, Electron lepton number conservation, Tau lepton number conservation, Mu lepton number conservation/
      DATA NLEPTONCSV/.FALSE.,.FALSE.,.FALSE.,.FALSE./
C--Let BH lose mass and angular momentum in production
      DATA MJLOST /.TRUE./
C--Use MINMSS for the minimum mass of the black holes generated 
      DATA USEMINMSSBH /.TRUE./
C--Mean mass fraction lost: <MLOST> = FMLOST*(ECM-MBH) (power law)
      DATA FMLOST /0.99D0/ 
C--Bias towards constant horizon angular velocity lines in the allowed region of the M,J plane
      DATA CVBIAS /.FALSE./
C--Stable remnant model
      DATA RMSTAB /.FALSE./
C--RECOIL option 1 (in the rest frame of initial BH)
      DATA RECOIL /2/
C--Use the estimate for <N> Nto terminate evaporation. Default is NBODYAVERAGE=.TRUE.
      DATA NBODYAVERAGE /.TRUE./
C--NBODYVAR ON
      DATA NBODYVAR /.FALSE./
C--NBODYPHASE OFF
      DATA NBODYPHASE /.FALSE./
C--Skip evaporation OFF
      DATA SKIP2REMNANT /.FALSE./
C--print out option (0=no printout, 1 =errors only, 2=errors+info)
      DATA IBHPRN /1/
C--Default les houches event file name and histories file name (must be 8 characters long)
      DATA LHEFILENAME /'lhouches'/
      DATA HISFILENAME /'histfile'/
C--Switch to decide whether to put the BH in the les houches file or not
      DATA BHLHOUCHES /.FALSE./
C--pdf's for the beams (use the generator default)
C--MUST BE THE SAME FOR BOTH BEAMS
#if defined(LHAPDF)
      DATA PDFGUP /0,0/
      DATA PDFSUP /10000,10000/
#elif defined(PDFLIB)
      DATA PDFGUP /3,3/
      DATA PDFSUP /89,89/
#else
      DATA PDFGUP /-1,-1/
      DATA PDFSUP /-1,-1/
#endif
      DATA LHAPDFSET /10000/
C--NMAXEV the maximum number of bh events to generate
      DATA CHNMAXEV /100/
C--Beam particles and energies (only proton and/or antiproton)
      DATA IDBMUP /2212,2212/
      DATA EBMUP /7000.0D0,7000.0D0/
#if defined(GINGRICH)
C--Brane tension and escape flag.
      DATA DGTENSION /1.0D3/
      DATA DGESCAPE /.FALSE./
C--String ball parameters
      DATA DGMS /1000.0D0/
      DATA DGGS /0.3/
      DATA DGSB /.FALSE./
#endif
      END

#if defined(CHARYBDIS)
CDECK  ID>, CHDATAEXTRA.
*CMZ :-        -21/08/06  10.33.53  by  Peter Richardson
*-- Author :    
C-----------------------------------------------------------------------
      BLOCK DATA CHDATAEXTRA
C-----------------------------------------------------------------------
C    BLOCK DATA TO SET INITIAL VALUES OF PARAMETERS IN STANDALONE MODE
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Seed for the random number generator (taken from herwig)
      INTEGER ISEED
      COMMON/HWSEED/ISEED(2)
      INTEGER NRN
      COMMON/CHRANDOM/NRN(2)
c--Set default seeds
      DATA ISEED/12345,67890/
      DATA NRN /245234,42542/
      END
#endif

CDECK  ID>, CHEVNT.
*CMZ :-        -24/01/08  16.02.10  by Marco Sampaio
*--   Author :    Chris Harris & Bryan Webber (for spin) 
*                              & Marco Sampaio (for spin)
C----------------------------------------------------------------------
      SUBROUTINE CHEVNT(GENEVT)
C----------------------------------------------------------------------
C     Black Hole Production and Decay
C----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
C--Les Houches run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &                IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &                XMAXUP(MAXPUP),LPRUP(MAXPUP)
C--event common block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
C--Local Variables
      DOUBLE PRECISION ECMS,ONE,TWO,PIFAC,THREE,ZERO,
     &     CHRGEN,RPOW,A0,A1,EMJ,FACT,QSQ,HCS,PROB,W1,W2,W3,W4,
     &     TEMPFACTR,FACTR,RCS,BHPOW,DISF1,DISF2,OMEGAD,DL2GT,PSTEST,
     $     CFACTOR(13,13),
#if !defined(GINGRICH)
     &     RHORSQ,PLPOW,MPFACT,PA(5),PB(5),PC(5),PCM,CHUPCM,PCMF(5),
#else
     &     RHORSQ,PLPOW,MPFACT,PA(12),PB(12),PC(5),PCM,CHUPCM,PCMF(5),
#endif
     &     XXMIN,XLMIN,EMSCA,XX(2),DISF(13,2),GEV2PB,CHFMAS,YRF(7),
     &     MBH,JBH,RH,OBLA,OMEG,THAW,CYR(0:7),TNOW,MSEM,PSEM,
     &     EEM,JEM,MEM,JNEW,MNEW,FREC,RTR,SSQ,SBH(3),RBH(3),
     &     ROT(3,3),PEM(5),QBH(4),FMAXEM,BHMFR,BHJFR,BFR,MBHFINALSQ,
     &     NAVERAGE,NFACTOR,TOTFLUX,JEMNBODY(21),MEMNBODY(21),
     &     JINTRINSIC,MINTRINSIC,BHJ,ASTARPREVIOUS,BSTNEW,BSTNOW
      INTEGER IBARYN,IQRK(MAXPUP),IGLU(MAXPUP),IQBR(MAXPUP),MTRY,
     &     NQRK,NGLU,NQBR,JQRK,JGLU,JQBR,J,ISTART,CHST,ILINE,
     &     IP,IQ,PT,I,BHCHRG,IDN(2),NTRY,NHTRY,SAVCHR,SAVBRN,ID(21),
     &     MHTRY,LTRY,LHTRY,CHFCHG,SPIN,FSPIN(21),NEX,NBODYNOW,
     &     STATNBODY,POLARIZ,FPOLARIZ(21),ILEPTON(0:3),SAVLEPTON(0:3)
      LOGICAL GENEVT,CHRLOG,FIRST,OKDEC,NBODYCONVERGENCE,LOWTH,HVYQONLY
      EXTERNAL CHFMAS,CHRGEN,CHUPCM,CHRLOG,CHFCHG,BSTNEW,BSTNOW
      PARAMETER(ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,ZERO=0.0D0,NHTRY=200)
      PARAMETER(MHTRY=20,LHTRY=200,GEV2PB=389379.D3)
      SAVE HCS,BHPOW,PLPOW,RPOW,A0,A1,FACT,MPFACT,PIFAC,NEX
      DATA FIRST/.TRUE./
      DATA HVYQONLY /.TRUE./
#if defined(GINGRICH)
      DOUBLE PRECISION THRESH,BRECOIL
      INTEGER BHBARY,CHST2
#endif
C--Yoshino-Rychkov factors for c-s and Jmax
C  from Phys.Rev.D71:104028,2005 [hep-th/0503171]
      DATA YRF/1.54,2.15,2.52,2.77,2.95,3.09,3.20/
      DATA CYR/0.42d0,0.62d0,0.73d0,0.79d0,0.83d0,0.86d0,0.88d0,0.89d0/

C-----------------------------------------------------------------------
C                             INITIALISATION
C-----------------------------------------------------------------------

C--Set FMAXEM, which controls whether or not to use the full spectrum without kinematical restrictions, for the evaporation step 
      IF(KINCUT)THEN
         FMAXEM=1D3
      ELSE
         FMAXEM=1D0
      END IF
C--Initial calculatioin of parameters and convention fixing
      IF(FIRST) THEN
         PIFAC = ACOS(-ONE)
         NEX=TOTDIM-4
C--In what follows it will be assumed the Giddings Thomas convention so the relevant conversions are made here if needed so that the internal planck mass (INTMPPL) is in this convention.
C--Convert MPLNCK to M_GT if above value is M_DL
         IF (MSSDEF.EQ.2) THEN
            DL2GT=(TWO**(TOTDIM-6)*PIFAC**(TOTDIM-5))**(1/(TOTDIM-TWO))
            INTMPL=DL2GT*MPLNCK
C--Convert MPLNCK to M_GT if above value is M_D
         ELSEIF (MSSDEF.EQ.3) THEN
            INTMPL=MPLNCK*TWO**(1/(TOTDIM-TWO))
         ELSE
            INTMPL=MPLNCK
         ENDIF
C--Calculate other parameters
C-------Some powers appearing in the expressions for the Horizon radius
         PLPOW=TWO/(TOTDIM-THREE)
         BHPOW=PLPOW-7.0D0
         RPOW=1/BHPOW
         IF (TOTDIM.EQ.11) THEN
            OMEGAD=(PIFAC**5)/12
         ELSEIF (TOTDIM.EQ.10) THEN
            OMEGAD=(32*PIFAC**4)/105
         ELSEIF (TOTDIM.EQ.9) THEN
            OMEGAD=(PIFAC**4)/3
         ELSEIF (TOTDIM.EQ.8) THEN
            OMEGAD=(16*PIFAC**3)/15
         ELSEIF (TOTDIM.EQ.7) THEN
            OMEGAD=PIFAC**3
         ELSEIF(TOTDIM.EQ.6) THEN
            OMEGAD=(8*PIFAC**2)/3
         ELSEIF(TOTDIM.EQ.5) THEN
            OMEGAD=2*PIFAC**2
         ENDIF
         MPFACT=INTMPL**(TOTDIM-2)
         RHFACT=(4*(2*PIFAC)**(TOTDIM-4))/((TOTDIM-2)*OMEGAD*MPFACT)
         FIRST = .FALSE.
      ENDIF

C-----------------------------------------------------------------------
C                              PRODUCTION
C-----------------------------------------------------------------------

C--Select parton-parton invariant mass
C---Centre of mass energy
      ECMS = TWO*SQRT(EBMUP(1)*EBMUP(2))
      A0=MINMSS**BHPOW
      A1=(MAXMSS**BHPOW-A0)
C---Select scale for the emission
      EMSCA=(A0+A1*CHRGEN(1))**RPOW
      BHMASS=EMSCA
      QSQ=EMSCA**2
      EMJ=EMSCA**(1-BHPOW)/BHPOW*A1

C--Calculate corresponding (r_h)**2
      RHORSQ=(RHFACT*BHMASS)**PLPOW
#if defined(GINGRICH)
      IF (DGSB) THEN
        IF (BHMASS .GE. DGMS/(DGGS*DGGS)) THEN
C---Black hole cross section
          RHORSQ=(RHFACT*BHMASS)**PLPOW
        ELSEIF (BHMASS .LT. DGMS/DGGS) THEN
C---String ball cross section
          RHORSQ = (DGGS*BHMASS/DGMS**2)**2 / PIFAC
        ELSE
C---Unitarity limted string ball cross section
          RHORSQ = 1.0 / (PIFAC*DGMS**2)
        ENDIF
      ENDIF
#endif
C--Select initial momentum fractions (i.e. limits in the integration over pdfs)
      XXMIN=QSQ/ECMS**2
      XLMIN=LOG(XXMIN)
      FACT=-GEV2PB*PIFAC*RHORSQ*EMJ*XLMIN*TWO/EMSCA
C--Include Yoshino-Rychkov factor & spin
      IF (YRCSEC) THEN
        FACT=FACT*YRF(TOTDIM-4)




C-----MODIFICATION----------------------------------------------------------------
      ELSE 
        FACT=FACT*4D0*(1D0/(1D0+((TOTDIM-2)/2D0)**2))**(2D0/(TOTDIM-3))
C--------------------------------------------------------------------------------





      ENDIF
C--Save C-factor values---------
      IF (.NOT.HVYQONLY) THEN
         CFACTOR(1,1)=8.14536D-5    
         CFACTOR(1,2)=8.06258D-5      
         CFACTOR(1,3)=0.000393687      
         CFACTOR(1,4)=0.000391965     
         CFACTOR(1,5)=0.000807906     
         CFACTOR(1,13)=0.00107013      
         CFACTOR(2,2)=7.98065D-5     
         CFACTOR(2,3)=0.000528489     
         CFACTOR(2,4)=0.000665558     
         CFACTOR(2,5)=0.00387594     
         CFACTOR(2,13)=0.00105928      
         CFACTOR(3,3)=0.001898     
         CFACTOR(3,4)=0.00256642     
         CFACTOR(3,5)=0.0038955     
         CFACTOR(3,13)=0.0051701      
         CFACTOR(4,4)=0.0034703     
         CFACTOR(4,5)=0.0193658
         CFACTOR(4,13)=0.00699052
         CFACTOR(5,5)=0.0138953
         CFACTOR(5,13)=0.0139894
         CFACTOR(13,13)=0.0140845
      ELSEIF (HVYQONLY) THEN
         CFACTOR(1,1)=0.    
         CFACTOR(1,2)=0.    
         CFACTOR(1,3)=0.     
         CFACTOR(1,4)=0.     
         CFACTOR(1,5)=0.    
         CFACTOR(1,13)=0.     
         CFACTOR(2,2)=0.     
         CFACTOR(2,3)=0.     
         CFACTOR(2,4)=0.     
         CFACTOR(2,5)=0.     
         CFACTOR(2,13)=0.      
         CFACTOR(3,3)=0.     
         CFACTOR(3,4)=0.     
         CFACTOR(3,5)=0.     
         CFACTOR(3,13)=0.      
         CFACTOR(4,4)=0.     
         CFACTOR(4,5)=0.
         CFACTOR(4,13)=0.
         CFACTOR(5,5)=0.0138953
         CFACTOR(5,13)=0.0139894
         CFACTOR(13,13)=0.0140845
      ENDIF
      DO IP=1,13
         DO IQ=1,13
            IF (IP.GT.IQ) THEN
               CFACTOR(IP,IQ)=CFACTOR(IQ,IP)
            ENDIF
            IF (IP.EQ.6.OR.IQ.EQ.6) THEN
               CFACTOR(IP,IQ)=0.
            ENDIF
         ENDDO
      ENDDO
      DO IP=1,13
         DO IQ=1,13
            IF(IP.GT.6.AND.IP.LT.13.AND.IQ.LE.6) THEN
               CFACTOR(IP,IQ)=CFACTOR(IP-6,IQ)
            ELSEIF(IP.LE.6.AND.IQ.GT.6.AND.IQ.LT.13) THEN
               CFACTOR(IP,IQ)=CFACTOR(IP,IQ-6)
            ELSEIF(IP.GT.6.AND.IP.LT.13.AND.IQ.GT.6.AND.IQ.LT.13) THEN
               CFACTOR(IP,IQ)=CFACTOR(IP-6,IQ-6)
            ELSEIF(IP.EQ.13.AND.IQ.GT.6.AND.IQ.LT.13) THEN
               CFACTOR(IP,IQ)=CFACTOR(IP,IQ-6)
            ELSEIF(IQ.EQ.13.AND.IP.GT.6.AND.IP.LT.13) THEN
               CFACTOR(IP,IQ)=CFACTOR(IP-6,IQ)
            ENDIF
         ENDDO
      ENDDO
C      DO IP=1,13
C         DO IQ=1,13
C            WRITE(*,*) IP, IQ, CFACTOR(IP,IQ)
C         ENDDO
C      ENDDO

C--Change momentum scale for calling pdfs (see comment above)
      IF (GTSCA) EMSCA=1/SQRT(RHORSQ)
      SCALUP = EMSCA
      CALL CHPDF(XXMIN,XLMIN,EMSCA,XX,DISF)
      HCS=ZERO
      FACTR=ZERO
      DO IP=1,13
         DO IQ=1,13
            TEMPFACTR=FACT*DISF(IP,1)*DISF(IQ,2)*CFACTOR(IP,IQ)
            FACTR=FACTR+TEMPFACTR
         ENDDO
      ENDDO

C-----MODIFICATION------------------------------------------------------
C      W1=FACT*DISF(13,1)*DISF(13,2)*0.01408
C      W2=W1+FACT*DISF(13,1)*DISF(5,2)*(0.02021/2.+0.00419/2.)
C      W3=W2+FACT*DISF(5,1)*DISF(13,2)*(0.02021/2.+0.00419/2.)
C      W4=W3+FACT*DISF(5,1)*DISF(5,2)*(0.00133/4.+0.03082/4.+0.00639/2.)
C      FACTR=W4
C
      HCS=FACTR
      XWGTUP = HCS
      IF(.NOT.GENEVT)THEN
c-----This block determines if a black hole with a mass above MINMSS would be generated for this weight and sets it to zero if not (when MJLOST=.TRUE. and USEMINMSSBH=.TRUE.) 
         IF(MJLOST.AND.USEMINMSSBH)THEN
            RCS=HCS*CHRGEN(1)
            HCS = ZERO
C            IF(RCS.GT.0.AND.RCS.LE.W1) THEN
C               IDN(1)=13
C               IDN(2)=13
C            ELSEIF(RCS.GT.W1.AND.RCS.LE.W2) THEN
C               IDN(1)=13
C               IDN(2)=5
C            ELSEIF(RCS.GT.W2.AND.RCS.LE.W3) THEN
C               IDN(1)=5
C               IDN(2)=13
C            ELSEIF(RCS.GT.W3.AND.RCS.LE.W4) THEN
C               IDN(1)=5
C               IDN(2)=5
C            ENDIF
C----------------------------------------------------------------------------------------------------
            DO 201 IP=1,13
               DO 101 IQ=1,13
                  FACTR=FACT*DISF(IP,1)*DISF(IQ,2)*CFACTOR(IP,IQ)
                  HCS=HCS+FACTR
                  IF (HCS.GT.RCS) GOTO 991
 101           CONTINUE
 201        CONTINUE
C--Select particle types just to use in the mass reduction routine to determine whether to accept or reject the weight.
 991        IDN(1)=IP
            IDN(2)=IQ
C--Use the les houches common block temporarily
            NUP = 2
            DO I=1,2
               IF(IDN(I).GT.6.AND.IDN(I).LE.12) THEN
                  IDN(I) = -MOD(IDN(I)-1,6)-1
               ELSEIF(IDN(I).EQ.13) THEN
                  IDN(I) = 21
               ENDIF
               IDUP(I) = IDN(I)
               PUP(5,I) = CHFMAS(IDUP(I))
               PUP(1,I) = ZERO
               PUP(2,I) = ZERO
               PUP(3,I) = XX(I)*EBMUP(I)
               PUP(4,I) = SQRT(PUP(3,I)**2+PUP(5,I)**2)
               ISTUP(I) = -1
               SPINUP(I) = 9D0
C               IF(PUP(4,I).GT.EBMUP(I))THEN
C                  WRITE(*,*)'CHEVNT WARNING:'
C                  WRITE(*,*)'   Bad momentum fraction. Trying again!'
C                  WRITE(*,*)'   If this error repeats too many times'//
C     &' please check your beam energies and MINMSS MAXMSS ranges.'
C                  GOTO 31
C               END IF
            ENDDO
            PUP(3,2) = -PUP(3,2)
C--   PCMF is the total momentum of the collision 
            CALL CHVSUM(4,PUP(1,1),PUP(1,2),PCMF)
            CALL CHUMAS(PCMF)
C--   Choose BH mass (and J) and compare with the MINMSS
            CALL CHYRMJLS(TOTDIM,BHJFR,BHMFR,BFR,CVBIAS,MJLOST)
            BHMASS=BHMFR*PCMF(5)
            IF(BHMASS.LT.MINMSS) XWGTUP=0D0
         ENDIF
C--------------end of decision block for the weight (see also below the rejection when generating events -- events below MINMSS must be discarded)
      ENDIF
      RCS=HCS*CHRGEN(1)
      HCS = ZERO
      DO IP=1,13
         DO IQ=1,13
            TEMPFACTR=FACT*DISF(IP,1)*DISF(IQ,2)*CFACTOR(IP,IQ)
            FACTR=FACTR+TEMPFACTR
         ENDDO
      ENDDO



C-----MODIFICATION-------------------------------------------------------------------------------
C      IF(RCS.GT.0.AND.RCS.LE.W1) THEN
C         IDN(1)=13
C         IDN(2)=13
C      ELSEIF(RCS.GT.W1.AND.RCS.LE.W2) THEN
C         IDN(1)=13
C         IDN(2)=5
C      ELSEIF(RCS.GT.W2.AND.RCS.LE.W3) THEN
C         IDN(1)=5
C         IDN(2)=13
C      ELSEIF(RCS.GT.W3.AND.RCS.LE.W4) THEN
C         IDN(1)=5
C         IDN(2)=5
C      ENDIF
C-----------------------------------------------------------------------------------------------







      DO 20 IP=1,13
      DO 10 IQ=1,13
      FACTR=FACT*DISF(IP,1)*DISF(IQ,2)*CFACTOR(IP,IQ)
      HCS=HCS+FACTR
      IF (HCS.GT.RCS) GOTO 99
 10   CONTINUE
 20   CONTINUE
C--Generate event
 99   IDN(1)=IP
      IDN(2)=IQ
C--put the incoming particles in the Les Houches common block
      NUP = 2
      DO I=1,2
         IF(IDN(I).GT.6.AND.IDN(I).LE.12) THEN
            IDN(I) = -MOD(IDN(I)-1,6)-1
         ELSEIF(IDN(I).EQ.13) THEN
            IDN(I) = 21
         ENDIF
         IDUP(I) = IDN(I)
         PUP(5,I) = CHFMAS(IDUP(I))
         PUP(1,I) = ZERO
         PUP(2,I) = ZERO
         PUP(3,I) = XX(I)*EBMUP(I)
         PUP(4,I) = SQRT(PUP(3,I)**2+PUP(5,I)**2)
         ISTUP(I) = -1
         SPINUP(I) = 9D0
C         IF(PUP(4,I).GT.EBMUP(I))THEN
C            WRITE(*,*)'CHEVNT WARNING:'
C            WRITE(*,*)'   Bad momentum fraction. Trying again!'
C            WRITE(*,*)'   If this error repeats too many times please'//
C     &' check your beam energies and MINMSS MAXMSS ranges.'
C            GOTO 31
C         END IF
      ENDDO
      PUP(3,2) = -PUP(3,2)
C--PCMF is the total momentum of the collision 
      CALL CHVSUM(4,PUP(1,1),PUP(1,2),PCMF)
      CALL CHUMAS(PCMF)
C--Calculate BH charge (actually 3*charge) assuming none lost
      BHCHRG=CHFCHG(IDN(1))+CHFCHG(IDN(2))
#if defined (GINGRICH)
C--Calculate BH baryon number
      BHBARY = 0
      DO I=1,2
        IF(ABS(IDN(I)).LE.6) BHBARY = BHBARY + SIGN(1,IDN(I))
      ENDDO
C--Initialize escape flag
      DGESCAPE = .FALSE.
#endif

C--Black Hole decay
      MTRY = 0
      CHST = BHCHRG
#if defined(GINGRICH)
      CHST2 = BHBARY
#endif

C--Choose initial BH mass and spin
      CALL CHYRMJLS(TOTDIM,BHJFR,BHMFR,BFR,CVBIAS,MJLOST)
      BHMASS=BHMFR*PCMF(5)
c----- Reject the event if the mass of the black hole is below the MINMSS threshold (when USEMINMSSBH=.TRUE.)
      IF((BHMASS.LT.MINMSS).AND.USEMINMSSBH) THEN
         XWGTUP=0D0
         RETURN
      END IF
c--- Start attepmts to evaporate
 35   MTRY = MTRY+1
      BHCHRG = CHST
#if defined(GINGRICH)
      BHBARY = CHST2
#endif
      NBHD=0
      NUP = 2
C---- Compute again the initial BHMASS in case the evaporation has failed and is being repeated from 35.
      BHMASS=BHMFR*PCMF(5)

c----- Select the spin
      SBH(3)=ZERO
      IF (BHSPIN) THEN
C--J corresponds to reduced mass at same impact parameter
         JBH=NINT(BHJFR*CYR(NEX)*PCMF(5)*SQRT(RHORSQ)*BFR)
C--Now combine this with intrinsic spin of colliding particles
C-----First select the helicities of the incoming particles and the wigner rotated (S,S_perpendicular) state. The probabilitites for this state are obtained from applying a wigner rotation of the state (S,S_z) obtained from combining the helicities of the incoming particles.
         CALL CHINTRINSICSPIN(IDN,JINTRINSIC,MINTRINSIC)
C-----Combine the previous wave with a (J,J_perpendicular)=(JBH,JBH) state
         CALL CHNEWJ(JBH,JINTRINSIC,-MINTRINSIC,JNEW,MNEW)
C-----Find the (tilted) total angular momentum probabilistically
         IF (JNEW-ABS(MNEW).LT.0.1D0) THEN
C----------Fixed orientation
            SBH(1)=MNEW
            SBH(2)=ZERO
            SBH(3)=ZERO
         ELSEIF (ABS(MNEW).GT.JNEW) THEN
            PRINT *,' CHEVNT: BAD MNEW,JNEW=',MNEW,JNEW
         ELSE
C----------Variable orientation
            RTR=SQRT(JNEW**2-MNEW**2)
            CALL CHRAZM(RTR,SBH(2),SBH(3))
            SBH(1)=MNEW
         ENDIF
         IF (JNEW.LT.0.1D0) THEN
C----------Choose Z axis at random
 461        SBH(1)=CHRGEN(1)
            SBH(2)=CHRGEN(2)
            SBH(3)=CHRGEN(3)
            SSQ=SBH(1)**2+SBH(2)**2+SBH(3)**2
            IF (SSQ.EQ.ZERO.OR.SSQ.GT.ONE) GOTO 461
         ENDIF
C----Reset JBH
         JBH=JNEW
C-----Finally rotate the total angular momentum vector obtained around the beam axis randomly to account for the axial symmetry.
         RTR=SQRT(SBH(1)**2+SBH(2)**2)
         CALL CHRAZM(RTR,SBH(1),SBH(2))
      ELSE
         JBH=ZERO
         SBH(1)=ZERO
         SBH(2)=ZERO
      ENDIF
C--Find initial BH 4-momentum and put
C  missing 4-mom ('graviton') in event record
      PA(5)=BHMASS    
      IF(MJLOST) THEN
         NUP=NUP+1
C--Generate mass lost power distributed with mean fraction FMLOST
         PUP(5,NUP)=(PCMF(5)-BHMASS)*CHRGEN(1)**(ONE/FMLOST-ONE) 
         PCM=CHUPCM(PCMF(5),PA(5),PUP(5,NUP))
         CALL CHDTWO(PCMF,PA,PUP(1,NUP),PCM,TWO,.TRUE.)
         IDUP(NUP)=39
         ISTUP(NUP)=1      
         MOTHUP(1,NUP) = 1
         MOTHUP(2,NUP) = 2
      ELSE
         CALL CHVEQU(5,PCMF,PA)
      END IF
#if defined(GINGRICH)
      DO I=6,12
        PA(I) = 0.0D0
      ENDDO
      DO I=1,12
        PB(I) = 0.0D0
      ENDDO
#endif

C-----------------------------------------------------------------------
C                           EVAPORATION PHASE
C-----------------------------------------------------------------------
      BHJ = JBH
C--Try again if parameters chosen are unphysical
      IF((NEX.EQ.1).AND.(BSTNOW(BHMASS,BHJ,INTMPL,NEX).GT.1)) THEN
         GOTO 35
      END IF
C--Compute BH parameters
      CALL CHPARA(BHMASS,BHJ,INTMPL,NEX,RH,OBLA,OMEG,TNOW)
C--Initialise the switch which controls if the Hawking temperature is too large
      LOWTH=.TRUE.
C---Start evaporation
 40   NUP = NUP+1
C--Set BHMASS to correct value (necessary in time varying case)
      IF (TIMVAR.AND.LOWTH) THEN 
         BHMASS=PA(5)
         BHJ=JBH
         ASTARPREVIOUS = OBLA/RH
      END IF
C--Compute BH parameters
      CALL CHPARA(BHMASS,BHJ,INTMPL,NEX,RH,OBLA,OMEG,TNOW)
C--Update temperature
      IF (TIMVAR.OR.NBHD.EQ.0) THAW=TNOW
C-----Freeze parameters if the computed value for the temperature is too large
      IF (THAW.GT.THWMAX) LOWTH=.FALSE.
      IF (.NOT.LOWTH) THEN
         THAW=THWMAX
         RH=(NEX+1+(NEX-1)*ASTARPREVIOUS**2)/
     &((4*PIFAC)*(1+ASTARPREVIOUS**2)*THAW)
         OBLA = ASTARPREVIOUS*RH
         OMEG = OBLA/(RH**2+OBLA**2)
      END IF
C--Select emission and compute estimate for the remaninig number of emissions
#if !defined(GINGRICH)
      IF(BHSPIN) THEN 
         CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX)
      ELSE
         CALL CHSELECT(PT,SPIN,BHCHRG,ZERO,RH,NEX,NFACTOR,TOTFLUX)
      END IF
#else
      IF(BHSPIN) THEN 
         CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
      ELSE
         CALL CHSELECT(PT,SPIN,BHCHRG,ZERO,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
      END IF
#endif
c--Estimate for the average number of particles to be emitted
      NAVERAGE=NFACTOR*RH*BHMASS
      MBH=PA(5)
C--Increase counter of number of emissions
      NBHD=NBHD+1
C--Save BH state into the history
      IF (NBHD.LE.MBHD) THEN
         CALL CHVEQU(5,PA,PBHD(1,NBHD))
#if defined(GINGRICH)
         CALL CHVEQU(NEX,PA(6),DGPSBD(1,NBHD))
#endif
         IF(JBH.LT.0.1D0)THEN
            SBHD(1,NBHD)=0D0
            SBHD(2,NBHD)=0D0
            SBHD(3,NBHD)=0D0
         ELSE
            CALL CHVEQU(3,SBH,SBHD(1,NBHD))
         END IF
         RHBHD(NBHD)=RH
         OBBHD(NBHD)=OBLA
         OMBHD(NBHD)=OMEG
         THBHD(NBHD)=THAW
         NFLUXBHD(NBHD)=TOTFLUX
      ELSE
         NBHD=NBHD-1
         WRITE(*,*)'CHEVNT: BH history truncated at NBHD=MBHD=',MBHD
         WRITE(*,*)'----- If you are running high multiplicity events'//
     &' you need to increase the parameter MBHD in charybdis*.inc and'//
     &' the parameter MAXNUP in the LES HOUCHES common block'//
     &' declarations in all *.f *.F files (including in the source'//
     &' file of any general purpose event generator you are'//
     &' interfacing to). It is recomended to set'//
     &' MAXNUP to the value of of at least MBHD+100, to allow for'//
     &' remnant particles.'
         WRITE(*,*)'NOTE: If you are interfacing CHARYBDIS to a general'
     &//' purpose event generator such as HERWIG or PYTHIA, you may'// 
     &' have to worry about array sizes for the internal event record'//
     &' (HEPEVT common block) of the general purpose event generator.'//
     &' For HERWIG this means that you will have to increase the'//
     &' value of the parameter NMXHEP in the herwig*.inc file. For'//
     &' PYTHIA you will have to replace all the PARAMETER (NMXHEP=...)'
     &//' statements in the source code.' 
         WRITE(*,*)' '
         STOP
      ENDIF
c-----Skip evaporation (when using nbodyaverage as criterion and naverage small enough), (or when KINCUT=.FALSE. and mass is low)
      IF((NBODYAVERAGE.AND.(NAVERAGE.LE.(NBODY-1))).OR.
     &((.NOT.KINCUT).AND.(RMBOIL).AND.(PA(5).LT.RMMINM)).OR.
     &((.NOT.KINCUT).AND.(.NOT.RMBOIL).AND.(PA(5).LT.INTMPL)).OR.
     &SKIP2REMNANT) GOTO 50
C--Select the decay product
      LTRY=0
      IF (BHSPIN) THEN
 45      MBH=PA(5)
         MSEM=CHFMAS(PT)
C--Check if there exists a non empty kinematically allowed region
         IF (MSEM.GT.MBH) THEN
            IF (LTRY.GT.LHTRY) THEN 
               WRITE(*,*)'CHEVNT WARNING: Too many tries!'
               WRITE(*,*)'---Could not find  MSEM < MBH*FMAXEM'//
     &'. Skipping to remnant phase.'
               NAVERAGE=NFACTOR*RH*MBH
               GOTO 50
            ELSE
               LTRY=LTRY+1
C--------------Try again. Reset BH charge and select a new emission.
               BHCHRG=BHCHRG+CHFCHG(PT)
#if defined (GINGRICH)
               IF(ABS(PT).LE.6) BHBARY=BHBARY+SIGN(1,PT)
               CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
#else
               CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX)
#endif
               GOTO 45
            ENDIF
         ENDIF
C--Generate a Hawking emission
         CALL CHEMIT(NEX,MBH,RH,OBLA,FMAXEM,MSEM,SPIN,EEM,JEM,MEM)
C--Stop the decay if the energy of the emitted object is too close to MBH or 
C  outside the kinematically allowed region (in the KINCUT=.TRUE. case)
         IF(RECOIL.EQ.1)THEN
            MBHFINALSQ=MBH**2-2*MBH*EEM+MSEM**2
         ELSE IF (RECOIL.EQ.2)THEN
            MBHFINALSQ=(MBH-EEM)*ABS(MBH-EEM)
         ELSE
            WRITE(*,*)'Invalid recoil option!'
            STOP
         END IF
C--Check if kinematicaly allowed
         IF (KINCUT.AND.MBHFINALSQ.LT.RMMINM**2) THEN
            NAVERAGE=NFACTOR*RH*MBH
C-----------Skip to remnant
            GOTO 50
         ELSE IF(.NOT.KINCUT.AND.(MBHFINALSQ.LT.0)) THEN
            IF (LTRY.GT.LHTRY) THEN 
               WRITE(*,*)'CHEVNT WARNING: Too many tries!'
               WRITE(*,*)'---Could not find kinematically allowed'//
     &' emission. Skipping to remnant phase.'
               NAVERAGE=NFACTOR*RH*MBH
               GOTO 50
            ELSE
C-----------Try again
               LTRY=LTRY+1
C--------------Try again. Reset BH charge and select a new emission.
               BHCHRG=BHCHRG+CHFCHG(PT)
#if defined (GINGRICH)
               IF(ABS(PT).LE.6) BHBARY=BHBARY+SIGN(1,PT)
               CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
#else
               CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX)
#endif
               GOTO 45
            ENDIF
         END IF
C--Generate momentum vector of the emission
         PSEM=SQRT(EEM**2-MSEM**2)
         IF (BHANIS) THEN
C----------Use spheroidal functions    
            CALL CHEMOM(PSEM,EEM,JEM,MEM,OBLA,PEM(1),SPIN,POLARIZ,PT)
         ELSE
C----------Use spherical harmonics
            CALL CHEMOM(PSEM,EEM,ZERO,ZERO,ZERO,PEM(1),SPIN,POLARIZ,PT)
         ENDIF
C--COMPUTE NEW ANGR MOM OF BH
         IF (BHJVAR.AND.THAW.LT.THWMAX) THEN
C----------Variable orientation (use Clebsch-Gordan coefficients)
            CALL CHNEWJ(JBH,JEM,MEM,JNEW,MNEW)
         ELSE
C----------Fixed orientation
            MNEW=JBH-MEM
            JNEW=ABS(MNEW)
         ENDIF
         IF((NEX.EQ.1).AND.(BSTNEW(MBH,JNEW,INTMPL,NEX,EEM,MSEM).GE.1)) 
     &THEN
            IF (LTRY.GT.LHTRY) THEN 
               WRITE(*,*)'CHEVNT WARNING: Too many tries!'
               WRITE(*,*)'---Could not find kinematically allowed'//
     &' emission. Skipping to remnant phase.'
               NAVERAGE=NFACTOR*RH*MBH
               GOTO 50
            ELSE
C-----------Try again
               LTRY=LTRY+1
C--------------Try again. Reset BH charge and select a new emission.
               BHCHRG=BHCHRG+CHFCHG(PT)
#if defined (GINGRICH)
               IF(ABS(PT).LE.6) BHBARY=BHBARY+SIGN(1,PT)
               CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
#else
               CALL CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX)
#endif
               GOTO 45
            ENDIF
         END IF

C--COMPUTE NEW ORIENTATION OF BH SPIN
         IF (JNEW-ABS(MNEW).LT.0.1D0) THEN
C----------Fixed orientation
            RBH(1)=ZERO
            RBH(2)=ZERO
            RBH(3)=MNEW
         ELSEIF (ABS(MNEW).GT.JNEW) THEN
            PRINT *,' CHEVNT: BAD MNEW,JNEW=',MNEW,JNEW
         ELSE
C----------Variable orientation
            RTR=SQRT(JNEW**2-MNEW**2)
            CALL CHRAZM(RTR,RBH(1),RBH(2))
            RBH(3)=MNEW
         ENDIF
         IF (JNEW.LT.0.1D0.AND.BHJVAR) THEN
C--CHOOSE Z AXIS AT RANDOM
 46         SBH(1)=CHRGEN(1)
            SBH(2)=CHRGEN(2)
            SBH(3)=CHRGEN(3)
            SSQ=SBH(1)**2+SBH(2)**2+SBH(3)**2
	    IF (SSQ.EQ.ZERO.OR.SSQ.GT.ONE) GOTO 46 
         END IF
         CALL CHUROT(SBH,ONE,ZERO,ROT)
         CALL CHUROB(ROT,PEM,PEM)
         CALL CHUROB(ROT,RBH,SBH)
C--COMPUTE RECOIL - Emission occurs in the rest frame of the initial BH for RECOIL=1. 
C  This should be valid since any relativistic wave (v=1) going outwards from the horizon
C  propagates on the initial background (any change in mass of the black hole can't propagate
C  ahead of the wave itself).
         PEM(5)=MSEM
         IF (RECOIL.EQ.1)THEN
C---------PEM is the momentum of the emission (in the rest frame of the initial BH)
            PEM(4)=EEM
C---------QBH is the final momentum of the BH (in the rest frane of the initial BH)
            QBH(4)=MBH-EEM
            DO IP=1,3
               QBH(IP)=-PEM(IP)
            ENDDO
C---------Boost to lab frame
            CALL CHULB4(PA,PEM,PC)
            CALL CHULB4(PA,QBH,PB)        
C---------Compute new BH mass and J
            MBH=SQRT(MBH**2-2*MBH*EEM+MSEM**2)
            JBH=JNEW
C---------Save masses
            PB(5)=MBH
            PC(5)=MSEM  
         ELSE IF (RECOIL.EQ.2)THEN
C---In this case the recoil is such that the Mass of the black hole is reduced by an amount which is the energy of the emission.
C---------PEM is the momentum of the emission (in the rest frame of the initial BH)
            PEM(4)=EEM-0.5D0*EEM**2/MBH+0.5D0*MSEM**2/MBH
C---------QBH is the final momentum of the BH (in the rest frane of the initial BH)
            QBH(4)=MBH-PEM(4)
            FREC=SQRT(PEM(4)**2-MSEM**2)/PSEM
            DO IP=1,3
               PEM(IP)=PEM(IP)*FREC
               QBH(IP)=-PEM(IP)
            ENDDO
C---------Boost to lab frame
            CALL CHULB4(PA,PEM,PC)
            CALL CHULB4(PA,QBH,PB)
C---------Compute new BH mass and J
            MBH=MBH-EEM
            JBH=JNEW
C---------Save masses
            PB(5)=MBH
            PC(5)=MSEM
         ELSE
            WRITE(*,*)'Invalid recoil option'
            STOP
         END IF
      ELSE
C--IF ROTATION OFF:
C--Obtain energy of emitted parton
 47      MBH=PA(5)
         PC(5)=CHFMAS(PT)
         IF(GRYBDY) THEN 
            CALL CHEMIT(NEX,MBH,RH,ZERO,FMAXEM,PC(5),SPIN,PC(4),JEM
     &,MEM)
         ELSE
            CALL CHNOGRYBDY(PC(4),SPIN,THAW*RH)
            PC(4)=PC(4)/RH
C--Check energy > mass for selected parton
            IF (PC(4).LT.PC(5)) THEN
               BHCHRG=BHCHRG+CHFCHG(PT)
#if defined (GINGRICH)
               IF(ABS(PT).LE.6) BHBARY=BHBARY+SIGN(1,PT)
#endif
               IF(IBHPRN.EQ.2) PRINT *,'CHEVNT: Not enough energy. Tryin
     &g again'
#if !defined(GINGRICH)
               CALL CHSELECT(PT,SPIN,BHCHRG,ZERO,RH,NEX,NFACTOR,TOTFLUX)
#else
               CALL CHSELECT(PT,SPIN,BHCHRG,ZERO,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
#endif
               GOTO 47
            ENDIF
         END IF
C--Check that emission is kinematically allowed
         IF (PC(4).GT.(PA(5)**2+PC(5)**2)/(TWO*PA(5))) THEN
            IF (KINCUT.OR.LTRY.GT.LHTRY) THEN 
               NAVERAGE=NFACTOR*RH*MBH
C-------------Skip to remnant
               GOTO 50
            ELSE
               LTRY=LTRY+1
C-------------Reset parameters and try again
               BHCHRG=BHCHRG+CHFCHG(PT)
#if defined (GINGRICH)
               IF(ABS(PT).LE.6) BHBARY=BHBARY+SIGN(1,PT)
               CALL CHSELECT(PT,SPIN,BHCHRG,ZERO,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
#else
               CALL CHSELECT(PT,SPIN,BHCHRG,ZERO,RH,NEX,NFACTOR,TOTFLUX)
#endif
               GOTO 47
            ENDIF
         ENDIF
C--Use 2-body phase space decay 
         PB(5)=SQRT(PA(5)**2+PC(5)**2-TWO*PA(5)*PC(4))
C--Calculate centre of mass momentum for decay
         PCM=SQRT(PC(4)**2-PC(5)**2)
C--PA is BH 5-momentum before, PB after decay
#if !defined(GINGRICH)
         CALL CHDTWO(PA,PB,PC,PCM,TWO,.TRUE.)
#else
         IF (PT .EQ. 39) THEN
           CALL CHDN(PA,PB,PC,PCM)
         ELSE
           CALL CHDTWO(PA,PB,PC,PCM,TWO,.TRUE.)
         ENDIF
#endif
      ENDIF
C--Add emitted particle to event record
      CALL CHVEQU(5,PC,PUP(1,NUP))
      IDUP(NUP)=PT
      ISTUP(NUP) = 1
      MOTHUP(1,NUP) = 1
      MOTHUP(2,NUP) = 2
C--Sort out polarizations
      IF(BHSPIN)THEN
         IF(SPIN.NE.0)THEN
            SPINUP(NUP)=POLARIZ
         ELSE IF(SPIN.EQ.0.AND.(PT.EQ.23.OR.ABS(PT).EQ.24))THEN
C--Longitudinal vector bosons
            SPINUP(NUP)=0D0          
         ELSE
            SPINUP(NUP)=9D0         
         END IF
      ELSE
         SPINUP(NUP)=9D0       
      END IF
C--Change BH 5-momentum back to PA
#if !defined(GINGRICH)
      CALL CHVEQU(5,PB,PA)
#else
      CALL CHVEQU(TOTDIM+1,PB,PA)
      IF (PT .EQ. 39) THEN
C--Energy threshold to leave brane (DMG).
        IF (MSSDEF .EQ. 2) THEN
C--DL definition
          THRESH = DGTENSION * MPLNCK 
     &           * (PA(5)/MPLNCK)**(3.0/(REAL(TOTDIM)-3.0))
        ELSEIF (MSSDEF .EQ. 3) THEN
C--PDG definition
          THRESH = DGTENSION * MPLNCK * ( (6.28**(TOTDIM-5)/2.0) 
     &           * (PA(5)/MPLNCK) )**(3.0/(REAL(TOTDIM)-3.0))
        ELSE
          WRITE(*,*) 'Invalid Planck scale definition.'
          STOP
        ENDIF
        THRESH = THRESH * (THRESH + 2.0*PA(5))
C--Square of momentum transverse to brane.
        BRECOIL = 0.0
        DO I=6,TOTDIM+1
          BRECOIL = BRECOIL + PA(I)**2
        ENDDO
        BRECOIL = BRECOIL / THRESH
        IF (BRECOIL .GT. 1.0) THEN
          DGESCAPE = .TRUE.
          GOTO 56
        ELSE
          DGESCAPE = .FALSE.
        ENDIF
      ELSE
        DGESCAPE = .FALSE.
      ENDIF
#endif
      GOTO 40

C-----------------------------------------------------------------------
C                             REMNANT DECAY
C-----------------------------------------------------------------------
 50   BHCHRG=BHCHRG+CHFCHG(PT)
#if defined(GINGRICH)
      IF (ABS(PT) .LE. 6) BHBARY = BHBARY + SIGN(1,PT)
#endif
      NUP = NUP-1
#if defined(GINGRICH)
 56   CONTINUE
#endif
C--Calculate the baryon number violation
C--Find the baryon number of the inital state
      IBARYN = 0
      DO I=1,2
        IF(ABS(IDUP(I)).LE.6) IBARYN = IBARYN-SIGN(1,IDUP(I))
      ENDDO
C--and the final state
      IF(NUP.GT.2)THEN
         DO I=3,NUP
            IF(ABS(IDUP(I)).LE.6) IBARYN = IBARYN+SIGN(1,IDUP(I))
         ENDDO
      END IF

C--Calculate the total lepton number violation
C--Find the lepton number of the inital state
      ILEPTON(0) = 0
      DO I=1,2
        IF(ABS(IDUP(I)).GE.11 .AND.ABS(IDUP(I)).LE.16) 
     &ILEPTON(0) = ILEPTON(0)-SIGN(1,IDUP(I))
      ENDDO
C--and the final state
      IF(NUP.GT.2)THEN
         DO I=3,NUP
            IF(ABS(IDUP(I)).GE.11 .AND.ABS(IDUP(I)).LE.16) 
     &ILEPTON(0) = ILEPTON(0)-SIGN(1,IDUP(I))
         ENDDO
      END IF  
C--Calculate the lepton number violation for each family
      DO J=0,2
         ILEPTON(J+1) = 0
C--Find the lepton number of the inital state
         DO I=1,2
            IF(ABS(IDUP(I)).GE.(11+2*J) .AND.ABS(IDUP(I)).LE.(12+2*J)) 
     &ILEPTON(J+1) = ILEPTON(J+1)-SIGN(1,IDUP(I))
         ENDDO
C--and the final state
         IF(NUP.GT.2)THEN
            DO I=3,NUP
            IF(ABS(IDUP(I)).GE.(11+2*J) .AND.ABS(IDUP(I)).LE.(12+2*J)) 
     &ILEPTON(J+1) = ILEPTON(J+1)-SIGN(1,IDUP(I))
            ENDDO
         END IF  
      END DO

      IF (RMSTAB) THEN
C--Stable remnant model -- only B=0 allowed
         IF (IBARYN.NE.0) GOTO 35
         NUP=NUP+1
         CALL CHVEQU(5,PA,PUP(1,NUP))
         IF (BHCHRG.EQ.0) THEN
            IDUP(NUP)= 50
         ELSEIF (BHCHRG.EQ.3) THEN
            IDUP(NUP)= 51
         ELSEIF (BHCHRG.EQ.-3) THEN
            IDUP(NUP)=-51
         ELSE
            PRINT *,'CHEVNT: BAD REMNANT CHARGE = ',BHCHRG
            GOTO 35
         ENDIF
         ISTUP(NUP)=2      
         MOTHUP(1,NUP) = 1
         MOTHUP(2,NUP) = 2
         NBODYNOW=0
         SPINUP(NUP)=9d0
         NUP=NUP-NBODY
         GOTO 110
      ENDIF
C--Remant decay --- Try to balance the baryon number
C-- N Variable, N-body decay
      IF(NBODYVAR) THEN
c---------Select NBODYNOW using Poisson distribution
         CALL CHVARNBODY(NAVERAGE,NBODYNOW)
      ELSE
         NBODYNOW=NBODY
      END IF
C---Save BH parameters
      NTRY = 0
      SAVCHR = BHCHRG
      SAVBRN = IBARYN
      SAVLEPTON(0)=ILEPTON(0)
      SAVLEPTON(1)=ILEPTON(1)
      SAVLEPTON(2)=ILEPTON(2)
      SAVLEPTON(3)=ILEPTON(3)

C---Try to find N particles
 60   NTRY = NTRY+1
      PCM = PA(5)
#if defined(GINGRICH)
C--Begin recoil.
      IF (DGESCAPE) THEN
C--BH has escaped from brane. Add parton(s) to event
C--record with 0 kinematics to conserve baryon number.
        NUP = NUP + 1
        PUP(1,NUP) = 0.0
        PUP(2,NUP) = 0.0
        PUP(3,NUP) = 0.0
        PUP(4,NUP) = 0.0
        PUP(5,NUP) = 0.0
        ISTUP(NUP) = 2
        MOTHUP(1,NUP) = 1
        MOTHUP(2,NUP) = 2
        IF (IBARYN .EQ. 0) THEN
          IDUP(NUP) = 21
        ELSEIF (IBARYN .EQ. 1) THEN
          IDUP(NUP) = -1
        BHCHRG = BHCHRG - CHFCHG(-1)
          BHBARY = BHBARY + IBARYN
          IBARYN = 0
        ELSEIF (IBARYN .EQ. -1) THEN
          IDUP(NUP) = 1
          BHCHRG = BHCHRG - CHFCHG(1)
          BHBARY = BHBARY + IBARYN
          IBARYN = 0
C--for |B|>1 
        ELSEIF (IBARYN .EQ. 2) THEN
          IDUP(NUP) = -1
          BHCHRG = BHCHRG - CHFCHG(-1)
          IBARYN = IBARYN - 1
          BHBARY = BHBARY + 1
          NUP = NUP + 1
          PUP(1,NUP) = 0.0
          PUP(2,NUP) = 0.0
          PUP(3,NUP) = 0.0
          PUP(4,NUP) = 0.0
          PUP(5,NUP) = 0.0
          IDUP(NUP) = -1
          BHCHRG = BHCHRG - CHFCHG(-1)
          IBARYN = IBARYN - 1
          BHBARY = BHBARY + 1
          ISTUP(NUP) = 2
          MOTHUP(1,NUP) = 1
          MOTHUP(2,NUP) = 2
        ELSEIF (IBARYN .EQ. -2) THEN
          IDUP(NUP) = 1
          BHCHRG = BHCHRG - CHFCHG(1)
          IBARYN = IBARYN + 1
          BHBARY = BHBARY - 1
          NUP = NUP + 1
          PUP(1,NUP) = 0.0
          PUP(2,NUP) = 0.0
          PUP(3,NUP) = 0.0
          PUP(4,NUP) = 0.0
          PUP(5,NUP) = 0.0
          IDUP(NUP) = 1
          BHCHRG = BHCHRG - CHFCHG(1)
          IBARYN = IBARYN + 1
          BHBARY = BHBARY - 1
          ISTUP(NUP) = 2
          MOTHUP(1,NUP) = 1
          MOTHUP(2,NUP) = 2
        ELSE
C--for |B|>2 
          NUP = NUP - 1 
          GOTO 35
        ENDIF
C--Skip N-body deay
        GOTO 111
      ENDIF
C--End recoil.
#endif
c-----Determine the statistics of the BH
      STATNBODY=(-1)**(INT(2*JBH))
      DO I=1,NBODYNOW
c--------Select decay products
#if !defined(GINGRICH)
         CALL CHSELECT(ID(I),FSPIN(I),BHCHRG,OBLA,RH,NEX,
     &NFACTOR,TOTFLUX)
#else
         CALL CHSELECT(ID(I),FSPIN(I),BHCHRG,OBLA,RH,NEX,
     &NFACTOR,TOTFLUX,BHBARY)
#endif
         IF(ABS(ID(I)).LE.6) THEN
            IBARYN = IBARYN+SIGN(1,ID(I))
         ENDIF
         PCM = PCM-CHFMAS(ID(I))
         STATNBODY=STATNBODY*(-1)**(INT(FSPIN(I)))
      ENDDO
c-----Accept only if baryon number and charge are conserved, and statistics is correct 
      IF(((IBARYN.NE.0.OR.BHCHRG.NE.0).OR.
     &    PCM.LT.ZERO.OR.STATNBODY.EQ.(-1)).AND.NTRY.LE.NHTRY) THEN
         BHCHRG = SAVCHR
         IBARYN = SAVBRN
C-------Try again
         GOTO 60
      ENDIF

      DO J=0,3
c-----Accept only if lepton number is conserved 
      IF(NLEPTONCSV(J).AND.(ILEPTON(J).NE.0).AND.NTRY.LE.NHTRY) THEN
         DO I=0,3
            ILEPTON(I) = SAVLEPTON(I)
         END DO
C-------Try again
         GOTO 60
      ENDIF
      END DO


C--Go back if number of tries for NBODY decay expires
      IF(NTRY.GT.NHTRY) THEN
         IF(MTRY.LE.MHTRY) THEN
C----------Try again whole decay
            IF(IBHPRN.EQ.2) THEN
               PRINT *,'CHEVNT: Attempt',MTRY,' failed.'
               PRINT *,'CHEVNT: Starting whole decay again'
            ENDIF
            GOTO 35
         ELSE
C----------Give up this event
            IF(IBHPRN.GT.1) 
     &           PRINT *,'CHEVNT: TOO MANY TRIES NEEDED'
            XWGTUP = ZERO
            RETURN
         ENDIF
      ENDIF
C--Perform the NBODY decay of the remnant
      PSTEST=PA(5)
      DO I=1,NBODYNOW
C------Put particle properties in the Les houches common block
         PUP(5,NUP+I) = CHFMAS(ID(I))
         PSTEST=PSTEST-PUP(5,NUP+I)
         ISTUP(NUP+I) = 2
         IDUP(NUP+I) = ID(I)
         MOTHUP(1,NUP+I) = 1
         MOTHUP(2,NUP+I) = 2
      ENDDO
C--check if decay is allowed (if not go to line 100 to force trying again the whole decay)
      IF (PSTEST.LE.ZERO) GOTO 100
      OKDEC=.TRUE.
c-----Phase space combined with Hawking spectra
      IF(.NOT.NBODYPHASE) THEN
C---Sets the oblateness to zero in the non rotating case
         IF(.NOT.BHSPIN) OBLA=0D0
c-----Select angular momenta
         CALL CHJSELECTNBODY(NEX,OBLA/RH,JBH,FSPIN(1),JEMNBODY(1),
     &MEMNBODY(1),NBODYNOW)
c-----Selects the momenta of the decay products and returns if it was successful
         CALL CHPSELECTNBODY(MBH,OBLA,RH,NEX,JEMNBODY(1),MEMNBODY(1),
     &FSPIN(1),PUP(1,NUP+1),NBODYNOW,NBODYCONVERGENCE,FPOLARIZ(1),ID(1))
c-----Try again if no success
         IF(.NOT.NBODYCONVERGENCE)THEN
            BHCHRG = SAVCHR
            IBARYN = SAVBRN
            GOTO 60
         END IF
c-----Boost all momenta back to lab frame
         DO I=1,NBODYNOW
            CALL CHULB4(PA,PUP(1,NUP+I),PUP(1,NUP+I))
         END DO
      ELSE
c-----Phase space decay
 105     IF(NBODYNOW.GT.5) THEN
c-----------Try again with NBODYNOW in the allowed range of 1-5
            BHCHRG = SAVCHR
            IBARYN = SAVBRN
            NBODYNOW=2+INT(CHRGEN(4)*4)
            GOTO 60
         END IF
         IF (NBODYNOW.EQ.2) THEN
            PCM = CHUPCM(PA(5),PUP(5,NUP+1),PUP(5,NUP+2))
#if !defined(GINGRICH)
            CALL CHDTWO(PA,PUP(1,NUP+1),PUP(1,NUP+2),PCM,TWO,.TRUE.)
#else
            IF ((ID(1) .EQ. 39) .AND. (ID(2) .EQ. 39)) THEN
              CALL CHDN(PA,PUP(1,NUP+1),PUP(1,NUP+2),PCM)
            ELSE
              CALL CHDTWO(PA,PUP(1,NUP+1),PUP(1,NUP+2),PCM,TWO,.TRUE.)
            ENDIF
#endif
         ELSEIF (NBODYNOW.EQ.3) THEN
            CALL CHDTHR(PA,PUP(1,NUP+1),PUP(1,NUP+2),PUP(1,NUP+3),OKDEC)
         ELSEIF (NBODYNOW.EQ.4) THEN
            CALL CHDFOR(PA,PUP(1,NUP+1),PUP(1,NUP+2),PUP(1,NUP+3),
     &           PUP(1,NUP+4),OKDEC)
         ELSEIF (NBODYNOW.EQ.5) THEN
            CALL CHDFIV(PA,PUP(1,NUP+1),PUP(1,NUP+2),PUP(1,NUP+3),
     &           PUP(1,NUP+4),PUP(1,NUP+5),OKDEC)
         ENDIF
         GOTO 106        
      END IF
 106  CONTINUE

      IF (OKDEC) GOTO 110
 100  IF(IBHPRN.GE.1) THEN
          PRINT *,'CHEVNT: FAILED',NBODYNOW,
     &     '-BODY REMNANT DECAY: M,Q=',PA(5),PSTEST
      END IF
C---Trying again whole decay if the decay in not OK
      GOTO 35
C---If the decay is OK, sort out polarizations
 110  DO I=1,NBODYNOW 
         IF(.NOT.NBODYPHASE)THEN
            IF(FSPIN(I).NE.0)THEN
               SPINUP(NUP+I)=FPOLARIZ(I)
            ELSE IF(FSPIN(I).EQ.0.AND.(ID(I).EQ.23.OR.ABS(ID(I)).EQ.24))
     &THEN
C---------Longitudinal vector bosons
               SPINUP(NUP+I)=0D0
            ELSE
               SPINUP(NUP+I)=9D0
            END IF
         ELSE
            SPINUP(NUP+I)=9D0
         END IF
      END DO
      NUP = NUP+NBODYNOW
      NBODYBHD=NBODYNOW
#if defined(GINGRICH)
 111  CONTINUE
#endif
C--Now sort out the colour connections
C--First the initial state
      NQRK = 0
      NGLU = 0
      NQBR = 0
      DO I=1,2
        IF(IDUP(I).GT.0.AND.IDUP(I).LE.6) THEN
           NQBR = NQBR+1
           IQBR(NQBR) = I
        ELSEIF(IDUP(I).LT.0.AND.IDUP(I).GE.-6) THEN
           NQRK = NQRK+1
           IQRK(NQRK) = I
        ELSEIF(IDUP(I).EQ.21) THEN
           NGLU = NGLU+1
           IGLU(NGLU) = I
        ENDIF
      ENDDO
C--Then the final state
      DO I=3,NUP
        IF(IDUP(I).GT.0.AND.IDUP(I).LE.6) THEN
           NQRK = NQRK+1
           IQRK(NQRK) = I
        ELSEIF(IDUP(I).LT.0.AND.IDUP(I).GE.-6) THEN
           NQBR = NQBR+1
           IQBR(NQBR) = I
        ELSEIF(IDUP(I).EQ.21) THEN
           NGLU = NGLU+1
           IGLU(NGLU) = I
        ENDIF
      ENDDO
C--zero the colour connections
      DO I=1,NUP
         DO J=1,2
            ICOLUP(J,I) = 0
         ENDDO
      ENDDO
C--Now make the colour connections
      JGLU = 0
      JQBR = 0
      JQRK = 0
      IF(NQRK.GT.0) THEN 
         I = IQRK(1)
         JQRK = 1
         J = 1
      ELSEIF(NGLU.GT.0) THEN
         I = IGLU(1)
         JGLU = 1
         J = 2
      ELSE
         GOTO 200
      ENDIF
      ISTART = I
      ILINE  = 500
C--Find an antiquark or gluon to pair this with
 150   IF(J.EQ.1.OR.J.EQ.2) THEN
         IF(JQBR.EQ.NQBR.AND.JGLU.EQ.NGLU) THEN
            GOTO 250
         ELSE
            PROB = DBLE(NQBR-JQBR)/DBLE(NQBR+NGLU-JQBR-JGLU)
C--Pair with an antiquark
            IF(CHRLOG(PROB)) THEN
               IF(JQBR.NE.NQBR) THEN
                 IF(NGLU.NE.JGLU.and.NQBR-JQBR.EQ.1) GOTO 150
                  JQBR = JQBR+1
                  IF(I.LE.2) THEN
                     ICOLUP(2,I         ) = ILINE
                  ELSE
                     ICOLUP(1,I         ) = ILINE
                  ENDIF
                  IF(IQBR(JQBR).LE.2) THEN
                     ICOLUP(1,IQBR(JQBR)) = ILINE
                  ELSE
                     ICOLUP(2,IQBR(JQBR)) = ILINE
                  ENDIF
                  ILINE = ILINE+1
                  I = IQBR(JQBR)
                  J = 3
                  GOTO 150
               ELSE
                  GOTO 250
               ENDIF
C--Pair with a gluon
            ELSE
               IF(JGLU.NE.NGLU) THEN
                  JGLU = JGLU+1
                  IF(I.LE.2) THEN
                     ICOLUP(2,I         ) = ILINE 
                  ELSE
                     ICOLUP(1,I         ) = ILINE 
                  ENDIF
                  IF(IGLU(JGLU).LE.2) THEN
                     ICOLUP(1,IGLU(JGLU)) = ILINE
                  ELSE
                     ICOLUP(2,IGLU(JGLU)) = ILINE
                  ENDIF
                  ILINE = ILINE+1
                  I = IGLU(JGLU)
                  J = 2
                  GOTO 150
               ELSE
                  GOTO 250
               ENDIF
            ENDIF
         ENDIF
C--Find a quark to pair this with
      ELSEIF(J.EQ.3) THEN
         IF(JQRK.NE.NQRK) THEN
            JQRK = JQRK+1
            IF(I.LE.2) THEN
               ICOLUP(2,I         ) = ILINE
            ELSE
               ICOLUP(1,I         ) = ILINE
            ENDIF
            IF(IQRK(JQRK).LE.2) THEN
               ICOLUP(1,IQRK(JQRK)) = ILINE
            ELSE
               ICOLUP(2,IQRK(JQRK)) = ILINE
            ENDIF
            ILINE = ILINE+1
            I = IQRK(JQRK)
            J = 1
            GOTO 150
         ELSE
            GOTO 250
         ENDIF
      ENDIF
 250  IF(I.LE.2) THEN
         ICOLUP(2,I)      = ILINE
      ELSE
         ICOLUP(1,I)      = ILINE
      ENDIF
      IF(ISTART.LE.2) THEN
         ICOLUP(1,ISTART) = ILINE
      ELSE
         ICOLUP(2,ISTART) = ILINE
      ENDIF
C--zero the flavours
 200  DO I=1,NUP
         IF(IDUP(I).GT.0.AND.IDUP(I).LE.6) THEN
            ICOLUP(2,I) = 0
         ELSEIF(IDUP(I).GE.-6.AND.IDUP(I).LT.0) THEN
            ICOLUP(1,I) = 0
         ELSEIF(IDUP(I).NE.21) THEN
            ICOLUP(1,I) = 0
            ICOLUP(2,I) = 0
         ENDIF
      ENDDO
      END


C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                          Production Routines
C----------------------------------------------------------------------
C----------------------------------------------------------------------


CDECK  ID>, CHPDF.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Peter Richardson
C----------------------------------------------------------------------
      SUBROUTINE CHPDF(XXMINB,XLMINB,EMSCAB,XXB,DISFB)
C----------------------------------------------------------------------
C     Subroutine to calculate the pdfs (generator dependent)
C----------------------------------------------------------------------
#if !defined(HERWIG)
      IMPLICIT NONE
#else
      INCLUDE 'HERWIG65.INC'
#endif
      INCLUDE 'charybdis2.inc'
C--Les Houches run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &                IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &                XMAXUP(MAXPUP),LPRUP(MAXPUP)
C--Local variables
      DOUBLE PRECISION EMSCAB,XXB(2),DISFB(13,2),XXMINB,
     &     XLMINB,CHRUNI,UPV,DNV,USEA,DSEA,STR,CHM,BTM,
     &     TOP,GLU
      INTEGER I
      EXTERNAL CHRUNI
      DOUBLE PRECISION VALUE(20)
      CHARACTER*20 PARM(20)
      LOGICAL FIRST
      DATA VALUE/20*0D0/,PARM/20*' '/
      DATA FIRST/.TRUE./
C--Some local variables use in the Pythia implementation
#if defined(PYTHIA)
      DOUBLE PRECISION XPQ(-25:25),Q2
#else
c--Some local variables used in the other implementations
      INTEGER IDB(2),J
      CHARACTER *8 DUMMY
#endif

C---Use generator internal pdfs
      IF(PDFGUP(1).LT.0) THEN
#if defined(PYTHIA)
C------pythia implementation
         XXB(1)=EXP(CHRUNI(0,0.0D0,XLMINB))
         XXB(2)=XXMINB/XXB(1)
         Q2 = EMSCAB**2
         CALL PYPDFU(IDBMUP(1),XXB(1),Q2,XPQ)
         DO I=1,6
            DISFB(I,1) = XPQ(I)
            DISFB(I+6,1) = XPQ(-I)
         ENDDO
         DISFB(13,1) = XPQ(21)
         CALL PYPDFU(IDBMUP(2),XXB(2),Q2,XPQ)
         DO I=1,6
            DISFB(I,2) = XPQ(I)
            DISFB(I+6,2) = XPQ(-I)
         ENDDO
         DISFB(13,2) = XPQ(21)
#elif defined(HERWIG)
C------herwig implementation
C---SET UP INITIAL STATE
         CALL HWUIDT(1,IDBMUP(1),IDB(1),DUMMY)
         PART1=DUMMY
         CALL HWUIDT(1,IDBMUP(2),IDB(2),DUMMY)
         PART2=DUMMY
         IPART1 = IDB(1)
         IPART2 = IDB(2)
         NHEP=1
         ISTHEP(NHEP)=101
         PHEP(1,NHEP)=0.0D0
         PHEP(2,NHEP)=0.0D0
         PHEP(3,NHEP)=EBMUP(1)
         PHEP(4,NHEP)=EBMUP(1)
         PHEP(5,NHEP)=RMASS(IDB(1))
         JMOHEP(1,NHEP)=0
         JMOHEP(2,NHEP)=0
         JDAHEP(1,NHEP)=0
         JDAHEP(2,NHEP)=0
         IDHW(NHEP)=IPART1
         IDHEP(NHEP)=IDPDG(IDB(1))
         NHEP=NHEP+1
         ISTHEP(NHEP)=102
         PHEP(1,NHEP)=0.0D0
         PHEP(2,NHEP)=0.0D0
         PHEP(3,NHEP)=-EBMUP(2)
         PHEP(4,NHEP)=EBMUP(2)
         PHEP(5,NHEP)=RMASS(IDB(2))
         JMOHEP(1,NHEP)=0
         JMOHEP(2,NHEP)=0
         JDAHEP(1,NHEP)=0
         JDAHEP(2,NHEP)=0
         IDHW(NHEP)=IPART2
         IDHEP(NHEP)=IDPDG(IPART2)
C---NEXT ENTRY IS OVERALL CM FRAME
         NHEP=NHEP+1
         IDHW(NHEP)=14
         IDHEP(NHEP)=0
         ISTHEP(NHEP)=103
         JMOHEP(1,NHEP)=NHEP-2
         JMOHEP(2,NHEP)=NHEP-1
         JDAHEP(1,NHEP)=0
         JDAHEP(2,NHEP)=0
         CALL CHVSUM(4,PHEP(1,NHEP-1),PHEP(1,NHEP-2),PHEP(1,NHEP))
         CALL CHUMAS(PHEP(1,NHEP))
         XLMIN = XLMINB
         EMSCA = EMSCAB
         XXMIN = XXMINB
         CALL HWSGEN(.TRUE.)
         DO I=1,2
            XXB(I) = XX(I)
            DO J=1,13
               DISFB(J,I) = DISF(J,I)
            ENDDO
         ENDDO
#endif
      ELSE
C--Otherwise use pdf libraries
C--Initialisation
         IF(FIRST) THEN
            IF(PDFGUP(1).NE.PDFGUP(2).OR.PDFSUP(1).NE.PDFSUP(2)) THEN
               print *,'MUST HAVE SAME PDF FOR BOTH BEAMS'
               STOP
            ENDIF
#if defined(LHAPDF)
            PARM(1)='DEFAULT'
C---------LHAPDF pdf set number (set this in main program)
            VALUE(1)=LHAPDFSET
C--------- The set number is stored in the relevant variables in the lhes houces common block so that the correct numbers appear in the lhes houches files 
            PDFSUP(1)=LHAPDFSET
            PDFSUP(2)=LHAPDFSET
#elif defined(PDFLIB)
C---------PDFLIB: Set PDFGUP and PDFSUP in main program according to the set you want to use (if not the default in CHDATA).
            PARM(1)='DEFAULT'
            VALUE(1)=1
#else
            PARM(1)='NPTYPE'
            VALUE(1)=1
#endif
            PARM(2)='NGROUP'
            VALUE(2)=PDFGUP(1)
            PARM(3)='NSET'
            VALUE(3)=PDFSUP(1)
            CALL PDFSET(PARM,VALUE)
            FIRST = .FALSE.
         ENDIF
C--------Generate momentum fractions
         XXB(1)=EXP(CHRUNI(0,0.0D0,XLMINB))
         XXB(2)=XXMINB/XXB(1)
C--------Evaluate distribution functions at the given momentum fractions
         DO I=1,2
            CALL STRUCTM(XXB(I),EMSCAB,
     &           UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)


            DISFB( 3,I) = STR
            DISFB( 4,I) = CHM
            DISFB( 5,I) = BTM
            DISFB( 6,I) = TOP
            DISFB( 9,I) = STR
            DISFB(10,I) = CHM
            DISFB(11,I) = BTM
            DISFB(12,I) = TOP
            DISFB(13,I) = GLU
            IF(IDBMUP(I).GT.0) THEN
               DISFB(1,I) = DNV+DSEA
               DISFB(2,I) = UPV+USEA
               DISFB(7,I) = DSEA
               DISFB(8,I) = USEA
            ELSE
               DISFB(1,I) = DSEA
               DISFB(2,I) = USEA
               DISFB(7,I) = DNV+DSEA
               DISFB(8,I) = UPV+USEA
            ENDIF
         ENDDO
      ENDIF

      END

#if defined(PDFLIB)
c---no routines to be added
#elif defined(LHAPDF)
c---no routines to be added
#else 
CDECK  ID>, PDFSET.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C----------------------------------------------------------------------
      SUBROUTINE PDFSET(PARM,VAL)
C----------------------------------------------------------------------
C     DUMMY SUBROUTINE: SET PDFGUP(I) and PDFSUP(I)
C     IN MAIN PROGRAM IF YOU USE PDFLIB CERN-LIBRARY
C     PACKAGE FOR NUCLEON STRUCTURE FUNCTIONS
C----------------------------------------------------------------------
      DOUBLE PRECISION VAL(20)
      CHARACTER*20 PARM(20)
      WRITE (6,10)
   10 FORMAT(/10X,'PDFSET CALLED BUT NOT LINKED')
      STOP
      END
CDECK  ID>, STRUCTM.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
C-----------------------------------------------------------------------
C     DUMMY SUBROUTINE
C-----------------------------------------------------------------------
      DOUBLE PRECISION X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU
      WRITE (6,10)
  10  FORMAT(/10X,'STRUCTM CALLED BUT NOT LINKED')
      STOP
      END
#endif

CDECK  ID>, .
*CMZ :-        -25/03/08  11.30.30  by  Jo Gaunt
*-- Author :   Jo Gaunt
C-----------------------------------------------------------------------
      SUBROUTINE CHYRMJLS(TOTDIM, JFRAC, MFRAC, BFRAC, CVBIASON,ON)
C-----------------------------------------------------------------------
C     This subroutine randomly generates an impact parameter B for the 
C     colliding partons, which is smaller than the Yoshino-Rychkov
C     calculation for BMAX (see Phys.Rev.D71:104028,2005 or 
C     hep-th/0503171). This is output as BFRAC = B / BMAX.
C     
C     The randomly generated B is then used to generate a mass and 
C     angular momentum of the black hole which remains after the 
C     production phase. These are generated such that they lie within 
C     the bounds calculated by Yoshino and Rychkov, but are otherwise 
C     randomly generated using a linear ramp distribution which favours
C     higher remaining J and M. The option to bias the probability 
C     distribution slightly such that it is more likely for the BH to 
C     end up with a similar angular velocity to that with which it 
C     started is included - to turn this on call the routine with 
C     CVBIAS set to .TRUE. 
C
C     The generated mass and angular momentum are output as fractions
C     of the largest mass and angular momentum that could end up in the
C     BH (2*mu and b*mu respectively) - as MFRAC and JFRAC.
C
C     N.B. All lengths in this subroutine are in units of the "R0" 
C     quantity defined in, for example, gr-qc/0209003 equation (17). 
C     All angular momenta are in units of b*mu, and all masses are in 
C     units of 2*mu.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL             ALLOW, CVBIASON,ON
      DOUBLE PRECISION    B, BMAX
      DOUBLE PRECISION    JFRAC, MFRAC, BFRAC
      DOUBLE PRECISION    MIRRED, MLB, CMPRND, CMPFIG
      DOUBLE PRECISION    OMEGAZ, MIRRZ, ASTARZ
      DOUBLE PRECISION    WASTAR, WOMEGA, ODLIM, L
      DOUBLE PRECISION    BLO, MLBLO, BHI, MLBHI, BDAT, MLBDAT
      DOUBLE PRECISION    CHRGEN
      DOUBLE PRECISION    ONE
      INTEGER             TOTDIM
      INTEGER             I, MLBFIL
      EXTERNAL            CHRGEN, CHMJLSPA
      CHARACTER*1         NCHAR
      CHARACTER*24        FNAME

      MLBFIL = 43
      ONE = 1.0D0

C     Parameters associated with the constant velocity bias. 
C     ODLIM is equal to the maximum fractional difference
C     between the angular velocity of a (J,M) point and the
C     angular velocity of JFRAC = MFRAC = 1 after which the
C     probability will not be enhanced.
      ODLIM = 0.2D0
C     L dictates by how much the probability will be 
C     enhanced around the line "ang velocity = initial ang velocity"
C     and how sharply peaked this enhancement is. If L is larger
C     the enhancement is smaller and less peaked.
      L = 0.4D0
      

C     Set BMAX equal to the appropriate Yoshino-Rychkov value.

      IF(TOTDIM==5) THEN
        BMAX = 1.145D0
      ELSEIF (TOTDIM==6) THEN
        BMAX = 1.333D0
      ELSEIF (TOTDIM==7) THEN
        BMAX = 1.441D0
      ELSEIF (TOTDIM==8) THEN
        BMAX = 1.515D0
      ELSEIF (TOTDIM==9) THEN
        BMAX = 1.570D0
      ELSEIF (TOTDIM==10) THEN
        BMAX = 1.613
      ELSEIF (TOTDIM==11) THEN
        BMAX = 1.648D0
      END IF

C     Randomly generate an impact parameter less than BMAX. 
C     The impact parameter is generated such that the second 
C     particle has a uniform probability of hitting anywhere
C     in a circle of radius BMAX around the first particle.

      BFRAC = CHRGEN(1)**0.5D0
      B = BMAX * BFRAC
      IF(B.GT.BMAX)B=BMAX

      IF(.NOT.ON) THEN
         MFRAC=1D0
         JFRAC=1D0
         GOTO 1
      END IF
     
C     The value of MLB associated with the values of TOTDIM and
C     B used is obtained by linear interpolation on the appropriate
C     Yoshino-Rychkov MLB data file.

      WRITE(NCHAR, '(I1)')  TOTDIM - 4

C     Select data file with the correct D.

      FNAME = './data_files/MLB_N'//NCHAR//'.data'

      OPEN(UNIT = MLBFIL, FILE = FNAME, STATUS = 'OLD')

      REWIND(MLBFIL)

      READ(MLBFIL, *) BDAT, MLBDAT

      DO
        BLO = BDAT
        MLBLO = MLBDAT
        READ(MLBFIL, *) BDAT, MLBDAT
        BHI = BDAT
        MLBHI = MLBDAT
        IF(BDAT>B) EXIT
      ENDDO

C     Linear interpolation to our value of B.

      MLB = MLBLO + (MLBHI - MLBLO) * (B-BLO)/(BHI - BLO)

      CLOSE(MLBFIL)

C     Calculate initial angular velocity of the black hole
C     (before production phase) - will be needed for
C     comparison if using constant velocity bias.

      CALL CHMJLSPA(ONE, ONE, B, TOTDIM, MIRRZ, OMEGAZ, ASTARZ)

      I = 2
      ALLOW = .FALSE.

      DO WHILE(.NOT. ALLOW)

C     Randomly generate an MFRAC and JFRAC. 
C     MFRAC is generated according to a linear ramp distribution 
C     which extends between MFRAC = MLB and MFRAC = 1 (one can 
C     show that the lowest mass allowed by the Yoshino-Rychkov
C     bounds is MLB). JFRAC is generated according to a similar
C     distribution which extends between 0 and 1.
      
        MFRAC = CHRGEN(I)
        MFRAC = MFRAC**(0.5)*(1-MLB) + MLB
        JFRAC = CHRGEN(I+1)
        JFRAC = JFRAC**(0.5)

C     Calculate angular velocity WOMEGA and irreducible mass MIRRED
C     associated with the (MFRAC, JFRAC) point generated.

        CALL CHMJLSPA(MFRAC, JFRAC, B, TOTDIM, MIRRED, WOMEGA, 
     &  WASTAR)


        IF(CVBIASON) THEN

C     CMPRND is a comparison random number.

          CMPRND = CHRGEN(I+2)
          I = I+3

C     Generate a comparison figure CMPFIG for the (MFRAC, JFRAC) point
C     according to its angular velocity. The CMPFIG distribution
C     is flat away from the "ang velocity = initial ang velocity" line,
C     but increases close to the line in a Breit-Wigner-type fashion. 

          IF(OMEGAZ.NE.0.0D0 .AND. ABS(WOMEGA/OMEGAZ - 1)<ODLIM .AND. 
     &    ABS(LOG(WASTAR/ASTARZ))<0.4) THEN

            CMPFIG = (L*L/4.0)/((1-WOMEGA/OMEGAZ)**2+L**2/4.0)

            ELSE 

            CMPFIG = (L*L/4.0)/((ODLIM)**2+L**2/4.0)

          END IF

        ELSE

          CMPFIG = 1.0D0
          CMPRND = 0.0D0
          I = I+2

        END IF

C     Accept generated point if it lies within the limits calculated by
C     Yoshino and Rychkov (i.e. if MIRRED >= MLB). If CVBIASON is .TRUE.
C     CMPFIG must also be greater than CMPRND for the point to be accepted
C     (where points close to the constant angular velocity line have higher
C     CMPFIGs and are more likely to be accepted).

        IF (MIRRED>=MLB .AND. ((.NOT. CVBIASON) .OR. 
     &      CMPFIG > CMPRND)) THEN

          ALLOW = .TRUE.
  
        END IF

      ENDDO

 1    CONTINUE
      END

      
CDECK  ID>, .
*CMZ :-        -25/03/08  11.30.30  by  Jo Gaunt
*-- Author :   Jo Gaunt
C-----------------------------------------------------------------------
      SUBROUTINE CHMJLSPA (M, J, B, D, MIRR, OMEGA, ASTAR)
C-----------------------------------------------------------------------
C     Subroutine associated with CHYRMJLS. Calculates the horizon radius,
C     irreducible mass, and horizon angular velocity of a black hole
C     with mass M*2*mu and angular momentum J*B*mu. Irreducible mass
C     is output in units of 2*mu, horizon radius in units of R0, and 
C     angular velocity in units of R0^-1.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION             M, J, B 
      INTEGER                      D
      DOUBLE PRECISION             K1, K2, F, FPRIME
      INTEGER                      DUMMYV
      DOUBLE PRECISION             RKDIFF, RKNEW
      DOUBLE PRECISION             PIFACT
      DOUBLE PRECISION             RHPOW, SPHDM3, SPHDM2
      DOUBLE PRECISION             MIRR, RK, OMEGA, ASTAR
      DOUBLE PRECISION             SPHARY(8)

      RK = 1 
      PIFACT = ACOS(-1.0)

C     Array of areas of the unit sphere in various numbers
C     of dimensions. SPHARY(N) is equal to the area of the 
C     unit (N+1)-sphere.

      SPHARY(1) = 4.0*PIFACT
      SPHARY(2) = 2.0*PIFACT**2
      SPHARY(3) = 8.0*PIFACT**2 / 3.0
      SPHARY(4) = PIFACT**3
      SPHARY(5) = 16.0*PIFACT**3 / 15.0
      SPHARY(6) = PIFACT**4 / 3.0
      SPHARY(7) = 32.0*PIFACT**4 / 105.0
      SPHARY(8) = PIFACT**5 / 12.0

C     Store the area of the unit (D-3)-sphere in SPHDM3
      SPHDM3 = SPHARY(D-4)
C     Store the area of the unit (D-2)-sphere in SPHDM2      
      SPHDM2 = SPHARY(D-3)

C     Schwarzchild radius of a D dimensional black hole of
C     mass M, raised to the power (D-3), in units of R0.
      RHPOW = (4*SPHDM3*M)/((D-2)*SPHDM2)

C     Calculate the horizon radius using a Newton-Raphson
C     procedure with 100 iterations.

      K1 = RHPOW
      K2 = (((D-2)*J)/(4*M))**2 * B**2

      DO DUMMYV = 1,100

        IF (D>=6) THEN
          F      = RK**(D-3) + K2*RK**(D-5)-K1
          FPRIME = (D-3)*RK**(D-4) + K2*(D-5)*RK**(D-6) - K1
        ELSE
          F = K1*RK**(5-D) - RK**2 - K2
          FPRIME = K1*(5-D)*RK**(4-D) - 2.0*RK
        END IF

        RKNEW = RK - F/FPRIME
        RKDIFF = RKNEW - RK
        RK = RKNEW
      END DO

      RK = ABS(RK)
      RKDIFF = ABS(RKDIFF)
      

C     Check that Newton Raphson iteration has converged. If so, use
C     Rk value found to calculate irreducible mass and horizon angular
C     velocity of the black hole. If not, display an error message,
C     and set irreducible mass to a large negative value such that
C     the (M,J) point is rejected in CHYRMJLS.

      IF (RKDIFF<=0.0001*RK) THEN
        
        MIRR = ((D-2)*SPHDM2)/(4*SPHDM3)
     &  * (RHPOW * RK)**((D-3.0)/(D-2.0))
        
        ASTAR = (D-2)*J*B/(4.0*M*RK)

        OMEGA = ASTAR/(RK*(1.0+ASTAR*ASTAR))
        
        ELSE
        MIRR = -1000.0
        OMEGA = 0.0
      END IF
 

      END

CDECK  ID>, CHINTRINSICSPIN.
*CMZ :-        -11/11/08  18.11.30  by Marco Sampaio
*-- Author :   Marco Sampaio 
C---------------------------------------------------------------------
      SUBROUTINE CHINTRINSICSPIN(IDN,J,M)
C-----------------------------------------------------------------------
C     Subroutine to select the partial wave from unpolarized collisions
C     of particles with IDs IDN, with respect to eigenstates of spin 
C     projected along a perpendicular axis to that of the initial
C     combined angular momentum states of the incident particles
C-----------------------------------------------------------------------
C  *INPUT:
C     --IDN: IDs of the incoming particles
C  *OUTPUT:
C     --J,M: Partial wave selected from the wigner rotated state
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--I/O
      INTEGER IDN(2)
      DOUBLE PRECISION J,M
C-- Local variables
      INTEGER ALIGN,MSIGN
      DOUBLE PRECISION CUMULATIVE,CHRGEN
      LOGICAL CHRLOG
      EXTERNAL CHRGEN,CHRLOG
      
      ALIGN=1
      MSIGN=1
C-- Select alignment of helicities of incoming particles
      IF(CHRLOG(0.5D0)) ALIGN=-1

      IF((IDN(1).EQ.21).AND.(IDN(2).EQ.21))THEN
C-- Gluon-Gluon collision
C------Total J from helicity combinations set by ALIGN
         J=ABS(1D0+ALIGN) 
C------Choose ABS(M)
         IF(J.EQ.0)THEN
            M=0D0
         ELSE
C------Flip sign of M half of the time
            IF(CHRLOG(0.5D0)) MSIGN=-1  
            CUMULATIVE=CHRGEN(1)*8D0
            IF(CUMULATIVE.LT.1)THEN
               M=2D0
            ELSE IF(CUMULATIVE.LT.5)THEN
               M=1D0
            ELSE
               M=0D0
            END IF
            M=MSIGN*M
         END IF
      ELSE IF ((IDN(1).EQ.21).OR.(IDN(2).EQ.21))THEN
C-- Gluon-Quark collision
C------Total J from helicity combinations set by ALIGN
         J=ABS(1D0+0.5D0*ALIGN)
C------Flip sign of M half of the time
         IF(CHRLOG(0.5D0)) MSIGN=-1   
C------Choose ABS(M)
         IF(J.EQ.0.5D0)THEN
            M=0.5D0
         ELSE
            IF(CHRGEN(1)*4D0.LT.1)THEN
               M=1.5D0
            ELSE
               M=0.5D0
            END IF
         END IF
         M=MSIGN*M
      ELSE
C-- Quark-Quark collision
C------Total J from helicity combinations set by ALIGN
         J=0.5D0*ABS(1D0+ALIGN)  
         M=0D0
C------Choose ABS(M)
C------Flip sign of M half of the time
         IF(CHRLOG(0.5D0)) MSIGN=-1 
         IF(CHRLOG(0.5D0).AND.(J.EQ.1D0)) M=1D0
         M=MSIGN*M
      END IF
      
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                          Evaporation Routines
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, CHSELECT.
*CMZ :-        -24/01/08  18.11.30  by Marco Sampaio
*-- Author :   Marco Sampaio (adapted from Chris Harris' routine in 
*              CHARYBDIS 1)
C---------------------------------------------------------------------
#if !defined(GINGRICH)
      SUBROUTINE CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX)
#else
      SUBROUTINE CHSELECT(PT,SPIN,BHCHRG,OBLA,RH,NEX,NFACTOR,TOTFLUX,
     &BHBARY)
#endif
C-----------------------------------------------------------------------
C     Subroutine to select the type of the next particle, by calculating
C     cumulants
C-----------------------------------------------------------------------
C  *INPUT:
C     --BHCHRG,OBLA,RH,NEX: Charge, oblateness and radius of the black 
C       hole; and number of extra dimensions
C  *OUTPUT:
C     --PT,SPIN,NFACTOR: Particle type, and spin, average multiplicity 
C       factor
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
C--Input/Output
      INTEGER PT,SPIN,BHCHRG,NEX
      DOUBLE PRECISION OBLA,RH,BMASS
C--Local variables


C-----MODIFICATION-----------------------------------------------------
      INTEGER CHFCHG,UFLUX,UEFLUX,MA,MSPIN,IASTAR,INDX
      DOUBLE PRECISION DSUM, KL
      DOUBLE PRECISION WGTUP,WGTCHM,WGTTOP,WGTDWN,WGTSTR,WGTBTM
      PARAMETER(MA=26,MSPIN=3)
      DOUBLE PRECISION  DFERMIONL(6), DFERMIONR(6)
      DATA DFERMIONL /0.0006851,0.0150148,0.0600244,0.0006851,
     &0.0150148,0.0600244/
      DATA DFERMIONR / 2.9D-8, 0.0021630, 0.4961677, 2.1D-7, 3.6D-6,
     &0.0028090/

      DOUBLE PRECISION PTYPE,CHRGEN,FLUX(MSPIN,MA),CFERMION(3),
     &DGAUGE,DHIGGS,PI
C----------------------------------------------------------------------




#if !defined (GINGRICH)
     &,CBOSON(5),ASTAR,ASTARMAX,F2S0,F2S1,F2S2,FRACDA,PTOT,NFACTOR,
     &TOTFLUX
      DOUBLE PRECISION EFLUX(MSPIN,MA),EF2S0,EF2S1,EF2S2,EPTOT
#else
     &,CBOSON(6),ASTAR,ASTARMAX,F2S0,F2S1,F2S2,F2S4,FRACDA,PTOT,NFACTOR,
     &TOTFLUX
      DOUBLE PRECISION EFLUX(MSPIN,MA),EF2S0,EF2S1,EF2S2,EF2S4,EPTOT
      INTEGER BHBARY
      DOUBLE PRECISION GFLUX(6)
      DATA GFLUX /0.6E0,0.91D0,1.9D0,2.5D0,5.1D0,7.6D0/ 
#endif

      PARAMETER(ASTARMAX=5D0)
      LOGICAL FIRST,CHRLOG
      DATA FIRST/.TRUE./
      EXTERNAL CHRGEN,CHFCHG,CHRLOG
      CHARACTER NEXTRA*1
      SAVE FLUX,EFLUX

      UFLUX=25
      UEFLUX=26
      NEXTRA=CHAR(ICHAR('0')+NEX)
C--Initialises the arrays with cumulative fluxes for each SPIN by 
C  reading the corresponding data file for a number of extra dimensions 
C  NEX. This occurs only once at FIRST call.
      IF (FIRST) THEN
         IF ((NEX.GT.0).AND.(NEX.LT.7))
     & THEN
            OPEN(UNIT=UFLUX,FILE='./data_files/flux_n'//NEXTRA//'.dat',
     &STATUS='OLD')
            OPEN(UNIT=UEFLUX,FILE='./data_files/eflux_n'//NEXTRA//'.dat'
     &,STATUS='OLD')
         ELSE
            PRINT *,' NO FLUX TABLE YET FOR NEX =',NEX
            STOP
         ENDIF

         READ(UFLUX,*)FLUX
         READ(UEFLUX,*)EFLUX         
         CLOSE(UNIT=UFLUX)
         CLOSE(UNIT=UEFLUX)

         FIRST=.FALSE.
      END IF
C--Select the correct IASTAR in the flux tables

      ASTAR=OBLA/RH
      IF (ASTAR.LE.ASTARMAX) THEN
         IASTAR=1+INT(ASTAR/ASTARMAX*(MA-1.))
      ELSE
         IASTAR=MA
      END IF
C--Interpolate the flux for each spin      
      IF(GRYBDY)THEN
         IF (IASTAR.LT.MA)THEN
            FRACDA=ASTAR*(MA-1D0)/ASTARMAX-(IASTAR-1) 
            F2S0=FLUX(1,IASTAR+1)*FRACDA+FLUX(1,IASTAR)*(1-FRACDA)
            F2S1=FLUX(2,IASTAR+1)*FRACDA+FLUX(2,IASTAR)*(1-FRACDA)
            F2S2=FLUX(3,IASTAR+1)*FRACDA+FLUX(3,IASTAR)*(1-FRACDA)
            EF2S0=EFLUX(1,IASTAR+1)*FRACDA+EFLUX(1,IASTAR)*(1-FRACDA)
            EF2S1=EFLUX(2,IASTAR+1)*FRACDA+EFLUX(2,IASTAR)*(1-FRACDA)
            EF2S2=EFLUX(3,IASTAR+1)*FRACDA+EFLUX(3,IASTAR)*(1-FRACDA)
         ELSE
            F2S0=FLUX(1,IASTAR)
            F2S1=FLUX(2,IASTAR)
            F2S2=FLUX(3,IASTAR)
            EF2S0=FLUX(1,IASTAR)
            EF2S1=FLUX(2,IASTAR)
            EF2S2=FLUX(3,IASTAR)
         END IF
#if defined (GINGRICH)
         F2S4  = F2S0*GFLUX(NEX-1)
         EF2S4 = F2S0*GFLUX(NEX-1)
#endif
      ELSE
C------Here the geometrical approximation was used
         F2S0=1D0
         F2S1=3D0/4D0
         F2S2=1D0
         EF2S0=1D0
         EF2S1=7D0/8D0
         EF2S2=1D0
#if defined (GINGRICH)
         F2S4=1D0
         EF2S4=1D0
#endif
      ENDIF
      IF((F2S0.LT.0).OR.(F2S1.LT.0).OR.(F2S2.LT.0).OR.(EF2S0.LT.0).OR.
     &(EF2S1.LT.0).OR.(EF2S2.LT.0))THEN 
         WRITE(*,*)'CH: WARNING - Negative fluxes not allowed'
         STOP
      END IF
C-----MODIFICATION----------------------------------------------------
      PI=ACOS(-1D0)
      KL=11.3D0*PI
C      WRITE(*,*) BMASS,MPLNCK

C-----Left handed-----------------------------------------------------
      DSUM=0D0
      DO 1 INDX=1,6
         DSUM=DSUM+DFERMIONL(INDX)
C         WRITE(*,*) INDX,DFERMION(INDX)
 1    CONTINUE

C-----Right handed---------------------------------------------------     
      DO 2 INDX=1,6
         DSUM=DSUM+DFERMIONR(INDX)
C         WRITE(*,*) INDX,DFERMION(INDX)
 2    CONTINUE      
      WGTUP=F2S1*6*(DFERMIONL(1)+DFERMIONR(1))
      WGTCHM=WGTUP+F2S1*6*(DFERMIONL(2)+DFERMIONR(2))
      WGTTOP=WGTCHM+F2S1*6*(DFERMIONL(3)+DFERMIONR(3))
      WGTDWN=WGTTOP+F2S1*6*(DFERMIONL(4)+DFERMIONR(4))
      WGTSTR=WGTDWN+F2S1*6*(DFERMIONL(5)+DFERMIONR(5))
      WGTBTM=WGTSTR+F2S1*6*(DFERMIONL(6)+DFERMIONR(6))

      DGAUGE=0.0281690
      DHIGGS=1D0     
C---------------------------------------------------------------------


C      WRITE(*,*) EPSLN,DFERMION(3)+DFERMION(9)
C--Calculate cumulants in case we are allowing light particles only
      IF(MSSDEC.EQ.1)THEN
C--Leptons (6 d.o.f. = 3*2*2 = #families * (particles + antiparticles)*2helicitystates.
         CFERMION(1)=F2S1*12.
C-- + Neutrinos (3*2 = #families * (particles + antiparticles) )
         CFERMION(2)=CFERMION(1)+F2S1*6.
C--   + 5 light Quarks (5*2*3*2 = #light quarks * (part. +antipart.) * #colours*2helicities)
         CFERMION(3)=CFERMION(2)+F2S1*60.
C-- + Gluons + Photon
C--------(8 gluons * 2 polarizations)
         CBOSON(1)=CFERMION(3)+F2S2*16.
C--------(2 polarizations)
         CBOSON(2)=CBOSON(1)+F2S2*2.
         CBOSON(3)=CBOSON(2)
         CBOSON(4)=CBOSON(2)
         CBOSON(5)=CBOSON(2)
         PTOT=CBOSON(2)
         EPTOT=EF2S1*12.+EF2S1*6.+EF2S1*60.+EF2S2*16.+EF2S2*2.
C--Light+W+Z
      ELSE IF(MSSDEC.EQ.2)THEN
C--Leptons (as before)
         CFERMION(1)=F2S1*12.
C-- + Neutrinos (as before)
         CFERMION(2)=CFERMION(1)+F2S1*6.
C-- + all quarks (as before with 5-->6)
         CFERMION(3)=CFERMION(2)+F2S1*72.
C-- + Gluons + Photon + W + Z 
C-------(as before)
         CBOSON(1)=CFERMION(3)+F2S2*16.
C-------(as before)
         CBOSON(2)=CBOSON(1)+F2S2*2.
C-------(same as for photon + scalar contribution from longitudinal mode of massive W and a factor of two accounting for W+ and W-)
         CBOSON(3)=CBOSON(2)+2*(F2S0+2.*F2S2)
C-------(same as for W but no factor of 2 because Z is not electrically charged)
         CBOSON(4)=CBOSON(3)+(F2S0+2.*F2S2)
         CBOSON(5)=CBOSON(4)
         PTOT=CBOSON(4)
         EPTOT=EF2S1*12.+EF2S1*6.+EF2S1*72.+EF2S2*16.+EF2S2*2.
     &+2*(EF2S0+2.*EF2S2)+(EF2S0+2.*EF2S2)
C--All the standard model particles
      ELSE IF(MSSDEC.EQ.5)THEN
         
         CFERMION(1)=0.
         CFERMION(2)=CFERMION(1)+0.
         CFERMION(3)=CFERMION(2)+F2S1*3*2*DSUM !TOP only 24.
         CBOSON(1)=CFERMION(3)+F2S2*16*DGAUGE !GLUON
         CBOSON(2)=CBOSON(1)+F2S2*2*DGAUGE
         CBOSON(3)=CBOSON(2)+2*(F2S0+2.*F2S2)*DGAUGE
         CBOSON(4)=CBOSON(3)+(F2S0+2.*F2S2)*DGAUGE
         CBOSON(5)=CBOSON(4)+F2S0*DHIGGS
         PTOT=CBOSON(5)
         EPTOT=EF2S1*6*DSUM+EF2S2*16*DGAUGE+EF2S2*2*DGAUGE
     &+2*(EF2S0+2.*EF2S2)*DGAUGE+(EF2S0+2.*EF2S2)*DGAUGE+EF2S0*DHIGGS

C         WRITE(*,*) CBOSON(5),CBOSON(6)
#if !defined (GINGRICH)
      ELSE
#else
      ELSE IF(MSSDEC.EQ.3)THEN
#endif
C--Leptons (as before)
         CFERMION(1)=F2S1*12.
C-- + Neutrinos (as before)
         CFERMION(2)=CFERMION(1)+F2S1*6.
C-- + all quarks (as before)
         CFERMION(3)=CFERMION(2)+F2S1*72.
C-- + Gluons + Photon + W + Z + Higgs (as before but with one extra real scalar d.o.f. for the Higgs)
         CBOSON(1)=CFERMION(3)+F2S2*16.
         CBOSON(2)=CBOSON(1)+F2S2*2.
         CBOSON(3)=CBOSON(2)+2*(F2S0+2.*F2S2)
         CBOSON(4)=CBOSON(3)+(F2S0+2.*F2S2)   
C--------Extra real scalar d.o.f for the Higgs
         CBOSON(5)=CBOSON(4)+F2S0 
         PTOT=CBOSON(5)
         EPTOT=EF2S1*12.+EF2S1*6.+EF2S1*72.+EF2S2*16.+EF2S2*2.
     &+2*(EF2S0+2.*EF2S2)+(EF2S0+2.*EF2S2)+EF2S0 
#if defined (GINGRICH)
C--All the standard model particles + graviton in D-dimensions.
      ELSE IF(MSSDEC.EQ.4)THEN
C--Leptons (as before)
         CFERMION(1)=F2S1*12.
C-- + Neutrinos (as before)
         CFERMION(2)=CFERMION(1)+F2S1*6.
C-- + all quarks (as before)
         CFERMION(3)=CFERMION(2)+F2S1*72.
C-- + Gluons + Photon + W + Z + Higgs (as before but with one extra real scalar d.o.f. for the Higgs)
         CBOSON(1)=CFERMION(3)+F2S2*16.
         CBOSON(2)=CBOSON(1)+F2S2*2.
         CBOSON(3)=CBOSON(2)+2*(F2S0+2.*F2S2)
         CBOSON(4)=CBOSON(3)+(F2S0+2.*F2S2)   
C--------Extra real scalar d.o.f for the Higgs
         CBOSON(5)=CBOSON(4)+F2S0 
C--------Graviton.
         CBOSON(6)=CBOSON(5)+F2S4         
         PTOT=CBOSON(6)
         EPTOT=EF2S1*12.+EF2S1*6.+EF2S1*72.+EF2S2*16.+EF2S2*2.
     &+2*(EF2S0+2.*EF2S2)+(EF2S0+2.*EF2S2)+EF2S0+EF2S4
     
#endif 
      END IF

C----Factor for the estimate of number of emissions in the remainder of the evaporation      
      IF(EPTOT.GT.0)THEN
         NFACTOR=PTOT/EPTOT
         TOTFLUX=PTOT/(2D0*ACOS(-1D0))/RH
      ELSE
         WRITE(*,*)'CHSELECT CRASH: integrated energy flux not positive'
         WRITE(*,*)'---Check if something is wrong with flux data files'
         STOP
      END IF

C--Assume a fermion (and change below when it is not. Note: SPIN is actually 2*spin) 
      SPIN=1
C--Choose decay product
      PTYPE=CHRGEN(2)*PTOT
C--ChangeS SPIN to 2 if we have a boson instead (except for one 
C  degree of freedom of massive gauge bosons)
      IF (PTYPE.GE.CFERMION(3)) THEN
         SPIN=2
      ENDIF
C--Choose particle type
      IF (PTYPE.LT.CFERMION(1)) THEN
C--Leptons
         PT=11+2*INT(3D0*PTYPE/CFERMION(1))
C--Ensure BH charge remains close to 0, and antiparticles half time
         IF (((CHFCHG(PT)*BHCHRG).LT.0).OR.
     &       ((BHCHRG.EQ.0).AND.CHRLOG(0.5D0))) THEN 
            PT=-PT
         ENDIF
      ELSEIF (PTYPE.LT.CFERMION(2)) THEN
C--Neutrinos
         PT=12+2*INT(3*(PTYPE-CFERMION(1))/(CFERMION(2)-CFERMION(1)))
         IF(CHRLOG(0.5D0)) PT =-PT
      ELSEIF (PTYPE.LT.CFERMION(3)) THEN
C--Quarksi---MODIFICATION------------------------------------------
         IF(MSSDEC.EQ.1)THEN
            PT=1+INT(5*(PTYPE-CFERMION(2))/(CFERMION(3)-CFERMION(2)))
         ELSE IF(MSSDEC.EQ.5)THEN
            IF(PTYPE.GT.CFERMION(2).AND.PTYPE.LE.WGTUP)THEN
               PT=2
            ELSEIF(PTYPE.GT.WGTUP.AND.PTYPE.LE.WGTCHM)THEN
               PT=4           
            ELSEIF(PTYPE.GT.WGTCHM.AND.PTYPE.LE.WGTTOP)THEN
               PT=6
            ELSEIF(PTYPE.GT.WGTTOP.AND.PTYPE.LE.WGTDWN)THEN
               PT=1
            ELSEIF(PTYPE.GT.WGTDWN.AND.PTYPE.LE.WGTSTR)THEN
               PT=3
            ELSEIF(PTYPE.GT.WGTSTR.AND.PTYPE.LE.WGTBTM)THEN
               PT=5
            ENDIF
C-----------------------------------------------------------------


         ELSE
            PT=1+INT(6*(PTYPE-CFERMION(2))/(CFERMION(3)-CFERMION(2))) 
         END IF
C--Ensure BH charge remains close to 0, and antiparticles half time
         IF (((CHFCHG(PT)*BHCHRG).LT.0).OR.
     &       ((BHCHRG.EQ.0).AND.CHRLOG(0.5D0))) THEN 
            PT=-PT
         ENDIF
      ELSEIF (PTYPE.LT.CBOSON(1)) THEN
C--Gluons
         PT=21
      ELSEIF (PTYPE.LT.CBOSON(2)) THEN
C--Photons
         PT=22
      ELSEIF (PTYPE.LT.CBOSON(3)) THEN
C--W+/W-
         IF(CHRGEN(4).GT.(F2S2*2D0)/(F2S2*2D0+F2S0))THEN
            SPIN=0
         END IF
         PT=24
C--Ensure BH charge remains close to 0, and antiparticles half time
         IF (((CHFCHG(PT)*BHCHRG).LT.0).OR.
     &       ((BHCHRG.EQ.0).AND.CHRLOG(0.5D0))) THEN 
            PT=-PT
         ENDIF
      ELSEIF (PTYPE.LT.CBOSON(4)) THEN
C--Z0
         IF(CHRGEN(4).GT.(F2S2*2D0)/(F2S2*2D0+F2S0))THEN
            SPIN=0
         END IF
         PT=23
#if !defined (GINGRICH)
      ELSE
#else
      ELSEIF (PTYPE.LE.CBOSON(5)) THEN
#endif

C--Higgs
         PT=25
         SPIN=0
#if defined (GINGRICH)
C--Graviton
      ELSEIF (PTYPE.LE.CBOSON(6)) THEN
         PT=39
C        There are no grey-body spectra for gravitons so we must
C        use the scalar, fermion, or vector boson spectra.
C        SPIN=0 Higgs, SPIN=1 fermions, SPIN=2 vector bosons.
         SPIN=1
      ELSE
        print *,'invalid'
        stop
#endif
      ENDIF
C--Calculate new charge of BH
      BHCHRG=BHCHRG-CHFCHG(PT)
#if defined(GINGRICH)
C-- Calculate new baryon number of BH
      BHBARY=BHBARY-SIGN(1,PT)
#endif
      END

CDECK  ID>, CHEMIT.
*CMZ :-        -24/01/08  15.58.30  by  Bryan Webber, Marco Sampaio
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHEMIT(NEX,MBH,RH,OBLA,FMAXEM,MSEM,SPIN,EEM,JEM,MEM)
C-----------------------------------------------------------------------
C     Generates an emission in the rest frame of the black hole by 
C     selecting the energy (EEM) and the mode (JEM,EM) within the 
C     kinematically allowed region in [MSEM,EMAX].
C-----------------------------------------------------------------------
C  *INPUT:
C     --NEX,MBH,RH,OBLA: are number of extra dimensions, mass of 
C       the black hole, horizon radius oblateness and Hawking 
C       temperature respectively.
C     --MSEM,SPIN: are the mass and spin of the particle emitted.
C  *OUTPUT:
C     --EEM,JEM,MEM: The energy and angular momentum of the particle 
C       emitted.
C----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
C--Input/Output:
      INTEGER NEX,SPIN
      DOUBLE PRECISION MBH,RH,OBLA,MSEM,EEM,JEM,MEM,FMAXEM
C--Local variables
      INTEGER NTRY,IASTAR,K,MA,ME,MK0,MK1,MK2,UCFLUXKE,ITRY
      PARAMETER (MA=26,MK0=169,MK1=156,MK2=120,ME=100)
      DOUBLE PRECISION ASTAR,CHRGEN, CFLUX2S0K(MK0+1,MA),
     &CFLUX2S1K(MK1+1,MA),CFLUX2S2K(MK2+1,MA),CFLUX2S0KE(ME+1,MK0,MA),
     &CFLUX2S1KE(ME+1,MK1,MA),CFLUX2S2KE(ME+1,MK2,MA),FRACDA,ASTARMAX,
     &EEMRH,EMAX
      PARAMETER(ASTARMAX=5D0)
      LOGICAL FIRST,CONVERGENCE
      EXTERNAL CHRGEN
      CHARACTER NEXTRA*1

      DATA NTRY,UCFLUXKE,FIRST/200000,23,.TRUE./
      SAVE CFLUX2S0K,CFLUX2S1K,CFLUX2S2K,CFLUX2S0KE,CFLUX2S1KE,
     &CFLUX2S2KE

C-----Compute maximum energy allowed kinematically for the recoil model in use (or if KINCUT then in any region)
      IF(RECOIL.EQ.1)THEN
         EMAX=(MBH**2+MSEM**2)/(2D0*MBH)*FMAXEM
      ELSE IF (RECOIL.EQ.2)THEN
         EMAX=MBH*FMAXEM
      ELSE
         WRITE(*,*)'Invalid recoil option!'
         STOP
      END IF
C--Initialises the arrays with cumulative fluxes by reading the 
C  corresponding data file for a number of extra dimensions NEX This 
C  occurs only once at FIRST call.

      IF (FIRST) THEN
         NEXTRA=CHAR(ICHAR('0')+NEX)
         IF ((NEX.GT.0).AND.(NEX.LT.7))
     &THEN
            OPEN(UNIT=UCFLUXKE,FILE='./data_files/cfluxke_n'//NEXTRA//
     &'.dat',STATUS='OLD')
         ELSE
            PRINT *,' NO CUMULATIVE TABLES YET FOR NEX =',NEX
            STOP
         ENDIF
         READ(UCFLUXKE,*)CFLUX2S0K,CFLUX2S1K,CFLUX2S2K,CFLUX2S0KE,
     &CFLUX2S1KE,CFLUX2S2KE
         CLOSE(UNIT=UCFLUXKE)

         FIRST=.FALSE.
      ENDIF

C--Select the correct IASTAR (the position of the astar imediately below 
C  the actual one) in the flux tables. Interpolation will be done between 
C  IASTAR and IASTAR+1.
      ASTAR=OBLA/RH
      IF (ASTAR.LE.ASTARMAX) THEN
         IASTAR=1+INT(ASTAR/ASTARMAX*(MA-1D0))
C------- FRACDA is the linear interpolation weight  at IASTAR+1 and 
C        1-FRACDA at IASTAR
         FRACDA=ASTAR*(MA-1.)/ASTARMAX-(IASTAR-1)  
      ELSE
         IASTAR=MA-1
         FRACDA=1D0  
      END IF
C--Select the mode (K) and then the energy (EEMRH)
      DO 5 ITRY=0,NTRY
         IF(SPIN.EQ.0)THEN
C-----Scalars
            CALL CHFINDK(CFLUX2S0K(1,IASTAR),CFLUX2S0K(1,IASTAR+1),MK0,
     &ASTAR,FRACDA,K)
            CALL CHFINDE(CFLUX2S0KE(1,K,IASTAR),CFLUX2S0KE(1,K,IASTAR+1)
     &,ME,ASTAR,FRACDA,MSEM*RH,EMAX*RH,CONVERGENCE,EEMRH)
            JEM=INT(SQRT(K-1D-2))
            MEM=K-JEM*(JEM+1)-1
         ELSE IF(SPIN.EQ.1)THEN
C-----Fermions
            CALL CHFINDK(CFLUX2S1K(1,IASTAR),CFLUX2S1K(1,IASTAR+1),MK1,
     &ASTAR,FRACDA,K)
            CALL CHFINDE(CFLUX2S1KE(1,K,IASTAR),CFLUX2S1KE(1,K,IASTAR+1)
     &,ME,ASTAR,FRACDA,MSEM*RH,EMAX*RH,CONVERGENCE,EEMRH)
            JEM=INT(SQRT(K-1D-2)-0.5D0)+0.5D0
            MEM=K-JEM*(JEM+1D0)-0.75D0
         ELSE IF(SPIN.EQ.2)THEN
C-----Vectors
            CALL CHFINDK(CFLUX2S2K(1,IASTAR),CFLUX2S2K(1,IASTAR+1),MK2,
     &ASTAR,FRACDA,K)   
            CALL CHFINDE(CFLUX2S2KE(1,K,IASTAR),CFLUX2S2KE(1,K,IASTAR+1)
     &,ME,ASTAR,FRACDA,MSEM*RH,EMAX*RH,CONVERGENCE,EEMRH)
            JEM=INT(SQRT(K+1D-2))
            MEM=K-JEM*(JEM+1D0)
         END IF
C------Get out of the loop if CONVERGENCE=.TRUE.
         IF(CONVERGENCE) GOTO 6
 5    CONTINUE
      WRITE(*,*)'WARNING CHEMIT: NO CONVERGENCE! BHMASS= ',EMAX
      WRITE(*,*)'Energy selected at random in the allowed region'
      EEMRH=RH*(MSEM+CHRGEN(1)*(EMAX-MSEM))
 6    CONTINUE
      EEM=EEMRH/RH

      END

CDECK  ID>, CHNOGRYBDY
*CMZ :-        
*-- Author :   Marco Sampaio (adapted from a part of Chris Harris' 
C                             routine CHUBHS in Charybdis 1)
C---------------------------------------------------------------------
      SUBROUTINE CHNOGRYBDY(ENERGY,SPIN,THAWRH)
C---------------------------------------------------------------------
C     Selects energy by using the geometrical approximations for the 
C     energy spectrum for Schwarzschild BH without taking greybody 
C     factors into account
C---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
      DOUBLE PRECISION ENERGY,SPCMAX,SPCVAL,THAWRH,ETRAT
      DOUBLE PRECISION CHRGEN
      INTEGER STAT,SPIN
      EXTERNAL CHRGEN

      STAT=(-1)**SPIN

C----These were calculated using the geometrical optics approximation
      IF (STAT.EQ.-1) THEN
         SPCMAX=0.5*THAWRH**2
      ELSE
         SPCMAX=2./3.*THAWRH**2
      ENDIF

C--MC for energy
 10   ENERGY=CHRGEN(1)*5.
      ETRAT=ENERGY/THAWRH
C-----Using the geometrical optics spectrum (sum over partial waves has been done before hence the omega**2 factor)
      SPCVAL=(ENERGY**2)/(EXP(ETRAT)-STAT)
      IF (SPCVAL.LT.SPCMAX*CHRGEN(1)) GOTO 10

      END

CDECK  ID>, CHEMOM.
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHEMOM(PSEM,EEM,JEM,MEM,OBLA,PEM,SPINTIMES2,POLARIZ,PT)
C-----------------------------------------------------------------------
C     Generates the momentum of the emission using spheroidal angular 
C     functions.
C-----------------------------------------------------------------------
C  *INPUT:
C     --PSEM,EEM,JEM,MEM: Momentum magnitude, Energy and angular 
C       momentum of the emitted particle
C     --OBLA: the oblateness of the BH
C     --SPINTIMES2, POLARIZ PT: Spin, helicity and particle type of the
C       emission 
C  *OUTPUT:
C     --PEM: The 3 momentum of the emission.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Input/Output
      DOUBLE PRECISION PSEM,EEM,JEM,MEM,PEM(3),OBLA
      INTEGER SPINTIMES2,POLARIZ,PT
C--Local variables
C----Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
C----Other variables
      INTEGER LTIMES2,MTIMES2,NLAST1,NLAST2,NYSIZE
      DOUBLE PRECISION PTRANSVERSE,CHRGEN,C,COSTHETA,
     &SPHRDLSERIES1(0:NMAX),CONSTANTS1(1:4),DELTAN1,
     &SPHRDLSERIES2(0:NMAX),CONSTANTS2(1:4),DELTAN2,A,Y(0:NDIVIDEMAX+1),
     &XOFY(1:NDIVIDEMAX+1)
      EXTERNAL CHRGEN

C---Some parameters to be used      
      LTIMES2=INT(2*JEM)
      MTIMES2=INT(2*MEM)
      C=OBLA*EEM
C---Initialise polarization. It will be flipped half of the time.
      POLARIZ=-1
      IF(CHRGEN(1).LT.0.5) POLARIZ=1
C-----Neutrinos (left handed) and antineutrinos (right handed)
      IF(PT.EQ.12.OR.PT.EQ.14.OR.PT.EQ.16) POLARIZ=-1
      IF(PT.EQ.-12.OR.PT.EQ.-14.OR.PT.EQ.-16) POLARIZ=1
C-----Flip to counterintuitive convention for the sign of SPINTIMES2 (it will be restored by the end of this routine). The reason for this is that the spheroidal functions were coded in first, using the predominant convention in the literature for which positive SPINTIMES2 neutrinos are emitted mostly antiparalelly to the spin axis of the black hole.
      SPINTIMES2=-SPINTIMES2*POLARIZ
c-----Initialise spheroidal function      
      CALL CHSPHRDLINI(LTIMES2,MTIMES2,SPINTIMES2,C,SPHRDLSERIES1,
     &CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2)
c-----Match two expansions for spheroidal function
      CALL CHSPHRDLMATCH(C,MTIMES2,SPINTIMES2,SPHRDLSERIES1,CONSTANTS1,
     &DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2,A)
c-----Compute the cumulative function
      CALL CHSPHRDLCUMULATIVE(C,LTIMES2,MTIMES2,SPINTIMES2,SPHRDLSERIES1
     &,CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2
     &,A,Y,XOFY,NYSIZE)
c-----Select a COSTHETA using the cumulative function
      CALL CHSPHRDLXRAND(Y,XOFY,NYSIZE,COSTHETA)
c-----Compute the momentum vector in rest frame of the BH
      PEM(3)=PSEM*COSTHETA
      PTRANSVERSE=SQRT(PSEM**2-PEM(3)**2)
      CALL CHRAZM(PTRANSVERSE,PEM(1),PEM(2))
C-----Restore convention
      SPINTIMES2=-SPINTIMES2*POLARIZ

      END

CDECK  ID>, CHNEWJ.
*CMZ :-        -22/08/06  15.58.30  by  Bryan Webber
*-- Author :    Bryan Webber
C----------------------------------------------------------------------
      SUBROUTINE CHNEWJ(JBH,JEM,MEM,JNEW,MNEW)
C----------------------------------------------------------------------
C     CHOOSE NEW BH J AND M_J USING CLEBSCH-GORDAN COEFFTS
C----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION JBH,JEM,MEM,JNEW,MNEW,SUMCG,RANCG,CHRGEN,CG
      EXTERNAL CHRGEN
C-----Set new J to minimum J allowed
      JNEW=JBH-JEM
C-----Set new M
      MNEW=JBH-MEM
C-----Initialise cumulative clebsch-gordan coefficients to zero
      SUMCG=0D0
C-----Generate the cumulant to be inverted
      RANCG=CHRGEN(1)
      DO WHILE (JNEW.LT.JBH+JEM+0.1D0)
C-----Calculate cumulants until larger than RANCG
         CALL CHSCGC(JBH,JEM,JNEW,JBH,-MEM,MNEW,CG)
         SUMCG=SUMCG+CG*CG
         IF (SUMCG.GT.RANCG) RETURN
         JNEW=JNEW+1D0
      ENDDO
C-----Use last one if (maybe because of numerical errors) the squares of the CG coefficients don't add up to 1.
      JNEW=JBH+JEM
      PRINT *,'CHNEWJ: NO CONVERGENCE, JBH,JNEW =',JBH,JNEW
      END

CDECK  ID>, CHFINDK.
*CMZ :-        -24/01/08  15.58.30 by  Marco Sampaio
*-- Author :    Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHFINDK(CFLUX2SK1,CFLUX2SK2,MK,ASTAR,FRACDA,K)
C-----------------------------------------------------------------------
C     Selects the mode K. K counts the modes ordered in JEM,MEM where 
C     we have them organised in sets of modes for JEM increasing (one 
C     set after the other) and within each set we have modes ordered in 
C     MEM from -JEM to JEM.
C     K is selected by looking at a cumulative table in K, where we have 
C     integrated over energy.
C-----------------------------------------------------------------------
C  *INPUT:
C     --CFLUX2SK1,CFLUX2SK2: arrays with cumulative fluxes from the mode 
C       K=1 up to MK (number of modes available), for IASTAR and 
C       IASTAR+1 respectively, so that we can interpolate in astar.
C     --ASTAR,FRACDA: the true astar and the interpolation weight 
C       respectively.
C  *OUTPUT:
C     --K: The selected mode.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Input/Output:
      INTEGER K,MK
      DOUBLE PRECISION CFLUX2SK1(MK+1),CFLUX2SK2(MK+1),FRACDA,ASTAR
C--Local variables:
      INTEGER ISEARCH,KLOW,KHIGH,KTRY
      DOUBLE PRECISION RNK,CINT,CHRGEN
      EXTERNAL CHRGEN

C--Generates the interpolated cumulative flux to be inverted 
      RNK=CHRGEN(1)*(CFLUX2SK1(MK)*(1-FRACDA)+CFLUX2SK2(MK)*FRACDA)
C--Set limits of K for its search
      KLOW=0
      KHIGH=MK
C--Bracketing
      DO 1 ISEARCH=1,MK
C-----Checks if the true point is in the first or second half of the 
C     interval, by calculating the cumulant in the middle point.
         KTRY=(KLOW+KHIGH)/2
         CINT=CFLUX2SK1(KTRY)*(1-FRACDA)+CFLUX2SK2(KTRY)*FRACDA
C-----Reduces the size of the interval by 1/2 
         IF(CINT.GT.RNK)THEN
            KHIGH=KTRY
         ELSE
            KLOW=KTRY
         END IF
         IF((KHIGH-KLOW).EQ.1)THEN
C--------K has been found!
            GOTO 2
         END IF  
 1    CONTINUE

 2    CONTINUE
C--End of search. KHIGH is the correct mode.
      K=KHIGH

      END

CDECK  ID>, CHFINDE.
*CMZ :-        -24/01/08  15.58.30 by  Marco Sampaio
*-- Author :    Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHFINDE(CFLUX2SKE1,CFLUX2SKE2,ME,ASTAR,FRACDA,EMIN,EMAX
     &,CONVERGENCE,EEM)
C-----------------------------------------------------------------------
C     Selects the energy (EEM) within the allowed region [EMIN:EMAX]
C-----------------------------------------------------------------------
C  *INPUT:
C     --CFLUX2SKE1,CFLUX2SKE2: arrays with cumulative fluxes for the 
C       selected mode for IASTAR and IASTAR+1 for variable energy.
C     --ME: Number of points in energy.
C     --ASTAR,FRACDA: the true astar and the interpolation weight 
C       respectively.
C     --CONVERGENCE: Deals with some extremely rare cases of no 
C       convergence when a mode with zero cumulants (in the 
C       kinematically allowed region) is selected. 
C  *OUTPUT:
C     --EEM: The selected energy.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Input/Output
      INTEGER ME
      DOUBLE PRECISION CFLUX2SKE1(ME+1),CFLUX2SKE2(ME+1),ASTAR,FRACDA,
     &EMIN,EMAX,EEM
      LOGICAL CONVERGENCE
C--Local variables
      INTEGER ISEARCH,NELOW,NEHIGH,NETRY,IEMIN,IEMAX
      DOUBLE PRECISION RN,RNC,CINT,MAXETABLE,FRACDC,CHRGEN,
     &DENOMDC,NUMDC,DIEMIN,DIEMAX,DI,DE
      PARAMETER(MAXETABLE=5D0,DE=5D-2)
      EXTERNAL CHRGEN
      CONVERGENCE=.TRUE.
C--Sets bounds for the position of the energy on the tables based on 
C  kinematical restrictions, i.e. smaller than MAXETABLE and larger than 
C  EMIN. Note that because we are taking the integer part IEMIN can 
C  correspond to an energy slightly below EMIN and IEMAX to an energy 
C  above EMAX.
      IEMIN=1+INT(EMIN/MAXETABLE*ME)
      IEMAX=2+INT(EMAX/MAXETABLE*ME)
C--Restricts the search to the maximum energy available if EMAX is larger 
C  than MAXETABLE. (Note: In CFLUX2SKE1(or 2) we have include the point at 
C  zero energy and that's why we have ME+1 points)
      IF(IEMAX.GT.(ME+1)) THEN
         IEMAX=ME+1
         EMAX=MAXETABLE
      END IF
C--Generates the interpolated cumulative flux to be inverted
      RN=CHRGEN(1)
      RNC=(CFLUX2SKE1(IEMIN)*(1-FRACDA)+CFLUX2SKE2(IEMIN)*FRACDA)*(1-RN)
     &+RN*(CFLUX2SKE1(IEMAX)*(1-FRACDA)+CFLUX2SKE2(IEMAX)*FRACDA)
C--Set limits for the search
      NELOW=IEMIN
      NEHIGH=IEMAX
C--Bracketing
      DO 1 ISEARCH=1,ME+1
C-----Checks if the true point is in the first or second half of the 
C     interval, by calculating the cumulant in the middle point.
         NETRY=(NELOW+NEHIGH)/2
         CINT=CFLUX2SKE1(NETRY)*(1-FRACDA)+CFLUX2SKE2(NETRY)*FRACDA
C--------Reduce the size of the interval by 1/2
         IF(CINT.GT.RNC)THEN
            NEHIGH=NETRY
         ELSE
            NELOW=NETRY
         END IF
         IF((NEHIGH-NELOW).EQ.1)THEN
C-------Bracketing complete
            GOTO 2
         END IF        
 1    CONTINUE

 2    CONTINUE
C--Calculates the quantities used to determine the interpolation weight 
C  (interpolation in energy). The numerator and denominator are separated to 
C  take into account division by zero in flat regions or regions of zero cumulant
      DENOMDC=(CFLUX2SKE1(NEHIGH)-CFLUX2SKE1(NEHIGH-1))*(1-FRACDA)+
     &(CFLUX2SKE2(NEHIGH)-CFLUX2SKE2(NEHIGH-1))*FRACDA
      NUMDC=RNC-(CFLUX2SKE1(NEHIGH-1)*(1-FRACDA)+CFLUX2SKE2(NEHIGH-1)*
     &FRACDA)

C--Calculate energy:
      IF((DENOMDC.NE.0).AND.(RNC.NE.0))THEN
C-----If everything OK (not a flat region or zero cumulant) then interpolates 
C     linearly and returns the energy
         FRACDC=NUMDC/DENOMDC
         EEM=(NEHIGH-1D0)*MAXETABLE/ME*FRACDC+(NEHIGH-2D0)*MAXETABLE/ME
     &*(1-FRACDC)
      ELSE IF((DENOMDC.EQ.0).AND.(RNC.NE.0)) THEN
C-----If the region selected is flat, returns some energy in the interval that
C     has been bracketed with uniform probability
         FRACDC=CHRGEN(2)
         EEM=(NEHIGH-1D0)*MAXETABLE/ME*FRACDC+(NEHIGH-2D0)*MAXETABLE/ME
     &*(1-FRACDC) 
      ELSE
C-----If region of zero probability (i.e. cumulant), return no convergence
         CONVERGENCE=.FALSE.
      END IF
C--If the energy selected is outside the kinematically allowed region 
C  (because the integers IEMAX, IEMIN don't exactly coincide with the
C  maximum and minimum energy) then return an energy with uniform probability
C  either around IEMAX or IEMIN if the kinematically allowed region is big 
C  enough. If the allowed region is smaller than the intervals around IEMAX or
C  IEMIN, then retun an energy with uniform probability in that region.    
      IF(EEM.GT.EMAX)THEN
C-----Energy slightly above EMAX case
         FRACDC=CHRGEN(4)
C-------Distance above EMAX:
         DIEMAX=IEMAX-1D0-1.*EMAX/MAXETABLE*ME
C-------Size of the kinematically allowed region
         DI=(EMAX-EMIN)*1./MAXETABLE*ME
         IF(DIEMAX.LT.DI)THEN
C----------Select energy in a small interval below EMAX
            EEM=EMAX-DIEMAX*FRACDC*DE
         ELSE
C----------Select energy in the narrow kinematically allowed region
            EEM=EMIN*(1-FRACDC)+EMAX*FRACDC        
         END IF
      ELSE IF(EEM.LT.EMIN)THEN
         FRACDC=CHRGEN(4)
C-------Distance below EMIN:
         DIEMIN=1D0+EMIN/MAXETABLE*ME-IEMIN
C-------Size of the kinematically allowed region
         DI=(EMAX-EMIN)*1./MAXETABLE*ME
         IF(DIEMIN.LT.DI)THEN
C----------Select energy in a small interval above EMIN
            EEM=EMIN+DIEMIN*FRACDC*DE 
         ELSE
C----------Select energy in the narrow kinematically allowed region
            EEM=EMIN*(1-FRACDC)+EMAX*FRACDC        
         END IF
      END IF

      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                     Remnant decay models Routines
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, CHDFIV.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Ian Knowles
C-----------------------------------------------------------------------
      SUBROUTINE CHDFIV(P0,P1,P2,P3,P4,P5,OKDEC)
C-----------------------------------------------------------------------
C     Generates 5-body decay 0->1+2+3+4+5 using pure phase space
C     (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRGEN,P0(5),P1(5),P2(5),P3(5),P4(5),P5(5),B,C,
     & AA,BB,CC,DD,EE,FF,TT,S1,RS1,GG,S2,RS2,HH,S3,PP,QQ,RR,SS,P1CM,
     & P2345(5),P2CM,P345(5),P3CM,P45(5),P4CM
      DOUBLE PRECISION TWO
      INTEGER MTRY,NTRY
      PARAMETER (MTRY=1000000)
      PARAMETER (TWO=2.D0)
      EXTERNAL CHRGEN
      LOGICAL OKDEC
      OKDEC=.FALSE.
      B=P0(5)-P1(5)
      C=P2(5)+P3(5)+P4(5)+P5(5)
      IF (B.LT.C) RETURN
      AA=(P0(5)+P1(5))**2
      BB=B**2
      CC=C**2
      DD=(P3(5)+P4(5)+P5(5))**2
      EE=(P4(5)+P5(5))**2
      FF=(P4(5)-P5(5))**2
      TT=(B-C)*P0(5)**11/729
      NTRY=0
C Select squared masses S1, S2 and S3 of 2345, 345 and 45 subsystems
  10  S1=BB+CHRGEN(1)*(CC-BB)
      NTRY=NTRY+1
      IF (NTRY.GT.MTRY) RETURN
      RS1=SQRT(S1)
      GG=(RS1-P2(5))**2
      S2=DD+CHRGEN(2)*(GG-DD)
      RS2=SQRT(S2)
      HH=(RS2-P3(5))**2
      S3=EE+CHRGEN(3)*(HH-EE)
      PP=(AA-S1)*(BB-S1)
      QQ=((RS1+P2(5))**2-S2)*(GG-S2)/S1
      RR=((RS2+P3(5))**2-S3)*(HH-S3)/S2
      SS=(S3-EE)*(S3-FF)/S3
      IF (PP*QQ*RR*SS*((GG-DD)*(HH-EE))**2.LT.TT*S1*S2*S3*CHRGEN(4)**2)
     & GOTO 10
C Do two body decays: 0-->1+2345, 2345-->2+345, 345-->3+45 and 45-->4+5
      P1CM=SQRT(PP/4)/P0(5)
      P2345(5)=RS1
      P2CM=SQRT(QQ/4)
      P345(5)=RS2
      P3CM=SQRT(RR/4)
      P45(5)=SQRT(S3)
      P4CM=SQRT(SS/4)
      CALL CHDTWO(P0   ,P1,P2345,P1CM,TWO,.TRUE.)
      CALL CHDTWO(P2345,P2,P345 ,P2CM,TWO,.TRUE.)
      CALL CHDTWO(P345 ,P3,P45  ,P3CM,TWO,.TRUE.)
      CALL CHDTWO(P45  ,P4,P5   ,P4CM,TWO,.TRUE.)
      OKDEC=.TRUE.
      END
CDECK  ID>, CHDFOR.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Ian Knowles
C-----------------------------------------------------------------------
      SUBROUTINE CHDFOR(P0,P1,P2,P3,P4,OKDEC)
C-----------------------------------------------------------------------
C     Generates 4-body decay 0->1+2+3+4 using pure phase space
C     (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRGEN,P0(5),P1(5),P2(5),P3(5),P4(5),B,C,AA,BB,
     & CC,DD,EE,TT,S1,RS1,FF,S2,PP,QQ,RR,P1CM,P234(5),P2CM,P34(5),P3CM
      DOUBLE PRECISION TWO
      INTEGER MTRY,NTRY
      PARAMETER (MTRY=1000000)
      PARAMETER (TWO=2.D0)
      EXTERNAL CHRGEN
      LOGICAL OKDEC
      OKDEC=.FALSE.
      B=P0(5)-P1(5)
      C=P2(5)+P3(5)+P4(5)
      IF (B.LT.C) RETURN
      AA=(P0(5)+P1(5))**2
      BB=B**2
      CC=C**2
      DD=(P3(5)+P4(5))**2
      EE=(P3(5)-P4(5))**2
      TT=(B-C)*P0(5)**7/16
      NTRY=0
C Select squared masses S1 and S2 of 234 and 34 subsystems
  10  S1=BB+CHRGEN(1)*(CC-BB)
      NTRY=NTRY+1
      IF (NTRY.GT.MTRY) RETURN
      RS1=SQRT(S1)
      FF=(RS1-P2(5))**2
      S2=DD+CHRGEN(2)*(FF-DD)
      PP=(AA-S1)*(BB-S1)
      QQ=((RS1+P2(5))**2-S2)*(FF-S2)/S1
      RR=(S2-DD)*(S2-EE)/S2
      IF (PP*QQ*RR*(FF-DD)**2.LT.TT*S1*S2*CHRGEN(3)**2) GOTO 10
C Do two body decays: 0-->1+234, 234-->2+34 and 34-->3+4
      P1CM=SQRT(PP/4)/P0(5)
      P234(5)=RS1
      P2CM=SQRT(QQ/4)
      P34(5)=SQRT(S2)
      P3CM=SQRT(RR/4)
      CALL CHDTWO(P0  ,P1,P234,P1CM,TWO,.TRUE.)
      CALL CHDTWO(P234,P2,P34 ,P2CM,TWO,.TRUE.)
      CALL CHDTWO(P34 ,P3,P4  ,P3CM,TWO,.TRUE.)
      OKDEC=.TRUE.
      END
CDECK  ID>, CHDTHR.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHDTHR(P0,P1,P2,P3,OKDEC)
C-----------------------------------------------------------------------
C     GENERATES THREE-BODY DECAY 0->1+2+3 DISTRIBUTED
C     ACCORDING TO PHASE SPACE * WEIGHT
C     (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRGEN,CHRUNI,A,B,C,D,AA,BB,CC,DD,EE,FF,PP,QQ,WW,
     & RR,PCM1,PC23,P0(5),P1(5),P2(5),P3(5),P23(5),TWO,ONE
      EXTERNAL CHRGEN,CHRUNI
      PARAMETER (ONE=1.0D0,TWO=2.D0)
      INTEGER MTRY,NTRY
      PARAMETER (MTRY=1000000)
      LOGICAL OKDEC
      OKDEC=.FALSE.
      A=P0(5)+P1(5)
      B=P0(5)-P1(5)
      C=P2(5)+P3(5)
      IF (B.LT.C) RETURN
      D=ABS(P2(5)-P3(5))
      AA=A*A
      BB=B*B
      CC=C*C
      DD=D*D
      EE=(B-C)*(A-D)
      A=0.5D0*(AA+BB)
      B=0.5D0*(CC+DD)
      C=4.0D0/(A-B)**2
      NTRY=0
C
C  CHOOSE MASS OF SUBSYSTEM 23 WITH PRESCRIBED DISTRIBUTION
C
   10 FF=CHRUNI(0,BB,CC)
      NTRY=NTRY+1
      IF (NTRY.GT.MTRY) RETURN
      PP=(AA-FF)*(BB-FF)
      QQ=(CC-FF)*(DD-FF)
      WW=ONE
      RR=EE*FF*CHRGEN(1)
      IF (PP*QQ*WW.LT.RR*RR) GOTO 10
C
C  FF IS MASS SQUARED OF SUBSYSTEM 23.
C
C  DO 2-BODY DECAYS 0->1+23, 23->2+3
C
      P23(5)=SQRT(FF)
      PCM1=SQRT(PP)*0.5D0/P0(5)
      PC23=SQRT(QQ)*0.5D0/P23(5)
      CALL CHDTWO(P0,P1,P23,PCM1,TWO,.TRUE.)
      CALL CHDTWO(P23,P2,P3,PC23,TWO,.TRUE.)
      OKDEC=.TRUE.
      END
CDECK  ID>, CHDTWO.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber & Mike Seymour
C-----------------------------------------------------------------------
      SUBROUTINE CHDTWO(P0,P1,P2,PCM,COSTH,ZAXIS)
C-----------------------------------------------------------------------
C     GENERATES DECAY 0 -> 1+2
C
C     PCM IS CM MOMENTUM
C
C     COSTH = COS THETA IN P0 REST FRAME (>1 FOR ISOTROPIC)
C     IF ZAXIS=.TRUE., COS THETA IS MEASURED FROM THE ZAXIS
C     IF .FALSE., IT IS MEASURED FROM P0'S DIRECTION 
C     (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRUNI,ONE,ZERO,PCM,COSTH,C,S,P0(5),P1(5),P2(5),
     & PP(5),R(9)
      LOGICAL ZAXIS
      EXTERNAL CHRUNI
      PARAMETER (ZERO=0.D0, ONE=1.D0)
C--CHOOSE C.M. ANGLES
      C=COSTH
      IF (C.GT.ONE) C=CHRUNI(0,-ONE,ONE)
      S=SQRT(ONE-C*C)
      CALL CHRAZM(PCM*S,PP(1),PP(2))
C--PP IS MOMENTUM OF 2 IN C.M.
      PP(3)=-PCM*C
      PP(4)=SQRT(P2(5)**2+PCM**2)
      PP(5)=P2(5)
C--ROTATE IF NECESSARY
      IF (COSTH.LE.ONE.AND..NOT.ZAXIS) THEN
        CALL CHUROT(P0,ONE,ZERO,R)
        CALL CHUROB(R,PP,PP)
      ENDIF
C--BOOST FROM C.M. TO LAB FRAME
      CALL CHULOB(P0,PP,P2)
      CALL CHVDIF(4,P0,P2,P1)
      END

CDECK  ID>, .
*CMZ :-        -27/04/08  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHVARNBODY(NAVERAGE,N)
C-----------------------------------------------------------------------
C     Returns a nonzero final multiplicity which is determined by the 
C     nonzero N region of a Poisson distribution
C-----------------------------------------------------------------------
C   I/O: 
C     NAVERAGE-> The average multiplicity
C     N-> Output multiplicity chosen according to the Poisson ditribution
C-----------------------------------------------------------------------
*Fixed parameters:
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER(NMAX=15)
*I/0:
      INTEGER N
      DOUBLE PRECISION NAVERAGE
*Local variables:
      INTEGER I
      DOUBLE PRECISION CHRGEN,RCUMULATIVE,CUMULATIVE,PROB
      EXTERNAL CHRGEN

      PROB=1D0
C--COMPUTES THE RANDOM NUMBER TIMES THE NORMALIZATION FACTOR
      RCUMULATIVE=(EXP(NAVERAGE)-1)*CHRGEN(1)
C--INITIALISES CUMULATIVE FUNCTION TO ZERO
      CUMULATIVE=0D0
      I=1
C--CALCULATES CUMULATIVE FUNTION AND WHENEVER IT IS LARGER THAN RCUMULATIVE IT SELECTS THAT I.
      DO WHILE((I.LT.NMAX).AND.(CUMULATIVE.LT.RCUMULATIVE))
         PROB=PROB*NAVERAGE/I
         CUMULATIVE=CUMULATIVE+PROB
         I=I+1
      END DO

C--RETURNS THE MULTIPLICITY + 1 BECAUSE 1 EMISSION MEANS REALLY 2-BODY FINAL DECAY
      N=I
 
      END

CDECK  ID>, CHEMIT.
*CMZ :-        -07/07/08  15.58.30  by  Marco Sampaio
*-- Author :    Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHJSELECTNBODY(NEX,ASTAR,JBH,SPIN,JEM,MEM,NBODYNOW)
C-----------------------------------------------------------------------
C     Selects a number (NBODYNOW) of partial wave numbers for the
C     emissions in the N-body final decay
C-----------------------------------------------------------------------
C  *INPUT:
C     --NEX,ASTAR,JBH: are number of extra dimensions, oblateness/rh and 
C     angular momentum of the black hole 
C     --SPIN: spin of the particle emitted.
C     ---NBODYNOW: Number of particles in the final decay
C  *OUTPUT:
C     --JEM,MEM: The Angular momenta of the particles 
C       emitted.
C----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
C--Input/Output:
      INTEGER NEX,NBODYNOW
      INTEGER SPIN(1:NBODYNOW)
      DOUBLE PRECISION ASTAR,JEM(1:NBODYNOW),MEM(1:NBODYNOW),JBH
C--Local variables
      INTEGER NBODYTRYMAX,IASTAR,K,MA,ME,MK0,MK1,MK2,UCFLUXKE,ITRY,
     &NTRY,I,IJ1,IJ2
      DOUBLE PRECISION JMK0,JMK1,JMK2
      PARAMETER (MA=26,MK0=169,MK1=156,MK2=120,ME=100,JMK0=12,
     &JMK1=11.5D0,JMK2=10,NBODYTRYMAX=50000)
      DOUBLE PRECISION CHRGEN,CFLUX2S0K(MK0+1,MA),
     &CFLUX2S1K(MK1+1,MA),CFLUX2S2K(MK2+1,MA),CFLUX2S0KE(ME+1,MK0,MA),
     &CFLUX2S1KE(ME+1,MK1,MA),CFLUX2S2KE(ME+1,MK2,MA),FRACDA,ASTARMAX,
     &JMIN,JMAX,MNOW
      PARAMETER(ASTARMAX=5D0)
      LOGICAL FIRST
      EXTERNAL CHRGEN
      CHARACTER NEXTRA*1

      DATA UCFLUXKE,FIRST/23,.TRUE./
      SAVE CFLUX2S0K,CFLUX2S1K,CFLUX2S2K,CFLUX2S0KE,CFLUX2S1KE,
     &CFLUX2S2KE

C--Initialises the arrays with cumulative fluxes by reading the 
C  corresponding data file for a number of extra dimensions NEX This 
C  occurs only once at FIRST call.
      IF (FIRST) THEN
         NEXTRA=CHAR(ICHAR('0')+NEX)
         IF ((NEX.GT.0).AND.(NEX.LT.7))
     &THEN
            OPEN(UNIT=UCFLUXKE,FILE='./data_files/cfluxke_n'//NEXTRA//
     &'.dat',STATUS='OLD')
         ELSE
            PRINT *,' NO CUMULATIVE TABLES YET FOR NEX =',NEX
            STOP
         ENDIF
         READ(UCFLUXKE,*)CFLUX2S0K,CFLUX2S1K,CFLUX2S2K,CFLUX2S0KE,
     &CFLUX2S1KE,CFLUX2S2KE
         CLOSE(UNIT=UCFLUXKE)

         FIRST=.FALSE.
      ENDIF

C--Select the correct IASTAR (the position of the astar imediately below 
C  the actual one) in the flux tables. Interpolation will be done between 
C  IASTAR and IASTAR+1.
      IF (ASTAR.LE.ASTARMAX) THEN
         IASTAR=1+INT(ASTAR/ASTARMAX*(MA-1D0))
C------- FRACDA is the linear interpolation weight  at IASTAR+1 and 
C        1-FRACDA at IASTAR
         FRACDA=ASTAR*(MA-1.)/ASTARMAX-(IASTAR-1)  
      ELSE
         IASTAR=MA-1
         FRACDA=1D0  
      END IF

C--Select the modes K for each of the emissions in the final decay
      DO 5 ITRY=0,NBODYTRYMAX
         JMIN=0D0
         JMAX=0D0
         MNOW=0D0
         DO 6 NTRY=1,NBODYNOW-1
C---------Find all K except the last one by generating cumulative functions distributed uniformly.
            IF(SPIN(NTRY).EQ.0)THEN
C-----Scalars
               CALL CHFINDK(CFLUX2S0K(1,IASTAR),CFLUX2S0K(1,IASTAR+1),
     &MK0,ASTAR,FRACDA,K)
               JEM(NTRY)=INT(SQRT(K-1D-2))
               MEM(NTRY)=K-JEM(NTRY)*(JEM(NTRY)+1)-1
c--------------keep track of the allowed ranges of quantum numbers for the combined angular momentum using the usual rules for the decomposition of the tensor product representation into irreducible representations.
               IF((JEM(NTRY).LT.JMAX).AND.(JEM(NTRY).GT.JMIN))THEN 
                  JMIN=(JEM(NTRY)-JMIN)-INT(JEM(NTRY)-JMIN)
               ELSE IF(JEM(NTRY).GE.JMAX)THEN
                  JMIN=JEM(NTRY)-JMAX
               ELSE
                  JMIN=JMIN-JEM(NTRY)
               END IF
               JMAX=JMAX+JEM(NTRY)
c--------------current total M for all these waves
               MNOW=MNOW+MEM(NTRY)
            ELSE IF(SPIN(NTRY).EQ.1)THEN
C-----Fermions
               CALL CHFINDK(CFLUX2S1K(1,IASTAR),CFLUX2S1K(1,IASTAR+1),
     &MK1,ASTAR,FRACDA,K)
               JEM(NTRY)=INT(SQRT(K-1D-2)-0.5D0)+0.5D0
               MEM(NTRY)=K-JEM(NTRY)*(JEM(NTRY)+1D0)-0.75D0
c--------------keep track of the allowed ranges of quantum numbers for the combined angular momentum using the usual rules for the decomposition of the tensor product representation into irreducible representations.
               IF((JEM(NTRY).LT.JMAX).AND.(JEM(NTRY).GT.JMIN))THEN 
                  JMIN=(JEM(NTRY)-JMIN)-INT(JEM(NTRY)-JMIN)
               ELSE IF(JEM(NTRY).GE.JMAX)THEN
                  JMIN=JEM(NTRY)-JMAX
               ELSE
                  JMIN=JMIN-JEM(NTRY)
               END IF
               JMAX=JMAX+JEM(NTRY)
c--------------current total M for all these waves
               MNOW=MNOW+MEM(NTRY)
            ELSE IF(SPIN(NTRY).EQ.2)THEN
C-----Vectors
               CALL CHFINDK(CFLUX2S2K(1,IASTAR),CFLUX2S2K(1,IASTAR+1),
     &MK2,ASTAR,FRACDA,K)
               JEM(NTRY)=INT(SQRT(K+1D-2))
               MEM(NTRY)=K-JEM(NTRY)*(JEM(NTRY)+1D0)
c--------------keep track of the allowed ranges of quantum numbers for the combined angular momentum using the usual rules for the decomposition of the tensor product representation into irreducible representations.               
               IF((JEM(NTRY).LT.JMAX).AND.(JEM(NTRY).GT.JMIN))THEN 
                  JMIN=(JEM(NTRY)-JMIN)-INT(JEM(NTRY)-JMIN)
               ELSE IF(JEM(NTRY).GE.JMAX)THEN
                  JMIN=JEM(NTRY)-JMAX
               ELSE
                  JMIN=JMIN-JEM(NTRY)
               END IF
               JMAX=JMAX+JEM(NTRY)

c--------------current total M for all these waves
               MNOW=MNOW+MEM(NTRY)
            END IF
 6       CONTINUE
C----Sets MNOW to the M that the last particle must have
         MNOW=JBH-MNOW
         IF(ABS(MNOW).GT.JMIN) JMIN=ABS(MNOW)
c----Try again if the wave is not allowed (either wrong statistics or JMAX too low) 
         IF((SPIN(NBODYNOW)/2D0-ABS(MNOW)-INT(SPIN(NBODYNOW)/2D0-
     &ABS(MNOW)).NE.0).OR.(JMAX.LT.SPIN(NBODYNOW)/2D0)
     &.OR.(JMAX.LT.JMIN)) THEN
            GOTO 4
         END IF
         MEM(NBODYNOW)=MNOW

c-----Find the last J by doing hit and miss
         IF(SPIN(NBODYNOW).EQ.0)THEN
C-----Scalars
            CALL CHFINDJNBODY(CFLUX2S0K(1,IASTAR),
     &CFLUX2S0K(1,IASTAR+1),MK0,ASTAR,FRACDA,JEM(NBODYNOW),MEM(NBODYNOW)
     &,JMIN,JMAX,SPIN(NBODYNOW),JMK0)
         ELSE IF(SPIN(NBODYNOW).EQ.1)THEN
C-----Fermions
            CALL CHFINDJNBODY(CFLUX2S1K(1,IASTAR),
     &CFLUX2S1K(1,IASTAR+1),MK1,ASTAR,FRACDA,JEM(NBODYNOW),MEM(NBODYNOW)
     &,JMIN,JMAX,SPIN(NBODYNOW),JMK1)
         ELSE IF(SPIN(NBODYNOW).EQ.2)THEN
C-----Vectors
            CALL CHFINDJNBODY(CFLUX2S2K(1,IASTAR),
     &CFLUX2S2K(1,IASTAR+1),MK2,ASTAR,FRACDA,JEM(NBODYNOW),MEM(NBODYNOW)
     &,JMIN,JMAX,SPIN(NBODYNOW),JMK2)
         END IF
c-----End search if successful, otherwise try again
         IF(JEM(NBODYNOW).NE.-1) GOTO 7

 4       CONTINUE

 5    CONTINUE

C---If the number of tries expires then we select the one of the lowest total J and lowest individual Js states at random
      JMAX=0D0
c-----compute total lowest J      
      DO 8 I=1,NBODYNOW
         JMAX=JMAX+SPIN(I)/2D0
 8    CONTINUE

c----- IF JMAX too small, compute the difference to JMAX and distribute it equally among all particles. 
      IF(JMAX.LT.JBH)THEN
c---------Difference to distribute among all particles
         IJ1=(JBH-JMAX)/NBODYNOW
c---------Number of the last particle to give the last amount+1
         IJ2=INT(JBH-JMAX)-IJ1*NBODYNOW
c---------Give the first IJ2 particles IJ1+1 units of angular momentum
         DO 80 I=1,IJ2
            JEM(I)=SPIN(I)/2D0+IJ1+1
            MEM(I)=JEM(I)
 80      CONTINUE
c---------Give all the other particles IJ1 units of angualar momentum
         DO 81 I=(IJ2+1),NBODYNOW
            JEM(I)=SPIN(I)/2D0+IJ1
            MEM(I)=JEM(I)
 81      CONTINUE  
      ELSE
         IJ1=INT(JMAX-JBH)
c---------set all the waves to M=J
         DO 90 I=1,NBODYNOW
            JEM(I)=SPIN(I)/2D0
            MEM(I)=JEM(I)
 90      CONTINUE
         IF(IJ1.NE.0)THEN
c---------reduce M by one for all particles with non-zero spin (fermion+vectors)until Mtotal = JBH
            DO 91 I=1,NBODYNOW
               IF((SPIN(I).GT.0D0).AND.(IJ1.NE.0)) THEN
                  JEM(I)=SPIN(I)/2D0
                  MEM(I)=JEM(I)-1D0
                  IJ1=IJ1-1
               END IF
 91         CONTINUE
         END IF
c---------reduce M by one for vectors until Mtotal = JBH if not yet there from the previous loop
         IF(IJ1.NE.0)THEN
            DO 92 I=1,NBODYNOW
               IF((SPIN(I).GT.1D0).AND.(IJ1.NE.0)) THEN
                  JEM(I)=SPIN(I)/2D0
                  MEM(I)=JEM(I)-1D0
                  IJ1=IJ1-1
               END IF
 92         CONTINUE
         END IF
      END IF
      
 7    CONTINUE

      END

CDECK  ID>, CHPSELECTNBODY.
*CMZ :-        -27/07/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHPSELECTNBODY(MBH,OBLA,RH,NEXTRA,JEM,MEM,SPINTIMES2,
     &PNBODY,NBODYNOW,CONVERGENCE,FPOLARIZ,PT)
C-----------------------------------------------------------------------
C     Selects the momenta of the N-body final decay particles
C-----------------------------------------------------------------------
C  *I/O:
C     --MBH,OBLA,RH,NEXTRA: BH parameters plus number of extra dimensions
C     --JEM,MEM,SPINTIMES2,PNBODY: Quantum numbers of emissions
C     --NBODYNOW: Number of particles in the N-body decay
C     --CONVERGENCE: Logical variable to determine is convergence occurred
C     --FPOLARIZ: Polarizations of the emission
C     --PT: Particle types for the emissions
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Input/Output
      INTEGER NBODYNOW,NEXTRA,FPOLARIZ(NBODYNOW)
      DOUBLE PRECISION MBH,OBLA,RH,JEM(NBODYNOW),MEM(NBODYNOW),
     &PNBODY(5,NBODYNOW)
      INTEGER SPINTIMES2(NBODYNOW),PT(NBODYNOW)
C--Local variables
C----Fixed parameters:
      DOUBLE PRECISION ASTARMAX
      INTEGER NMAX,NDIVIDEMAX,MAXSEARCH,MAXTRY,I,J,MA
      PARAMETER(NMAX=250,NDIVIDEMAX=150,MAXSEARCH=500,MAXTRY=5000,MA=26,
     &ASTARMAX=5D0)
C----Other variables
      INTEGER LTIMES2(NBODYNOW),MTIMES2(NBODYNOW),NLAST1(NBODYNOW),
     &NLAST2(NBODYNOW),NYSIZE(NBODYNOW),IASTAR,K(NBODYNOW),
     &HTIMES2(NBODYNOW)
      DOUBLE PRECISION CHRGEN,COSTHETA,SPHRDLSERIES1(0:NMAX,NBODYNOW),
     &CONSTANTS1(1:4,NBODYNOW),DELTAN1(NBODYNOW),
     &SPHRDLSERIES2(0:NMAX,NBODYNOW),CONSTANTS2(1:4,NBODYNOW),
     &DELTAN2(NBODYNOW),A(NBODYNOW),Y(0:NDIVIDEMAX+1),
     &XOFY(1:NDIVIDEMAX+1),DISTMAX,PI,CNBODY(NBODYNOW),MU,MODPEM,
     &PTRANSVERSE,ASTAR,PQTRANSFER(5,NBODYNOW),FRACDA,ERH(NBODYNOW),
     &DISTVAL,MUI(0:NBODYNOW),RANDOM(NBODYNOW),
     &Y2(0:NDIVIDEMAX+1,NBODYNOW),
     &NORM(NBODYNOW),XNOW(NBODYNOW),ENBODYQ(NBODYNOW)
      LOGICAL CONVERGENCE
      EXTERNAL CHRGEN
      PI=ACOS(-1D0)
      CONVERGENCE=.TRUE.

      DO J=1,NBODYNOW
         LTIMES2(J)=INT(2*JEM(J))
         MTIMES2(J)=INT(2*MEM(J))
      END DO
c-----Use cumulative function to select orientation of the momenta in the cm frame if 2-body decay
      IF(NBODYNOW.EQ.2)THEN
C-----POLARIZATIONS
         DO J=1,2
            FPOLARIZ(J)=1
            IF(CHRGEN(1).LT.0.5)FPOLARIZ(J)=-1
C-----Neutrinos (left handed) and antineutrinos (right handed)
          IF(PT(J).EQ.12.OR.PT(J).EQ.14.OR.PT(J).EQ.16) FPOLARIZ(J)=-1
          IF(PT(J).EQ.-12.OR.PT(J).EQ.-14.OR.PT(J).EQ.-16) FPOLARIZ(J)=1
         END DO
c-----Magnitude of momenta fixed by conservation
         PNBODY(4,1)=MBH/2D0+(PNBODY(5,1)**2-PNBODY(5,2)**2)/(2D0*MBH)
         PNBODY(4,2)=MBH/2D0+(PNBODY(5,2)**2-PNBODY(5,1)**2)/(2D0*MBH)
         MODPEM=SQRT(PNBODY(4,1)**2-PNBODY(5,1)**2)
c-----Compute spheroidal series for each particle in the n-body decay
         DO 3 J=1,2
C---------Compute HTIMES2 (basically SPINTIMES2 with a sign according to polarization) in a counterintuitive convention (because it is the one used in the spheroidal functions).
            HTIMES2(J)=-SPINTIMES2(J)*FPOLARIZ(J)
            CNBODY(J)=OBLA*PNBODY(4,J)
C---------Avoid numerically unstable areas of phase space (which should be very rare anyway) by imposing a maximum C
            IF(CNBODY(J).GT.50D0)CNBODY(J)=50D0
            CALL CHSPHRDLINI(LTIMES2(J),MTIMES2(J),HTIMES2(J),
     &CNBODY(J),SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),NLAST1(J),
     &SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J))
            CALL CHSPHRDLMATCH(CNBODY(J),MTIMES2(J),HTIMES2(J),
     &SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),NLAST1(J),
     &SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J),A(J))
 3       CONTINUE
c------Compute cumulative function for this case
         CALL CHSPHRDLCUMULATIVE2(CNBODY,LTIMES2,MTIMES2,
     &HTIMES2,SPHRDLSERIES1,CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,
     &CONSTANTS2,DELTAN2,NLAST2,A,Y,XOFY,NYSIZE)
c------Select orientation
         CALL CHSPHRDLXRAND(Y,XOFY,NYSIZE,COSTHETA)
c------Determine momentum of particle 1
         PNBODY(3,1)=MODPEM*COSTHETA
         PTRANSVERSE=SQRT(MODPEM**2-PNBODY(3,1)**2)
         CALL CHRAZM(PTRANSVERSE,PNBODY(1,1),PNBODY(2,1))
c------Determine the momentum of particle 2 (same magnitude opposite direction)
         PNBODY(1,2)=-PNBODY(1,1)
         PNBODY(2,2)=-PNBODY(2,1)
         PNBODY(3,2)=-PNBODY(3,1)    
      ELSE
c------Use hit and miss for a chain with N-particles with momentum conservation
         ASTAR=OBLA/RH
         MUI(0)=0D0
         DO 4 J=1,NBODYNOW
c------Sum of masses up to J
            MUI(J)=MUI(J-1)+PNBODY(5,J)
 4       CONTINUE
c------Availlable energy MU
         MU=MBH-MUI(NBODYNOW)
c------Initialise K, IASTAR and FRACDA
         DO 5 J=1,NBODYNOW
            K(J)=INT((JEM(J)+SPINTIMES2(J)/2D0)*(JEM(J)+1-
     &SPINTIMES2(J)/2D0)+MEM(J)+1-SPINTIMES2(J)/2D0)
            IF (ASTAR.LE.5D0) THEN
               IASTAR=1+INT(ASTAR/ASTARMAX*(MA-1D0))
C------- FRACDA is the linear interpolation weight  at IASTAR+1 and 
C        1-FRACDA at IASTAR
               FRACDA=ASTAR*(MA-1D0)/ASTARMAX-(IASTAR-1)
            ELSE
               IASTAR=MA-1
               FRACDA=1D0  
            END IF       
 5       CONTINUE

c-----Search for the maximum
C-------Initialise last momentum transfer to BH momentum
         PQTRANSFER(5,NBODYNOW)=MBH
         PQTRANSFER(4,NBODYNOW)=MBH
         PQTRANSFER(3,NBODYNOW)=0D0
         PQTRANSFER(2,NBODYNOW)=0D0
         PQTRANSFER(1,NBODYNOW)=0D0
         DISTMAX=0D0
c--------do search
         DO 1 I=0,MAXSEARCH
            RANDOM(1)=0D0
            RANDOM(NBODYNOW)=1D0
c-----------Generate NODYNOW-2 random numbers
            DO J=2,NBODYNOW-1
               RANDOM(J)=CHRGEN(2)
            END DO
c-----------Sort them in increasing order
            CALL CHSPHRDLQCKSRT(RANDOM,NBODYNOW)
            DO J=NBODYNOW,2,-1
c--------------Calculate invariant mass of the next momentum transfer
               PQTRANSFER(5,J-1)=MUI(J-1)+RANDOM(J-1)*MU
c--------------Energy of the Jth particle in the frame of the momentum transfer  
               ENBODYQ(J)=(PQTRANSFER(5,J)**2-PQTRANSFER(5,J-1)**2+
     &PNBODY(5,J)**2)/(2D0*PQTRANSFER(5,J))
c--------------Two body decay of the Jth momentum transfer into J-1th momentum transfer and the Jth particle
               CALL CHDTWO(PQTRANSFER(1,J),PNBODY(1,J),PQTRANSFER(1,J-1)
     &,SQRT(ENBODYQ(J)**2-PNBODY(5,J)**2),2D0,.TRUE.)
            END DO
c-----------Set last momentum transfer to the last emission
            CALL CHVEQU(4,PQTRANSFER(1,1),PNBODY(1,1))

            DO J=1,NBODYNOW
c-----------Computes input C and X for the momenta we are trying and choose polarisations at random
C-----POLARIZATIONS
               FPOLARIZ(J)=1
               IF(CHRGEN(3).LT.0.5)FPOLARIZ(J)=-1
C-----Neutrinos (left handed) and antineutrinos (right handed)
          IF(PT(J).EQ.12.OR.PT(J).EQ.14.OR.PT(J).EQ.16) FPOLARIZ(J)=-1
          IF(PT(J).EQ.-12.OR.PT(J).EQ.-14.OR.PT(J).EQ.-16) FPOLARIZ(J)=1
C---------Compute HTIMES2 (basically SPINTIMES2 with a sign according to polarization) in a counterintuitive convention (because it is the one used in the spheroidal functions).
               HTIMES2(J)=-SPINTIMES2(J)*FPOLARIZ(J)
               ERH(J)=PNBODY(4,J)*RH
               CNBODY(J)=ASTAR*ERH(J)
               IF(CNBODY(J).GT.50D0)CNBODY(J)=50D0
               XNOW(J)=PNBODY(3,J)/SQRT(PNBODY(1,J)**2+PNBODY(2,J)**2+
     &PNBODY(3,J)**2)
c-----------Initialises spheroidal function and calculates normalisation using an alternative integrator (chsphrdlcumulative3) with a bit less precision (to make it faster)
               CALL CHSPHRDLINI(LTIMES2(J),MTIMES2(J),HTIMES2(J),
     &CNBODY(J),SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),NLAST1(J),
     &SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J))
               CALL CHSPHRDLMATCH(CNBODY(J),MTIMES2(J),HTIMES2(J),
     &SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),NLAST1(J),
     &SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J),A(J))
               CALL CHSPHRDLCUMULATIVE3(CNBODY(J),LTIMES2(J),MTIMES2(J),
     &HTIMES2(J),SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),
     &NLAST1(J),SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J),
     &A(J),Y2(0,J),NYSIZE(J))
c------------Returns normalisation
               NORM(J)=Y2(NYSIZE(J),J)
            END DO
c------------Evaluate prefactor in the weight
            DISTVAL=SQRT(ENBODYQ(2)**2-PNBODY(5,2)**2)/
     &(PNBODY(4,1)**2*PNBODY(4,2)**2*ENBODYQ(2)**2)
            DO J=3,NBODYNOW
               DISTVAL=DISTVAL*SQRT((ENBODYQ(J)**2-PNBODY(5,J)**2)/
     &(PNBODY(4,J)**2-PNBODY(5,J)**2))
            END DO
c-----------Evaluate spheroidal and Hawking factors of the weight
            CALL CHNBODYHITMISS(NEXTRA,K,ERH,IASTAR,FRACDA,CNBODY,
     &LTIMES2,MTIMES2,HTIMES2,SPHRDLSERIES1,CONSTANTS1,DELTAN1,
     &NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2,A,NORM,1D0,
     &XNOW,DISTVAL,NBODYNOW)
c-----------Replace by new maximum
            IF(DISTVAL.GT.DISTMAX) THEN
               DISTMAX=DISTVAL
            END IF
 1      CONTINUE

c--------Start hit and miss
         DO 7 I=0,MAXTRY
c-----------Select a point in phase space with uniform probability as in the search
            DO J=2,NBODYNOW-1
               RANDOM(J)=CHRGEN(5)
            END DO
            CALL CHSPHRDLQCKSRT(RANDOM,NBODYNOW)
            DO J=NBODYNOW,2,-1
               PQTRANSFER(5,J-1)=MUI(J-1)+RANDOM(J-1)*MU
               ENBODYQ(J)=(PQTRANSFER(5,J)**2-PQTRANSFER(5,J-1)**2+
     &PNBODY(5,J)**2)/(2D0*PQTRANSFER(5,J))
               CALL CHDTWO(PQTRANSFER(1,J),PNBODY(1,J),PQTRANSFER(1,J-1)
     &,SQRT(ENBODYQ(J)**2-PNBODY(5,J)**2),2D0,.TRUE.)
            END DO
            CALL CHVEQU(4,PQTRANSFER(1,1),PNBODY(1,1))
c-----------Compute the normalisation as in the search above
            DO J=1,NBODYNOW
C-----POLARIZATIONS
               FPOLARIZ(J)=1
               IF(CHRGEN(6).LT.0.5)FPOLARIZ(J)=-1
C-----Neutrinos (left handed) and antineutrinos (right handed)
          IF(PT(J).EQ.12.OR.PT(J).EQ.14.OR.PT(J).EQ.16) FPOLARIZ(J)=-1
          IF(PT(J).EQ.-12.OR.PT(J).EQ.-14.OR.PT(J).EQ.-16) FPOLARIZ(J)=1
C---------Compute HTIMES2 (basically SPINTIMES2 with a sign according to polarization) in a counterintuitive convention (because it is the one used in the spheroidal functions).
               HTIMES2(J)=-SPINTIMES2(J)*FPOLARIZ(J)
               ERH(J)=PNBODY(4,J)*RH
               CNBODY(J)=ASTAR*ERH(J)
               IF(CNBODY(J).GT.50D0)CNBODY(J)=50D0
               XNOW(J)=PNBODY(3,J)/SQRT(PNBODY(1,J)**2+PNBODY(2,J)**2+
     &PNBODY(3,J)**2)
               CALL CHSPHRDLINI(LTIMES2(J),MTIMES2(J),HTIMES2(J),
     &CNBODY(J),SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),NLAST1(J),
     &SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J))
               CALL CHSPHRDLMATCH(CNBODY(J),MTIMES2(J),HTIMES2(J),
     &SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),NLAST1(J),
     &SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J),A(J))
               CALL CHSPHRDLCUMULATIVE3(CNBODY(J),LTIMES2(J),MTIMES2(J),
     &HTIMES2(J),SPHRDLSERIES1(0,J),CONSTANTS1(1,J),DELTAN1(J),
     &NLAST1(J),SPHRDLSERIES2(0,J),CONSTANTS2(1,J),DELTAN2(J),NLAST2(J),
     &A(J),Y2(0,J),NYSIZE(J))
               NORM(J)=Y2(NYSIZE(J),J)
            END DO
c-----------Compute weight as above... 
            DISTVAL=SQRT(ENBODYQ(2)**2-PNBODY(5,2)**2)/
     &(PNBODY(4,1)**2*PNBODY(4,2)**2*ENBODYQ(2)**2)
            DO J=3,NBODYNOW
               DISTVAL=DISTVAL*SQRT((ENBODYQ(J)**2-PNBODY(5,J)**2)/
     &(PNBODY(4,J)**2-PNBODY(5,J)**2))
            END DO
c-----------...and perform hit and miss
            CALL CHNBODYHITMISS(NEXTRA,K,ERH,IASTAR,FRACDA,CNBODY,
     &LTIMES2,MTIMES2,HTIMES2,SPHRDLSERIES1,CONSTANTS1,DELTAN1,
     &NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2,A,NORM,DISTMAX,
     &XNOW,DISTVAL,NBODYNOW)
            IF(DISTVAL.GE.CHRGEN(7)) GOTO 72
c-----------Try again
 7       CONTINUE
         CONVERGENCE=.FALSE.
 72      CONTINUE
      END IF
         
      END

CDECK  ID>, CHFINDJNBODY.
*CMZ :-        -07/07/08  15.58.30 by  Marco Sampaio
*-- Author :    Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHFINDJNBODY(CFLUX2SK1,CFLUX2SK2,MK,ASTAR,FRACDA,JNOW,
     &MNOW,JMIN,JMAX,SPINTIMES2,JMK)
C-----------------------------------------------------------------------
C     Performs hit and miss for the J selection using the weight function
c     for the last wave (after integrating out a kroenecker delta)
C-----------------------------------------------------------------------
C  *I/O:
C     --CFLUX2SK1,CFLUX2SK2: arrays with cumulative fluxes from the mode 
C       K=1 up to MK (number of modes available), for IASTAR and 
C       IASTAR+1 respectively, so that we can interpolate in astar.
C     --ASTAR,FRACDA: the true astar and the interpolation weight 
C       respectively.
C     --JNOW,MNOW,JMIN,JMAX,SPINTIMES2,JMK: J,M of the partial wave, 
C       minimum J, maximum J, twice the spin of the particle and j of 
C       the last allowed mode
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Input/Output:
      INTEGER MK,SPINTIMES2
      DOUBLE PRECISION CFLUX2SK1(MK+1),CFLUX2SK2(MK+1),FRACDA,ASTAR,JMIN
     &,JMAX,MNOW,JMK,JNOW
C--Local variables:
      INTEGER KNOW
      DOUBLE PRECISION CHRGEN,CHRUNI,JRANDOM,PNBODY,SPIN
      EXTERNAL CHRGEN,CHRUNI
      SPIN=SPINTIMES2/2D0
c------selects a J in the range we have access to
      JRANDOM=INT(CHRUNI(1,0D0,JMK+1-SPIN))+SPIN
      PNBODY=0D0
c------If in the allowed range, determine K and the weight (PNBODY)
      IF((JMIN.LE.JRANDOM).AND.(JRANDOM.LE.JMAX)) THEN
         KNOW=(JRANDOM+SPIN)*(JRANDOM+1-SPIN)+MNOW+1-SPIN
         IF(KNOW.GT.1) THEN
            PNBODY=(CFLUX2SK1(KNOW)-CFLUX2SK1(KNOW-1))*
     &(1-FRACDA)+(CFLUX2SK2(KNOW)-CFLUX2SK2(KNOW-1))*FRACDA
         ELSE
            PNBODY=CFLUX2SK1(KNOW)*(1-FRACDA)+
     &CFLUX2SK2(KNOW)*FRACDA
         END IF
      END IF
c------Do hit and miss by comparing with the upper bound on the maximum given from normalisation
      IF(PNBODY.GT.CHRGEN(2)*(CFLUX2SK1(MK+1)*(1-FRACDA)+CFLUX2SK2(MK+1)
     &*FRACDA)) THEN
c-------accept
         JNOW=JRANDOM
      ELSE
c-------reject
         JNOW=-1D0
      ENDIF

      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                           Special functions
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLINI(LTIMES2,MTIMES2,HTIMES2,C,SPHRDLSERIES1,
     &CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2)
C-----------------------------------------------------------------------
C     Generates the first estimate of the series expansion
C-----------------------------------------------------------------------
C   INPUT: 
C     LTIMES2-> 2L 
C     MTIMES2-> 2M 
C     HTIMES2-> "Helicity" component (times 2)
C     C-> a_\star \times z
C-----------------------------------------------------------------------
C   OUTPUT:
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so far
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER LTIMES2,MTIMES2,HTIMES2,NLAST1,NLAST2
      DOUBLE PRECISION C,SPHRDLSERIES1(0:NMAX),CONSTANTS1(1:4),
     &DELTAN1,SPHRDLSERIES2(0:NMAX),CONSTANTS2(1:4),DELTAN2
c*Local variables:
      INTEGER N,NMAXTHISMODE,NEIGEN,N1,NNN
      DOUBLE PRECISION ALPHA,BETA,DIAG(1:NMAX),DIAGBACK(1:NMAX),
     &SUBDIAG(1:NMAX),SUBDIAGBACK(1:NMAX),ERROREIGEN,EIGENPREVIOUS,
     &UNDFLOW,MAXRELERROR

      UNDFLOW=1E-200
      MAXRELERROR=1E-12

      ALPHA=ABS(MTIMES2-HTIMES2)/4D0
      BETA=ABS(MTIMES2+HTIMES2)/4D0
      NMAXTHISMODE=NMAX-INT(ABS(MTIMES2-ABS(HTIMES2))/2.)

*Estimates the order of truncation of the series
      NLAST1=10+INT(SQRT(LTIMES2*(LTIMES2+2)/4.+
     &(C+ABS(MTIMES2/2.)+ABS(HTIMES2)+3)**2))
      IF(NLAST1.GT.NMAXTHISMODE)THEN
         NLAST1=NMAXTHISMODE 
      END IF
      NLAST2=NLAST1

*Initialises constants used in the first recurrence relation (around x=-1) except for the eigenvalue.
      CONSTANTS1(1)=2*ALPHA
      CONSTANTS1(2)=(2*BETA-1+4*C)/2D0
      CONSTANTS1(4)=ALPHA+BETA-1-HTIMES2/2D0

      DO 1500 N=1,NLAST1
         DIAG(N)=N*((N+CONSTANTS1(1))/2D0+CONSTANTS1(2))
         SUBDIAG(N)=-SQRT(C*N*(N+CONSTANTS1(1))*(N+1+CONSTANTS1(4)))
         DIAGBACK(N)=DIAG(N)
         SUBDIAGBACK(N)=SUBDIAG(N)
 1500 CONTINUE

*Calculates the eigenvalue for the first expansion
      NEIGEN=NLAST1
      IF(C.GE.1E-3)THEN
         DO 1550 N=0,NMAX
            EIGENPREVIOUS=CONSTANTS1(3)
            CALL CHSPHRDLTQLI(DIAG,SUBDIAG,NEIGEN)
            CALL CHSPHRDLQCKSRT(DIAG,NEIGEN)
            IF((ABS(HTIMES2).EQ.2).AND.(ABS(MTIMES2).EQ.0))THEN
               CONSTANTS1(3)=DIAG((LTIMES2-ABS(MTIMES2))/2)
            ELSE
               CONSTANTS1(3)=DIAG(1+(LTIMES2-ABS(MTIMES2))/2)
            END IF
            NEIGEN=NEIGEN+10
            IF(ABS(CONSTANTS1(3)).GE.UNDFLOW)THEN
               ERROREIGEN=ABS((CONSTANTS1(3)
     &-EIGENPREVIOUS)/CONSTANTS1(3))
            ELSE
               ERROREIGEN=MAXRELERROR+1.
            END IF
            IF(((ERROREIGEN).LT.MAXRELERROR).OR.(NEIGEN.GT.NMAX))THEN
               GOTO 1555
            END IF
*Because the previous routines destroy DIAG, it is reinitialised by looking in the backup arrays BACK
            DO 1551 N1=1,NEIGEN-10
               DIAG(N1)=DIAGBACK(N1)
               SUBDIAG(N1)=SUBDIAGBACK(N1) 
 1551       CONTINUE

*Calculates more terms in the diagonal and subdiagonal (bigger matrix)
            DO 1552 N1=NEIGEN-9,NEIGEN
               DIAG(N1)=N1*((N1+CONSTANTS1(1))/2.+CONSTANTS1(2))
               SUBDIAG(N1)=-SQRT(C*N1*(N1+CONSTANTS1(1))
     &*(N1+1+CONSTANTS1(4)))
               DIAGBACK(N1)=DIAG(N1)
               SUBDIAGBACK(N1)=SUBDIAG(N1)
 1552       CONTINUE
            
 1550    CONTINUE
      ELSE
         IF((ABS(HTIMES2).EQ.2).AND.(ABS(MTIMES2).EQ.0))THEN
            NNN=(LTIMES2-ABS(MTIMES2))/2
            CONSTANTS1(3)=NNN*(CONSTANTS1(2)+0.5D0*(NNN+CONSTANTS1(1)))
         ELSE
            NNN=(LTIMES2-ABS(MTIMES2))/2+1
            CONSTANTS1(3)=NNN*(CONSTANTS1(2)+0.5D0*(NNN+CONSTANTS1(1)))
         END IF
      END IF
 1555 CONTINUE
         

*Initialises remaining constants for the first expansion      
      DELTAN1=-C/(2.+CONSTANTS1(1))
      SPHRDLSERIES1(0)=1D0
      SPHRDLSERIES1(1)=1D0/2D0+CONSTANTS1(2)/(1D0+CONSTANTS1(1))-
     &CONSTANTS1(3)/(1D0+CONSTANTS1(1))

*Initialises constants used in the second recurrence relation (around x=1)
      CONSTANTS2(1)=2*BETA
      CONSTANTS2(2)=(2*ALPHA-1+4*C)/2D0
      CONSTANTS2(3)=CONSTANTS1(3)+2*C*(ALPHA-BETA-HTIMES2/2D0)
      CONSTANTS2(4)=ALPHA+BETA-1+HTIMES2/2D0
     
      DELTAN2=-C/(2D0+CONSTANTS2(1))
      SPHRDLSERIES2(0)=1D0
      SPHRDLSERIES2(1)=1D0/2D0+CONSTANTS2(2)/(1D0+CONSTANTS2(1))-
     &CONSTANTS2(3)/(1D0+CONSTANTS2(1)) 

*Reccurence relation    
      DO 1003 N=2,NLAST1
         SPHRDLSERIES1(N)=(1D0/2D0+CONSTANTS1(2)/(N+CONSTANTS1(1))-
     &CONSTANTS1(3)/(N+CONSTANTS1(1))/N)*SPHRDLSERIES1(N-1)+
     &(1D0+CONSTANTS1(4)/N)*DELTAN1*SPHRDLSERIES1(N-2)
         DELTAN1=DELTAN1*(1D0-1D0/(N+1D0+CONSTANTS1(1)))

         SPHRDLSERIES2(N)=(1D0/2D0+CONSTANTS2(2)/(N+CONSTANTS2(1))-
     &CONSTANTS2(3)/(N+CONSTANTS2(1))/N)*SPHRDLSERIES2(N-1)+
     &(1D0+CONSTANTS2(4)/N)*DELTAN2*SPHRDLSERIES2(N-2)
         DELTAN2=DELTAN2*(1D0-1D0/(N+1D0+CONSTANTS2(1)))
 1003 CONTINUE
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLMATCH(C,MTIMES2,HTIMES2,SPHRDLSERIES1,
     &CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2
     &,A)
C-----------------------------------------------------------------------
C     Matches the two distributions in an intermediate region, checks
C     errors and finds the normalizing factor
C-----------------------------------------------------------------------
C   INPUT: 
C     C-> a_\star \times z
C     MTIMES1-> Twice the angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C-----------------------------------------------------------------------
C   OUTPUT:
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so far
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C     A-> The relative normalisation (S1=A*S2)
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER MTIMES2,HTIMES2,NLAST1,NLAST2
      DOUBLE PRECISION C,SPHRDLSERIES1(0:NMAX),CONSTANTS1(1:4),
     &DELTAN1,SPHRDLSERIES2(0:NMAX),CONSTANTS2(1:4),DELTAN2,A

*Points generation variables
      INTEGER NMATCH,NMATCHCOUNT,NMAXTHISMODE
      DOUBLE PRECISION LOGAAVERAGE,LOG2AAVERAGE,DLOGAAVERAGE,
     &UNDFLOW,MAXRELERROR,X,RELERR21,RELERR22,RELERR1,RELERR2,S1,S2

      MAXRELERROR=1E-3
      UNDFLOW=1E-200
      NMAXTHISMODE=NMAX-INT(ABS(MTIMES2-ABS(HTIMES2))/4.)
 
C--Calculates some points in the region of matching and tests the error. 
C  If the error is small calculates the relative normalisation,otherwise 
C  calculates more coefficients.

 1002 NMATCHCOUNT=0
      RELERR21=0
      RELERR22=0  
      LOGAAVERAGE=0
      LOG2AAVERAGE=0   
C--Calculates a first estimate for the relative normalisation by looking 
C  at 20 points in the matching region with both expansions
      DO 2002 NMATCH=0,20
         X=-0.25+NMATCH*0.5/20.
         CALL CHSERIES1ERR(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,S1,
     &RELERR1)
         CALL CHSERIES2ERR(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,S2,
     &RELERR2)

         IF(RELERR21.LT.RELERR1)THEN
            RELERR21=RELERR1
         ELSE IF(RELERR21.LT.RELERR2)THEN
            RELERR22=RELERR2
         END IF

         IF ((ABS(S1).GT.UNDFLOW).AND.(ABS(S2).GT.UNDFLOW)) THEN
            NMATCHCOUNT=NMATCHCOUNT+1
            LOGAAVERAGE=(LOGAAVERAGE*(NMATCHCOUNT-1)+LOG(ABS(S1))
     &-LOG(ABS(S2)))/NMATCHCOUNT
            LOG2AAVERAGE=(LOG2AAVERAGE*(NMATCHCOUNT-1)+(LOG(ABS(S1))
     &-LOG(ABS(S2)))**2)/NMATCHCOUNT   
         END IF
 2002 CONTINUE

*calculates the relative error in A
      DLOGAAVERAGE=SQRT(ABS(LOG2AAVERAGE-LOGAAVERAGE*LOGAAVERAGE))
      IF(((DLOGAAVERAGE.GT.MAXRELERROR).OR.(RELERR21.GT.MAXRELERROR)
     &.OR.(RELERR22.GT.MAXRELERROR)))THEN

C--Uses all the errors calculated before and calculates more coefficients
C  for the expansion if needed by going back to the begining
         IF((RELERR21.GT.MAXRELERROR).AND.(NLAST1.LT.NMAXTHISMODE)) THEN
            CALL CHSPHRDLMORE1(NMAXTHISMODE,SPHRDLSERIES1,CONSTANTS1,
     &DELTAN1,NLAST1)
         ELSE IF((DLOGAAVERAGE.GT.MAXRELERROR).AND.
     &(NLAST1.LT.NMAXTHISMODE))THEN
            CALL CHSPHRDLMORE1(NMAXTHISMODE,SPHRDLSERIES1,CONSTANTS1,
     &DELTAN1,NLAST1)
         END IF

         IF((RELERR22.GT.MAXRELERROR).AND.(NLAST2.LT.NMAXTHISMODE)) THEN
            CALL CHSPHRDLMORE2(NMAXTHISMODE,SPHRDLSERIES2,CONSTANTS2,
     &DELTAN2,NLAST2)
         ELSE IF((DLOGAAVERAGE.GT.MAXRELERROR).AND.
     &(NLAST2.LT.NMAXTHISMODE))THEN
            CALL CHSPHRDLMORE2(NMAXTHISMODE,SPHRDLSERIES2,CONSTANTS2,
     &DELTAN2,NLAST2)
         END IF
         IF ((NLAST1.LT.NMAXTHISMODE).OR.(NLAST2.LT.NMAXTHISMODE))THEN
            GOTO 1002
         END IF
            
      END IF

      A=EXP(LOGAAVERAGE)
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLCUMULATIVE(C,LTIMES2,MTIMES2,HTIMES2,
     &SPHRDLSERIES1,CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,
     &DELTAN2,NLAST2,A,Y,XOFY,NYSIZE)
C-----------------------------------------------------------------------
C     Calculates an array of points for the integral of S^2(x)+S^2(-x) 
C     and its inverse
C-----------------------------------------------------------------------
C   INPUT: 
C     C-> a_\star \times z
C     LTIMES2-> Twice angular momentum quantum number
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so fa
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C     A-> The relative normalisation (S1=A*S2)
C-----------------------------------------------------------------------
C   OUTPUT:
C     Y-> Cumulative function
C     XOFY-> Inverse of the cumulative function
C     NYSIZE-> Size of the array Y
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER LTIMES2,MTIMES2,HTIMES2,NLAST1,NLAST2,NYSIZE
      DOUBLE PRECISION C,SPHRDLSERIES1(0:NMAX),CONSTANTS1(1:4),
     &DELTAN1,SPHRDLSERIES2(0:NMAX),CONSTANTS2(1:4),DELTAN2,A

*Local variables:
      INTEGER N,NDIVIDE,NDIVIDEMIN
      DOUBLE PRECISION DELTAX,DELTAXPREVIOUS,DISTRIBUTION, DISTRIBUTION1
     &,DISTRIBUTION2,DISTPREVIOUS,DERIVATIVE1,DERIVATIVE2,X,
     &S1,S2,Y(0:NDIVIDEMAX+1),XOFY(1:NDIVIDEMAX+1)  


C--Choose an estimate for the minimum number of points in the integration
C  and use it as first guess. 
      NDIVIDEMIN=10+5*INT((LTIMES2-ABS(HTIMES2))/2.)
      NDIVIDE=NDIVIDEMIN
C--If larger than the maximum number of points use the NDIVIDEMAX
      IF (NDIVIDE.GT.NDIVIDEMAX) THEN
         NDIVIDE=NDIVIDEMAX
         NDIVIDEMIN=NDIVIDEMAX
      END IF
      
C-- Calculates the first point used in the integration and stores its 
C   symmetric times the DX in Y(0) so that the first point Y(1), which
C   corresponds to X(1)=-1 will be zero  
      X=-1D0
      DELTAX=2D0/NDIVIDE
      DELTAXPREVIOUS=DELTAX

      CALL CHSERIES1(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,S1)
      DISTPREVIOUS=(S1/A)**2
      Y(0)=-DISTPREVIOUS*DELTAXPREVIOUS

C--Integrates and checks for precision. The step of integration is decided 
C  for the next integration in the current loop, by evaluating the spheroidal 
C  function at the current point, then at the next point X+DX and a point at 
C  X+DX/2. The step DX is reduced if the plot is not linear enough.
      DO 1004 N=0,NDIVIDEMAX
         IF ((X.GE.1.).AND.(X-DELTAX.LE.0.9999D0))THEN
C--------If the previous point calculated was at x<1 we are calculating the
C        final point so X is set to 1, and we go out of the loop to calculate it
            X=1D0
            NYSIZE=N+1
            GOTO 2004
         ELSE IF (X.GE.1)THEN
C--------If the second condition in the previous if doesn't hold, then the final 
C        point has been already calculated
            NYSIZE=N
            GOTO 2005
         END IF
C--Calculates the current point
      IF(X.LT.0)THEN
         CALL CHSERIES1(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,S1)
         DISTRIBUTION=(S1/A)**2
      ELSE
         CALL CHSERIES2(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,S2)
         DISTRIBUTION=S2**2
      END IF
C--Calculates the possible next points
 1504 IF(X.LT.0)THEN
         CALL CHSERIES1(X+DELTAX,C,MTIMES2,HTIMES2,SPHRDLSERIES1,
     &NLAST1,S1)
         DISTRIBUTION1=(S1/A)**2
         CALL CHSERIES1(X+DELTAX/2D0,C,MTIMES2,HTIMES2,SPHRDLSERIES1,
     &NLAST1,S1)
         DISTRIBUTION2=(S1/A)**2
      ELSE
         CALL CHSERIES2(X+DELTAX,C,MTIMES2,HTIMES2,SPHRDLSERIES2,
     &NLAST2,S2)
         DISTRIBUTION1=S2**2
         CALL CHSERIES2(X+DELTAX/2D0,C,MTIMES2,HTIMES2,SPHRDLSERIES2,
     &NLAST2,S2)
         DISTRIBUTION2=S2**2
      END IF
C------Estimates derivative with the two next points
         DERIVATIVE1=(DISTRIBUTION1-DISTRIBUTION)/DELTAX
         DERIVATIVE2=2.*(DISTRIBUTION2-DISTRIBUTION)/DELTAX
C------Tests for error in slope of the curve for interpolations where the curve is not flat 
C      (the derivative in units of the height of the plot is large enough)
         IF(ABS(DERIVATIVE1/DISTRIBUTION1).GE.0.1)THEN
C---------Tests for error in slope and decreasis step by 1/2 (only down to a minimum)
            IF((ABS((DERIVATIVE2-DERIVATIVE1)/DERIVATIVE1).GE.0.1).AND.
     &((NDIVIDEMAX-N).GT.(2*(NDIVIDE-N)))) THEN
               NDIVIDE=2*NDIVIDE
               DELTAX=2D0/NDIVIDE
               GOTO 1504
            END IF
C------Otherwise if the curve is flat enough, increase the step (only up to a maximum)
         ELSE IF(NDIVIDE.GT.NDIVIDEMIN)THEN
            NDIVIDE=NDIVIDE/2
            DELTAX=2D0/NDIVIDE
         END IF
C--Store the point calculated
         XOFY(N+1)=X
         Y(N+1)=Y(N)+(DISTRIBUTION+DISTPREVIOUS)*DELTAXPREVIOUS/2D0
C--Store the current stepsize to be used in the next integration and the current height of 
C  the distribution         
         DELTAXPREVIOUS=DELTAX
         DISTPREVIOUS=DISTRIBUTION
C--Go to the next X by adding DX
         X=X+DELTAX  
C--Check for the maximum of the distribution
 1004 CONTINUE

C--Calculate the last point (if not yet done)
 2004 CALL CHSERIES2(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,S2)
      DISTRIBUTION=S2**2
      XOFY(NYSIZE)=1D0
      Y(NYSIZE)=Y(NYSIZE-1)+(DISTRIBUTION+DISTPREVIOUS)*DELTAXPREVIOUS
     &/2D0

 2005 CONTINUE

      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLCUMULATIVE2(C,LTIMES2,MTIMES2,
     &HTIMES2,SPHRDLSERIES1,CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,
     &CONSTANTS2,DELTAN2,NLAST2,A,Y,XOFY,NYSIZE)
C-----------------------------------------------------------------------
C     Computes cumulative function for product of spheroidal functions
C-----------------------------------------------------------------------
C   INPUT: 
C     C-> a_\star \times z
C     LTIMES2-> Twice angular momentum quantum number
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so fa
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C     A-> The relative normalisation (S1=A*S2)
C-----------------------------------------------------------------------
C   OUTPUT:
C     Y-> Cumulative function
C     XOFY-> Inverse of the cumulative function
C     NYSIZE-> Size of the array Y
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER LTIMES2(2),MTIMES2(2),HTIMES2(2),NLAST1(2),NLAST2(2),
     &NYSIZE
      DOUBLE PRECISION C(2),SPHRDLSERIES1(0:NMAX,2),
     &CONSTANTS1(1:4,2),DELTAN1(2),SPHRDLSERIES2(0:NMAX,2),
     &CONSTANTS2(1:4,2),DELTAN2(2),A(2)

*Local variables:
      INTEGER N,NDIVIDE,NDIVIDEMIN
      DOUBLE PRECISION DELTAX,DELTAXPREVIOUS,DISTRIBUTION, DISTRIBUTION1
     &,DISTRIBUTION2,DISTPREVIOUS,DERIVATIVE1,DERIVATIVE2,X,
     &S1(2),S2(2),Y(0:NDIVIDEMAX+1),XOFY(1:NDIVIDEMAX+1)  


C--Choose an estimate for the minimum number of points in the integration
C  and use it as first guess. 
      NDIVIDEMIN=10+5*INT((LTIMES2(1)-ABS(HTIMES2(1)))/2.)
      NDIVIDE=NDIVIDEMIN
C--If larger than the maximum number of points use the NDIVIDEMAX
      IF (NDIVIDE.GT.NDIVIDEMAX) THEN
         NDIVIDE=NDIVIDEMAX
         NDIVIDEMIN=NDIVIDEMAX
      END IF
      
C-- Calculates the first point used in the integration and stores its 
C   symmetric times the DX in Y(0) so that the first point Y(1), which
C   corresponds to X(1)=-1 will be zero  
      X=-1D0
      DELTAX=2D0/NDIVIDE
      DELTAXPREVIOUS=DELTAX
      CALL CHSERIES1(X,C(1),MTIMES2(1),HTIMES2(1),SPHRDLSERIES1(0,1),
     &NLAST1(1),S1(1))
      CALL CHSERIES1(X,C(2),MTIMES2(2),HTIMES2(2),SPHRDLSERIES1(0,2),
     &NLAST1(2),S1(2))
      DISTPREVIOUS=(S1(1)/A(1))**2*(S1(2)/A(2))**2
      Y(0)=-DISTPREVIOUS*DELTAXPREVIOUS
C--Integrates and checks for precision. The step of integration is decided 
C  for the next integration in the current loop, by evaluating the spheroidal 
C  function at the current point, then at the next point X+DX and a point at 
C  X+DX/2. The step DX is reduced if the plot is not linear enough.
      DO 1004 N=0,NDIVIDEMAX
         IF ((X.GE.1.).AND.(X-DELTAX.LE.0.9999D0))THEN
C--------If the previous point calculated was at x<1 we are calculating the
C        final point so X is set to 1, and we go out of the loop to calculate it
            X=1D0
            NYSIZE=N+1
            GOTO 2004
         ELSE IF (X.GE.1)THEN
C--------If the second condition in the previous "if" doesn't hold, then the final 
C        point has been already calculated
            NYSIZE=N
            GOTO 2005
         END IF
C--Calculates the current point
      IF(X.LT.0)THEN
        CALL CHSERIES1(X,C(1),MTIMES2(1),HTIMES2(1),SPHRDLSERIES1(0,1),
     &NLAST1(1),S1(1))
        CALL CHSERIES1(X,C(2),MTIMES2(2),HTIMES2(2),SPHRDLSERIES1(0,2),
     &NLAST1(2),S1(2))
      DISTRIBUTION=(S1(1)/A(1))**2*(S1(2)/A(2))**2
      ELSE
        CALL CHSERIES2(X,C(1),MTIMES2(1),HTIMES2(1),SPHRDLSERIES2(0,1),
     &NLAST2(1),S2(1))
        CALL CHSERIES2(X,C(2),MTIMES2(2),HTIMES2(2),SPHRDLSERIES2(0,2),
     &NLAST2(2),S2(2))
        DISTRIBUTION=S2(1)**2*S2(2)**2
      END IF
C--Calculates the possible next points
 1504 IF(X.LT.0)THEN
      CALL CHSERIES1(X+DELTAX,C(1),MTIMES2(1),HTIMES2(1),
     &SPHRDLSERIES1(0,1),NLAST1(1),S1(1))
      CALL CHSERIES1(X+DELTAX,C(2),MTIMES2(2),HTIMES2(2),
     &SPHRDLSERIES1(0,2),NLAST1(2),S1(2))
      DISTRIBUTION1=(S1(1)/A(1))**2*(S1(2)/A(2))**2
      CALL CHSERIES1(X+DELTAX/2D0,C(1),MTIMES2(1),HTIMES2(1),
     &SPHRDLSERIES1(0,1),NLAST1(1),S1(1))
      CALL CHSERIES1(X+DELTAX/2D0,C(2),MTIMES2(2),HTIMES2(2),
     &SPHRDLSERIES1(0,2),NLAST1(2),S1(2))
      DISTRIBUTION2=(S1(1)/A(1))**2*(S1(2)/A(2))**2

      ELSE
      CALL CHSERIES2(X+DELTAX,C(1),MTIMES2(1),HTIMES2(1),
     &SPHRDLSERIES2(0,1),NLAST2(1),S2(1))
      CALL CHSERIES2(X+DELTAX,C(2),MTIMES2(2),HTIMES2(2),
     &SPHRDLSERIES2(0,2),NLAST2(2),S2(2))
      DISTRIBUTION1=S2(1)**2*S2(2)**2
      CALL CHSERIES2(X+DELTAX/2D0,C(1),MTIMES2(1),HTIMES2(1),
     &SPHRDLSERIES2(0,1),NLAST2(1),S2(1))
      CALL CHSERIES2(X+DELTAX/2D0,C(2),MTIMES2(2),HTIMES2(2),
     &SPHRDLSERIES2(0,2),NLAST2(2),S2(2))
      DISTRIBUTION2=S2(1)**2*S2(2)**2

      END IF
C------Estimates derivative with the two next points
         DERIVATIVE1=(DISTRIBUTION1-DISTRIBUTION)/DELTAX
         DERIVATIVE2=2.*(DISTRIBUTION2-DISTRIBUTION)/DELTAX
C------Tests for error in slope of the curve for imterpolations where the curve is not flat 
C      (the derivative in units of the height of the plot is large enough)
         IF(ABS(DERIVATIVE1/DISTRIBUTION1).GE.0.1)THEN
C---------Tests for error in slope and decreasis step by 1/2 (only down to a minimum)
            IF((ABS((DERIVATIVE2-DERIVATIVE1)/DERIVATIVE1).GE.0.1).AND.
     &((NDIVIDEMAX-N).GT.(2*(NDIVIDE-N)))) THEN
               NDIVIDE=2*NDIVIDE
               DELTAX=2./NDIVIDE
               GOTO 1504
            END IF
C------Otherwise if the curve is flat enough, increase the step (only up to a maximum)
         ELSE IF(NDIVIDE.GT.NDIVIDEMIN)THEN
            NDIVIDE=NDIVIDE/2
            DELTAX=2D0/NDIVIDE
         END IF
C--Store the point calculated
         XOFY(N+1)=X
         Y(N+1)=Y(N)+(DISTRIBUTION+DISTPREVIOUS)*DELTAXPREVIOUS/2.
C--Store the current stepsize to be used in the next integration and the current height of 
C  the distribution         
         DELTAXPREVIOUS=DELTAX
         DISTPREVIOUS=DISTRIBUTION
C--Go to the next X by adding DX
         X=X+DELTAX  
 1004 CONTINUE
      NYSIZE=NDIVIDEMAX
C--Calculate the last point (if not yet done)
 2004 CALL CHSERIES2(X,C(1),MTIMES2(1),HTIMES2(1),SPHRDLSERIES2(0,1),
     &NLAST2(1),S2(1))
      CALL CHSERIES2(X,C(2),MTIMES2(2),HTIMES2(2),SPHRDLSERIES2(0,2),
     &NLAST2(2),S2(2))
      DISTRIBUTION=S2(1)**2*S2(2)**2
      XOFY(NYSIZE)=1D0
      Y(NYSIZE)=Y(NYSIZE-1)+(DISTRIBUTION+DISTPREVIOUS)*DELTAXPREVIOUS
     &/2D0

 2005 CONTINUE
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLCUMULATIVE3(C,LTIMES2,MTIMES2,HTIMES2,
     &SPHRDLSERIES1,CONSTANTS1,DELTAN1,NLAST1,SPHRDLSERIES2,CONSTANTS2,
     &DELTAN2,NLAST2,A,Y,NYSIZE)
C-----------------------------------------------------------------------
C     Calculates an array of points for the integral of S^2(x)+S^2(-x) 
C     and its inverse. This is a lower precision routine to speed up hit 
C     or miss MC
C-----------------------------------------------------------------------
C   INPUT: 
C     C-> a_\star \times z
C     LTIMES2-> Twice angular momentum quantum number
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so fa
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C     A-> The relative normalisation (S1=A*S2)
C-----------------------------------------------------------------------
C   OUTPUT:
C     Y-> Cumulative function
C     NYSIZE-> Size of the array Y
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=50)
*INPUT and OUTPUT:
      INTEGER LTIMES2,MTIMES2,HTIMES2,NLAST1,NLAST2,NYSIZE
      DOUBLE PRECISION C,SPHRDLSERIES1(0:NMAX),CONSTANTS1(1:4),
     &DELTAN1,SPHRDLSERIES2(0:NMAX),CONSTANTS2(1:4),DELTAN2,A

*Local variables:
      INTEGER N
      DOUBLE PRECISION DELTAX,DISTRIBUTION,DISTPREVIOUS,X,S1,S2,
     &Y(0:NDIVIDEMAX+1)

C-- Calculates the first point used in the integration and stores its 
C   symmetric times the DX in Y(0) so that the first point Y(1), which
C   corresponds to X(1)=-1 will be zero.  
      X=-1D0
      DELTAX=4D-2

      CALL CHSERIES1(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,S1)
      DISTPREVIOUS=(S1/A)**2
      Y(0)=-DISTPREVIOUS*DELTAX

C--Integrates
      DO 1004 N=0,50
         IF ((X.GE.1.).AND.(X-DELTAX.LE.0.9999D0))THEN
C--------If the previous point calculated was at x<1 we are calculating the
C        final point so X is set to 1, and we go out of the loop to calculate it
            X=1D0
            NYSIZE=N+1
            GOTO 2004
         ELSE IF (X.GE.1)THEN
C--------If the second condition in the previous "if" doesn't hold, then the final 
C        point has been already calculated
            NYSIZE=N
            GOTO 2005
         END IF
C--Calculates the current point
      IF(X.LT.0)THEN
         CALL CHSERIES1(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,S1)
         DISTRIBUTION=(S1/A)**2
      ELSE
         CALL CHSERIES2(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,S2)
         DISTRIBUTION=S2**2
      END IF

         Y(N+1)=Y(N)+(DISTRIBUTION+DISTPREVIOUS)*DELTAX/2D0
C--Store the current height of the distribution for the next step of integration         
         DISTPREVIOUS=DISTRIBUTION
C--Go to the next X by adding DX
         X=X+DELTAX  
 1004 CONTINUE
      NYSIZE=50
C--Calculate the last point (if not yet done)
 2004 CALL CHSERIES2(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,S2)
      DISTRIBUTION=S2**2
      Y(NYSIZE)=Y(NYSIZE-1)+(DISTRIBUTION+DISTPREVIOUS)*DELTAX/2D0

 2005 CONTINUE

      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHNBODYHITMISS(NEX,K,ERH,IASTAR,FRACDA,C,LTIMES2,
     &MTIMES2,HTIMES2,SPHRDLSERIES1,CONSTANTS1,DELTAN1,NLAST1,
     &SPHRDLSERIES2,CONSTANTS2,DELTAN2,NLAST2,A,NORM,DISTMAX,
     &X,DISTVAL,NBODYNOW)
C-----------------------------------------------------------------------
C     Calculates the weight for the hit ans miss of the NBODY decay
C-----------------------------------------------------------------------
C   INPUT: 
C     NEX-> number of extra dimensions
C     K-> partial waves
C     ERH-> Energies times RH
C     IASTAR-> Location of current astar in the tables
C     FRACDA-> Interpolaion fraction
C     C-> a_\star \times z
C     LTIMES2-> Twice angular momentum quantum number
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so fa
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C     A-> The relative normalisation (S1=A*S2)
C     NORM-> Normalising factors
C     DISTMAX-> Maximum of the probability distribution
C     X-> cos(theta) for each particle
C     NBODYNOW-> current multiplicity of the decay
C-----------------------------------------------------------------------
C   OUTPUT:
C     DISTVAL-> Value of the probability distribution at this phase 
C               space point
C-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX,MA,ME,MK0,MK1,MK2,UFLUXKE,MAXETABLE,
     &NBODYNOW
      PARAMETER(NMAX=250,NDIVIDEMAX=150,MA=26,MK0=169,MK1=156,MK2=120,
     &ME=100,MAXETABLE=5D0)
*INPUT and OUTPUT:
      INTEGER LTIMES2(NBODYNOW),MTIMES2(NBODYNOW),HTIMES2(NBODYNOW),
     &NLAST1(NBODYNOW)
     &,NLAST2(NBODYNOW),NEX,K(NBODYNOW),IASTAR
      DOUBLE PRECISION C(NBODYNOW),SPHRDLSERIES1(0:NMAX,NBODYNOW),
     &CONSTANTS1(1:4,NBODYNOW),DELTAN1(NBODYNOW),
     &SPHRDLSERIES2(0:NMAX,NBODYNOW),
     &CONSTANTS2(1:4,NBODYNOW),DELTAN2(NBODYNOW),A(NBODYNOW),
     &ERH(NBODYNOW),DISTMAX,
     &X(NBODYNOW),FRACDA,DISTVAL,NORM(NBODYNOW)
*Local variables:
      INTEGER J
      DOUBLE PRECISION DISTNOW,S1,S2,CHRGEN,FLUX2S0KE(ME+3,MK0,MA),
     &FLUX2S1KE(ME+3,MK1,MA),FLUX2S2KE(ME+3,MK2,MA),INOW,DI
      LOGICAL FIRSTNBODY
      EXTERNAL CHRGEN
      CHARACTER NEXTRA*1

      DATA UFLUXKE,FIRSTNBODY/24,.TRUE./
      SAVE FLUX2S0KE,FLUX2S1KE,FLUX2S2KE
      
C--Initialises the arrays with fluxes by reading the 
C  corresponding data file for a number of extra dimensions NEX This 
C  occurs only once at FIRSTNBODY call.
      IF (FIRSTNBODY) THEN
         NEXTRA=CHAR(ICHAR('0')+NEX)
         IF ((NEX.GT.0).AND.(NEX.LT.7))
     &THEN
            OPEN(UNIT=UFLUXKE,FILE='./data_files/fluxke_n'//NEXTRA//
     &'.dat',STATUS='OLD')
         ELSE
            PRINT *,' NO CUMULATIVE TABLES YET FOR NEX =',NEX
            STOP
         ENDIF
         READ(UFLUXKE,*)FLUX2S0KE,FLUX2S1KE,FLUX2S2KE
         CLOSE(UNIT=UFLUXKE)
         FIRSTNBODY=.FALSE.
      ENDIF

C--Computes distribution through the product of all factors
      DO 1 J=1,NBODYNOW
         INOW=INT(ERH(J)/MAXETABLE*ME)+1
         DI=ERH(J)/MAXETABLE*ME-INOW+1
         IF(HTIMES2(J).EQ.0)THEN
C-----Scalars
            DISTNOW=FLUX2S0KE(INOW,K(J),IASTAR)*(1-FRACDA)+
     &FLUX2S0KE(INOW,K(J),IASTAR+1)*FRACDA
            DISTNOW=DISTNOW*(1-DI)+DI*(FLUX2S0KE(INOW+1,K(J),IASTAR)
     &*(1-FRACDA)+FLUX2S0KE(INOW+1,K(J),IASTAR+1)*FRACDA)
         ELSE IF(HTIMES2(J).EQ.1) THEN
C-----Fermions
            DISTNOW=FLUX2S1KE(INOW,K(J),IASTAR)*(1-FRACDA)+
     &FLUX2S1KE(INOW,K(J),IASTAR+1)*FRACDA
            DISTNOW=DISTNOW*(1-DI)+DI*(FLUX2S1KE(INOW+1,K(J),IASTAR)
     &*(1-FRACDA)+FLUX2S1KE(INOW+1,K(J),IASTAR+1)*FRACDA)
         ELSE IF (HTIMES2(J).EQ.2) THEN
C-----Vectors
            DISTNOW=FLUX2S2KE(INOW,K(J),IASTAR)*(1-FRACDA)+
     &FLUX2S2KE(INOW,K(J),IASTAR+1)*FRACDA
            DISTNOW=DISTNOW*(1-DI)+DI*(FLUX2S2KE(INOW+1,K(J),IASTAR)
     &*(1-FRACDA)+FLUX2S2KE(INOW+1,K(J),IASTAR+1)*FRACDA)
         END IF
c-----Spheroidal factor
         IF(X(J).LT.0)THEN
            CALL CHSERIES1(X(J),C(J),MTIMES2(J),HTIMES2(J),
     &SPHRDLSERIES1(0,J),NLAST1(J),S1)
            DISTVAL=DISTNOW*DISTVAL*(S1/A(J))**2/NORM(J)
         ELSE
            CALL CHSERIES2(X(J),C(J),MTIMES2(J),HTIMES2(J),
     &SPHRDLSERIES2(0,J),NLAST2(J),S2)
c-----All the factors
            DISTVAL=DISTNOW*DISTVAL*S2**2/NORM(J)
         END IF
 1    CONTINUE
c-----Divide by the maximum to get the weight (if not searching for the maximum, i.e if DISTMAX different from 1)
      DISTVAL=DISTVAL/DISTMAX

      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLXRAND(Y,XOFY,NYSIZE,XRAND)
C-----------------------------------------------------------------------
C     Returns an X with correct distribution given by S^2(x)+S^2(-x)
C-----------------------------------------------------------------------
C   INPUT: 
C     Y-> Cumulative function
C     XOFY-> Inverse of the cumulative function
C     NYSIZE-> Size of the array Y
C-----------------------------------------------------------------------
C   OUTPUT:
C     XRAND-> COSTHETA generated according to the angular distribution
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER NYSIZE
      DOUBLE PRECISION XRAND,Y(0:NDIVIDEMAX+1),XOFY(1:NDIVIDEMAX+1) 
*Local variables:
      INTEGER N,NTEST1,NTEST2,NTEST
      DOUBLE PRECISION YRAND,CHRGEN 
      EXTERNAL CHRGEN 
         
C--Generates a random number between 0 and 1 i.e. the cumulant to be inverted to get X
      NTEST1=1
      NTEST2=NYSIZE
      YRAND=CHRGEN(1)
C-- Bracketing routine to find the inverse of YRAND
      DO 3000 N=0,NYSIZE
      NTEST=NTEST1+(NTEST2-NTEST1)/2
         IF(Y(NTEST)/Y(NYSIZE).LE.YRAND)THEN
            NTEST1=NTEST
            IF((NTEST2-NTEST1).EQ.1) GOTO 3200
         ELSE 
            IF ((NTEST-NTEST1).LE.1) THEN
               NTEST2=NTEST1+1
               GOTO 3200
            ELSE
               NTEST2=NTEST
            END IF
         END IF
 3000 CONTINUE

 3200 CONTINUE
C--Bracketing complete. Now interpolate in X.
      XRAND=XOFY(NTEST1)+(XOFY(NTEST2)-XOFY(NTEST1))
     &*(YRAND*Y(NYSIZE)-Y(NTEST1))/(Y(NTEST2)-Y(NTEST1))
C--Renormalise if there is a small numerical error that pushes abs(xrand) outside allowed range
      IF(ABS(XRAND).GT.1) XRAND=XRAND/ABS(XRAND)

      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLMORE1(NMAXTHISMODE,SPHRDLSERIES1,CONSTANTS1,
     &DELTAN1,NLAST1)
C-----------------------------------------------------------------------
C     Generates the first estimate of the series expansion around x=-1
C-----------------------------------------------------------------------
C   I/O:
C     MTIMES2-> Twice angular momentum along z quantum number
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     CONSTANTS1-> The fixed constants used in the reccurence
C     DELTAN1-> Another coefficient calculated recursively
C     NLAST1-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER NLAST1
      DOUBLE PRECISION SPHRDLSERIES1(0:NMAX),CONSTANTS1(1:4),DELTAN1
*Local variables:
      INTEGER N,NHIGH,NMAXTHISMODE

C--Calculates 10 more terms (if possible) otherwise calculates as many 
C  terms as allowed from the size of the arrays. This upper limit 
C  is set in such a way that prevents overflow when we evaluate 
C  the series which will contain power of the type ~2**N then multiplied
C  by coefficients of order ~2**(-N)small numbers.          
      NHIGH=NLAST1+10      
      IF(NHIGH.GT.NMAXTHISMODE)THEN
         NHIGH=NMAXTHISMODE
      END IF
C--Calculates more coefficients using the recursion relation
      DO 2003 N=(NLAST1+1),NHIGH
         SPHRDLSERIES1(N)=(0.5D0+CONSTANTS1(2)/(N+CONSTANTS1(1))-
     &CONSTANTS1(3)/(N+CONSTANTS1(1))/N)*SPHRDLSERIES1(N-1)+
     &(1D0+CONSTANTS1(4)/N)*DELTAN1*SPHRDLSERIES1(N-2)
         DELTAN1=DELTAN1*(1D0-1D0/(N+1D0+CONSTANTS1(1)))
 2003 CONTINUE
C--New size for the array of coefficients for the expansion
      NLAST1=NHIGH
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLMORE2(NMAXTHISMODE,SPHRDLSERIES2,CONSTANTS2,
     &DELTAN2,NLAST2)
C-----------------------------------------------------------------------
C     Generates the first estimate of the series expansion around x=1
C-----------------------------------------------------------------------
C   I/O:
C     MTIMES2-> Twice angular momentum along z quantum number
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     CONSTANTS2-> The fixed constants used in the reccurence
C     DELTAN2-> Another coefficient calculated recursively
C     NLAST2-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER NLAST2
      DOUBLE PRECISION SPHRDLSERIES2(0:NMAX),CONSTANTS2(1:4),DELTAN2
*Local variables:
      INTEGER N,NHIGH,NMAXTHISMODE

C--Everything similar to the function CHSPHRDLMORE1 (see comments there)  
      NHIGH=NLAST2+10      
      IF(NHIGH.GT.NMAXTHISMODE)THEN
         NHIGH=NMAXTHISMODE
      END IF

      DO 3003 N=(NLAST2+1),NHIGH
         SPHRDLSERIES2(N)=(1D0/2D0+CONSTANTS2(2)/(N+CONSTANTS2(1))-
     &CONSTANTS2(3)/(N+CONSTANTS2(1))/N)*SPHRDLSERIES2(N-1)+
     &(1D0+CONSTANTS2(4)/N)*DELTAN2*SPHRDLSERIES2(N-2)
         DELTAN2=DELTAN2*(1D0-1D0/(N+1D0+CONSTANTS2(1)))
 3003 CONTINUE

      NLAST2=NHIGH
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSERIES1ERR(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,
     &S1,RELERR1)
C-----------------------------------------------------------------------
C     Calculates the series at a specific point X with the expansion
C     around x=-1 and estimates the error
C-----------------------------------------------------------------------
C   INPUT: 
C     X-> COS(THETA)
C     C-> a times omega
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     NLAST1-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
C   OUTPUT:
C     S1-> The angular function evaluated at x=X
C     RELERR1-> The relative error evaluated by the first neglected term
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER MTIMES2,HTIMES2,NLAST1
      DOUBLE PRECISION C,SPHRDLSERIES1(0:NMAX),S1,RELERR1,X

*Local variables:
      INTEGER N
      DOUBLE PRECISION FIRSTNEGLECT

      S1=0.
      RELERR1=0.
C--Calculates the value of the spheroidal function using the series
      DO 4003 N=0,NLAST1
         S1=S1+SPHRDLSERIES1(N)*
     &SQRT(SQRT((1+X)**(4*N+ABS(MTIMES2-HTIMES2)))) 
 4003 CONTINUE
      N=N-1
C--Estimates error by looking at the last contribution for the series
      FIRSTNEGLECT=SPHRDLSERIES1(N)*
     &SQRT(SQRT((1+X)**(4*N+ABS(MTIMES2-HTIMES2))))
      IF (S1.NE.0) THEN
         RELERR1=ABS(FIRSTNEGLECT/S1)
      END IF
C--Computes the result for the spheroidal function
      S1=EXP(-C*X)*SQRT(SQRT(((1-X)**ABS(MTIMES2+HTIMES2))))*S1
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSERIES1(X,C,MTIMES2,HTIMES2,SPHRDLSERIES1,NLAST1,S1)
C-----------------------------------------------------------------------
C     Calculates the series at a specific point X with the expansion
C     around x=-1
C-----------------------------------------------------------------------
C   INPUT: 
C     X-> COS(THETA)
C     C-> a times omega
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES1-> array of coefficients of the expansion
C     NLAST1-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
C   OUTPUT:
C     S1-> The angular function evaluated at x=X
C     RELERR1-> The relative error evaluated by the first neglected term
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER MTIMES2,HTIMES2,NLAST1
      DOUBLE PRECISION C,SPHRDLSERIES1(0:NMAX),S1,X
*Local variables:
      INTEGER N

      S1=0D0
      DO 4003 N=0,NLAST1
         S1=S1+SPHRDLSERIES1(N)*
     &SQRT(SQRT((1+X)**(4*N+ABS(MTIMES2-HTIMES2)))) 
 4003 CONTINUE

      S1=EXP(-C*X)*SQRT(SQRT(((1-X)**ABS(MTIMES2+HTIMES2))))*S1
      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSERIES2ERR(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,
     &S2,RELERR2)
C-----------------------------------------------------------------------
C     Calculates the series at a specific point X with the expansion
C     around x=1
C-----------------------------------------------------------------------
C   INPUT: 
C     X-> COS(THETA)
C     C-> a times omega
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     NLAST2-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
C   OUTPUT:
C     S2-> The angular function evaluated at x=X
C     RELERR2-> The relative error evaluated through the first neglected 
C               term
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER MTIMES2,HTIMES2,NLAST2
      DOUBLE PRECISION C,SPHRDLSERIES2(0:NMAX),S2,RELERR2,X

*Local variables:
      INTEGER N
      DOUBLE PRECISION FIRSTNEGLECT
C--Same comments as CHSERIES1ERR
      S2=0.
      RELERR2=0.
      DO 5003 N=0,NLAST2
         S2=S2+SPHRDLSERIES2(N)*
     &SQRT(SQRT((1-X)**(4*N+ABS(MTIMES2+HTIMES2))))  
 5003 CONTINUE
      N=N-1
      FIRSTNEGLECT=SPHRDLSERIES2(N)*
     &SQRT(SQRT((1-X)**(4*N+ABS(MTIMES2+HTIMES2))))

      IF (S2.NE.0) THEN
         RELERR2=ABS(FIRSTNEGLECT/S2)
      END IF

      S2=EXP(C*X)*SQRT(SQRT(((1.+X)**ABS(MTIMES2-HTIMES2))))*S2     

      END

CDECK  ID>, .
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSERIES2(X,C,MTIMES2,HTIMES2,SPHRDLSERIES2,NLAST2,S2)
C-----------------------------------------------------------------------
C     Calculates the series at a specific point X with the expansion
C     around x=1
C-----------------------------------------------------------------------
C   INPUT: 
C     X-> COS(THETA)
C     C-> a times omega
C     MTIMES2-> Twice angular momentum along z quantum number
C     HTIMES2-> "Helicity" component (times 2)
C     SPHRDLSERIES2-> array of coefficients of the expansion
C     NLAST2-> The order of the expansion calculated so far
C-----------------------------------------------------------------------
C   OUTPUT:
C     S2-> The angular function evaluated at x=X
C     RELERR2-> The relative error evaluated through the first neglected 
C               term
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*INPUT and OUTPUT:
      INTEGER MTIMES2,HTIMES2,NLAST2
      DOUBLE PRECISION C,SPHRDLSERIES2(0:NMAX),S2,X
*Local variables:
      INTEGER N

      S2=0D0
      DO 5003 N=0,NLAST2
         S2=S2+SPHRDLSERIES2(N)*
     &SQRT(SQRT((1-X)**(4*N+ABS(MTIMES2+HTIMES2))))  
 5003 CONTINUE

      S2=EXP(C*X)*SQRT(SQRT(((1.+X)**ABS(MTIMES2-HTIMES2))))*S2

      END

CDECK  ID>, CHSPHRDLTQLI.
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C----------------------------------------------------------------------- 
      SUBROUTINE CHSPHRDLTQLI(DIAG,SUBDIAG,N)
C-----------------------------------------------------------------------
C     Calculates the eigenvalues of the symmetric NxN tridiagonal matrix
C     defined by DIAG and SUBDIAG. This is an adaptation of the routine 
C     TQLI by the end of section 11.3 of "Numerical recipes-The art of 
C     scientific computing" by Press, W.H. et all. 
C-----------------------------------------------------------------------
C   I/O: 
C     DIAG-> Array of size NMAX+1 with the diagonal elements
C     SUBDIAG-> Array of size NMAX with the subdiagonal elements
C     N-> Size of the truncated matrix to be diagonalized
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*I/O 
      INTEGER N
      DOUBLE PRECISION DIAG(1:NMAX),SUBDIAG(1:NMAX)
*Local variables
      INTEGER L,M,ITER,I
      DOUBLE PRECISION SMALLSUBDIAG,G,R,S,C,P,F,B

      IF(N.GT.1) THEN
         SUBDIAG(N)=0.
         DO 15 L=1,N
            ITER=0
 1          DO 12 M=L,N-1
               SMALLSUBDIAG=ABS(DIAG(M))+ABS(DIAG(M+1))
               IF((ABS(SUBDIAG(M))+SMALLSUBDIAG).EQ.SMALLSUBDIAG) THEN
                  GOTO 2
               END IF
 12         CONTINUE
            M=N
 2          IF(M.NE.L)THEN
               IF(ITER.EQ.50)THEN
                  WRITE(*,*)'ERROR (In SPHRDLTQLI): TOO MANY ITERATIONS'
                  STOP
               END IF
               ITER=ITER+1
               G=(DIAG(L+1)-DIAG(L))/(2.*SUBDIAG(L))
               R=SQRT(G**2+1.)
               G=DIAG(M)-DIAG(L)+SUBDIAG(L)/(G+SIGN(R,G))
               S=1.
               C=1.
               P=0.
               DO 14 I=M-1,L,-1
                  F=S*SUBDIAG(I)
                  B=C*SUBDIAG(I)
                  IF(ABS(F).GE.ABS(G))THEN
                     C=G/F
                     R=SQRT(C**2+1.)
                     SUBDIAG(I+1)=F*R
                     S=1./R
                     C=C*S
                  ELSE
                     S=F/G
                     R=SQRT(S**2+1.)
                     SUBDIAG(I+1)=G*R
                     C=1./R
                     S=S*C
                  END IF
                  G=DIAG(I+1)-P
                  R=(DIAG(I)-G)*S+2.*C*B
                  P=S*R
                  DIAG(I+1)=G+P
                  G=C*R-B
 14            CONTINUE
               DIAG(L)=DIAG(L)-P
               SUBDIAG(L)=G
               SUBDIAG(M)=0.
               GOTO 1
            END IF
 15      CONTINUE
      END IF
      END

CDECK  ID>, CHSPHRDLQCKSRT.
*CMZ :-        -20/06/07  18.11.30  by  Marco Sampaio
*-- Author :   Marco Sampaio
C-----------------------------------------------------------------------
      SUBROUTINE CHSPHRDLQCKSRT(DIAG,N)
C-----------------------------------------------------------------------
C     Sorts out the elements of DIAG in increasing order, using an 
C     adaptation of the quicksort routine given by the end of section 
C     8.4 of "Numerical recipes-The art of scientific computing" 
C     by Press, W.H. et all. 
C-----------------------------------------------------------------------
C   I/O: 
C     DIAG-> The array of numbers to be sorted
C     N-> The size of the array
C-----------------------------------------------------------------------
*Fixed parameters:
      INTEGER NMAX,NDIVIDEMAX
      PARAMETER(NMAX=250,NDIVIDEMAX=150)
*I/0:
      INTEGER N
      DOUBLE PRECISION DIAG(1:N)
*Local variables:
      INTEGER M,NSTACK
      PARAMETER(M=7)
      INTEGER ISTACK(1:NMAX),JSTACK,L,IR,I,J,IQ
      DOUBLE PRECISION A
      NSTACK=NMAX      

      JSTACK=0
      L=1
      IR=N
 10   IF((IR-L).LT.M)THEN
         DO 13 J=L+1,IR
            A=DIAG(J)
            DO 11 I=J-1,1,-1
               IF(DIAG(I).LE.A) THEN
                  GOTO 12
               END IF
               DIAG(I+1)=DIAG(I)
 11         CONTINUE
            I=0
 12         DIAG(I+1)=A
 13      CONTINUE
         IF(JSTACK.EQ.0)THEN
            GOTO 40
         END IF
         IR=ISTACK(JSTACK)
         L=ISTACK(JSTACK-1)
         JSTACK=JSTACK-2
      ELSE
         I=L
         J=IR
         IQ=L+INT((IR-L)*RAND())
         A=DIAG(IQ)
         DIAG(IQ)=DIAG(L)
 20      CONTINUE
 21      IF(J.GT.0)THEN
            IF(A.LT.DIAG(J))THEN
               J=J-1
               GOTO 21
            END IF
         END IF
         IF(J.LE.I)THEN
            DIAG(I)=A
            GOTO 30
         END IF
         DIAG(I)=DIAG(J)
         I=I+1
 22      IF(I.LE.N) THEN
            IF(A.GT.DIAG(I))THEN
               I=I+1
               GOTO 22
            END IF
         END IF
         IF(J.LE.I)THEN
            DIAG(J)=A
            I=J
            GOTO 30
         END IF
         DIAG(J)=DIAG(I)
         J=J-1
         GOTO 20
 30      JSTACK=JSTACK+2
         IF(JSTACK.GT.NSTACK) THEN
            WRITE(*,*)'ERROR (in SPHRDLQCKSRT): NSTACK must be larger'
            STOP
         END IF
         IF((IR-1).GE.(I-L))THEN
            ISTACK(JSTACK)=IR
            ISTACK(JSTACK-1)=I+1
            IR=I-1
         ELSE
            ISTACK(JSTACK)=I-1
            ISTACK(JSTACK-1)=L
            L=I+1
         END IF
      END IF
      GOTO 10

 40   CONTINUE
      END

CDECK  ID>, CHSCGC.
*CMZ :-        -22/08/06  15.58.30  by  Bryan Webber
*-- Author :    Arturo Quirantes Sierra
C------------------------------------------------------------------
      SUBROUTINE CHSCGC(AJ,BJ,CJ,AM,BM,CM,CG)
C	SUBROUTINE NED(AJ,BJ,CJ,AM,BM,CM,CG)
C	ARTURO QUIRANTES SIERRA
C	Department of Applied Physics, Faculty of Sciences
C	University of Granada, 18071 Granada (SPAIN)
C	http://www.ugr.es/local/aquiran/codigos.htm
C	aquiran@ugr.es
C
C	Last update: 20 May 2.003	
C
C	Subroutine NED
C	to calculate Clebsch-Gordan coefficients
C
C	You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine
C	Input:
C		AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices)
C	Output:
C		CG=C-G(AJ,BJ,CJ,AM,BM,CM)
C
C------------------------------------------------------------------
C--BW fix 20/08/06: replaced obsolete arithmetic if statements
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER QSIZE
      PARAMETER (QSIZE=500)
      DIMENSION Q(QSIZE,QSIZE)
      DOUBLE PRECISION CG
      DOUBLE PRECISION AJ,BJ,CJ,AM,BM,CM
      INTEGER ZZ
      ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
      IF (ZZ.GT.QSIZE) THEN
         WRITE(*,*)'CHSCGC: NED returns CG=0: needs an array size'//
     &' of at least',ZZ
         WRITE(*,*)'If you are generating events with very'//
     &' large angular momentum you will need to increase the QSIZE'//
     &' parameter in the subroutine CHSGC of the charybdis*.F file '
     &//'to a larger value.'
         WRITE(*,*)''
           CG=0D0
           RETURN
        ENDIF
	DO 2 I=1,ZZ
           Q(I,1)=1.D0
           Q(I,I)=1.D0
 2      CONTINUE
	DO 3 I=2,ZZ-1
        DO 3 K=2,I
           Q(I+1,K)=Q(I,K-1)+Q(I,K)
 3      CONTINUE
	CG=0.D0
	JA=AJ+AM+1.01D0
	MA=AJ-AM+1.01D0
	JB=BJ+BM+1.01D0
	MB=BJ-BM+1.01D0
	JC=CJ+CM+1.01D0
	MC=CJ-CM+1.01D0
	LA=BJ+CJ-AJ+1.01D0
	LB=CJ+AJ-BJ+1.01D0
	LC=AJ+BJ-CJ+1.01D0
	LT=AJ+BJ+CJ+1.01D0
	D=DABS(AM+BM-CM)-0.01D0
	IF (D.LE.0.D0) THEN
           LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
           IF (LD.GT.0) THEN
              JA2=AJ+AJ+AM+AM
              JB2=BJ+BJ+BM+BM
              JC2=CJ+CJ-CM-CM
              I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
              IF (I2.EQ.0) THEN
                 FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
                 FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
                 FN=FN/Q(JA+MA-1,JA)
                 FN=FN/Q(JB+MB-1,JB)
                 FN=FN/Q(JC+MC-1,JC)
                 K0=MAX(0,LC-JA,LC-MB)+1
                 K1=MIN(LC,MA,JB)
                 X=0.D0
                 DO 50 K=K0,K1
                    X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
 50              CONTINUE
                 IP=K1+LB+JC
                 P=1-2*(IP-IP/2*2)
                 CG=P*X*DSQRT(FN)
C	What we've calculated is a Wigner 3-j coefficient
C	Next, we'll turn it into a Clebsch-Gordan coefficient
                 CG=CG*DSQRT(2*CJ+1)*(-1)**IDNINT(AJ-BJ-CM)
              ENDIF
           ENDIF
        ENDIF
	END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C          Parameter calculation and particle properties Routines
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, CHPARA.
*CMZ :-        -22/08/06  15.58.30  by  Bryan Webber
*-- Author :    Bryan Webber
C------------------------------------------------------------------
      SUBROUTINE CHPARA(MBH,JBH,MPL,NEX,RH,OBLA,OMEG,THAW)
C------------------------------------------------------------------
C--COMPUTE BH PARAMETERS FROM ITS MASS AND ANGR MOM
C------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
      DOUBLE PRECISION MBH,JBH,MPL,PIFAC,RTPI,CRS(0:7),PEX,AST,BST,
     & CST,DST,FST,SST,TINY,RS,RH,OBLA,OMEG,THAW
      INTEGER NEX,NTRY,MTRY
      LOGICAL FIRST
      DATA TINY,MTRY,FIRST/1D-8,100,.TRUE./
      SAVE FIRST,PIFAC,RTPI,CRS
      IF (FIRST) THEN
         PIFAC=ACOS(-1D0)
         RTPI=SQRT(PIFAC)
C----Prefactors for the expression of the radius dependent on the number of extra dimensions
         CRS(0)=1D0/(2*PIFAC)
         CRS(1)=4D0/(3*PIFAC)
         CRS(2)=3D0/2D0
         CRS(3)=32D0/5D0
         CRS(4)=10D0*PIFAC
         CRS(5)=384D0/7D0*PIFAC
         CRS(6)=105D0*PIFAC**2
         CRS(7)=2048D0*PIFAC**2/3D0
         FIRST=.FALSE.
      ENDIF
      PEX=1D0/(NEX+1D0)
C--SCHWARZSCHILD RADIUS
      RS=(CRS(NEX)*MBH/MPL)**PEX/MPL

C--ANGULAR MOMENTUM PARAMETER
      BST=(NEX+2D0)*JBH/(2D0*RS*MBH)
      IF ((NEX.EQ.0.AND.BST.GT.0.5D0)
     &.OR.(NEX.EQ.1.AND.BST.GT.1.0D0)) THEN
         PRINT *,'CHPARA: JBH TOO LARGE (ABOVE UPPER BOUND IN TOTDIM=',
     &NEX+4,'), BST=',BST
	 PRINT *,'NEX,JBH,RS,MBH =',NEX,JBH,RS,MBH
         STOP
      ENDIF
C--SOLVING AST/(1+AST**2)**PEX=BST (first guess in AST <<1 approximation)
      AST=BST
      IF (BST.EQ.0D0) GOTO 20
      DO NTRY=1,MTRY
         FST=AST/(1D0+AST**2)**PEX
         DST=FST*(1D0/AST-2D0*AST*PEX/(1D0+AST**2))
         SST=(BST-FST)/DST
         AST=AST+SST
         IF (SST.LT.TINY) GOTO 20
      ENDDO
C--Trying again with first guess in AST >> 1 approximation 
      AST=BST**(-1/(2*PEX-1))
      DO NTRY=1,MTRY
         FST=AST/(1D0+AST**2)**PEX
         DST=FST*(1D0/AST-2D0*AST*PEX/(1D0+AST**2))
         SST=(BST-FST)/DST
         AST=AST+SST
         IF (SST.LT.TINY) GOTO 20
      ENDDO
      IF(MBH.LT.RMMINM) THEN
C--Tiny masses give large astar values whose exact covergence are beyond numerical precision. The remnant phase which will follow immediately does not use it. So we can safely use the last best guess.
         GOTO 20
      END IF
      PRINT *,'CHPARA: NO CONVERGENCE, AST,BST=',AST,BST
      STOP
 20   CST=1D0+AST**2
C--HORIZON RADIUS
      RH=RS/CST**PEX
C--OBLATENESS
      OBLA=AST*RH
C--ANGULAR VELOCITY
      OMEG=AST/(CST*RH)
C--HAWKING TEMPERATURE
#if !defined(GINGRICH)      
      THAW=(NEX+1D0+(NEX-1D0)*AST**2)/(4D0*PIFAC*CST*RH)
#else
      IF (MBH .GE. DGMS/(DGGS*DGGS) .OR. (.NOT.DGSB)) THEN
        THAW=(NEX+1D0+(NEX-1D0)*AST**2)/(4D0*PIFAC*CST*RH)
      ELSE
        THAW=DGMS/(2.0D0*DSQRT(2.0D0)*PIFAC)
      ENDIF
#endif
      END

CDECK  ID>, BSTNOW
*CMZ :-        -19/11/09  15.58.30  by  Marco Sampaio
*-- Author :    Marco Sampaio
C------------------------------------------------------------------
      FUNCTION BSTNOW(MBH,JBH,MPL,NEX)
C------------------------------------------------------------------
C--COMPUTE BH PARAMETERS FROM ITS MASS AND ANGR MOM
C------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
      DOUBLE PRECISION BSTNOW,MBH,JBH,MPL,PIFAC1,RTPI1,CRS1(0:7),PEX,RS
      INTEGER NEX
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST,PIFAC1,RTPI1,CRS1
      IF (FIRST) THEN
         PIFAC1=ACOS(-1D0)
         RTPI1=SQRT(PIFAC1)
C----Prefactors for the expression of the radius dependent on the number of extra dimensions
         CRS1(0)=1D0/(2*PIFAC1)
         CRS1(1)=4D0/(3*PIFAC1)
         CRS1(2)=3D0/2D0
         CRS1(3)=32D0/5D0
         CRS1(4)=10D0*PIFAC1
         CRS1(5)=384D0/7D0*PIFAC1
         CRS1(6)=105D0*PIFAC1**2
         CRS1(7)=2048D0*PIFAC1**2/3D0
         FIRST=.FALSE.
      ENDIF
      PEX=1D0/(NEX+1D0)
C--SCHWARZSCHILD RADIUS
      RS=(CRS1(NEX)*MBH/MPL)**PEX/MPL

C--ANGULAR MOMENTUM PARAMETER
      BSTNOW=(NEX+2D0)*JBH/(2D0*RS*MBH)

      END

CDECK  ID>, BSTNEW
*CMZ :-        -19/11/09  15.58.30  by  Marco Sampaio
*-- Author :    Marco Sampaio
C------------------------------------------------------------------
      FUNCTION BSTNEW(MBH,JBH,MPL,NEX,EEM,MSEM)
C------------------------------------------------------------------
C--COMPUTE BH PARAMETERS FROM ITS MASS AND ANGR MOM
C------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
      DOUBLE PRECISION BSTNEW,MBH,JBH,MPL,PIFAC1,RTPI1,CRS1(0:7),PEX,RS,
     &EEM,MSEM,MBHLOCAL
      INTEGER NEX
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST,PIFAC1,RTPI1,CRS1
      IF (FIRST) THEN
         PIFAC1=ACOS(-1D0)
         RTPI1=SQRT(PIFAC1)
C----Prefactors for the expression of the radius dependent on the number of extra dimensions
         CRS1(0)=1D0/(2*PIFAC1)
         CRS1(1)=4D0/(3*PIFAC1)
         CRS1(2)=3D0/2D0
         CRS1(3)=32D0/5D0
         CRS1(4)=10D0*PIFAC1
         CRS1(5)=384D0/7D0*PIFAC1
         CRS1(6)=105D0*PIFAC1**2
         CRS1(7)=2048D0*PIFAC1**2/3D0
         FIRST=.FALSE.
      ENDIF
      PEX=1D0/(NEX+1D0)
C--SCHWARZSCHILD RADIUS
      IF(RECOIL.EQ.1)THEN
         MBHLOCAL=SQRT(MBH**2-2*MBH*EEM+MSEM**2)
      ELSE IF(RECOIL.EQ.2)THEN
         MBHLOCAL=MBH-EEM
      ELSE
         WRITE(*,*)"INVALID RECOIL OPTION"
         STOP
      END IF
      RS=(CRS1(NEX)*MBHLOCAL/MPL)**PEX/MPL

C--ANGULAR MOMENTUM PARAMETER
      BSTNEW=(NEX+2D0)*JBH/(2D0*RS*MBHLOCAL)

      END

CDECK  ID>, CHFCHG.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Peter Richardson
C----------------------------------------------------------------------
      FUNCTION CHFCHG(ID)
C----------------------------------------------------------------------
C     Function to return the charge of a SM particle input is PDG code
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER CHFCHG,ID
      IF(ABS(ID).LE.6) THEN
         IF(MOD(ID,2).EQ.0) THEN
            CHFCHG = 2
         ELSE
            CHFCHG = -1
         ENDIF
      ELSEIF(ABS(ID).GE.11.AND.ABS(ID).LE.16) THEN
         IF(MOD(ID,2).EQ.0) THEN
            CHFCHG = 0
         ELSE
            CHFCHG = -3
         ENDIF
      ELSEIF((ID.GE.21.AND.ID.LE.23).OR.ID.EQ.25) THEN
         CHFCHG = 0
      ELSEIF(ABS(ID).EQ.24) THEN
         CHFCHG = 3
#if defined (GINGRICH)
      ELSEIF(ID.EQ.39) THEN
         CHFCHG = 0
#endif
      ELSE
         WRITE(*,*)'CHFCHG WARNING UNKNOWN PARTICLE',ID
         STOP
      ENDIF
      IF(ID.LT.0) CHFCHG = -CHFCHG
      END

CDECK  ID>, CHFMAS.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Peter Richardson
C----------------------------------------------------------------------
      FUNCTION CHFMAS(ID)
C----------------------------------------------------------------------
C     Function to return the mass of a SM particle input is PDG code
C     This is generator dependent
C----------------------------------------------------------------------
#if defined(PYTHIA)
C--get the masses from the PYTHIA common block
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      INTEGER PYCOMP
      EXTERNAL PYCOMP
#elif defined(HERWIG)
      INCLUDE 'HERWIG65.INC'
#endif
      INTEGER ID
      DOUBLE PRECISION CHFMAS
#if defined(PYTHIA)
      CHFMAS  = PMAS(PYCOMP(ID),1)
#elif defined(HERWIG)
C--get the masses from the HERWIG common block
      IF(ABS(ID).LE.6) THEN
         CHFMAS = RMASS(ABS(ID))
      ELSEIF(ABS(ID).GE.11.AND.ABS(ID).LE.16) THEN
         CHFMAS = RMASS(110+ABS(ID))
      ELSEIF(ID.EQ.21) THEN
         CHFMAS = RMASS(13)
      ELSEIF(ID.EQ.22) THEN
         CHFMAS = 0.0D0
      ELSEIF(ID.EQ.23) THEN
         CHFMAS = RMASS(200)
      ELSEIF(ABS(ID).EQ.24) THEN
         CHFMAS = RMASS(198)
      ELSEIF(ID.EQ.25) THEN
         CHFMAS = RMASS(201)
#if defined(GINGRICH)
      ELSEIF(ID.EQ.39) THEN
         CHFMAS = 0.0D0
#endif
      ELSE
         PRINT *,'CHFMAS: UNKNOWN PARTICLE',ID
         STOP
      ENDIF
#else
C--Masses put in by hand
      CHFMAS=0D0
      IF(ABS(ID).EQ.1)THEN
         CHFMAS=0.3200D0
      ELSEIF(ABS(ID).EQ.2)THEN
         CHFMAS=0.3200D0
      ELSEIF(ABS(ID).EQ.3)THEN
         CHFMAS=0.5000D0
      ELSEIF(ABS(ID).EQ.4)THEN 
         CHFMAS=1.5500D0
      ELSEIF(ABS(ID).EQ.5)THEN 
         CHFMAS=4.9500D0
      ELSEIF(ABS(ID).EQ.6)THEN 
         CHFMAS=174.30D0
      ELSEIF(ABS(ID).EQ.11)THEN 
         CHFMAS=5.11D-04
      ELSEIF(ABS(ID).EQ.13)THEN 
         CHFMAS=.10566D0
      ELSEIF(ABS(ID).EQ.15)THEN 
         CHFMAS=1.7771D0
      ELSEIF(ABS(ID).EQ.23)THEN 
         CHFMAS=91.188D0
      ELSEIF(ABS(ID).EQ.24)THEN 
         CHFMAS=80.420D0
      ELSEIF(ABS(ID).EQ.25)THEN 
         CHFMAS=115.00D0
#if !defined(GINGRICH)
      ELSE IF (ABS(ID).GT.25)THEN
#else
      ELSEIF(ID.EQ.39) THEN
         CHFMAS = 0.0D0
      ELSE IF (ABS(ID).GT.25)THEN
#endif
         PRINT *,'CHFMAS: UNKNOWN PARTICLE',ID
         STOP
      ENDIF
#endif
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                          Kinematics Routines
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, CHULB4.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHULB4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(1)*PS(1)+PI(2)*PS(2)
     &         +PI(3)*PS(3)+PI(4)*PS(4))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) + FN*PS(1)
        PF(2)= PI(2) + FN*PS(2)
        PF(3)= PI(3) + FN*PS(3)
        PF(4)= PF4
      END IF
      END

CDECK  ID>, CHULOB.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHULOB(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4,5) = (PX,PY,PZ,E,M)
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION PS(5),PI(5),PF(5)
      CALL CHULB4(PS,PI,PF)
      PF(5)= PI(5)
      END
CDECK  ID>, CHUMAS.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHUMAS(P)
C-----------------------------------------------------------------------
C     PUTS INVARIANT MASS IN 5TH COMPONENT OF VECTOR
C     (NEGATIVE SIGN IF SPACELIKE) (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHUSQR,P(5)
      EXTERNAL CHUSQR
      P(5)=CHUSQR((P(4)+P(3))*(P(4)-P(3))-P(1)**2-P(2)**2)
      END
CDECK  ID>, CHUPCM.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      FUNCTION CHUPCM(EM0,EM1,EM2)
C-----------------------------------------------------------------------
C     C.M. MOMENTUM FOR DECAY MASSES EM0 -> EM1 + EM2
C     SET TO -1 BELOW THRESHOLD (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHUPCM,EM0,EM1,EM2,EMS,EMD
      EMS=ABS(EM1+EM2)
      EMD=ABS(EM1-EM2)
      IF (EM0.LT.EMS.OR.EM0.LT.EMD) THEN
        CHUPCM=-1.0D0
      ELSEIF (EM0.EQ.EMS.OR.EM0.EQ.EMD) THEN
        CHUPCM=0.0D0
      ELSE
        CHUPCM=SQRT((EM0+EMD)*(EM0-EMD)*
     &              (EM0+EMS)*(EM0-EMS))*0.5D0/EM0
      ENDIF
      END
CDECK  ID>, CHUROB.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHUROB(R,P,Q)
C-----------------------------------------------------------------------
C     ROTATES VECTORS BY INVERSE OF ROTATION MATRIX R 
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION S1,S2,S3,R(3,3),P(3),Q(3)
      S1=P(1)*R(1,1)+P(2)*R(2,1)+P(3)*R(3,1)
      S2=P(1)*R(1,2)+P(2)*R(2,2)+P(3)*R(3,2)
      S3=P(1)*R(1,3)+P(2)*R(2,3)+P(3)*R(3,3)
      Q(1)=S1
      Q(2)=S2
      Q(3)=S3
      END
CDECK  ID>, CHUROT.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHUROT(P,CP,SP,R)
C-----------------------------------------------------------------------
C     R IS ROTATION MATRIX TO GET FROM VECTOR P TO Z AXIS, FOLLOWED BY
C     A ROTATION BY PSI ABOUT Z AXIS, WHERE CP = COS-PSI, SP = SIN-PSI
C     (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION WN,CP,SP,PTCUT,PP,PT,CT,ST,CF,SF,P(3),R(3,3)
      DATA WN,PTCUT/1.D0,1.D-20/
      PT=P(1)**2+P(2)**2
      PP=P(3)**2+PT
      IF (PT.LE.PP*PTCUT) THEN
         CT=SIGN(WN,P(3))
         ST=0.0D0
         CF=1.0D0
         SF=0.0D0
      ELSE
         PP=SQRT(PP)
         PT=SQRT(PT)
         CT=P(3)/PP
         ST=PT/PP
         CF=P(1)/PT
         SF=P(2)/PT
      END IF
      R(1,1)= CP*CF*CT+SP*SF
      R(1,2)= CP*SF*CT-SP*CF
      R(1,3)=-CP*ST
      R(2,1)=-CP*SF+SP*CF*CT
      R(2,2)= CP*CF+SP*SF*CT
      R(2,3)=-SP*ST
      R(3,1)= CF*ST
      R(3,2)= SF*ST
      R(3,3)= CT
      END
CDECK  ID>, CHUSQR.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      FUNCTION CHUSQR(X)
C-----------------------------------------------------------------------
C     SQUARE ROOT WITH SIGN RETENTION (extracted from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHUSQR,X
      CHUSQR=SIGN(SQRT(ABS(X)),X)
      END

CDECK  ID>, CHVDIF.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHVDIF(N,P,Q,R)
C-----------------------------------------------------------------------
C     VECTOR DIFFERENCE (extracted from HERWIG)
C-----------------------------------------------------------------------
      DOUBLE PRECISION P(N),Q(N),R(N)
      INTEGER N,I
      DO 10 I=1,N
   10 R(I)=P(I)-Q(I)
      END
CDECK  ID>, CHVEQU.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHVEQU(N,P,Q)
C-----------------------------------------------------------------------
C     VECTOR EQUALITY (extracted from HERWIG)
C-----------------------------------------------------------------------
      DOUBLE PRECISION P(N),Q(N)
      INTEGER N,I
      DO 10 I=1,N
   10 Q(I)=P(I)
      END
CDECK  ID>, CHVSUM.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHVSUM(N,P,Q,R)
C-----------------------------------------------------------------------
C    VECTOR SUM (extracted from HERWIG)
C-----------------------------------------------------------------------
      DOUBLE PRECISION P(N),Q(N),R(N)
      DO 10 I=1,N
   10 R(I)=P(I)+Q(I)
      END

CDECK  ID>, CHRAZM
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE CHRAZM(PT,PX,PY)
C-----------------------------------------------------------------------
C     RANDOMLY ROTATED 2-VECTOR (PX,PY) OF LENGTH PT
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRGEN,PT,PX,PY,C,S,CS,QT,ONE,ZERO
      PARAMETER(ONE=1.0D0, ZERO=0.0D0)
      EXTERNAL CHRGEN
   10 C=2.0D0*CHRGEN(1)-1.0D0
      S=2.0D0*CHRGEN(2)-1.0D0
      CS=C*C+S*S
      IF (CS.GT.ONE .OR. CS.EQ.ZERO) GOTO 10
      QT=PT/CS
      PX=(C*C-S*S)*QT
      PY=2.0D0*C*S*QT
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                        Random number generators
C----------------------------------------------------------------------
C----------------------------------------------------------------------

CDECK  ID>, CHRGEN.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      FUNCTION CHRGEN(I)
C-----------------------------------------------------------------------
C     Random number generator (just calls HERWIG or PYTHIA one)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRGEN
      INTEGER I
#if defined(PYTHIA)
C--use the PYTHIA one
      DOUBLE PRECISION PYR
      EXTERNAL PYR
      CHRGEN =  PYR(I)
#else
C--use HERWIG random number generator
      DOUBLE PRECISION HWRGEN
      EXTERNAL HWRGEN
      CHRGEN = HWRGEN(I)
#endif
      END

CDECK  ID>, CHRLOG.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      FUNCTION CHRLOG(A)
C-----------------------------------------------------------------------
C     Returns .TRUE. with probability A
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CHRGEN,A,R
      LOGICAL CHRLOG
      EXTERNAL CHRGEN
      CHRLOG=.TRUE.
      R=CHRGEN(1)
      IF(R.GT.A) CHRLOG=.FALSE.
      END
CDECK  ID>, CHRUNI.
*CMZ :-        -17/07/03  18.11.30  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      FUNCTION CHRUNI(I,A,B)
C-----------------------------------------------------------------------
C     Uniform random random number in range [A,B]
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      DOUBLE PRECISION CHRUNI,CHRGEN,A,B,RN
      INTEGER I
      EXTERNAL CHRGEN
      RN=CHRGEN(I)
      CHRUNI=A+RN*(B-A)
      END

C---Random number generator taken from herwig for standalone charybdis
#if defined(CHARYBDIS)
CDECK  ID>, HWRSET.
*CMZ :-        -26/04/91  12.42.30  by  Federico Carminati
*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION HWRSET(JSEED)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     SETTING SEEDS (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION HWRSET
      COMMON/HWSEED/ISEED(2)
      INTEGER ISEED
      INTEGER JSEED(2)
      HWRSET=0.0D0
      IF (JSEED(1).EQ.0.OR.JSEED(2).EQ.0) THEN
        WRITE(*,*)'INVALID SEED'
        STOP
      ENDIF
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
      END

CDECK  ID>, HWRGEN.
*CMZ :-        -26/04/91  12.42.30  by  Federico Carminati
*-- Author :    F. James, modified by Mike Seymour
*- Split in 3 files by M. Kirsanov. Initial seeds ISEED set in HWUDAT
C-----------------------------------------------------------------------
      FUNCTION HWRGEN(I)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C     (taken from HERWIG)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION HWRGEN
      COMMON/HWSEED/ISEED(2)
      INTEGER ISEED
      INTEGER I,K,IZ

      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      HWRGEN=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      END
#endif
#if defined(GINGRICH)
C-----------------------------------------------------------------------
      SUBROUTINE CHDN(P0,P1,P2,PCM)
C-----------------------------------------------------------------------
C     GENERATES DECAY 0 -> 1+2 in N-DIMENSIONS
C     Doug Gingrich
C     PCM IS CM MOMENTUM
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'charybdis2.inc'
      INTEGER I,J

      DOUBLE PRECISION PCM,P0(12),P1(12),P2(5),PP(12)
      DOUBLE PRECISION CHRGEN,PI,R,X1,X2
      DOUBLE PRECISION PF4,FN,X(10),PN2(12)
      EXTERNAL CHRGEN
      PARAMETER (PI=3.1416)
C
C.....Generate random point on TOTDIM-1 sphere.
      R = 0.0
      DO I=1,TOTDIM/2
        J = 2*I-1
C.......Get 2 random numbers at a time.
        X1 = CHRGEN(1)
        X2 = CHRGEN(1)
        X(J) = SQRT(-2.0*LOG(X1))*COS(2.0*PI*X2)
        R = R + X(J)*X(J)
        IF (J .NE. TOTDIM-1) THEN
          X(J+1) = SQRT(-2.0*LOG(X1))*SIN(2.0*PI*X2)
          R = R + X(J+1)*X(J+1)
        ENDIF
      ENDDO
C
C.....Normalize to unit sphere.
      R = SQRT(R)
      DO I=1,TOTDIM-1
        X(I) = X(I) / R
      ENDDO
C
C.....Kinematics in CMS system of graviton (particle 2).
      PP(1) = PCM*X(1)
      PP(2) = PCM*X(2)
      PP(3) = -PCM*X(3)
      PP(4) = SQRT(P2(5)**2 + PCM**2)
      PP(5) = P2(5)
      DO I=1,TOTDIM-4
        PP(5+I) = PCM*X(3+I)
      ENDDO
C
C.....Boost from CMS to lab frame.
      IF (P0(4) .EQ. P0(5)) THEN
        DO I=1,TOTDIM+1
          IF (I .NE. 5) PN2(I) = PP(I)
        ENDDO
      ELSE
        PF4 = 0
        DO I=1,TOTDIM+1
          IF (I .NE. 5) PF4 = PF4 + PP(I)*P0(I)
        ENDDO
        PF4 = PF4/P0(5)
        FN   = (PF4 + PP(4)) / (P0(4) + P0(5))
        DO I=1,TOTDIM+1
          IF (I .EQ. 4) THEN
            PN2(4) = PF4
          ELSEIF (I .EQ. 5) THEN
            PN2(5)= PP(5)
          ELSE
            PN2(I) = PP(I) + FN*P0(I)
          ENDIF
       ENDDO
      ENDIF
C
C.....Kinematics in lab of black hole (particle 1).
      DO I=1,TOTDIM+1
        IF (I .NE. 5) P1(I) = P0(I) - PN2(I)
      ENDDO
      DO I=1,5
        P2(I) = PN2(I)
      ENDDO
C
      RETURN
      END
#endif
