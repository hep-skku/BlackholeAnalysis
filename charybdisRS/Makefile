#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#---------------------Makefile for black hole generator----------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#------------User defined options: Please adapt to your needs----------------
#----------------------------------------------------------------------------
# Choice of implementation: Here you can choose to run in CHARYBDIS, HERWIG or PYTHIA mode:
#	GENERATOR = CHARYBDIS compiles standalone parton level charybdis,
#	GENERATOR = HERWIG interfaces to herwig  
# 	GENERATOR = PYTHIA interfaces to pythia. 
GENERATOR = CHARYBDIS
#----------------------------------------------------------------------------
# Choice to use extra features such as string ball model and D-dimensional graviton, by setting to GINGRICH. Leave blank otherwise. 
EXTRA =
#----------------------------------------------------------------------------
# Name of the CHARYBDIS source file
#CHSRC = original/charybdis2-1.0.3.F   
CHSRC = charybdis2-1.0.3.F   
#----------------------------------------------------------------------------
# Name of the HERWIG source file
HERWIGSRC = herwig6510.f
#----------------------------------------------------------------------------
# Name of the PYTHIA source file
PYTHIASRC = pythia-6.4.27.f   
#----------------------------------------------------------------------------
# Option to use external pdfs. You can choose:
#	LHAPDF (for standalone mode only); 
#	PDFLIB 
#	or leave blank to use herwig or pythia internal pdfs
PDFCHOICE = LHAPDF
#----------------------------------------------------------------------------
# Location of your pdf library (Uncomment one of the lines as appropriate and
# adapt the path name)
#
##PDFLIB  examples:
#LPDF   = -L/cern/pro/lib -lpdflib -lpacklib -lmathlib -lkernlib 
# with some installations of cernlib this may need to be changed to
#LPDF   = -L/cern/pro/lib -lpdflib804 -lpacklib -lmathlib -lkernlib
#
##LHAPDF example: 
#LPDF   = -L/lha_pdf_instalation_path/lib/ -lLHAPDF
LPDF   = -L/store/tools/sw/slc6_x86_64_gcc491/lib -lLHAPDF
#----------------------------------------------------------------------------
# Set compiler
f77 = gfortran -I.
#----------------------------------------------------------------------------
#------END of user defined options (no need to change anything below)--------
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
#----------------------Compiler flags and root names-------------------------
#----------------------------------------------------------------------------
# GENFLAG defines the macro used to identify which compilation of the source code to perform, SRC is the source code filename for the general purpose event generator (herwig or pythia), OBJS the corresponding object code filename, MAIN the mainprogram file name.
#----------------------------------------------------------------------------
# HERWIG implementation
ifeq ("$(GENERATOR)","HERWIG")
GENFLAG = -DHERWIG
SRC   = $(HERWIGSRC)
OBJS  = $(patsubst %.f, %.o, $(HERWIGSRC))
MAIN  = mainherwig.o
endif
#----------------------------------------------------------------------------
# PYTHIA implementation
ifeq ("$(GENERATOR)","PYTHIA")
GENFLAG = -DPYTHIA
SRC   = $(PYTHIASRC)
OBJS  = $(patsubst %.f, %.o, $(PYTHIASRC))
MAIN  = mainpythia.o
endif
#----------------------------------------------------------------------------
# Standalone parton level implementation
ifeq ("$(GENERATOR)","CHARYBDIS")
GENFLAG = -DCHARYBDIS
SRC   = 
OBJS  = 
MAIN  = maincharybdis.o
endif	
#----------------------------------------------------------------------------
# Macros for the executable (CHEXE) and object file (CHOBJ)
#----------------------------------------------------------------------------
CHOBJ = $(patsubst %.F, %.o, $(CHSRC))
CHEXE = $(patsubst %.F, %.exe, $(CHSRC))
#----------------------------------------------------------------------------
# PDF library flags for linking PDFLIB or LHAPDF
ifeq ("$(PDFCHOICE)","PDFLIB")
PDFFLAG = -D$(PDFCHOICE)
LIB = $(LPDF)
else
PDFFLAG = 
LIB   =
endif 
ifeq ("$(PDFCHOICE)","LHAPDF")
PDFFLAG = -D$(PDFCHOICE)
LIB = $(LPDF)
endif
#----------------------------------------------------------------------------
# Define flag for extra features
ifeq ("$(EXTRA)","GINGRICH")
EXTRAFLAG = -DGINGRICH
else
EXTRAFLAG = 
endif	
#----------------------------------------------------------------------------

# Rules to create executable
default: $(CHEXE)
$(CHEXE): $(OBJS) $(CHOBJ) $(MAIN)
	$(f77) -o $(CHEXE) $(CHOBJ) $(OBJS) $(MAIN) $(LIB)
$(OBJS): $(SRC)
	$(f77) -c $(SRC)	
$(CHOBJ): $(CHSRC)
	$(f77) -c $(CHSRC) $(GENFLAG) $(PDFFLAG) $(EXTRAFLAG)
maincharybdis.o: maincharybdis.f
	$(f77) -c maincharybdis.f $(GENFLAG)
mainherwig.o: mainherwig.f
	$(f77) -c mainherwig.f $(GENFLAG)
mainpythia.o: mainpythia.f
	$(f77) -c mainpythia.f $(GENFLAG)

# Remove executable and object files
clean:
	rm charybdis2*.exe *.o
