all: tracking

CC	         = mpicxx
CXX	         = mpicxx

DEPDIR           := deps
LDFLAGS		 = 
FFLAGS		 = 

CPPFLAGS         = $(CXXFLAGS_EXT) $(PETSC_INC) -DPETSC_VER=$(PETSC_VER) -std=c++11 -Wno-unused-result

FPPFLAGS         =
LOCDIR		 = 

SOURCE0          = $(wildcard *.cpp)
# Remove dataPartExtr.cpp from the source files.

DEP_DIR          = deps
DEPENDS          = $(addprefix $(DEP_DIR)/, $(patsubst %.cpp, %.d, $(SOURCE0)))
-include           $(DEPENDS)

# OBJSC            = $(patsubst %.cpp, %.o, $(SOURCE0))
# $(OBJSC): compile makefile.seawulf

SOURCE           = $(filter-out dataPartExtr.cpp reyTest.cpp, $(SOURCE0))
OBJSC            = $(patsubst %.cpp, %.o, $(SOURCE))

LIBDIR           = $(PETSC_DIR) $(HYPREDIR) -L $(BLASDIR) -L $(LAPACK_DIR)

LIBFLAG	 = -lpthread -l rt -l dl $(PETSC_LIB) $(HYPRELIB) $(LIB_EXTR) $(LAPACK_LIB) -l $(BLASLIB) -l gfortran

tracking: ${OBJSC}
	$(CXX) -o tracking $(CPPFLAGS) $(OBJSC) $(LIBDIR) $(LIBFLAG)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<
	if [ ! -d $(DEP_DIR) ]; then mkdir $(DEP_DIR); fi
	$(CXX) -MM $(CPPFLAGS) $< > $(DEP_DIR)/$*.d


clean:
	rm -f deps/*.d
	rm -f *.o
	rm -f tracking


dataPartExtr: dataPartExtr.o ${OBJSC}
	$(CXX) -o dataPartExtr $(CPPFLAGS) $(filter-out main.o, $(OBJSC)) dataPartExtr.o tecio64.a $(LIBDIR) $(LIBFLAG)
