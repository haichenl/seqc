MEX    := /home/haichen/program/bin/MATLAB/R2013a/bin/mex

MEXFLAGS := -O2 -fpic

MEXSUFFIX=$(shell mexext)

CNDOSRC1=Cpp_Cndo2CalcCartesianMatrixElementsByGTOExpansionLoop
MNDOSRC1=Cpp_MndoGetNddoRepulsionIntegral

build: $(CNDOSRC1).($MEXSUFFIX) $(MNDOSRC1).($MEXSUFFIX)
	cp $(CNDOSRC1).$(MEXSUFFIX) ../+SEQC/@Cndo2
	cp $(MNDOSRC1).$(MEXSUFFIX) ../+SEQC/@Mndo

%.($MEXSUFFIX): %.cpp
	${MEX} CPPFLAGS="${MEXFLAGS}" $<

clean: 
	rm *.$(MEXSUFFIX) ../+SEQC/@Cndo2/$(CNDOSRC1).$(MEXSUFFIX) ../+SEQC/@Mndo/$(MNDOSRC1).$(MEXSUFFIX)

