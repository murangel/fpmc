# 2013/03 On Lxplus, one needs 64bit libraries

#export CERN=/afs/cern.ch/sw/lcg/external/cernlib/2006
#export CERN_LEVEL=slc4_amd64_gcc4
#export PATH=$CERN/$CERN_LEVEL/bin:$PATH

setenv LHAPDFLIB  /afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.8/x86_64-slc5-gcc46-opt/lib/
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.8/x86_64-slc5-gcc46-opt/lib/
