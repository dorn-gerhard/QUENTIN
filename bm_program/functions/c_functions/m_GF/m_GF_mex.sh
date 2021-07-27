MATLAB="/afs/itp.tugraz.at/opt/matlab/R2018b"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/afs/itp.tugraz.at/user/dorn/.matlab/R2018b"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for m_GF" > m_GF_mex.mki
echo "CC=$CC" >> m_GF_mex.mki
echo "CFLAGS=$CFLAGS" >> m_GF_mex.mki
echo "CLIBS=$CLIBS" >> m_GF_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> m_GF_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> m_GF_mex.mki
echo "CXX=$CXX" >> m_GF_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> m_GF_mex.mki
echo "CXXLIBS=$CXXLIBS" >> m_GF_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> m_GF_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> m_GF_mex.mki
echo "LDFLAGS=$LDFLAGS" >> m_GF_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> m_GF_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> m_GF_mex.mki
echo "Arch=$Arch" >> m_GF_mex.mki
echo "LD=$LD" >> m_GF_mex.mki
echo OMPFLAGS= >> m_GF_mex.mki
echo OMPLINKFLAGS= >> m_GF_mex.mki
echo "EMC_COMPILER=gcc" >> m_GF_mex.mki
echo "EMC_CONFIG=optim" >> m_GF_mex.mki
"/afs/itp.tugraz.at/opt/matlab/R2018b/bin/glnxa64/gmake" -j 1 -B -f m_GF_mex.mk
