CC="/usr/bin/gcc"
                CXX="g++"
                CFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -D_GNU_SOURCE -DMATLAB_MEX_FILE "
                CXXFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -D_GNU_SOURCE -DMATLAB_MEX_FILE "
                COPTIMFLAGS="-O2 -fwrapv -DNDEBUG"
                CXXOPTIMFLAGS="-O2 -fwrapv -DNDEBUG"
                CDEBUGFLAGS="-g"
                CXXDEBUGFLAGS="-g"
                LD="/usr/bin/gcc"
                LDXX="g++"
                LDFLAGS="-pthread -Wl,--no-undefined -Wl,-rpath-link,/afs/itp.tugraz.at/opt/matlab/R2018b/bin/glnxa64 -shared  -L"/afs/itp.tugraz.at/opt/matlab/R2018b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -Wl,--version-script,"/afs/itp.tugraz.at/opt/matlab/R2018b/extern/lib/glnxa64/mexFunction.map""
                LDDEBUGFLAGS="-g"
