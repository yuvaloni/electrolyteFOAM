the solver was developed and tested on the ESI version

compilation instructions:
extract the contest of electrolytesimpleFoam into path/to/openfoam/applications/solvers
open bash and run the following commands:
```console
cd path/to/openfoam
WM_PROJECT_DIR="$PWD"
source etc/bashrc
./applications/Allwmake
```

to run the test case:
extract the contents of test case folder
run the following commands:
```console
blockmesh -case path/to/testcase
electrolyteFoam -case path/to/testcase
```
