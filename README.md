# microBUDE

An even more reduced version of miniBUDE with just a single source file.
Currently, only the OpenMP model is implemented.

Both Makefile and CMakeList is included:

```shell
> make
g++ -std=c++17 -Wall -Wno-sign-compare -Ofast -march=native -fopenmp -g3 -c main.cpp -o main.o
g++ -std=c++17 -Wall -Wno-sign-compare -Ofast -march=native -fopenmp -g3 main.o -o microbude
```

Set `CXX` and `CXXFLAGS` to customise the compiler and flags.
Define `PPWI`, `ITERS`, and `DECK_PATH` using -D flags (i.e `-DPPWI=4`) to adjust problem sizes.

The kernel `fasten_main` is annotated with `extern "C"` so that it's not completely mangled.