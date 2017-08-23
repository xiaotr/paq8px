set PATH=c:\mingw\mingw64\bin;c:\mingw\msys\1.0\bin;%PATH%
g++ misc.cpp main.cpp filter.cpp encoder.cpp filter/*.cpp -DWINDOWS -lz -Wall -Wextra -O3 -static -static-libgcc -opaq8px.exe 
pause