# dominoes

Tool for calculation of dynamics of a chain of dominoes  
Theoretical part of the SOWAS project, summer term 2018  

---
**Dependencies:**  
 * GNU Scientific Library  
   Version with CMake support can be found [here](https://github.com/ampl/gsl)  

---
**How to compile:**

    clang++ .\src\main.cpp .\src\DominoChain.cpp -Wall -Wextra -DHAVE_INLINE
    -isystem <path\to\gsl\> -I .\include\ -L <path\to\gsl\libs> -l gsl
    -l gslcblas -o main.exe

**Compile with video support**

    clang++ .\tests\video_test.cpp .\src\DominoChain.cpp -I .\include\
    -I D:\source\gsl\build-dir\ -I D:\source\opencv\build\include\
    -L D:\source\gsl\build-dir\Debug\ -L D:\source\opencv\build\x64\vc15\lib\
    -l gsl -l gslcblas -l opencv_world341 -l opencv_world341d -o
    .\tests\video_test.exe

---
*Instructions to be added*
