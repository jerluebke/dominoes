# dominoes

Tool for calculation of dynamics of a chain of dominoes  
Theoretical part of the SOWAS project, summer term 2018  

---
**Dependencies:**  
 * GNU Scientific Library  
   Version with CMake support can be found [here](https://github.com/ampl/gsl)  

---
**How to compile:**

    clang++ .\src\main.cpp .\src\DominoChain.cpp -Wall -Wextra -fexceptions
    -DHAVE_INLINE -isystem <path\to\gsl\> -I .\include\ -L <path\to\gsl\libs>
    -l gsl -gslcblas -o main.exe

---
*Instructions to be added*
