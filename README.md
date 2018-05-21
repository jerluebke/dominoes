# dominoes

Tool for calculation of dynamics of a chain of dominoes  
Theoretical part of the SOWAS project, summer term 2018  

---
**Dependencies:**  
 * GNU Scientific Library  
   Version with CMake support can be found [here](https://github.com/ampl/gsl)  
 * OpenCV  
   make sure to have the OpenCV `.dll`s in your directory  

---
**How to compile:**

    clang++ .\test\main.cpp .\src\DominoChain.cpp -Wall -Wextra -DHAVE_INLINE
    -isystem <path\to\gsl\> -I .\include\ -L <path\to\gsl\libs> -l gsl
    -l gslcblas -o main.exe

**Compile with video support**

    clang++ .\tests\video_test.cpp .\src\DominoChain.cpp -Wall -Wextra
    -DHAVE_INLINE -DVIDEO -isystem <path\to\gsl> -I .\include\
    -I <path\to\opencv> -L <path\to\gsl\libs> -l gsl -l gslcblas
    -L <path\to\opencv\libs> -l opencv_world341 -l opencv_world341d
    -o .\tests\video_test.exe

---
*Instructions to be added*
