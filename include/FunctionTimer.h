#pragma once

/*!
 * Erlaubt es zu messen, wie lange (in Sekunden) die Ausführung einer Funktion/Methode oder 
 * eines Blockes benötigt.

 * Usage for global functions:
 * timeFunctionCall(&functionFoo, elapsedTime).function();
 *
 * For Member-Functions internal:
 * (this->*(timeFunctionCall(&ClassFoo::functionBar, elapsedTime).function))();
 *
 * For Member-Function external:
 * (objectFoo->*(timeFunctionCall(&ClassBar::functionBar, elapsedTime).function))(parameter); bzw.
 * (objectFoo.*(timeFunctionCall(&ClassBar::functionBar, elapsedTime.function))(parameter); wenn objectFoo auf dem Stack liegt.
 *
 * For Overloaded Member-Functions:
 * (objectFoo.*(timeFunctionCall<returnType (ClassBar::*)(parameter)>(&ClassBar::functionBar, elapsedTime).function))(parameter).
 * (timeFunctionCall has to be specialized with the Signature of the Method that should be timed.
 *
 * For BlockTiming:
 * double elapsedTime;
 * {
 *    BlockTimer timer(elapsedTime);
 *    // do what you want to time
 * }
 * std::cout << elapsedTime << "s" << std::endl;
 *
 *
 * Die Klasse misst, wie lange eine Objektinstanz existiert.
 * Da das Objekt nach Beenden der timeFunctionCall Methode zerstört wird,
 * entspricht dies der Ausführungszeit der übergebenen Methode.
 * Wenn eine Membermethode getimed werden soll, muss explizit angegeben werden,
 * zu welchem Objekt der Funktionspointer gehört, auf den zugegriffen werden muss.
 * Dies geschieht über die Ausdrücke this->* bzw. objectFoo->* / objectFoo.* .
 *
 * timeFunctionCall gibt ein FunctionTimer Objekt zurück. Auf diesem kann die zu
 * timende Methode über den Member "function" aufgerufen werden. Somit kann normal
 * mit Parametern oder Rückgabewerten gearbeitet werden.
 *
 * Wenn FunctionTimer mit <void *> spezialisiert wird (oder über den typedef BlockTimer),
 * kann auch Code getimed werden, welcher sich nicht in einer Methode befindet.
 * Dazu muss dieser Code lediglich von einem Block umschlossen werden, in welchem 
 * als erste Anweisung ein FunctionTimer-Objekt angelegt wird. Sobald der Block verlassen
 * wird, befindet sich in elapsedTime die Zeit, die zur Ausführung benötigt wurde.
 //*/

#include <omp.h>

template <typename T>
class FunctionTimer
{
public:
	FunctionTimer(T function, double& time) : function(function), time(time), stop(0)
	{
		start = omp_get_wtime();
	}

	~FunctionTimer(){
		stop = omp_get_wtime();
		time = stop - start;
	}

	T function;

private:
	double &time, start, stop;
};

template<>
class FunctionTimer<void *>
{
public:
	FunctionTimer(double& time) : time(time), stop(0)
	{
		start = omp_get_wtime();
	}

	~FunctionTimer()
	{
		stop = omp_get_wtime();
		time = stop - start;
	}

private:
	double &time, start, stop;

};

typedef FunctionTimer<void *> BlockTimer;

template<typename T>
FunctionTimer<T> timeFunctionCall(T function, double& time)
{
	return FunctionTimer<T>(function, time);
}

#define TIME_THIS(x, elapsedTime) { \
	BlockTimer t(elapsedTime);      \
	do { x; } while(0);             \
}
