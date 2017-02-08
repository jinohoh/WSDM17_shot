#ifndef _UTIL_
#define _UTIL_

#include <cstdlib> // srand() and drand()
#include <random>  //mt19937 stuffs
#include <time.h> // used for seed srand();

#ifndef _WIN32
#ifndef _TICK_
#define _TICK_

#ifdef __APPLE__
#define CLOCK_MONOTONIC 0
#include <sys/time.h>
//clock_gettime is not implemented on OSX
inline int clock_gettime(int /*clk_id*/, struct timespec* t) {
	struct timeval now;
	int rv = gettimeofday(&now, NULL);
	if (rv) 
		return rv;
	t->tv_sec  = now.tv_sec;
	t->tv_nsec = now.tv_usec * 1000;
	return 0;
}
#else
#include <time.h> // Unix based system specific
#endif // #ifdef __APPLE__ 

inline double GetTickCount(void)
{
	struct timespec now;
	if (clock_gettime(CLOCK_MONOTONIC, &now))
		return 0;
	return now.tv_sec * 1000.0 + now.tv_nsec / 1000000.0;
};

#endif // #ifndef _TICK_
#endif // #ifndef _WIN32


std::mt19937 generator;

void inline _srand(time_t x){
#ifdef _WIN32
	srand(x);
#else   
	generator.seed(x);
	//srand48(x);
#endif  
};       
        
std::uniform_real_distribution<double> dis(0.0, 1.0);

double inline _drand(){
#ifdef _WIN32   
return double(rand())/RAND_MAX;
#else
return dis(generator);
//return drand48();
#endif
};

#endif
