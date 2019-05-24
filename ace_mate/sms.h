#ifndef _SMS_H
#define _SMS_H

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(WIN64) || defined(__WIN64__) || defined(_WIN64)
#define SMSWINDOWS
#define _CRT_SECURE_NO_DEPRECATE
#endif

#if defined (__linux) || defined (linux) || defined(__linux__)
#define SMSLINUX
#endif

#if defined (__MACH) || defined (MACH) || defined(__MACH__)
#define SMSMAC
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef SMSWINDOWS
#include "direct.h"
#endif 

#ifdef SMSLINUX
#include <unistd.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <dlfcn.h>
#endif 

#ifdef SMSMAC
#include <dlfcn.h>
#endif

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b)) 
#endif

#ifdef SMSWINDOWS
#define TIMEB _timeb
#define ISNAN _isnan
#define FINITE _finite
#define DLLEXPORT __declspec(dllexport)
#define FTIME _ftime
#define GETCWD _getcwd
#endif

#ifdef SMSLINUX
#define GETCWD getcwd
#define _copysign copysign
#define FTIME ftime
#define TIMEB timeb
#define _MAX_PATH 4096
#define _isnan isnan
#define _finite finite
#define ISNAN isnan
#define FINITE finite
#define CALLBACK __attribute__((__stdcall__))
#define DLLEXPORT
#define _inline
#endif 

#ifdef SMSMAC
#define GETCWD getcwd
#define _copysign copysign
#define FTIME ftime
#define TIMEB timeb
#define _MAX_PATH 4096
#define _isnan isnan
#define _finite finite
#define ISNAN isnan
#define FINITE finite
#define CALLBACK __attribute__((__stdcall__))
#define DLLEXPORT
#define _inline
#endif 


double SMSDot(double a[],double b[],int n);

double SMSSum(double a[],int n);

double SMSKDelta(int i,int j);

double SMSDeltaPart(double *a,int i, int j,int k);

void SMSMove(double a[],double b[],int n);

void SMSZero(double a[],int n);

#define Power(x, y) (pow((double)(x), (double)(y)))
#define Sqrt(x)        (sqrt((double)(x)))
#define Cbrt(x) ( x<0 ? -pow(-(double)(x),1./3.) : pow((double)(x),1./3.) )
#define Abs(x)        (fabs((double)(x)))
#define Less(a,b,c)        (a<b && b<c)
#define LessEqual(a,b,c)        (a<=b && b<=c)
#define Greater(a,b,c)        (a>b && b>c)
#define GreaterEqual(a,b,c)        (a>=b && b>=c)

#endif