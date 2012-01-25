/* file: math.h */

/*
** Copyright  (C) 1985-2007 Intel Corporation. All rights reserved.
**
** The information and source code contained herein is the exclusive property
** of Intel Corporation and may not be disclosed, examined, or reproduced in
** whole or in part without explicit written authorization from the Company.
**
*/


/* Include_next should be before guard macros in order to at last reach system header */
#if defined(__PURE_SYS_C99_HEADERS__)
# include_next <math.h> /* utilize sys header */
#else

#if !defined(__PURE_INTEL_C99_HEADERS__)
# include_next <math.h> /* utilize and expand sys header */
#endif

#ifndef __MATH_H_INCLUDED
#define __MATH_H_INCLUDED

/* Check usage correctness */
#if !defined(__INTEL_COMPILER)
# error "This Intel <math.h> is for use with only the Intel compilers!"
#endif

#if defined(__linux__) || defined (__APPLE__) || defined(__FreeBSD__) || defined(__QNX__)
# if defined(fpclassify)
#  undef fpclassify
# endif
# if defined(finite)
#  undef finite
# endif
# if defined(isnan)
#  undef isnan
# endif
# if defined(isinf)
#  undef isinf
# endif
# if defined(isnormal)
#  undef isnormal
# endif
# if defined(isfinite)
#  undef isfinite
# endif
# if defined(signbit)
#  undef signbit
# endif
#endif /* __linux__ */

#if defined(__linux__) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__QNX__)
# if defined(isgreater)
#  undef isgreater
# endif
# if defined(isless)
#  undef isless
# endif
# if defined(isgreaterequal)
#  undef isgreaterequal
# endif
# if defined(islessequal)
#  undef islessequal
# endif
# if defined(islessgreater)
#  undef islessgreater
# endif
# if defined(isunordered)
#  undef isunordered
# endif
#endif /* __linux__ || __APPLE__ || __FreeBSD__ || __QNX__ */

#if (!defined (__linux__) && !defined(__APPLE__) && !defined(__FreeBSD__) && !defined(__QNX__)) || defined (__PURE_INTEL_C99_HEADERS__)
#if !(defined (DECIMAL_DIG))
#define DECIMAL_DIG 21
#endif
#if defined (FLT_EVAL_METHOD)
#if (FLT_EVAL_METHOD == 0)
typedef float  float_t;
typedef double double_t;
#elif (FLT_EVAL_METHOD == 1)
typedef double float_t;
typedef double double_t;
#elif (FLT_EVAL_METHOD == 2)
typedef long double float_t;
typedef long double double_t;
#else
typedef float   float_t;
typedef double  double_t;
#endif
#else
#define FLT_EVAL_METHOD -1
typedef float   float_t;
typedef double  double_t;
#endif /* FLT_EVAL_METHOD */
#endif /* (!__linux__ && !__APPLE__ && !__FreeBSD__ && !__QNX__) || __PURE_INTEL_C99_HEADERS__ */


#if defined(__APPLE__)
#undef NAN
#endif
#if (!defined(__linux__) && !defined(__FreeBSD__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
static unsigned int __qnan[] = {0x7fc00001};
#define NAN (*((float *)__qnan))
#endif


#if defined(__PURE_INTEL_C99_HEADERS__) || !(defined(__linux__) || defined(__unix__) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__QNX__)) /* We need to define FP_ILOGB0, FP_ILOGBNAN */
#if (defined(__linux__) || defined(__unix__) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__QNX__)) && !(defined(__ECL) || defined(__ECC))   /* Linux/unix or MacOS */           
#define FP_ILOGB0	(-2147483647 - 1)
#define FP_ILOGBNAN (-2147483647 - 1)
#else /* Windows and Intel(R) Itanium(R) architecture */
#define FP_ILOGB0	(-2147483647 - 1)
#define FP_ILOGBNAN	2147483647
#endif
#endif

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__FreeBSD__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
static unsigned int __huge_valf[] = {0x7f800000};
#define HUGE_VALF (*((double *)__huge_valf))
static unsigned int __huge_vall[] = {0, 0, 0x00007f80, 0};
#define HUGE_VALL (*((long double *)__huge_vall))

static unsigned int __huge_val[] = {0, 0x7ff00000};
#define HUGE_VAL (*((double *)__huge_val))
static unsigned int __infinity[] = {0x7f800000};
#define INFINITY (*((float *)__infinity))
#endif /* (!__linux__ && !__APPLE__) || __PURE_INTEL_C99_HEADERS__ */

#if defined(__LONG_DOUBLE_SIZE__)	/* Compiler-predefined macros. If defined, should be 128|80|64 */
# define __IMFLONGDOUBLE (__LONG_DOUBLE_SIZE__)
#else
# define __IMFLONGDOUBLE 64
#endif

#if defined(__cplusplus)    /* C linkage */
extern "C" {
#endif

#define _LIBIMF_EXT extern

#if defined(_DLL) && (defined(_WIN32) || defined(_WIN64)) && !defined(_LIBMMDS)   /* Windows DLL */
# define _LIBIMF_PUBAPI __declspec(dllimport) __cdecl
# define _LIBIMF_PUBVAR __declspec(dllimport)
#elif defined(__unix__) || defined(__APPLE__) || defined(__QNX__) /* Linux, MacOS or QNX */
# define _LIBIMF_PUBAPI /* do not change this line! */
# define _LIBIMF_PUBVAR
#else                                             /* Windows static */
# define _LIBIMF_PUBAPI __cdecl
# define _LIBIMF_PUBVAR
#endif

/* for Microsoft compatibility */
#if !defined (__linux__) && !defined(__APPLE__) && !defined(__FreeBSD__) && !defined(__QNX__) && !defined(__INTEL_COMPILER)
#define _LIBIMF_LONGLONG __int64
#else
#define _LIBIMF_LONGLONG long long int
#endif
#define _LIBIMF_LONGINT  long int
#define _LIBIMF_INT      int
#define _LIBIMF_DOUBLE   double
#define _LIBIMF_FLOAT    float
#define _LIBIMF_XDOUBLE  long double
#define _LIBIMF_VOID     void
#define _LIBIMF_CCHAR    const char


/*-- Classification macros --*/

#if defined (__IWMMXT__) || defined(__PURE_INTEL_C99_HEADERS__) || defined(_WIN32) || defined(_WIN64)
#if defined(__APPLE__)
#define FP_NAN          (1)
#define FP_INFINITE     (2)
#define FP_ZERO         (3)
#define FP_NORMAL       (4)
#define FP_SUBNORMAL    (5)
#elif defined (__FreeBSD__)
# define FP_NAN         (2)
# define FP_INFINITE    (1)
# define FP_ZERO        (16)
# define FP_SUBNORMAL   (8)
# define FP_NORMAL      (4)
#else
# define FP_NAN       (0)
# define FP_INFINITE  (1)
# define FP_ZERO      (2)
# define FP_SUBNORMAL (3)
# define FP_NORMAL    (4)
#endif
#endif

#if !defined FP_NAN
# define FP_NAN       (0)
#endif
#if !defined FP_INFINITE
# define FP_INFINITE  (1)
#endif
#if !defined FP_ZERO
# define FP_ZERO      (2)
#endif
#if !defined FP_SUBNORMAL
# define FP_SUBNORMAL (3)
#endif
#if !defined FP_NORMAL
# define FP_NORMAL    (4)
#endif

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_INT   _LIBIMF_PUBAPI abs( _LIBIMF_INT __x );
#endif

#if defined (__APPLE__)
#define _LIBIMF_DBL_XDBL    _LIBIMF_XDOUBLE
#else
#define _LIBIMF_DBL_XDBL    _LIBIMF_DOUBLE
#endif

_LIBIMF_EXT _LIBIMF_INT fpclassifyf    ( _LIBIMF_FLOAT              __x );
_LIBIMF_EXT _LIBIMF_INT fpclassify     ( _LIBIMF_DBL_XDBL               __x );
_LIBIMF_EXT _LIBIMF_INT fpclassifyd    ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT fpclassifyl    ( _LIBIMF_XDOUBLE            __x );

_LIBIMF_EXT _LIBIMF_INT __fpclassifyf  ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __fpclassify   ( _LIBIMF_DBL_XDBL               __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __fpclassifyd  ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __fpclassifyl  ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT isinff         ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT isinf          ( _LIBIMF_DBL_XDBL               __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT isinfd         ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT isinfl         ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT __isinff       ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __isinf        ( _LIBIMF_DBL_XDBL               __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __isinfd       ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __isinfl       ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT isnanf         ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT isnan          ( _LIBIMF_DBL_XDBL               __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT isnand         ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT isnanl         ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT __isnanf       ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __isnan        ( _LIBIMF_DBL_XDBL               __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __isnand       ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __isnanl       ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT isnormalf      ( _LIBIMF_FLOAT              __x );
_LIBIMF_EXT _LIBIMF_INT isnormal       ( _LIBIMF_DBL_XDBL               __x );
_LIBIMF_EXT _LIBIMF_INT isnormald      ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT isnormall      ( _LIBIMF_XDOUBLE            __x );

_LIBIMF_EXT _LIBIMF_INT __isnormalf    ( _LIBIMF_FLOAT              __x );
_LIBIMF_EXT _LIBIMF_INT __isnormal     ( _LIBIMF_DBL_XDBL               __x );
_LIBIMF_EXT _LIBIMF_INT __isnormald    ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __isnormall    ( _LIBIMF_XDOUBLE            __x );

_LIBIMF_EXT _LIBIMF_INT isfinitef      ( _LIBIMF_FLOAT              __x );
_LIBIMF_EXT _LIBIMF_INT isfinite       ( _LIBIMF_DBL_XDBL               __x );
_LIBIMF_EXT _LIBIMF_INT isfinited      ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT isfinitel      ( _LIBIMF_XDOUBLE            __x );

_LIBIMF_EXT _LIBIMF_INT __isfinitef    ( _LIBIMF_FLOAT              __x );
_LIBIMF_EXT _LIBIMF_INT __isfinite     ( _LIBIMF_DBL_XDBL               __x );
_LIBIMF_EXT _LIBIMF_INT __isfinited    ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __isfinitel    ( _LIBIMF_XDOUBLE            __x );

_LIBIMF_EXT _LIBIMF_INT finitef        ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT finite         ( _LIBIMF_DOUBLE             __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT finited        ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT finitel        ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT __finitef      ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __finite       ( _LIBIMF_DOUBLE             __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __finited      ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __finitel      ( _LIBIMF_XDOUBLE            __x ) __THROW;

_LIBIMF_EXT _LIBIMF_INT signbitf       ( _LIBIMF_FLOAT              __x );
_LIBIMF_EXT _LIBIMF_INT signbit        ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT signbitd       ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT signbitl       ( _LIBIMF_XDOUBLE            __x );

_LIBIMF_EXT _LIBIMF_INT __signbitf     ( _LIBIMF_FLOAT              __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __signbit      ( _LIBIMF_DOUBLE             __x ) __THROW;
_LIBIMF_EXT _LIBIMF_INT __signbitd     ( _LIBIMF_DOUBLE             __x );
_LIBIMF_EXT _LIBIMF_INT __signbitl     ( _LIBIMF_XDOUBLE            __x ) __THROW;


#define __IMFC99MACRO1ARG_ALL( __x__, __func__, __fprefix__, __fsuffix__, \
                                                __dprefix__, __dsuffix__, \
                                                __lprefix__, __lsuffix__) \
    (( sizeof( __x__ ) > sizeof( double ))                                \
     ? __lprefix__##__func__##__lsuffix__( (long double)(__x__) )         \
     : (( sizeof( __x__ ) == sizeof( float ))                             \
        ? __fprefix__##__func__##__fsuffix__(    (float)(__x__) )         \
        : __dprefix__##__func__##__dsuffix__(   (double)(__x__) )         \
       )                                                                  \
    )


#if defined (__FreeBSD__)

#define fpclassify( __x__ ) __IMFC99MACRO1ARG_ALL( __x__, fpclassify, __, f, __, d, __, l)
#define isinf( __x__ )      __IMFC99MACRO1ARG_ALL( __x__, isinf,      __, f,   ,  , __, l)
#define isnan( __x__ )      __IMFC99MACRO1ARG_ALL( __x__, isnan,        , f,   ,  , __, l)
#define isnormal( __x__ )   __IMFC99MACRO1ARG_ALL( __x__, isnormal,   __, f, __,  , __, l)
#define isfinite( __x__ )   __IMFC99MACRO1ARG_ALL( __x__, isfinite,   __, f, __,  , __, l)
#define signbit( __x__ )    __IMFC99MACRO1ARG_ALL( __x__, signbit,    __, f, __,  , __, l)

#elif defined (__APPLE__)

#define fpclassify( __x__ ) __IMFC99MACRO1ARG_ALL( __x__, fpclassify, __, f, __, d, __,  )
#define isinf( __x__ )      __IMFC99MACRO1ARG_ALL( __x__, isinf,      __, f, __, d, __,  )
#define isnan( __x__ )      __IMFC99MACRO1ARG_ALL( __x__, isnan,      __, f, __, d, __,  )
#define isnormal( __x__ )   __IMFC99MACRO1ARG_ALL( __x__, isnormal,   __, f, __, d, __,  )
#define isfinite( __x__ )   __IMFC99MACRO1ARG_ALL( __x__, finite,     __, f, __,  , __, l)
#define signbit( __x__ )    __IMFC99MACRO1ARG_ALL( __x__, signbit,    __, f, __, d, __, l)

#else

#define fpclassify( __x__ ) __IMFC99MACRO1ARG_ALL( __x__, fpclassify, __, f, __,  , __, l)
#define isinf( __x__ )      __IMFC99MACRO1ARG_ALL( __x__, isinf,      __, f, __,  , __, l)
#define isnan( __x__ )      __IMFC99MACRO1ARG_ALL( __x__, isnan,      __, f, __,  , __, l)
#define isnormal( __x__ )   __IMFC99MACRO1ARG_ALL( __x__, isnormal,   __, f, __,  , __, l)
#define isfinite( __x__ )   __IMFC99MACRO1ARG_ALL( __x__, finite,     __, f, __,  , __, l)
#define signbit( __x__ )    __IMFC99MACRO1ARG_ALL( __x__, signbit,    __, f, __,  , __, l)

#endif


/* Comparison macros */

_LIBIMF_EXT _LIBIMF_INT isgreaterf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT isgreater( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT isgreaterl( long double xl, long double yl );
_LIBIMF_EXT _LIBIMF_INT __isgreaterf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT __isgreater( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT __isgreaterl( long double xl, long double yl );

_LIBIMF_EXT _LIBIMF_INT isgreaterequalf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT isgreaterequal( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT isgreaterequall( long double xl, long double yl );
_LIBIMF_EXT _LIBIMF_INT __isgreaterequalf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT __isgreaterequal( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT __isgreaterequall( long double xl, long double yl );

_LIBIMF_EXT _LIBIMF_INT islessf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT isless( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT islessl( long double xl, long double yl );
_LIBIMF_EXT _LIBIMF_INT __islessf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT __isless( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT __islessl( long double xl, long double yl );

_LIBIMF_EXT _LIBIMF_INT islessequalf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT islessequal( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT islessequall( long double xl, long double yl );
_LIBIMF_EXT _LIBIMF_INT __islessequalf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT __islessequal( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT __islessequall( long double xl, long double yl );

_LIBIMF_EXT _LIBIMF_INT islessgreaterf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT islessgreater( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT islessgreaterl( long double xl, long double yl );
_LIBIMF_EXT _LIBIMF_INT __islessgreaterf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT __islessgreater( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT __islessgreaterl( long double xl, long double yl );

_LIBIMF_EXT _LIBIMF_INT isunorderedf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT isunordered( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT isunorderedl( long double xl, long double yl );
_LIBIMF_EXT _LIBIMF_INT __isunorderedf( float xf, float yf );
_LIBIMF_EXT _LIBIMF_INT __isunordered( double xd, double yd );
_LIBIMF_EXT _LIBIMF_INT __isunorderedl( long double xl, long double yl );

#define __IMFC99MACRO2ARG( __x__, __y__, __func__ ) \
	((( sizeof( __x__ ) > sizeof( double )) || ( sizeof( __y__ ) > sizeof( double ))) \
	 ? __func__##l( (long double)(__x__), (long double)(__y__) ) \
	 : ((( sizeof( __x__ ) + sizeof( __y__ )) == (2*sizeof( float ))) \
		? __func__##f( (float)(__x__), (float)(__y__) ) \
		: __func__( (double)(__x__), (double)(__y__) )))

#define isgreater( __x__, __y__ )       __IMFC99MACRO2ARG( __x__, __y__, __isgreater )
#define isgreaterequal( __x__, __y__ )  __IMFC99MACRO2ARG( __x__, __y__, __isgreaterequal )
#define isless( __x__, __y__ )          __IMFC99MACRO2ARG( __x__, __y__, __isless )
#define islessequal( __x__, __y__ )     __IMFC99MACRO2ARG( __x__, __y__, __islessequal )
#define islessgreater( __x__, __y__ )   __IMFC99MACRO2ARG( __x__, __y__, __islessgreater )
#define isunordered( __x__, __y__ )     __IMFC99MACRO2ARG( __x__, __y__, __isunordered )

/*-- Real functions --*/

/* Radian argument trigonometric functions */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI acos( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI asin( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI atan( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI atan2( _LIBIMF_DOUBLE __y, _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI cos( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI sin( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI tan( _LIBIMF_DOUBLE __x );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI acosf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI asinf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI atanf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI atan2f( _LIBIMF_FLOAT __y, _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI cosf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI sinf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI tanf( _LIBIMF_FLOAT __x );
#endif

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI cot( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI cotf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI cotl( _LIBIMF_XDOUBLE __x );

#if (defined(__APPLE__) || (!defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99)))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sincos( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE *__psin, _LIBIMF_DOUBLE *__pcos );
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sincosf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT *__psin, _LIBIMF_FLOAT *__pcos );
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sincosl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE *__psin, _LIBIMF_XDOUBLE *__pcos );
#endif /* __USE_ISOC99 */

/* Degree argument trigonometric functions */

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI cosd( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI cosdf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI cosdl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI sind( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI sindf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sindl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI tand( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI tandf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI tandl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI cotd( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI cotdf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI cotdl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI acosd( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI acosdf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI acosdl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI asind( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI asindf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI asindl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI atand( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI atandf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atandl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI atand2( _LIBIMF_DOUBLE __y, _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI atand2f( _LIBIMF_FLOAT __y, _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atand2l( _LIBIMF_XDOUBLE __y, _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI atan2d( _LIBIMF_DOUBLE __y, _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI atan2df( _LIBIMF_FLOAT __y, _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atan2dl( _LIBIMF_XDOUBLE __y, _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sincosd( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE *__psin, _LIBIMF_DOUBLE *__pcos );
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sincosdf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT *__psin, _LIBIMF_FLOAT *__pcos );
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sincosdl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE *__psin, _LIBIMF_XDOUBLE *__pcos );

/* Hyperbolic functions */

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI acosh( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI acoshf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI asinh( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI asinhf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI atanh( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI atanhf( _LIBIMF_FLOAT __x );
#endif  /* (!__linux__ && !__APPLE__ && !__QNX__) || __PURE_INTEL_C99_HEADERS__ */

_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sinhcosh( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE *__psinh, _LIBIMF_DOUBLE *__pcosh );
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sinhcoshf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT *__psinh, _LIBIMF_FLOAT *__pcosh );
_LIBIMF_EXT _LIBIMF_VOID     _LIBIMF_PUBAPI sinhcoshl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE *__psinh, _LIBIMF_XDOUBLE *__pcosh );

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI cosh( _LIBIMF_DOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI sinh( _LIBIMF_DOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI tanh( _LIBIMF_DOUBLE __x );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI coshf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI sinhf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI tanhf( _LIBIMF_FLOAT __x );
#endif

/* Exponential functions */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI exp( _LIBIMF_DOUBLE __x );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI expf( _LIBIMF_FLOAT __x );
#endif

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI expm1( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI expm1f( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI exp2( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI exp2f( _LIBIMF_FLOAT __x );

#endif /* __USE_ISOC99 */

#if (!defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI exp10( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI exp10f( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI exp10l( _LIBIMF_XDOUBLE __x );
#endif /* __USE_ISOC99 */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI frexp( _LIBIMF_DOUBLE __x, _LIBIMF_INT *__exp );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI ldexp( _LIBIMF_DOUBLE __x, _LIBIMF_INT __exp );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI frexpf( _LIBIMF_FLOAT __x, _LIBIMF_INT *__exp );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI ldexpf( _LIBIMF_FLOAT __x, _LIBIMF_INT __exp );
#endif

#if (!defined(__APPLE__) && !defined(__QNX__) && (!defined(__linux__) || !defined (__USE_ISOC99) || (defined (__USE_ISOC99) && !defined (__USE_MISC) && !defined(__USE_XOPEN_EXTENDED)))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI scalb( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
#endif

#if !defined(__QNX__) && (!defined(__linux__) || !defined (__USE_ISOC99) || (defined (__USE_ISOC99) && !defined (__USE_MISC) && !defined (__USE_XOPEN_EXTENDED)) || defined(__PURE_INTEL_C99_HEADERS__))
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI scalbf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI scalbl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
#endif

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI scalbn( _LIBIMF_DOUBLE __x, _LIBIMF_INT __n );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI scalbnf( _LIBIMF_FLOAT __x, _LIBIMF_INT __n );

#endif  /* !__linux__ && !__APPLE__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI scalbln( _LIBIMF_DOUBLE __x, _LIBIMF_LONGINT __n );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI scalblnf( _LIBIMF_FLOAT __x, _LIBIMF_LONGINT __n );
#endif /* __USE_ISOC99 */

/* Logarithmic functions */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI log( _LIBIMF_DOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI log10( _LIBIMF_DOUBLE __x );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI logf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI log10f( _LIBIMF_FLOAT __x );
#endif

#if (!defined (__APPLE__) && !defined (__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI log2( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI log2f( _LIBIMF_FLOAT __x );
#endif /* __USE_ISOC99 */

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI log1p( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI log1pf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI logb( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI logbf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_INT      _LIBIMF_PUBAPI ilogb( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_INT      _LIBIMF_PUBAPI ilogbf( _LIBIMF_FLOAT __x );

#endif  /* !__linux__ && !__APPLE__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI modf( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE *__iptr );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI modff( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT *__iptr );
#endif

/* Power/root/abs functions */

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI cbrt( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI cbrtf( _LIBIMF_FLOAT __x );

#endif  /* !__linux__ && !__APPLE__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI fabs( _LIBIMF_DOUBLE __x );


_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI pow( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI sqrt( _LIBIMF_DOUBLE __x );

#endif  /*__PURE_INTEL_C99_HEADERS__*/


#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI fabsf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI hypot( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI hypotf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI powf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI sqrtf( _LIBIMF_FLOAT __x );
#endif

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI invsqrt( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI invsqrtf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI invsqrtl( _LIBIMF_XDOUBLE __x );

/* Error and gamma functions */

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI erf( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI erff( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI erfc( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI erfcf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI lgamma( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI lgammaf( _LIBIMF_FLOAT __x );

#endif  /* !__linux__ && !__APPLE__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI erfinv( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI erfinvf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI erfinvl( _LIBIMF_XDOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI gamma_r( _LIBIMF_DOUBLE __x, _LIBIMF_INT *__signgam );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI gammaf_r( _LIBIMF_FLOAT __x, _LIBIMF_INT *__signgam );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI gammal_r( _LIBIMF_XDOUBLE __x, _LIBIMF_INT *__signgam );

#if (!defined(__linux__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI gamma( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI gammal( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI gammaf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI lgammaf_r( _LIBIMF_FLOAT __x, _LIBIMF_INT *__signgam );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI lgamma_r( _LIBIMF_DOUBLE __x, _LIBIMF_INT *__signgam );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI lgammal_r( _LIBIMF_XDOUBLE __x, _LIBIMF_INT *__signgam );
#endif  /* !__linux__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */


#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI tgamma( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI tgammaf( _LIBIMF_FLOAT __x );
#endif /* __USE_ISOC99 */

/* Nearest integer functions */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI ceil( _LIBIMF_DOUBLE __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI floor( _LIBIMF_DOUBLE __x );

#endif  /*__PURE_INTEL_C99_HEADERS__*/


#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI ceilf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI floorf( _LIBIMF_FLOAT __x );
#endif

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI nearbyint( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI nearbyintf( _LIBIMF_FLOAT __x );
#endif /* __USE_ISOC99 */

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI rint( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI rintf( _LIBIMF_FLOAT __x );

#endif  /* !__linux__ && !__APPLE__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_LONGINT  _LIBIMF_PUBAPI lrint( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGINT  _LIBIMF_PUBAPI lrintf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_LONGLONG _LIBIMF_PUBAPI llrint( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGLONG _LIBIMF_PUBAPI llrintf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI round( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI roundf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_LONGINT  _LIBIMF_PUBAPI lround( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGINT  _LIBIMF_PUBAPI lroundf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_LONGLONG _LIBIMF_PUBAPI llround( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGLONG _LIBIMF_PUBAPI llroundf( _LIBIMF_FLOAT __x );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI trunc( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI truncf( _LIBIMF_FLOAT __x );
#endif /* __USE_ISOC99 */

/* Remainder functions */

#if defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI fmod( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );

#endif  /*__PURE_INTEL_C99_HEADERS__*/

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && (!defined(_M_IA64) || (_MSC_VER < 1300))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI fmodf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );
#endif

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI remainder( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI remainderf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );

#endif  /* !__linux__ || __PURE_INTEL_C99_HEADERS__ */

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI remquo( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y, _LIBIMF_INT *__quo );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI remquof( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y, _LIBIMF_INT *__quo );
#endif /* __USE_ISOC99 */

/* Manipulation functions */

#if (!defined(__linux__) && !defined(__QNX__) && !defined(__APPLE__)) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI significand( _LIBIMF_DOUBLE __x );
#endif
#if !defined(__linux__) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI significandf( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI significandl( _LIBIMF_XDOUBLE __x );
#endif

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI copysign( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI copysignf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI nextafter( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI nextafterf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );

#endif  /* !__linux__ || __PURE_INTEL_C99_HEADERS__ */

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI nexttoward( _LIBIMF_DOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI nexttowardf( _LIBIMF_FLOAT __x, _LIBIMF_XDOUBLE __y );
#endif /* __USE_ISOC99 */

/* Maximum, minimum, and positive difference functions */

#if (!defined (__APPLE__) && !defined(__QNX__) && !(defined(__linux__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI fdim( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI fdimf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI fmax( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI fmaxf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI fmin( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI fminf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );

/* Floating multiply-add */

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI fma( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y, _LIBIMF_DOUBLE __z );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI fmaf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y, _LIBIMF_FLOAT __z );
#endif /* __USE_ISOC99 */

/* Other - special functions */

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI annuity( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI compound( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );

_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI annuityf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI annuityl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );

_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI compoundf( _LIBIMF_FLOAT __x, _LIBIMF_FLOAT __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI compoundl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );

#if (!defined(__linux__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI j0f( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI j1f( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI jnf( _LIBIMF_INT __n, _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI y0f( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI y1f( _LIBIMF_FLOAT __x );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI ynf( _LIBIMF_INT __n, _LIBIMF_FLOAT __x );
#endif  /* !__linux__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

#if (!(defined(__APPLE__) && !defined(_ANSI_SOURCE)) && !defined(__linux__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)

_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI j0( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI j1( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI jn( _LIBIMF_INT __n, _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI y0( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI y1( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI yn( _LIBIMF_INT __n, _LIBIMF_DOUBLE __x );

#endif  /* !__linux__ && !_ANSI_SOURCE && !__APPLE__ && !__QNX__ || __PURE_INTEL_C99_HEADERS__ */

#if !defined (__PURE_INTEL_C99_HEADERS__) && defined(__APPLE__)
#if defined(_ANSI_SOURCE) || defined(_POSIX_C_SOURCE)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI gamma( _LIBIMF_DOUBLE __x );
#if defined(_XOPEN_SOURCE)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI significand( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI gamma_r( _LIBIMF_DOUBLE __x, _LIBIMF_INT *__signgam );
_LIBIMF_EXT _LIBIMF_FLOAT    _LIBIMF_PUBAPI lgammaf_r( _LIBIMF_FLOAT __x, _LIBIMF_INT *__signgam );
#endif /* _XOPEN_SOURCE */
#endif /* _ANSI_SOURCE || _POSIX_C_SOURCE */
#endif /* __PURE_INTEL_C99_HEADERS__ */

#if (!defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__)) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI acosl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI asinl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atanl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atan2l( _LIBIMF_XDOUBLE __y, _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI cosl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sinl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI tanl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI acoshl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI asinhl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atanhl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI coshl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sinhl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI tanhl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI expl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI frexpl( _LIBIMF_XDOUBLE __x, _LIBIMF_INT *__exp );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI ldexpl( _LIBIMF_XDOUBLE __x, _LIBIMF_INT __exp );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI scalbnl( _LIBIMF_XDOUBLE __x, _LIBIMF_INT __n );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI logl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI log10l( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI log1pl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI logbl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_INT      _LIBIMF_PUBAPI ilogbl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI modfl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE *__iptr );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI cbrtl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fabsl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI hypotl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI powl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sqrtl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI erfl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI erfcl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI lgammal( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI ceill( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI floorl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI rintl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fmodl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI remainderl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI copysignl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI nextafterl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
#endif

#if (!defined(__QNX__) && !(defined(__unix__) && defined(__USE_ISOC99))) || defined(__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI exp2l( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI expm1l( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI scalblnl( _LIBIMF_XDOUBLE __x, _LIBIMF_LONGINT __n );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI log2l( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI remquol( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y, _LIBIMF_INT *__quo );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI roundl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGINT  _LIBIMF_PUBAPI lroundl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGLONG _LIBIMF_PUBAPI llroundl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGINT  _LIBIMF_PUBAPI lrintl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_LONGLONG _LIBIMF_PUBAPI llrintl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI truncl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI nearbyintl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI tgammal( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI nexttowardl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fdiml( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fmaxl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fminl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fmal( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y, _LIBIMF_XDOUBLE __z );
#endif

#if (__IMFLONGDOUBLE == 64)	/* MS compatibility */
# define acosdl       acosd
# define acoshl       acosh
# define acosl        acos
# define annuityl     annuity
# define asindl       asind
# define asinhl       asinh
# define asinl        asin
# define atan2dl      atan2d
# define atan2l       atan2
# define atand2l      atand2
# define atandl       atand
# define atanhl       atanh
# define atanl        atan
# define cbrtl        cbrt
# define ceill        ceil
# define compoundl    compound
# define copysignl    copysign
# define cosdl        cosd
# define coshl        cosh
# define cosl         cos
# define cotdl        cotd
# define cotl         cot
# define erfcl        erfc
# define erfl         erf
# define exp10l       exp10
# define exp2l        exp2
# define expl         exp
# define expm1l       expm1
# define fabsl        fabs
# define fdiml        fdim
# define floorl       floor
# define fmal         fma
# define fmaxl        fmax
# define fminl        fmin
# define fmodl        fmod
# define frexpl       frexp
# define gammal       gamma
# define gammal_r     gamma_r
# define hypotl       hypot
# define ilogbl       ilogb
# define invsqrtl     invsqrt
# define  fpclassifyl  __fpclassify
# define  isfinitel    __finite
# define  isinfl       __isinf
# define  isnanl       __isnan
# define  isnormall    __isnormal
# define  signbitl     __signbit

# define __fpclassifyl __fpclassify
# define __finitel     __finite
# define __isinfl      __isinf
# define __isnanl      __isnan
# define __isnormall   __isnormal
# define __signbitl    __signbit

# define isgreaterl         __isgreater
# define isgreaterequall    __isgreaterequal
# define islessl            __isless
# define islessequall       __islessequal
# define islessgreaterl     __islessgreater
# define isunorderedl       __isunordered
# define __isgreaterl       __isgreater
# define __isgreaterequall  __isgreaterequal
# define __islessl          __isless
# define __islessequall     __islessequal
# define __islessgreaterl   __islessgreater
# define __isunorderedl     __isunordered
# define ldexpl       ldexp
# define lgammal      lgamma
# define lgammal_r    lgamma_r
# define llrintl      llrint
# define llroundl     llround
# define log10l       log10
# define log1pl       log1p
# define log2l        log2
# define logbl        logb
# define logl         log
# define lrintl       lrint
# define lroundl      lround
# define modfl        modf
# define nearbyintl   nearbyint
# define nextafterl   nextafter
# define nexttowardl  nexttoward
# define powl         pow
# define remainderl   remainder
# define remquol      remquo
# define rintl        rint
# define roundl       round
# define scalbl       scalb
# define scalblnl     scalbln
# define scalbnl      scalbn
# define significandl significand
# define sincosdl     sincosd
# define sincosl      sincos
# define sindl        sind
# define sinhcoshl    sinhcosh
# define sinhl        sinh
# define sinl         sin
# define sqrtl        sqrt
# define tandl        tand
# define tanhl        tanh
# define tanl         tan
# define tgammal      tgamma
# define truncl       trunc
#elif (defined (_WIN32) || defined(_WIN64)) && !defined(__PURE_INTEL_C99_HEADERS__)
# undef sinl
# undef cosl
# undef tanl
# undef asinl
# undef acosl
# undef atanl
# undef atan2l
# undef hypotl
# undef sqrtl
# undef expl
# undef logl
# undef log10l
# undef sinhl
# undef coshl
# undef tanhl
# undef ceill
# undef floorl
# undef fabsl
# undef fmodl
# undef modfl
# undef frexpl
# undef ldexpl
# undef powl

_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sinl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI cosl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI tanl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI asinl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI acosl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atanl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI atan2l( _LIBIMF_XDOUBLE __y, _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI hypotl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sinhl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI coshl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI tanhl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI expl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI logl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI log10l( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI frexpl( _LIBIMF_XDOUBLE __x, _LIBIMF_INT *__exp );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI ldexpl( _LIBIMF_XDOUBLE __x, _LIBIMF_INT __exp );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI modfl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE *__iptr );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fabsl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI powl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI sqrtl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI ceill( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI floorl( _LIBIMF_XDOUBLE __x );
_LIBIMF_EXT _LIBIMF_XDOUBLE  _LIBIMF_PUBAPI fmodl( _LIBIMF_XDOUBLE __x, _LIBIMF_XDOUBLE __y );
#endif /*(__IMFLONGDOUBLE == 64)*/


/* MS compatible exception handling */

/* Exception type passed in the type field of exception struct */

#define _DOMAIN    1 /* argument domain error */
#define _SING      2 /* argument singularity */
#define _OVERFLOW  3 /* overflow range error */
#define _UNDERFLOW 4 /* underflow range error */
#if !defined(__linux__) || defined(__PURE_INTEL_C99_HEADERS__)
#define _TLOSS      5   /* total loss of precision */
#define _PLOSS      6   /* partial loss of precision */
#endif

typedef struct ____exception {
    _LIBIMF_INT     type;
    _LIBIMF_CCHAR  *name;
    _LIBIMF_DOUBLE  arg1;
    _LIBIMF_DOUBLE  arg2;
    _LIBIMF_DOUBLE  retval;
} ___exception;

#if defined(__linux__) || defined(__APPLE__) || defined(__QNX__)
# if defined(__cplusplus)
#  define __exception ____exception   /* map 'struct __exception'  to 'struct ____exception'  */
# else     /*__cplusplus*/
#  define   exceptionf ____exceptionf /* map 'struct   exceptionf' to 'struct ____exceptionf' */
#  define   exceptionl ____exceptionl /* map 'struct   exceptionl' to 'struct ____exceptionl' */
# endif    /*__cplusplus*/
#else /* Win32 or Win64 */
# define _exception ____exception
# define _exceptionf ____exceptionf
# define _exceptionl ____exceptionl
#endif

typedef struct ____exceptionf {
    _LIBIMF_INT    type;
    _LIBIMF_CCHAR *name;
    _LIBIMF_FLOAT  arg1;
    _LIBIMF_FLOAT  arg2;
    _LIBIMF_FLOAT  retval;
} ___exceptionf;

typedef struct ____exceptionl {
    _LIBIMF_INT      type;
    _LIBIMF_CCHAR   *name;
    _LIBIMF_XDOUBLE  arg1;
    _LIBIMF_XDOUBLE  arg2;
    _LIBIMF_XDOUBLE  retval;
} ___exceptionl;

#if !defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && defined (__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_INT _LIBIMF_PUBAPI _matherr( struct ____exception  *__e );
#endif /*!__linux__ && !__APPLE__*/
_LIBIMF_EXT _LIBIMF_INT _LIBIMF_PUBAPI matherrf( struct ____exceptionf *__e );
_LIBIMF_EXT _LIBIMF_INT _LIBIMF_PUBAPI matherrl( struct ____exceptionl *__e );

#if !defined(__linux__) && !defined(__APPLE__) && !defined(__QNX__) && defined (__PURE_INTEL_C99_HEADERS__)
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _hypot( _LIBIMF_DOUBLE __x, _LIBIMF_DOUBLE __y );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _j0( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _j1( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _jn( _LIBIMF_INT __n, _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _y0( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _y1( _LIBIMF_DOUBLE __x );
_LIBIMF_EXT _LIBIMF_DOUBLE   _LIBIMF_PUBAPI _yn( _LIBIMF_INT __n, _LIBIMF_DOUBLE __x );
#endif
/* User-installable exception handlers
 *
 *  Static redefinition of matherr() is useful only for statically linked
 *  libraries. When Libm is built as a DLL, the Libm's matherr() is already
 *  loaded into the DLL and (statically) linked.  In this case, the only way
 *  to replace the library default matherr() with your matherr() is to use
 *  the matherr() exchange functions (see description below).
 *
 *  1. In user code, implement your own substitute matherr() function.
 *  2. To install it, call __libm_setusermatherr(), with your
 *     function as an argument. Note that the __libm_setusermatherr()
 *     returns the address of previously defined matherr. If you save
 *     the address, you can use it later to restore the original matherr().
 *  3. Your matherr() will now be installed! Your matherr() will be called
 *     instead of the default matherr().
 */

typedef _LIBIMF_INT ( _LIBIMF_PUBAPI *___pmatherr )( struct ____exception  *__e );
typedef _LIBIMF_INT ( _LIBIMF_PUBAPI *___pmatherrf )( struct ____exceptionf *__e );
typedef _LIBIMF_INT ( _LIBIMF_PUBAPI *___pmatherrl )( struct ____exceptionl *__e );

_LIBIMF_EXT ___pmatherr  _LIBIMF_PUBAPI __libm_setusermatherr( ___pmatherr  __user_matherr );
_LIBIMF_EXT ___pmatherrf _LIBIMF_PUBAPI __libm_setusermatherrf( ___pmatherrf __user_matherrf );
_LIBIMF_EXT ___pmatherrl _LIBIMF_PUBAPI __libm_setusermatherrl( ___pmatherrl __user_matherrl );

/* Standard conformance support */

#if (!defined(__linux__) || !defined(__USE_MISC)) || defined (__PURE_INTEL_C99_HEADERS__)
typedef enum ___LIB_VERSIONIMF_TYPE {
     _IEEE_ = -1    /* IEEE-like behavior    */
    ,_SVID_         /* SysV, Rel. 4 behavior */
    ,_XOPEN_        /* Unix98                */
    ,_POSIX_        /* Posix                 */
    ,_ISOC_         /* ISO C9X               */
} _LIB_VERSIONIMF_TYPE;
#else
# define _LIB_VERSIONIMF_TYPE _LIB_VERSION_TYPE
#endif

_LIBIMF_EXT _LIB_VERSIONIMF_TYPE _LIBIMF_PUBVAR _LIB_VERSIONIMF;

#undef _LIBIMF_EXT
#undef _LIBIMF_PUBAPI
#undef _LIBIMF_PUBVAR
#undef _LIBIMF_LONGLONG
#undef _LIBIMF_LONGINT
#undef _LIBIMF_INT
#undef _LIBIMF_DOUBLE
#undef _LIBIMF_FLOAT
#undef _LIBIMF_XDOUBLE
#undef _LIBIMF_VOID
#undef _LIBIMF_CCHAR

#ifdef __STDC_WANT_DEC_FP__ 
extern _Decimal32   acosd32  (_Decimal32 x);
extern _Decimal64   acosd64  (_Decimal64 x);
extern _Decimal128  acosd128 (_Decimal128 x);
extern _Decimal32   asind32  (_Decimal32 x);
extern _Decimal64   asind64  (_Decimal64 x);
extern _Decimal128  asind128 (_Decimal128 x);
extern _Decimal32   atand32  (_Decimal32 x);
extern _Decimal64   atand64  (_Decimal64 x);
extern _Decimal128  atand128 (_Decimal128 x);
extern _Decimal32   atan2d32  (_Decimal32 y, _Decimal32 x);
extern _Decimal64   atan2d64  (_Decimal64 y, _Decimal64 x);
extern _Decimal128  atan2d128 (_Decimal128 y, _Decimal128 x);
extern _Decimal32   cosd32  (_Decimal32 x);
extern _Decimal64   cosd64  (_Decimal64 x);
extern _Decimal128  cosd128 (_Decimal128 x);
extern _Decimal32   sind32  (_Decimal32 x);
extern _Decimal64   sind64  (_Decimal64 x);
extern _Decimal128  sind128 (_Decimal128 x);
extern _Decimal32   tand32  (_Decimal32 x);
extern _Decimal64   tand64  (_Decimal64 x);
extern _Decimal128  tand128 (_Decimal128 x);
extern _Decimal32   acoshd32  (_Decimal32 x);
extern _Decimal64   acoshd64  (_Decimal64 x);
extern _Decimal128  acoshd128 (_Decimal128 x);
extern _Decimal32   asinhd32  (_Decimal32 x);
extern _Decimal64   asinhd64  (_Decimal64 x);
extern _Decimal128  asinhd128 (_Decimal128 x);
extern _Decimal32   atanhd32  (_Decimal32 x);
extern _Decimal64   atanhd64  (_Decimal64 x);
extern _Decimal128  atanhd128 (_Decimal128 x);
extern _Decimal32   coshd32  (_Decimal32 x);
extern _Decimal64   coshd64  (_Decimal64 x);
extern _Decimal128  coshd128 (_Decimal128 x);
extern _Decimal32   sinhd32  (_Decimal32 x);
extern _Decimal64   sinhd64  (_Decimal64 x);
extern _Decimal128  sinhd128 (_Decimal128 x);
extern _Decimal32   tanhd32  (_Decimal32 x);
extern _Decimal64   tanhd64  (_Decimal64 x);
extern _Decimal128  tanhd128 (_Decimal128 x);
extern _Decimal32   expd32   (_Decimal32 x);
extern _Decimal64   expd64   (_Decimal64 x);
extern _Decimal128  expd128  (_Decimal128 x);
extern _Decimal32   exp2d32  (_Decimal32 x);
extern _Decimal64   exp2d64  (_Decimal64 x);
extern _Decimal128  exp2d128 (_Decimal128 x);
extern _Decimal32   expm1d32  (_Decimal32 x);
extern _Decimal64   expm1d64  (_Decimal64 x);
extern _Decimal128  expm1d128 (_Decimal128 x);
extern _Decimal32   frexpd32  (_Decimal32 x, int *i);
extern _Decimal64   frexpd64  (_Decimal64 x, int *i);
extern _Decimal128  frexpd128 (_Decimal128 x, int *i);
extern int   ilogbd32  (_Decimal32 x);
extern int   ilogbd64  (_Decimal64 x);
extern int   ilogbd128 (_Decimal128 x);
extern _Decimal32   ldexpd32  (_Decimal32 x, int n);
extern _Decimal64   ldexpd64  (_Decimal64 x, int n);
extern _Decimal128  ldexpd128 (_Decimal128 x, int n);
extern _Decimal32   logd32  (_Decimal32 x);
extern _Decimal64   logd64  (_Decimal64 x);
extern _Decimal128  logd128 (_Decimal128 x);
extern _Decimal32   log10d32  (_Decimal32 x);
extern _Decimal64   log10d64  (_Decimal64 x);
extern _Decimal128  log10d128 (_Decimal128 x);
extern _Decimal32   log2d32  (_Decimal32 x);
extern _Decimal64   log2d64  (_Decimal64 x);
extern _Decimal128  log2d128 (_Decimal128 x);
extern _Decimal32   log1pd32  (_Decimal32 x);
extern _Decimal64   log1pd64  (_Decimal64 x);
extern _Decimal128  log1pd128 (_Decimal128 x);
extern _Decimal32   logbd32  (_Decimal32 x);
extern _Decimal64   logbd64  (_Decimal64 x);
extern _Decimal128  logbd128 (_Decimal128 x);
extern _Decimal32   modfd32  (_Decimal32 x, _Decimal32 *iptr);
extern _Decimal64   modfd64  (_Decimal64 x, _Decimal64 *iptr);
extern _Decimal128  modfd128 (_Decimal128 x, _Decimal128 *iptr);
extern _Decimal32   scalbnd32  (_Decimal32 x, int n);
extern _Decimal64   scalbnd64  (_Decimal64 x, int n);
extern _Decimal128  scalbnd128 (_Decimal128 x, int n);
extern _Decimal32   scalblnd32  (_Decimal32 x, long int n);
extern _Decimal64   scalblnd64  (_Decimal64 x, long int n);
extern _Decimal128  scalblnd128 (_Decimal128 x, long int n);
extern _Decimal32   cbrtd32  (_Decimal32 x);
extern _Decimal64   cbrtd64  (_Decimal64 x);
extern _Decimal128  cbrtd128 (_Decimal128 x);
extern _Decimal32   fabsd32  (_Decimal32 x);
extern _Decimal64   fabsd64  (_Decimal64 x);
extern _Decimal128  fabsd128 (_Decimal128 x);
extern _Decimal32   hypotd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   hypotd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  hypotd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   powd32   (_Decimal32 x, _Decimal32 y);
extern _Decimal64   powd64   (_Decimal64 x, _Decimal64 y);
extern _Decimal128  powd128  (_Decimal128 x, _Decimal128 y);
extern _Decimal32   sqrtd32  (_Decimal32 x);
extern _Decimal64   sqrtd64  (_Decimal64 x);
extern _Decimal128  sqrtd128 (_Decimal128 x);
extern _Decimal32   erfd32   (_Decimal32 x);
extern _Decimal64   erfd64   (_Decimal64 x);
extern _Decimal128  erfd128  (_Decimal128 x);
extern _Decimal32   erfcd32  (_Decimal32 x);
extern _Decimal64   erfcd64  (_Decimal64 x);
extern _Decimal128  erfcd128 (_Decimal128 x);
extern _Decimal32   lgammad32  (_Decimal32 x);
extern _Decimal64   lgammad64  (_Decimal64 x);
extern _Decimal128  lgammad128 (_Decimal128 x);
extern _Decimal32   tgammad32  (_Decimal32 x);
extern _Decimal64   tgammad64  (_Decimal64 x);
extern _Decimal128  tgammad128 (_Decimal128 x);
extern _Decimal32   ceild32   (_Decimal32 x);
extern _Decimal64   ceild64   (_Decimal64 x);
extern _Decimal128  ceild128  (_Decimal128 x);
extern _Decimal32   floord32  (_Decimal32 x);
extern _Decimal64   floord64  (_Decimal64 x);
extern _Decimal128  floord128 (_Decimal128 x);
extern _Decimal32   nearbyintd32  (_Decimal32 x);
extern _Decimal64   nearbyintd64  (_Decimal64 x);
extern _Decimal128  nearbyintd128 (_Decimal128 x);
extern _Decimal32      rintd32    (_Decimal32 x);
extern _Decimal64      rintd64    (_Decimal64 x);
extern _Decimal128     rintd128   (_Decimal128 x);
extern long int        lrintd32   (_Decimal32 x);
extern long int        lrintd64   (_Decimal64 x);
extern long int        lrintd128  (_Decimal128 x);
extern long long int   llrintd32  (_Decimal32 x);
extern long long int   llrintd64  (_Decimal64 x);
extern long long int   llrintd128 (_Decimal128 x);
extern _Decimal32      roundd32   (_Decimal32 x);
extern _Decimal64      roundd64   (_Decimal64 x);
extern _Decimal128     roundd128  (_Decimal128 x);
extern long int        lroundd32  (_Decimal32 x);
extern long int        lroundd64  (_Decimal64 x);
extern long int        lroundd128 (_Decimal128 x);
extern long long int   llroundd32  (_Decimal32 x);
extern long long int   llroundd64  (_Decimal64 x);
extern long long int   llroundd128 (_Decimal128 x);
extern _Decimal32   truncd32  (_Decimal32 x);
extern _Decimal64   truncd64  (_Decimal64 x);
extern _Decimal128  truncd128 (_Decimal128 x);
extern _Decimal32   fmodd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   fmodd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  fmodd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   remainderd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   remainderd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  remainderd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   copysignd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   copysignd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  copysignd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   nand32  (char *tagp);
extern _Decimal64   nand64  (char *tagp);
extern _Decimal128  nand128 (char *tagp);
extern _Decimal32   nextafterd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   nextafterd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  nextafterd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   nexttowardd32  (_Decimal32 x, _Decimal128 y);
extern _Decimal64   nexttowardd64  (_Decimal64 x, _Decimal128 y);
extern _Decimal128  nexttowardd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   fdimd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   fdimd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  fdimd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   fmaxd32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   fmaxd64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  fmaxd128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   fmind32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   fmind64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  fmind128 (_Decimal128 x, _Decimal128 y);
extern _Decimal32   fmad32  (_Decimal32 x,  _Decimal32 y, _Decimal32 z);
extern _Decimal64   fmad64  (_Decimal64 x,  _Decimal64 y, _Decimal64 z);
extern _Decimal128  fmad128 (_Decimal128 x, _Decimal128 y, _Decimal128 z);
extern _Decimal32   quantized32  (_Decimal32 x, _Decimal32 y);
extern _Decimal64   quantized64  (_Decimal64 x, _Decimal64 y);
extern _Decimal128  quantized128 (_Decimal128 x, _Decimal128 y);
extern int   samequantumd32  (_Decimal32 x, _Decimal32 y);
extern int   samequantumd64  (_Decimal64 x, _Decimal64 y);
extern int   samequantumd128 (_Decimal128 x, _Decimal128 y);
extern int   quantexpd32  (_Decimal32 x);
extern int   quantexpd64  (_Decimal64 x);
extern int   quantexpd128 (_Decimal128 x);
extern int   isnand32  (_Decimal32 x);
extern int   isnand64  (_Decimal64 x);
extern int   isnand128 (_Decimal128 x);
extern int   isinfd32  (_Decimal32 x);
extern int   isinfd64  (_Decimal64 x);
extern int   isinfd128 (_Decimal128 x);
extern int   isfinited32  (_Decimal32 x);
extern int   isfinited64  (_Decimal64 x);
extern int   isfinited128 (_Decimal128 x);
extern int   isnormald32  (_Decimal32 x);
extern int   isnormald64  (_Decimal64 x);
extern int   isnormald128 (_Decimal128 x);
extern int   signbitd32  (_Decimal32 x);
extern int   signbitd64  (_Decimal64 x);
extern int   signbitd128 (_Decimal128 x);
extern int   fpclassifyd32  (_Decimal32 x);
extern int   fpclassifyd64  (_Decimal64 x);
extern int   fpclassifyd128 (_Decimal128 x);
extern int   isunorderedd32  (_Decimal32 x, _Decimal32 y);
extern int   isunorderedd64  (_Decimal64 x, _Decimal64 y);
extern int   isunorderedd128 (_Decimal128 x, _Decimal128 y);
extern int   isgreaterd32  (_Decimal32 x, _Decimal32 y);
extern int   isgreaterd64  (_Decimal64 x, _Decimal64 y);
extern int   isgreaterd128 (_Decimal128 x, _Decimal128 y);
extern int   isgreaterequald32  (_Decimal32 x, _Decimal32 y);
extern int   isgreaterequald64  (_Decimal64 x, _Decimal64 y);
extern int   isgreaterequald128 (_Decimal128 x, _Decimal128 y);
extern int   islessd32  (_Decimal32 x, _Decimal32 y);
extern int   islessd64  (_Decimal64 x, _Decimal64 y);
extern int   islessd128 (_Decimal128 x, _Decimal128 y);
extern int   islessequald32  (_Decimal32 x, _Decimal32 y);
extern int   islessequald64  (_Decimal64 x, _Decimal64 y);
extern int   islessequald128 (_Decimal128 x, _Decimal128 y);
#endif  /*__STDC_WANT_DEC_FP__*/

#if defined(__cplusplus)    /* end C linkage */
}
#endif

#endif  /*__MATH_H_INCLUDED*/
#endif  /* usage of sys headers */


