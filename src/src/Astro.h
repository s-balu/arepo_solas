#pragma once

#define MaxCharactersInLine (256)
#define NONE    (-1)

#define	SQ(x)		((x)*(x))
#define	CUBE(x)		((x)*(x)*(x))
#define	SWAP(type,x,y)	{type t = x; x = y; y = t; }
#define	swap(x,y)	(tmp=(x),(x)=(y),(y)=tmp)
#define	MAX(x,y)	(((x)>(y))?(x):(y))
#define	MIN(x,y)	(((x)<(y))?(x):(y))
#define	ABS(x)		(((x)>0)?(x):(-(x)))
#define	SGN(x)		(((x)>0)?(+1):(((x)<0)?(-1):(0)))
#define BitSign(x)      (((x)==1)?(+1):(-1))
#define Sign(x)         (((x)>=0)?(+1):(-1))
#define Parity(x)   (((x)%2==0)?(+1):(-1))

#define	NORM(x)		( sqrt( SQ((x)[0]) + SQ((x)[1]) + SQ((x)[2]) ) )
#define	NORM2(x)	( SQ((x)[0]) + SQ((x)[1]) + SQ((x)[2]) )
#define	DISTANCE(x,y)	( sqrt(SQ((x)[0]-(y)[0])+SQ((x)[1]-(y)[1])+SQ((x)[2]-(y)[2])) )
#define	DISTANCE2(x,y)	(SQ((x)[0]-(y)[0])+SQ((x)[1]-(y)[1])+SQ((x)[2]-(y)[2]))
#define	DOT_PRODUCT(x,y)	(((x)[0]*(y)[0])+((x)[1]*(y)[1])+((x)[2]*(y)[2]))

#define dprint(a)       printf(" ### " #a " = %d:%s:line%d:%s()\n", a,__FILE__,__LINE__,__FUNCTION__)
#define dlprint(a)      printf(" ### " #a " = %ld:%s:line%d:%s()\n", a,__FILE__,__LINE__,__FUNCTION__)
#define gprint(a)       printf(" ### " #a " = %g:%s:line%d:%s()\n", a,__FILE__,__LINE__,__FUNCTION__)
#define sprint(a)       printf(" ### " #a " %s:line%d:%s()\n", a, __FILE__,__LINE__,__FUNCTION__)
#define pprintl(a)      printf(" ### " #a " = %p:%s:line%d:%s()\n", a,__FILE__,__LINE__,__FUNCTION__)

#define FileOpen(fp, fname, mode)				\
	(fp) = fopen((fname),(mode));				\
	if ((fp) == NULL){					\
		fprintf(stderr,"Can't open file %s\n",(fname));	\
		exit(1);					\
	}

#define Snprintf(buf,expression, ...)                       \
    if( snprintf((buf),MaxCharactersInLine,expression, __VA_ARGS__ ) >= (MaxCharactersInLine) ){ \
        fprintf(stderr,"buffer overflow: %s:line%d:%s\n",       \
                __FILE__,__LINE__,__FUNCTION__);              \
        exit(1);                                            \
    }
