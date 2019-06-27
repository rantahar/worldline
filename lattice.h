#ifndef LATTICE_H
#define LATTICE_H

#ifndef MAIN
#define EXTERN extern
#else
#define EXTERN
#endif


#define ND 2
#define NDIRS (2*ND)

#define TUP 0
#define XUP 1
#define TDN 2
#define XDN 3


/* Neighbour index arrays, to be filled at the beginning
 */
EXTERN int *tup,*xup,*tdn,*xdn;

/* Utilities */
/* Functions for fetching neighboring coordinates */
static inline int tdir(int t, int dir){
  if( dir == TUP ) return tup[t];
  if( dir == TDN ) return tdn[t];
  return(t);
}

static inline int xdir(int x, int dir){
  if( dir == XUP ) return xup[x];
  if( dir == XDN ) return xdn[x];
  return(x);
}

/* Opposite of a direction */
static inline int opp_dir(int dir){
  return ( dir + ND ) % NDIRS;
}


/* Print error and exit */
static inline void errormessage( char * message ){
  fprintf( stderr, "%s", message );
  exit(1);
}


/* Ask for parameter */
static inline void get_int( char * name, int * dest ){
  printf(" %s :", name);
  if( scanf("%d", dest) == 0 ){
    char message[60];
    sprintf(message, "Missing parameter %s\n", name);
    errormessage(message);
  }
}

static inline void get_double( char * name, double * dest ){
  printf(" %s :", name);
  if( scanf("%lf", dest) == 0 ){
    char message[60];
    sprintf(message, "Missing parameter %s\n", name);
    errormessage(message);
  }
}

static inline void get_long( char * name, long * dest ){
  printf(" %s :", name);
  if( scanf("%ld", dest) == 0 ){
    char message[60];
    sprintf(message, "Missing parameter %s\n", name);
    errormessage(message);
  }
}

static inline void get_char( char * name, char * dest ){
  printf(" %s :", name);
  if( scanf("%s", dest) == 0 ){
    char message[60];
    sprintf(message, "Missing parameter %s\n", name);
    errormessage(message);
  }
}


#endif