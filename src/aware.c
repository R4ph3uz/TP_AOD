/**
 * \file aware.c
 * \brief cache aware implementation that computes the distance between two genetic sequences 
 *
 * Documentation: see aware.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in iteratif.h
 */

#include "aware.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> 

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
    
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 
#define Z 4096
#define L 4

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
} ;

long EditDistance_Aware(char* A, size_t lengthA, char* B, size_t lengthB)
{
    _init_base_match() ;
    struct MemoContext ctx;
    if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
    {  ctx.X = A ;
       ctx.M = lengthA ;
       ctx.Y = B ;
       ctx.N = lengthB ;
    }
    else
    {  ctx.X = B ;
       ctx.M = lengthB ;
       ctx.Y = A ;
       ctx.N = lengthA ;
    }
    size_t M = ctx.M ;
    size_t N = ctx.N ;
    
    long K_c = 1;
    long K_l = (Z - 2*L - 5*L) / (L + 1);

    /* Allocation and initialization of tab_phi_N and tab_phi_M to NOT_YET_COMPUTED*/
    long* tab_phi_N = (long *) malloc ( (N+1) * sizeof(long));
    long* tab_phi_M = (long *) malloc ( (M+1) * sizeof(long));
    if (tab_phi_N == NULL) { perror("EditDistance_Iterativ: malloc of tab_phi_N" ); exit(EXIT_FAILURE); }
    for (int j=0; j<=N; ++j) tab_phi_N[j] = NOT_YET_COMPUTED ;
    if (tab_phi_M == NULL) { perror("EditDistance_Iterativ: malloc of tab_phi_M" ); exit(EXIT_FAILURE); }
    for (int i=0; i<=M; ++i) tab_phi_M[i] = NOT_YET_COMPUTED ;

    /*On va parcourir notre table de phi(i,j) par blocs*/

    /*On initialise phi(M,N)*/
    tab_phi_N[N] = 0;
    tab_phi_M[M] = 0;

    for (long j = N-1; j >= 0; j--) tab_phi_N[j] = (isBase(ctx.Y[j]) ? INSERTION_COST : 0) + tab_phi_N[j+1] ;
    for (long i = M-1; i >= 0; i--) tab_phi_M[i] = (isBase(ctx.X[i]) ? INSERTION_COST : 0) + tab_phi_M[i+1] ;

    for (long I = M; I >0; I -= K_l){
        long i_end = (I - K_l > 0) ? (I - K_l) : 0;

        for (long J = N; J >0; J -= K_c){
            long j_end = (J - K_c > 0) ? (J - K_c) : 0;

            for (long i = I-1; i >= i_end; i--){
                long temp1 = tab_phi_M[i] ;

                for (long j = J-1; j >= j_end; j--){
                    long temp2;

                    if (! isBase(ctx.X[i])){ /* skip ccharacter in Xi that is not a base */
                        ManageBaseError( ctx.X[i] ) ;
                        temp2 = tab_phi_N[j] ;
                    }
                    else if (! isBase(ctx.Y[j])){  /* skip ccharacter in Yj that is not a base */
                        ManageBaseError( ctx.Y[j] ) ;
                        temp2 = temp1 ;
                    }
                    else{ 
                        long min = /* initialization  with cas 1*/
                            ( isUnknownBase(ctx.X[i]) ?  SUBSTITUTION_UNKNOWN_COST 
                                    : ( isSameBase(ctx.X[i], ctx.Y[j]) ? 0 : SUBSTITUTION_COST ) 
                            )
                            + ((j==(J-1) && i==(I-1))? tab_phi_M[i+1] : tab_phi_N[j+1]) ; 

                        long cas2 = INSERTION_COST + tab_phi_N[j] ;      
                        if (cas2 < min) min = cas2 ;

                        long cas3 = INSERTION_COST + temp1 ;      
                        if (cas3 < min) min = cas3 ; 

                        temp2 = min ;
                    }
                    
                    tab_phi_N[j+1] = temp1;
                    temp1 = temp2;
                }
                
                tab_phi_M[i+1]=tab_phi_N[j_end];
                tab_phi_M[i]=temp1;
                tab_phi_N[j_end]= temp1;
            }
        } 
        tab_phi_M[i_end] = tab_phi_N [N];
    }

    long res = tab_phi_N[0]; 
    
    /* Deallocation of tab_phi_N */
    free( tab_phi_N ) ;
    free( tab_phi_M ) ;

    return res ;
}

