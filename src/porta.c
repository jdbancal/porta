/*******************************************************************************

Copyright (C) 1997-2009 Thomas Christof and Andreas Loebel
 
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.
 
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA 02111-1307 USA
 

FILENAME: porta.c

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/
// Last edit on 4.4.2013 by J-D Bancal
// File modified to cope with larger integers in the denominator
/*  LAST EDIT: Fri Sep 20 13:36:44 2002 by Andreas Loebel (opt0.zib.de)  */
/* $Id: porta.c,v 1.5 2009/09/21 07:46:48 bzfloebe Exp $ */


#define _CRT_SECURE_NO_WARNINGS 1


#include "porta.h"
#include "arith.h"
#include "common.h"
#include "log.h"
#include "inout.h"
#include "mp.h"
#include "four_mot.h"
#include "portsort.h"


// The following variables were previously instanciated in the header file. Now they have been externalized here. J-D B
int mp_state;

void (*RAT_add)(),
	 (*RAT_sub)(),
	 (*RAT_mul)(),
	 (*RAT_row_prim)(RAT *, RAT *, RAT *, int), // Line completed on 26.2.2013 by J-D B to be compatible with the c++ library interface portalib.cpp
	 (*RAT_assign)(),
	 (*writeline)();

RAT RAT_const[2],var[4];

RAT *ar1,*ar2,*ar3,*ar4,*ar5,*ar6;
long nel_ar1,nel_ar2,nel_ar3,nel_ar4,nel_ar5,nel_ar6;

int  maxlist, total_size;

int  dim,
     equa,    /* number of equalities */
     ineq,    /* number of inequalities */
     conv, 
     cone, 
     points,
     blocks;

FILE *fp,*prt;
char * allo();

int option, allowed_options;



//FILE *logfile;


/* On 28.2.2013, by J-D B.
 * The main function was moved to file xporta.c so that porta.h can be
 * seen as a library by other programs. */

//int main( int argc, char *argv[] )
//{
		  //// A simple polytope :
		  ///*long int data[] = {1, 0, -4, 0, 0, 1, 0, 3, 0, 0, 1, 0, 0, 0, -2, 1, 0, 0, 0, 2};
		  //int dimension = 4;
		  //int nbPoints = 4;*/
		  
		  //// Another simple polytope :
		  ///*long int data[] = {1, 2, 3, 4, 5, 6, 7, 8};
		  //int dimension = 3;
		  //int nbPoints = 2;*/
		  
		  //// The CHSH polytope :
		  ////long int data[] = {1, -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		  ////int dimension = 8;
		  ////int nbPoints = 16;
		  
		  ////quick_and_dirty_poi_conv_call(&data, dimension, nbPoints);
	
		  ////exit(0);
	
	
    //int i, ieq_file, start;
    //char outfname[20];
    //char fname[20];
    //int   poi_file;
    //int   rowl_inar, ierl;
    //int  *indx = (int *)0;      
    //int equa_in,ineq_in, ineq_out;
    //FILE *outfp;

    

    //printf("\nPORTA - a POlyhedron Representation Transformation Algorithm\n");
    //printf(  "Version %s\n\n", VERSION );

    //printf( "Written by Thomas Christof (Uni Heidelberg)\n" );
    //printf( "Revised by Andreas Loebel (ZIB Berlin)\n\n" );

    //printf( "PORTA is free software and comes with ABSOLUTELY NO WARRENTY! You are welcome\n" );
    //printf( "to use, modify, and redistribute it under the GNU General Public Lincese.\n\n" ); 
    
    //printf( "This is the program XPORTA from the PORTA package.\n\n" );

    //if( argc <= 2 )
    //{
        //printf( "For more information read the manpages about porta.\n\n" );
        //exit(-1);
    //}
            
    ///* 17.01.1994: include logging on file porta.log */
    //logfile = fopen( "porta.log", "a" );
    //if( !logfile )
        //fprintf( stderr, "can't open logfile porta.log\n" );
    //else
    //{
        //porta_log( "\n\n\nlog for " );
        //for( i = 0; i < argc; i++ )
            //porta_log( "%s ", argv[i] );
        //porta_log( "\n\n" );
    //}
            
    //init_total_time();
    
    //initialize();

    //prt = stdout;
    //get_options(&argc,&argv);
    
    //if (option & Protocol_to_file) 
    //{
        //strcat(*argv,".prt");
        //prt = fopen(*argv,"w");
        //(*argv)[strlen(*argv)-4] = '\0';
    //}
    //setbuf(prt,CP 0);
    
    //set_I_functions();
    //SET_MP_not_ready;
    //ieq_file = !strcmp(*argv+strlen(*argv)-4,".ieq");
    //poi_file = !strcmp(*argv+strlen(*argv)-4,".poi");
    //printf("%i, %i",ieq_file, poi_file);
    //if (!poi_file && !ieq_file)
        //msg( "invalid format of command line", "", 0 );
    
    ///*
     //* change by M.S. 5.6.92:
     //* read_input_file writes to the output file, if is_set(Sort).
     //*/
    //outfp = 0;
    //strcpy(outfname,*argv);
    //if (is_set(Sort) && poi_file) 
    //{
        //strcat(outfname,".poi");
        //outfp = wfopen(outfname);
    //}
    //if (is_set(Sort) && ieq_file) 
    //{
        //strcat(outfname,".ieq");
        //fprintf(prt,"outfname = %s\n",outfname);
        //fflush(stdout);

        ///* 17.01.1994: include logging on file porta.log */
        //porta_log( "outfname = %s\n",outfname );
        //fflush(logfile);

        //outfp = wfopen(outfname);
    //}

    //if (is_set(Fmel) && ieq_file) 
    //{  
        ///* ONLY FM-ELIMINATON */
        
        //int *elim_ord,nel;
            //char *cp1;
        //elim_ord = 0;
        //cp1=strdup("ELIMINATION_ORDER");
        //if(!cp1)
            //msg( "allocation of new space failed", "", 0 );   
        //ineq = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,cp1,
                               //&elim_ord,"\0",(int **)&i,"\0",(RAT **)&i);
        //free(cp1);
        //sort_eqie_cvce(ar1,ineq,dim+2,&equa_in,&ineq_in);
        //ineq = ineq_in+equa_in;
        ///*     elim_ord = check_and_reorder_elim_ord(elim_ord,&nel); */
        //reorder_var(ineq,ar1,&ar2,(int *)&nel_ar2,&nel,&elim_ord,&indx);
        ///*     indx = elim_ord+nel; */
        ///* 
         //* Transform all inequalities into "<= 1" or "<=-1" inequalities,
         //* if possible.
         //* If the right-hand side is 0, divide the numerators by their gcd
         //* and the denominators by their gcd.
         //* (This is not really necessary).
         //*/
        ///*
           //for (i = 0; i < ineq; i++) 
           //(* RAT_row_prim)(ar2+(dim+2)*i,ar2+(dim+1)*i,ar2+(dim+2)*i+dim,dim+1);
        //*/
        //if(is_set(Long_arithmetic))
        //{
            //RAT a, b;
            //SET_MP_ready;
            //memset( &a, 0, sizeof(a) );
            //memset( &b, 0, sizeof(b) );
            //arith_overflow_func(0,0,a,b,0);
        //}
        //for (i = 0; i < ineq; i++) 
            //(* RAT_row_prim)(ar2+(dim+1)*i,ar2+(dim+1)*i,
                             //ar2+(dim+1)*i+dim,dim+1);
        ///*     nel_ar2 = ineq*(dim+1);  
         //*/
        //equa = 0;
        //ineq_out = ineq;
        //gauss(0, dim+1, dim+1, dim-nel, equa_in, &ineq_out, &equa, indx);
        //for (; (*indx) < 0; indx++);
        //ierl = dim-nel-equa +1;
        ///* row-length of inequalities after fourier_motzkin
         //* =  number of noneliminated variables+1 */
        //nel = nel - (ineq - ineq_out);  /* number of variables to be elim. */
        //ineq = ineq_out;
        //fourier_motzkin(0,ineq-equa,dim+1-equa_in,nel,poi_file,indx,elim_ord);
        //if ((MP_realised && return_from_mp()) || !MP_realised) 
            //sort(no_denom(ierl, 0, ineq,1), ierl, 0, ineq);
        //write_ieq_file(*argv,outfp,equa,ineq,dim+1,0,
                       //ineq,0,ierl,indx);
        
    //}
    //else if (is_set(Sort)) 
    //{
        //points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,"\0",(int **)&i,
                                 //"\0",(int **)&i,"\0",(RAT **)&i);
        //if (ieq_file)
            //sort_eqie_cvce(ar1,points,dim+2,&equa,&ineq);
        //listptoar(ar1,points,ieq_file?dim+2:dim+1,0); 
        //if (ieq_file) 
        //{ 
            //if (equa) sort(1,dim+1,0,equa);
            //if (ineq) sort(1,dim+1,equa,points); 
            //write_ieq_file(*argv,outfp,equa,0,dim+1,0,ineq,equa,dim+1,0);
        //}
        //else 
        //{
            //sort(1,dim+1,0,points);
            //for (cone = 0; !(porta_list[cone]->sys[dim].num); cone++);
            //write_poi_file(*argv,outfp,dim,0,0,cone,0,points-cone,cone);
        //}
    //}  
    //else if ((is_set(Traf) || is_set(Dim)) && poi_file) 
    //{
        //points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
                                 //"\0",(int **)&i,"\0",(int **)&i,"\0",
                                 //(RAT **)&i);
        //gentableau(ar1,1,&rowl_inar,&indx); 
        //if(is_set(Long_arithmetic))
        //{
            //RAT a, b;
            //SET_MP_ready;
            //memset( &a, 0, sizeof(a) );
            //memset( &b, 0, sizeof(b) );
            //arith_overflow_func(0,0,a,b,0);
        //}
        //ineq = (cone == points) ? dim : dim + 1;
        //ineq_out = ineq;  /*not used further */
        //gauss(1, points+dim+1,dim+1,dim,ineq,&ineq_out, &equa, indx);
        ///* make indx point to the system-variable section */
        //for (; (*indx) < 0; indx++);
        //if (is_set(Dim)) 
        //{
            //char str[100];
            //fprintf (prt,"\nDIMENSION OF THE POLYHEDRON : %i\n\n",dim-equa);
            
            ///* 17.01.1994: include logging on file porta.log */
            //porta_log( "\nDIMENSION OF THE POLYHEDRON : %i\n\n", dim-equa );
            
            //sprintf (str,"echo 'DIMENSION OF THE POLYHEDRON : %i' | cat >> %s",
                     //dim-equa,argv[0]);
            //system(str);
            //if (equa) 
            //{
                //fprintf(prt,"equations :\n");

                ///* 17.01.1994: include logging on file porta.log */
                //porta_log( "equations :\n");
                
                //listptoar(ar4,equa,dim+1,0); 
                //if ((MP_realised && return_from_mp()) || !MP_realised) 
                    //sort(no_denom(dim+1,0,equa,1), dim+1,0,equa);
                //start = 1;
                ///* last argument of writesys was lost? 
                   //writesys(prt,0,equa,dim+1,0,0,'=');
                   //*/
                //writesys(prt,0,equa,dim+1,0,0,'=', &start);

                ///* log equation system */
                //start = 1;
                //writesys(logfile,0,equa,dim+1,0,0,'=', &start);
            //}
        //}
        //else  
        //{
            ///* POINTS TO INEQUALITIES */
            //sprintf(fname,"%s.ieq",*argv);
            //fourier_motzkin(fname,ineq-equa,points+dim+1-ineq,
                            //points-ineq+equa,poi_file,indx,0);
            //if (is_set(Validity_table_out)) 
                //red_test(indx,ar1,&rowl_inar);
            //if ((MP_realised && return_from_mp()) || !MP_realised) 
            //{
                //if (equa) sort(no_denom(dim+1, ineq, ineq+equa,1), 
                               //dim+1, ineq, ineq+equa);
                //sort(no_denom(dim+1-equa, 0, ineq,1), 
                     //dim+1-equa, 0, ineq);
            //}
            //write_ieq_file(*argv,outfp,equa,ineq,
                           //dim+1,0,ineq,0,dim+1-equa,indx);
        //}
    //}
    //else if (is_set(Traf) && ieq_file) 
    //{
        ///* INEQUALITIES TO POINTS */
        //RAT *inner,*iep;
        //char *cp1;
        //cp1=strdup("VALID");
        //inner = 0;
        //if(!cp1)
            //msg( "allocation of new space failed", "", 0 );   
            
        //points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
                                 //"\0",(int **)&i,
                                 //"\0",(int **)&i,cp1,&inner);
        //free(cp1);
        //ar6 = inner; if (inner) nel_ar6 = dim;
        //sort_eqie_cvce(ar1,points,dim+2,&equa_in,&ineq_in);
        //iep = ar1+equa_in*(dim+2);
        ///* first equations then inequalities */
        //points = ineq_in;
        //polarformat(iep,&equa_in,ineq_in,inner);
        //gentableau(iep,0,&rowl_inar,&indx);
        //if(is_set(Long_arithmetic)) 
        //{
            //RAT a, b;
            //SET_MP_ready;
            //memset( &a, 0, sizeof(a) );
            //memset( &b, 0, sizeof(b) );
            //arith_overflow_func(0,0,a,b,0);
        //}
        //ineq = (cone == points) ? dim : dim + 1;
        //ineq_out = ineq; /* not used further */
        //gauss(1, points+dim+1,dim+1,dim,ineq,&ineq_out, &equa, indx);
        ///* make indx point to the x-variable section */
        //for (; (*indx) < 0; indx++);
        //fourier_motzkin(0,ineq-equa,points+dim+1-ineq,
                        //points-ineq+equa,poi_file,indx,0);
        //if (is_set(Validity_table_out)) 
            //red_test(indx,iep,&rowl_inar);
        //if (cone >= dim-equa)
            //origin_add(rowl_inar,iep); 
        //resubst(inner,equa_in,indx);
        //if ((MP_realised && return_from_mp()) || !MP_realised) 
        //{
            //if (equa)
                //sort(no_denom(dim+1, ineq, ineq+equa,1), 
                     //dim+1,ineq,ineq+equa);
            //sort(0, dim+1, 0, ineq);
        //}
        //for (cone = 0; !(porta_list[cone]->sys[dim].num); cone++);
        //conv = ineq - cone;
        //if (!MP_realised) no_denom(dim+1, 0, cone,1);
        //write_poi_file(*argv,outfp,dim,equa,ineq,cone,0,conv,cone);
    //}
    //else 
        //msg( "invalid format of command line", "", 0 );



    ///* 17.01.1994: include logging on file porta.log */
    //fclose( logfile );
    //exit(0);
//}


            







int *check_and_reorder_elim_ord( int *eo, int *nel )
/*
 * This routine is not needed any more.
 *
 * invert the elimination order:
 * The "invertion" h satisfies:
 * h[5-1] = 3-1     if 3 is the 5th variable to be eliminated
 * etc. up to h[*nel-1], where *nel is the number of variables to be eliminated
 * h[*nel]   = number of first variable not to be elim.
 * h[*nel+1] = number of second variable not to be eliminated,
 * etc.
 */
{
    int *h,i,j;
    
    if (!eo)
        msg ( "Need 'ELIMINATION_ORDER' to eliminate variables", "", 0 );
    
    h  = (int *) allo(CP 0,0,(dim+1)*sizeof(int));
    
    j = h[dim] = dim;
    for (i = dim-1; i >= 0; i--) 
        if (!(eo)[i])
            h[--j] = i;
    
    *nel = 0;
    do 
    {
        for (i = 0; i < dim; i++)
            if ((eo)[i] == *nel+1) 
            {
                h[(*nel)++] = i;
                (eo)[i] = 0;
                break;
            }
        if (i == dim) 
            break;
    } 
    while (1);
    
    for (i = 0; i < dim; i++)
        if ((eo)[i])
            msg ("Invalid format of 'ELIMINATION_ORDER' line", "", 0 );
    
    allo(CP eo, dim*sizeof(int),0);
    return(h);
    
}










int intcompare( int *i, int *j )
{
    return(*i - *j);
}







int *elim_in;

void polarformat( RAT *inieq, int *equa_in, int ineq_in, RAT *inner )
{
    
    int i,j,col,row,itr,pivcol,sysrow;
    RAT *o,*n,*pivot;
    
    sysrow = dim + 2;
    
    if (inner != 0) 
    {
        for (i = 0; i < *equa_in; i++) 
        {
            vecpr(inner,ar1+i*sysrow,var+3,dim); 
            (*RAT_sub)(*(ar1+(i+1)*sysrow-2),var[3],var+3);
            if (var[3].num)
                msg("%sinput point not valid (equalnum : %i)","",i+1);
        }
        for (i = 0; i < ineq_in; i++) 
        {
            vecpr(inner,inieq+i*sysrow,var+3,dim); 
            (*RAT_sub)(*(inieq+(i+1)*sysrow-2),var[3],var+3);
            if (var[3].num < 0)
                msg("%sinput point not valid (inequalnum : %i)","",i+1);
        }
        
    }
    
    if (*equa_in) 
    {  
        
        fprintf(prt,"Gauss - elimination (input - equations) :\n" );
        
        /* 17.01.1994: include logging on file porta.log */
        porta_log( "Gauss - elimination (input - equations) :\n" );
        
        for (i = 0; i < *equa_in+ineq_in; i++) 
        {
            allo_list(i,0,0);
            porta_list[i]->sys = ar1+i*(dim+2);
        }
        
        /* ELIMINATE INPUT - EQUALITIES */
        
        elim_in = (int *) allo(CP elim_in,1,U (dim+1)*sizeof(int));
        for (col = 0; col < dim; col++)
            *(elim_in+col) =  col ;
        
        for (itr = 0; itr < *equa_in; itr++) 
        {
            
            pivot = porta_list[itr]->sys;
            for(pivcol = 0; pivcol < dim && !pivot->num; pivcol++,pivot++);
            
            if (pivcol == dim) 
            {
                if (!(pivot->num)) 
                {
                    
                    (*RAT_row_prim)(porta_list[*equa_in-1]->sys,
                                    porta_list[itr]->sys, pivot,dim+1);
                    (*equa_in)--;
                    itr--;
                    continue;
                }
                else
                    msg( "input equality system has no solution", "", 0 );
            }

            (*RAT_row_prim)(porta_list[itr]->sys,porta_list[itr]->sys, 
                            pivot,dim+1);
            
            for (row = itr+1; row < *equa_in+ineq_in; row++)  
                if (row != itr) {
                gauss_calcnewrow(porta_list[itr]->sys,porta_list[row]->sys,
                                 pivcol,porta_list[row]->sys,0,sysrow-1);
                (*RAT_row_prim)(porta_list[row]->sys,porta_list[row]->sys,
                                porta_list[row]->sys+dim,dim+1);
            }
            
            fprintf(prt," elimination of variable %d\n", pivcol+1);

            /* 17.01.1994: include logging on file porta.log */
            porta_log( " elimination of variable %d\n", pivcol+1);
            
            for (col = 0; *(elim_in+col) != pivcol ; col++);
            for (; col < dim-itr-1 ; col++)
                *(elim_in+col) = *(elim_in+col+1);
            for (j = dim-itr-1; j < dim-1  ; j++)
                *(elim_in+j) = *(elim_in+j+1);
            *(elim_in+dim-1) = pivcol;
            
        }
        fprintf(prt,"\n" );

        /* 17.01.1994: include logging on file porta.log */
        porta_log( "\n" );
    } 
    
    if (inner) 
        for (i = 1; i <= ineq_in; i++) 
        {
            vecpr(inner,inieq+(i-1)*sysrow,var+3,dim);
            (*RAT_sub)(*(inieq+i*sysrow-2),var[3],inieq+i*sysrow-2);
        }
    
    if (*equa_in > 0) 
    {
        dim -= *equa_in;
        n = o = inieq;
        for (row = 0; row < ineq_in; row++) 
        {
            for (col = 0; col < dim; col++)
                (*RAT_assign)(n++,inieq+row*sysrow+*(elim_in +col));  
            (*RAT_assign)(n++,inieq+(row+1)*sysrow-2);  
            (*RAT_assign)(n++,inieq+(row+1)*sysrow-1);  
        }
        sysrow -= *equa_in;
    }
    
    if (!inner) 
    {
        for (i = 1; i <= ineq_in; i++)
            if ((inieq+i*sysrow-2)->num < 0) 
                msg( "no valid point known ", "", 0 );
        
    }
    
    
    /* SCALE LEFT SIDE TO ONE */
    
    for (o = n = inieq,i = 0; i < ineq_in; o += dim+2, n += dim+1, i++) 
        (*RAT_row_prim)(o,n,o+dim,dim+1);
    
}








void resubst( RAT *inner, int equa_in, int indx[] )
{
    int sysrow,ie,j;
    RAT *new_sys;

    
    new_sys = ar3 + (equa+ineq)*(dim+equa_in+1);
    while (new_sys > ar3+nel_ar3-1) 
    {
        reallocate(ineq, &new_sys);
        new_sys = ar3 + (equa+ineq)*(dim+equa_in+1);
    }
    
    if (equa) 
    {
        /* inequality dimension to dim */
        new_sys = ar3 + ineq*(dim+1);
        for (ie = ineq-1; ie >= 0; ie--, new_sys -= dim+1) 
            blow_up(new_sys,ie,indx,dim-equa,dim);
    } 
    
    if (equa_in) 
    {
        /* append equations on sys */
        new_sys = ar3 + (ineq+equa)*(dim+equa_in+1);
        for (ie = equa+ineq-1; ie >= 0; ie--, new_sys -= (equa_in+dim+1)) 
            blow_up(new_sys,ie,elim_in,dim,dim+equa_in);
    } 
    
    if (inner != 0) 
    {   
        for (ie = 0; ie < ineq; ie++)
            if ((porta_list[ie]->sys+dim+equa_in)->num)
                for (j = 0; j < dim+equa_in; j++)
                    (*RAT_add)(*(porta_list[ie]->sys+j),*(inner+j),
                               porta_list[ie]->sys+j);
    } 
    
    if (equa_in > 0) 
    { 
        
        fprintf(prt,"solving linear equality system ");
        fflush(prt);
                
        /* 17.01.1994: include logging on file porta.log */
        porta_log( "solving linear equality system ");
        fflush(logfile);
        
        sysrow = equa_in + dim + 2;
        for (ie = 0; ie < ineq+equa; ie++) 
        {
            if (ie%50 == 0) 
            {
                fprintf(prt,".");
                fflush(prt);
           
                /* 17.01.1994: include logging on file porta.log */
                porta_log( ".");
                fflush(logfile);
            }
            backwsubst(porta_list[ie]->sys,sysrow,equa_in);
        }
        
        dim += equa_in;
        fprintf(prt,"\n");

        /* 17.01.1994: include logging on file porta.log */
        porta_log( "\n");
    }
    
}








void blow_up( RAT *nptr, int ie, int *el, int eldim, int fdim )
/* blows up to get space for the eliminated components */
{
    int elp,not_elp,col,i;
    
    elp = fdim-1;
    not_elp = eldim-1;
    
    (*RAT_assign)(nptr,(porta_list[ie]->sys+eldim));
    nptr--;  
    for (col = fdim-1; col >= 0 ; col--) 
    {
        for (i = eldim; i < fdim; i++)
            if (el[i] == col)
                break;
        if (i < fdim) 
        {
            (*RAT_assign)(nptr,RAT_const);
            nptr--;
            if (elp > eldim) 
                elp--;
        }
        else  
        {
            (*RAT_assign)(nptr,porta_list[ie]->sys+not_elp);
            nptr--;
            not_elp--;
        }
    }
    porta_list[ie]->sys = nptr+1;
    
}







void backwsubst( RAT *ieptr, int sysrow, int equa_in )
/* computed the values of the components eliminated by the input equals */
{
    int j,col,i,ind;
    RAT *rptr;
    
    for (i = 0; i < equa_in; i++) 
    {
        rptr = ar1+i*sysrow;
        (*RAT_assign)(rptr+sysrow-1,RAT_const);
        for (j = 0; j < dim; j++) 
        { 
            col = *(elim_in+j);
            (*RAT_mul)(*(ieptr+col), *(rptr+col), var+3);
            (*RAT_add)(rptr[sysrow-1], var[3], rptr+sysrow-1);
        }
        if (ieptr[dim+equa_in].num)
            (*RAT_sub)(rptr[sysrow-2],rptr[sysrow-1], rptr+sysrow-1);
        else
            (rptr+sysrow-1)->num *= -1;
    }
    
    for (i = equa_in-1; i >= 0; i--) 
    {
        rptr = ar1+i*sysrow;
        (*RAT_assign)(var+2,RAT_const);
        for (j = i+1; j < equa_in; j++) 
        { 
            col = *(elim_in+dim+j);
            (*RAT_mul)( rptr[col],*(ieptr+col), var+3);
            (*RAT_add)( var[2], var[3], var+2);
        }
        ind = *(elim_in+dim+i);
        (*RAT_sub)(rptr[sysrow-1],var[2],ieptr+ind );
        (*RAT_mul)(*(ieptr+ind),*(rptr+ind),(ieptr+ind));
    }
    
}








void writemat( RAT *ptr, int row, int col )
{ 
    int i,j;
    
    for(i = 0; i<row; i++) 
    {
        for(j=0; j<col; j++,ptr++)
        {
            fprintf(prt,"%ld/%i ",ptr->num,ptr->den.i);
            
            /* 17.01.1994: include logging on file porta.log */
            porta_log( "%ld/%i ",ptr->num,ptr->den.i);
        }
        fprintf(prt,"\n");

        /* 17.01.1994: include logging on file porta.log */
        porta_log( "\n");
    }
    fprintf(prt,"\n");
    
    /* 17.01.1994: include logging on file porta.log */
    porta_log( "\n");
}


// This function computes the GCD of two integers
// code taken from rosettacode.org, added by J-D B on 14.4.2013.
long int longgcd(long int u, long int v)
{
	long int t;
	while (v) {
		t = u; 
        u = v; 
        v = t % v;
	}
    return labs(u);
};






//// OLD VERSION :
//int no_denom( int sysrow, int first, int last, int outmsg )
///*
 //* Make a fractional point integer by multiplying it with a positive number
 //* that is as small as possible.
 //*/
//{
    //int ret=1,scm,old_scm,i,j,ie,*denom;
    
    ////char *unptr;
    
    ////printf("Ok, WHAT ABOUT HERE1?\n");
    ////unptr = calloc(10,1);
    ////printf("Ok, WHAT ABOUT HERE2?\n");
    
    ////free(unptr);
    ////printf("Ok, WHAT ABOUT HERE3?\n");
    

    ////denom = (int *) allo(CP 0,0,U sysrow*sizeof(int));
    
    ////printf("address of denom = %p\n", ( void * ) denom);
    
    ////printf("DIDN'T EVEN REACH HERE!");

	////printf("right now...");
//////#if defined WIN32 || defined __CYGWIN32__
    ////free(denom);
//////#else // WIN32
//////    cfree(denom);
//////#endif // WIN32
	////printf("OK LAH");


    //denom = (int *) allo(CP 0,0,U sysrow*sizeof(int));


    //if(outmsg)
    //{ 
        //fprintf(prt,"transformation to integer values ");

        ///* 17.01.1994: include logging on file porta.log */
        //porta_log( "transformation to integer values ");
    //}
    
    //for (ie = first; ie < last; ie++) 
    //{
        
        //for (i = 0; i < sysrow; i++)
            //denom[i] = (porta_list[ie]->sys+i)->den.i;
        
        //qsort(CP denom,sysrow,sizeof(int),
              //(int(*)(const void*,const void*))intcompare);
        //old_scm = scm = denom[0];
        //for (i = 0,j = 0; i < sysrow; i++)
            //if (i > 0 && denom[i] != denom[i-1]) 
            //{
                //scm = old_scm*denom[i];
                //if (scm/denom[i] != old_scm) 
                //{
                    //ret = 0;
                    //break;
                //}
                //else
                    //old_scm = scm;
                //denom[j++] = denom[i];
            //}
        //if (i != sysrow)
            //continue;
        //for (i = 0; i < sysrow; i++)
            //denom[i] = scm/denom[i];
        //scm = scm/gcdrow(denom,j);
        
        //(porta_list[ie]->sys+sysrow-1)->num *= scm;
        //for (i = 0; i < sysrow-1; i++) 
        //{
            //(porta_list[ie]->sys+i)->num = 
                //(scm/(porta_list[ie]->sys+i)->den.i)
                //*(porta_list[ie]->sys+i)->num;
            //(porta_list[ie]->sys+i)->den.i = 1;
        //}
        
    //}
    
    //if(outmsg)
    //{
        //fprintf(prt,"\n");

        ///* 17.01.1994: include logging on file porta.log */
        //porta_log( "\n");
    //}   
    
////    printf("address of denom = %p\n", ( void * ) denom);
//#if defined WIN32 || defined __CYGWIN32__ // Comment to link with c++ rather than g++ by J-D B : cfree gives the following error :... *** glibc detected *** ./test: free(): invalid pointer: 0x00007fe897c5dfdc ***
    //free(denom);
//#else // WIN32
    //cfree(denom);
//#endif // WIN32

    //return(ret);
    
//}




// This function was modified by J-D Bancal on 4.4.2013 to make it able
// to remove more fractions.
int no_denom( int sysrow, int first, int last, int outmsg )
/*
 * Make a fractional point integer by multiplying it with a positive number
 * that is as small as possible.
 */
{
	fprintf(prt, "Arrived here with sysrow = %i, first = %i, last = %i, outmsg = %i\n", sysrow, first, last, outmsg);
    int ret=1,i,j,ie;//,*denom;
    int t;
    
//    denom = (int *) allo(CP 0,0,U sysrow*sizeof(int));
    typedef long int intDenom;
    // Here we use long int for the denominator, because product of denominators can become very big eventually...
    intDenom scm,old_scm,*denom, gicidi,scmtmp;
    denom = (intDenom *) allo(CP 0,0,U sysrow*sizeof(intDenom));
    
    if(outmsg)
    { 
        fprintf(prt,"transformation to integer values ");

        /* 17.01.1994: include logging on file porta.log */
        porta_log( "transformation to integer values ");
    }
    
    for (ie = first; ie < last; ie++) 
    {
        
        for (i = 0; i < sysrow; i++)
            denom[i] = (porta_list[ie]->sys+i)->den.i;
        
        qsort(CP denom,sysrow,sizeof(intDenom),
              (intDenom(*)(const void*,const void*))intcompare);
        old_scm = scm = denom[0];
        
        for (i = 0,j = 0; i < sysrow; i++)
            if (i > 0 && denom[i] != denom[i-1]) 
            {
                // Here we compute the LCM (least common multiple) of the two different numbers, using the fact that lcm(a,b)=|a*b|/gcd(a,b)
                // First we compute the GCD of the two numbers
                //gicidi = lgcd(old_scm, denom[i]);
                /*intDenom tt, u, v;
                u = old_scm;
                v = denom[i];
                while (v) {
                    tt = u; 
                    u = v; 
                    v = tt % v;
                }
                gicidi = labs(u);*/
                gicidi = longgcd(old_scm, denom[i]);
				// Now we continue with the LCM
				scmtmp = labs(old_scm)/gicidi;
                scm = scmtmp*labs(denom[i]);
                if (scm/labs(denom[i]) != scmtmp) // To check whether we multiplied too big numbers...
                {
					fprintf(prt, "int version of the error:\n");
					fprintf(prt, "Error for i=%i, j=%i\n", i, j);
					fprintf(prt, "old_scm = %i, scm = %i, denom[i] = %i, ratio = %i, scmtmp = %li\n", old_scm, scm, denom[i], scm/denom[i], scmtmp);
					fprintf(prt, "gicidi = %i\n, ", gicidi);

					int tmpi;
					for (tmpi = 0; tmpi < i; tmpi++)
						fprintf(prt, "denom[%i] = %i\n", tmpi, denom[tmpi]);

					fprintf(prt, "long int version of the error:\n");
					fprintf(prt, "Error for i=%li, j=%li\n", i, j);
					fprintf(prt, "old_scm = %li, scm = %li, denom[i] = %li, ratio = %li, scmtmp = %li\n", old_scm, scm, denom[i], scm/denom[i], scmtmp);
					fprintf(prt, "gicidi = %li\n, ", gicidi);

					for (tmpi = 0; tmpi < i; tmpi++)
					{
						fprintf(prt, "num[%li] = %li\n", tmpi, (porta_list[ie]->sys+i)->num);
						fprintf(prt, "denom[%li] = %li\n", tmpi, denom[tmpi]);
					}

					fprintf(prt, "The full line is : ");
					for (tmpi = 0; tmpi < sysrow; tmpi++)
					{
						fprintf(prt, "%li/%li ", (porta_list[ie]->sys+i)->num, denom[tmpi]);
					}
					fprintf(prt, "\n\n");


                    ret = 0;
                    break;
                }
                else
                    old_scm = scm;
                denom[j++] = denom[i];
            }
        if (i != sysrow)
            continue;
		for (i = 0; i < sysrow; i++)
            denom[i] = scm/denom[i];
        scm = scm/longgcdrow(denom,j);
        
        //printf("%d\n\n",scm);

        (porta_list[ie]->sys+sysrow-1)->num *= scm; // Careful that the produced number here could be larger than the largest representable integer... (This is not checked for at the moment...)
        for (i = 0; i < sysrow-1; i++) 
        {
            (porta_list[ie]->sys+i)->num = 
                (scm/(porta_list[ie]->sys+i)->den.i)
                *(porta_list[ie]->sys+i)->num;
            (porta_list[ie]->sys+i)->den.i = 1;
        }
        
    }
    
    if(outmsg)
    {
        fprintf(prt,"\n");

        /* 17.01.1994: include logging on file porta.log */
        porta_log( "\n");
    }   
    
    free(denom);

    return(ret);
    
}







void origin_add( int rowl_inar, RAT *inieq )
{
    int i,j;
    unsigned *m;
    listp lp;
    
    m = 0;
    for (i = 0; i < ineq; i++) 
    {
        for (j = 0; j < dim-equa; j++)
            if (( (porta_list[i]->sys)+j )->num)
                break;
        if (j == dim-equa && ( (porta_list[i]->sys)+j )->num) 
        {     
            m = porta_list[i]->mark;
            break;
        }
    }
    
    if (!m && cone >= dim-equa) 
    {
        m = (unsigned *) 1;
        allo_list(ineq+equa,&m,points/32+2);
        
        lp = porta_list[ineq+equa];
        for (j = ineq+equa; j != ineq; j--)
            porta_list[j] = porta_list[j-1];
        porta_list[ineq] = lp;
        /* 
         * Change by M.S. 5.6.92:
         * A bus error occurred in the following line with the .ieq file:
         * DIM=1
         * INEQUALITIES_SECTION
         * x1=0
         * with or without VALID points.
         * ineq is 0 in this case.
         if (porta_list[ineq-1]->sys+2*(dim+1-equa) > ar3+nel_ar3-1)
         reallocate(ineq, &(RAT *) i);
         
         porta_list[ineq]->sys = porta_list[ineq-1]->sys+dim+1-equa;
         */
        if (ineq) 
        {
            if (porta_list[ineq-1]->sys+2*(dim+1-equa) > ar3+nel_ar3-1)
                reallocate(ineq, (RAT **)&i);
            porta_list[ineq]->sys = porta_list[ineq-1]->sys+dim+1-equa;
        }
        else 
        {
            if (!ar3) 
            {
                nel_ar3 = FIRST_SYS_EL+ dim+1-equa;
                ar3 = (RAT *) RATallo(CP ar3,0,U nel_ar3);
            }
            else if (ar3+dim+1-equa > ar3+nel_ar3-1)
                reallocate(ineq, (RAT**)&i);
            porta_list[ineq]->sys = ar3;
        }

        for (j = 0; j < dim-equa; j++) 
            (*RAT_assign)(porta_list[ineq]->sys+j,RAT_const);
        (*RAT_assign)(porta_list[ineq]->sys+j,RAT_const+1);
        ineq++;
    }
    
    if (m) 
    {
        for (i=0; i< points/32+2; i++)
            m[i] = 0;
        for (i=1; i <= points; i++)
            if (!inieq[i*rowl_inar-1].num)
                domark(m,i-1);
    }
    
    
} 









void gentableau( RAT *ar1p, int poi_to_ieq, int *rowl_inar, int **indx )
/*****************************************************************/
/*
 * Copy the (points * dim+1) matrix ar1 (where the rows are the points) 
 * to the (dim+1 * points+dim+1) matrix ar2 
 * transposed (so that the points are columns), 
 * append a negative identity matrix of dimension "dim",
 * and a right-hand side vector of (0,1) of size: (dim,1).
 * The last row of ar2 corresponds to an equation bx=1,
 * where b_p = 1 of the corresponding point is a conv-point,
 * and b_p = 0, if it is a cone_point.
 * ar2 actually contains space for dim+1 more rows 
 * for temporary usage in gauss.
 * Compute cone and conv, set indx[] and porta_list[]->sys.
 */
{
    int sysrow,i,j;
    RAT *r_sp;

    if (poi_to_ieq == 0) 
    {
        /* KH(0 (!) ,....) */
        /* Append a conv-point 0 to array ar1 */
        for (r_sp = ar1p+(dim+1)*points, i = 0; i < dim; r_sp++, i++) 
            (*RAT_assign)(r_sp,RAT_const);
        (*RAT_assign)(r_sp,RAT_const+1);
        points++;
    }
    
    sysrow = dim+points+1;
    *rowl_inar = dim+1;

    nel_ar2 = 2*(dim+1)*sysrow;
    ar2 = (RAT *) RATallo(CP ar2,1, nel_ar2);
    
    for (j = 1; j <= points; j++) 
    {
        for (i = 1; i<= dim; i++)
            ar2[(i-1)*sysrow+j-1] = ar1p[(j-1)*(dim+1)+i-1];
        ar2[dim*sysrow+j-1] = ar1p[(j-1)*(dim+1)+dim];
        if (ar1p[(j-1)*(dim+1)+dim].num) conv++; else cone++;
    }
    for (j = points+1; j <= sysrow; j++) 
        for (i = 1; i<= dim+1; i++) 
        {
            ar2[(i-1)*sysrow+j-1].num = (j-points == i) ? -1 : 0;
            ar2[(i-1)*sysrow+j-1].den.i = 1;
        }
    ar2[sysrow*(dim+1)-1].num = 1;
    ar2[sysrow*(dim+1)-1].den.i = 1;
    
    for (r_sp = ar2,i = 0; i <= 2*dim+1; r_sp += sysrow,i++) 
    {
        allo_list(i,0,0);
        porta_list[i]->sys =  r_sp;
    }
    
    *indx = (int *) allo(CP (*indx),0,U (points+dim+1)*sizeof(int));
    for (i = 0; i < points; i++)
        (*indx)[i] =  -i-1 ;
    for (i = 0; i < dim; i++)
        (*indx)[points+i] =  i ;
    (*indx)[points+dim] = 0;
    
    equa=0;
    ineq=0;
}


/*****************************************************************/

void reorder_var( int ineq, RAT *ar1, RAT **ar2, int *nel_ar2, int *nel,
                  int **elim_ord, int **indx )
/*****************************************************************/
/*
 * Copy the (ineq * dim+1) matrix ar1 into the (2ineq * dim+1) matrix ar2,
 * so that the first variable to be eliminated is the first column of ar2,
 * the second variable is the second column, etc.
 * The extra "ineq" rows in ar2 are used for temporary storage in gauss.
 * Upon input, elim_ord[j-1] = i > 0 means that variable j
 * is the i-th variable to be eliminated.
 * elim_ord[j-1] = 0 means that variable j is not eliminated.
 * Upon output, elim_ord[0,1,...,dim] is set to 0,1,...,dim.
 * indx[] gives the original names of the variables (minus 1),
 *    which are made negative for the variables to be eliminated
 *    and nonnegative for the others.
 *    The negativity is used in gauss().
 * nel is the number of variables to be eliminated.
 */
{
    int sysrow,i,j,col,max;
    RAT *r1ptr, *r2ptr;
    
    if (!(*elim_ord))
        msg( "Need 'ELIMINATION_ORDER' to eliminate variables", "", 0 );
    
    sysrow = dim+1;
    
    *nel_ar2 = 2*ineq*sysrow;
    *ar2 = (RAT *) RATallo( CP (*ar2), 0, *nel_ar2);
    
    /* copy columns belonging to variables that will be eliminated */
    
    *nel = 0;   /* number of variables to be eliminated          */
    max = 0;    /* highest number occurring in the elimination order */
    for (i = 0; i < dim; i++) 
    {
        col = (*elim_ord)[i];
        if (!col) 
            continue;
        if (col > max) 
            max = col;
        if (col < 0) 
            msg( "Invalid format of 'ELIMINATION_ORDER' line", "", 0 );
        col = col-1;
        r1ptr = ar1+i;
        r2ptr = *ar2+col;
        for (j = 0; j < ineq; j++ ) 
        {
            (*RAT_assign)(r2ptr,r1ptr);
            r1ptr += dim + 2;
            r2ptr += dim + 1;
        }
        (*nel)++;
    }
    if (max != *nel) 
        msg( "Invalid format of 'ELIMINATION_ORDER' line", "", 0 );
    
    /* copy columns belonging to variables that will not be eliminated */
    col = *nel;     /* number of variables to eliminate */
    for (i = 0; i < dim; i++) 
    {
        if ((*elim_ord)[i]) 
            continue;
        r1ptr = ar1+i;
        r2ptr = *ar2+col;
        for (j = 0; j < ineq; j++ ) 
        {
            (*RAT_assign)(r2ptr,r1ptr);
            r1ptr += dim + 2;
            r2ptr += dim + 1;
        }
        col++;
    }
    
    /* Copy the right-hand side vector */
    r1ptr = ar1+dim;
    r2ptr = *ar2+dim;
    for (j = 0; j < ineq; j++ ) 
    {
        (*RAT_assign)(r2ptr,r1ptr);
        r1ptr += dim + 2;
        r2ptr += dim + 1;
    }
    
    /* Assign porta_list[]->sys */
    r2ptr = *ar2;
    for (i = 0; i <= 2*ineq; i++) 
    {
        allo_list(i,0,0);
        porta_list[i]->sys = r2ptr;
        r2ptr += sysrow;
    }

    /* indx */
    *indx = (int *) allo(CP (*indx),1,U (dim+1)*sizeof(int));
    for (i = 0; i < dim+1; i++) 
        (*indx)[i] = 0;
    j = *nel;
    for (i = 0; i < dim; i++) 
    {
        col = (*elim_ord)[i];
        if (col) 
            (*indx)[col-1] = -i-1;
        else 
            (*indx)[j++] = i;
    }
    (*indx)[dim] = 0;
    /* no "holes" in indx for the first "nel" elements? */
    for (i = 0; i < *nel; i++)
        if ((*indx)[i] == 0)
            msg( "Invalid format of 'ELIMINATION_ORDER' line", "", 0 );
        
    /* elim_ord */
    *elim_ord = (int *) allo( CP (*elim_ord),U (dim*sizeof(int)),
                             U ((dim+1)*sizeof(int)));
    for (i = 0; i < dim+1; i++)
        (*elim_ord)[i] = i;
}










void reallocate( int ie, RAT **sysptr )
/*
 * enlarge array "ar3",
 * and make the pointers porta_list[]->sys, porta_list[]->ptr, and sysptr
 * point to their respective new positions.
 * "ie" is the current number of inequalities in Fourier_Motzkin,
 * "sysptr" contains a pointer to some place in "ar3".
 */
{
    int i,*save_sys,*save_ptr,save_sp;
    
    save_sys = (int *) allo(CP 0,0,U ie*sizeof(int));
    save_ptr = (int *) allo(CP 0,0,U ie*sizeof(int));
    
    save_sp = (*sysptr == 0) ? -1 : (*sysptr - ar3);
    for (i = 0; i < ie; i++) 
    {
        if( porta_list[i]->sys == (RAT *)0 )
            save_sys[i] = -1;
        else
            save_sys[i] = porta_list[i]->sys - ar3;
        
        if( porta_list[i]->ptr == (RAT *)0 )
            save_ptr[i] = -1;
        else
            save_ptr[i] = porta_list[i]->ptr - ar3;
    }
    
    nel_ar3 += INCR_SYS_EL;
    ar3 = (RAT *) RATallo(CP ar3,nel_ar3-INCR_SYS_EL,U nel_ar3);
    /* fprintf(prt,"New space allocated \n"); */
    
    /* 
     * "ar3" has maybe changed its value, 
     * therefore sysptr has to be recomputed as ar3[save_sp],
     * and porta_list[i]->sys and porta_list[i]->ptr also have to be recomputed,
     * so they have the same offset in ar3 as before.
     */
    *sysptr = (save_sp == -1) ? 0 : (ar3 + save_sp);
    for (i = 0; i < ie; i++) 
    {
        porta_list[i]->sys =  (save_sys[i] == -1) ? 0 : (ar3 + save_sys[i]);
        porta_list[i]->ptr =  (save_ptr[i] == -1) ? 0 : (ar3 + save_ptr[i]);
    }
    
    /*
     * Change by M.S. 3.6.1992:
     * Free save_sys and save_ptr by using allo(), not cfree().
     *
     cfree(save_sys);
     cfree(save_ptr);
     */
    allo(CP save_sys,U ie*sizeof(int), 0);
    allo(CP save_ptr,U ie*sizeof(int), 0);
}


#if !defined WIN32
#include <unistd.h>
#include <sys/times.h>
#include <sys/time.h>
#endif
#include <time.h>







float time_used()
{
#if defined WIN32

    time_t now;
    time(&now);
    return (float) now;

#else // WIN32

    struct tms now;
    int sec, hun;
    float total;

#ifndef CLK_TCK
    int clocks_per_second = 0;
    
    clocks_per_second = sysconf( _SC_CLK_TCK );
    times( &now );
    
    hun  = ( ( now.tms_utime % clocks_per_second ) * 100 ) / clocks_per_second;
    sec = ( now.tms_utime / clocks_per_second );
#else
    times( &now );

    hun  = ( ( now.tms_utime % CLK_TCK ) * 100 ) / CLK_TCK;
    sec = ( now.tms_utime / CLK_TCK );
#endif

    total = (float)sec + ( (float)(hun) / 100.0 );
    
    return( total );

#endif // WIN32
}

     





static double initial_time;


void init_total_time()
{
#if defined WIN32

    time_t now;
    time(&now);

    initial_time = (double)now;

#else // WIN32

    struct timeval tp;
    struct timezone tzp;
    
    gettimeofday( &tp, &tzp );
    initial_time = (double)(tp.tv_sec) 
        + 0.000001 * (double)(tp.tv_usec);

#endif // WIN32
}








float total_time()
{
#if defined WIN32

    time_t now;
    time(&now);

    return (float)(now - initial_time);

#else // WIN32

    struct timeval tp;
    struct timezone tzp;
    
    gettimeofday( &tp, &tzp );
    return (double)(tp.tv_sec) 
        + 0.000001 * (double)(tp.tv_usec)
        - initial_time;

#endif // WIN32
}






///* The following added by J-D Bancal on 22.2.2012 to call porta as a c library.
	//The data variable is expected to be of the standard *.ext format, i.e.
	
	//[1 -1
	//1  2]   ,  which is given as data=[1 -1 1 2] (i.e. line by line)
	
	//for a 1-dimensional polytope with extreme points x=-1, and x=+2 (dimension = 1, nbPoints = 2).

//*/
//int quick_and_dirty_poi_conv_call(long int *data, int dimension, int nbPoints)
//{
    //int ieq_file, start;
    //char outfname[20];
    //char fname[20];
    //int   poi_file;
    //int   rowl_inar, ierl;
    //int  *indx = (int *)0;      
    //int equa_in,ineq_in, ineq_out;
    //FILE *outfp;

    //int i,j; // Added by J-D B, used as summation indices later.

    

    //printf("\nPORTA - a POlyhedron Representation Transformation Algorithm\n");
    //printf(  "Version %s\n\n", VERSION );

    //printf( "Written by Thomas Christof (Uni Heidelberg)\n" );
    //printf( "Revised by Andreas Loebel (ZIB Berlin)\n\n" );

    //printf( "PORTA is free software and comes with ABSOLUTELY NO WARRENTY! You are welcome\n" );
    //printf( "to use, modify, and redistribute it under the GNU General Public Lincese.\n\n" ); 
    
    //printf( "This is the program XPORTA from the PORTA package.\n\n" );

///*    if( argc <= 2 )
    //{
        //printf( "For more information read the manpages about porta.\n\n" );
        //exit(-1);
    //}*/
            
    ///* 17.01.1994: include logging on file porta.log */
    //logfile = fopen( "porta.log", "a" );
    //if( !logfile )
        //fprintf( stderr, "can't open logfile porta.log\n" );
    //else
    //{
///*        porta_log( "\n\n\nlog for " );
        //for( i = 0; i < argc; i++ )
            //porta_log( "%s ", argv[i] );
        //porta_log( "\n\n" );*/
        //porta_log("\n\n\nlog for quick_and_dirty_poi_conv_call\n\n");
    //}
            
    //init_total_time();
    
    //initialize();

    //prt = stdout;

    ///*get_options(&argc,&argv);*/
    ///*We bypass the previous line which gives errors, and just use the following,
     //* which corresponds to having the option "-T" :*/
    //option |= Traf;
                //allowed_options = Traf|
                    //Chernikov_rule_off|Validity_table_out|
                    //Redundance_check|Statistic_of_coefficients|
                    //Protocol_to_file|Opt_elim|Long_arithmetic;
    
    //if (option & Protocol_to_file) 
    //{
        //printf("\n\n\n\n Reached here!!!!! (Should not).\n\n\n\n");
///*        strcat(*argv,".prt");
        //prt = fopen(*argv,"w");
        //(*argv)[strlen(*argv)-4] = '\0';*/
    //}
    //setbuf(prt,CP 0);
    
    //set_I_functions();
    //SET_MP_not_ready;
    //ieq_file = 0; /*!strcmp(*argv+strlen(*argv)-4,".ieq");*/
    //poi_file = 1; /*!strcmp(*argv+strlen(*argv)-4,".poi");*/
    
    //if (!poi_file && !ieq_file)
        //msg( "invalid format of command line", "", 0 );
    
    ///*
     //* change by M.S. 5.6.92:
     //* read_input_file writes to the output file, if is_set(Sort).
     //*/
    //outfp = 0;
    ///*strcpy(outfname,*argv);
    //if (is_set(Sort) && poi_file) 
    //{
        //strcat(outfname,".poi");
        ///*outfp = wfopen(outfname);*/ /* Actually, we don't want to open any file since the data has been passedpassed to us in the argument. */
    ///*}
    //if (is_set(Sort) && ieq_file) 
    //{
        //strcat(outfname,".ieq");
        //fprintf(prt,"outfname = %s\n",outfname);
        //fflush(stdout);

        ///* 17.01.1994: include logging on file porta.log */
        ///*porta_log( "outfname = %s\n",outfname );
        //fflush(logfile);

        ///*outfp = wfopen(outfname);*/ /* Actually, we don't want to open any file since the data has been passed to us in the argument. */
    ////}

    //if (is_set(Fmel) && ieq_file) 
    //{  
        ///* ONLY FM-ELIMINATON */
        
        //int *elim_ord,nel;
            //char *cp1;
        //elim_ord = 0;
        //cp1=strdup("ELIMINATION_ORDER");
        //if(!cp1)
            //msg( "allocation of new space failed", "", 0 );
            
        ///*ineq = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,cp1,
                               //&elim_ord,"\0",(int **)&i,"\0",(RAT **)&i);*/
        ///* We need to replace the above line with something more adapted.... */
        //ineq = 0;
        
        
        //free(cp1);
        //sort_eqie_cvce(ar1,ineq,dim+2,&equa_in,&ineq_in);
        //ineq = ineq_in+equa_in;
        ///*     elim_ord = check_and_reorder_elim_ord(elim_ord,&nel); */
        //reorder_var(ineq,ar1,&ar2,(int *)&nel_ar2,&nel,&elim_ord,&indx);
        ///*     indx = elim_ord+nel; */
        ///* 
         //* Transform all inequalities into "<= 1" or "<=-1" inequalities,
         //* if possible.
         //* If the right-hand side is 0, divide the numerators by their gcd
         //* and the denominators by their gcd.
         //* (This is not really necessary).
         //*/
        ///*
           //for (i = 0; i < ineq; i++) 
           //(* RAT_row_prim)(ar2+(dim+2)*i,ar2+(dim+1)*i,ar2+(dim+2)*i+dim,dim+1);
        //*/
        //if(is_set(Long_arithmetic))
        //{
            //RAT a, b;
            //SET_MP_ready;
            //memset( &a, 0, sizeof(a) );
            //memset( &b, 0, sizeof(b) );
            //arith_overflow_func(0,0,a,b,0);
        //}
        //for (i = 0; i < ineq; i++) 
            //(* RAT_row_prim)(ar2+(dim+1)*i,ar2+(dim+1)*i,
                             //ar2+(dim+1)*i+dim,dim+1);
        ///*     nel_ar2 = ineq*(dim+1);  
         //*/
        //equa = 0;
        //ineq_out = ineq;
        //gauss(0, dim+1, dim+1, dim-nel, equa_in, &ineq_out, &equa, indx);
        //for (; (*indx) < 0; indx++);
        //ierl = dim-nel-equa +1;
        ///* row-length of inequalities after fourier_motzkin
         //* =  number of noneliminated variables+1 */
        //nel = nel - (ineq - ineq_out);  /* number of variables to be elim. */
        //ineq = ineq_out;
        //fourier_motzkin(0,ineq-equa,dim+1-equa_in,nel,poi_file,indx,elim_ord);
        //if ((MP_realised && return_from_mp()) || !MP_realised) 
            //sort(no_denom(ierl, 0, ineq,1), ierl, 0, ineq);
        ///*write_ieq_file(*argv,outfp,equa,ineq,dim+1,0,
                       //ineq,0,ierl,indx);*/
        ///* Here we still need to return the output inequalities computed ... */
        
    //}
    ///*else if (is_set(Sort)) 
    //{
        ///*points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,"\0",(int **)&i,
                                 //"\0",(int **)&i,"\0",(RAT **)&i);*/
        ///* Line above to be replaced ... */
        ///*if (ieq_file)
            //sort_eqie_cvce(ar1,points,dim+2,&equa,&ineq);
        //listptoar(ar1,points,ieq_file?dim+2:dim+1,0); 
        //if (ieq_file) 
        //{ 
            //if (equa) sort(1,dim+1,0,equa);
            //if (ineq) sort(1,dim+1,equa,points); 
            //write_ieq_file(*argv,outfp,equa,0,dim+1,0,ineq,equa,dim+1,0);
        //}
        //else 
        //{
            //sort(1,dim+1,0,points);
            //for (cone = 0; !(porta_list[cone]->sys[dim].num); cone++);
            //write_poi_file(*argv,outfp,dim,0,0,cone,0,points-cone,cone);
        //}
    //}  */
    //else if ((is_set(Traf) || is_set(Dim)) && poi_file) 
    //{
///*        points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
                                 //"\0",(int **)&i,"\0",(int **)&i,"\0",
                                 //(RAT **)&i);*/
                                 
        ///* Here we replace the previous line by the variables given in argument 
           //section appended by J-D Bancal, 22.2.2012. */
        //points = nbPoints;
        //dim = dimension;

        ///* We allocate the memory needed and initialize the table */
        //ar1 = (RAT *) RATallo(ar1, 0, (points+1)*(dim+1)); // Note that for some reason porta seems to allocate memory for one more extremal point... so we also do that here.
        //for (i = 0; i < points; i++)
        //{
            //for (j = 0; j < dim+1; j++)
            //{
            	    //if (j != dim) /* We put the constant term at the end of the column. */
            	    //{
////            	    	printf("Element (i,j)=(%i,%i) vaut %li\n", i, j, data[i*(dim+1)+j+1]);
            	        //ar1[i*(dim+1)+j].num = data[i*(dim+1)+j+1];
            	    //}
            	    //else 
            	    //{
////            	    	printf("Element (i,j)=(%i,%i) vaut %li\n", i, j, data[i*(dim+1)]);
            	        //ar1[i*(dim+1)+j].num = data[i*(dim+1)];
            	    //}
            //}
        //}
        ///* We also prepare the valid point */
        //ar6 = (RAT *) RATallo(ar6, 0, (points+1)*dim);
        //for (j = 0; j < dim; j++)
   	        //ar6[j] = ar1[j];

        
		///* Here we examine what the previous function produced. */   
		///*printf("points = %i\n", points);
		//printf("address of ar1 = %p\n", ( void * ) ar1);
		//printf("first few values of ar1 = ...\n");
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[0].num, ar1[0].den.i, ar1[1].num, ar1[1].den.i, ar1[2].num, ar1[2].den.i, ar1[3].num, ar1[3].den.i, ar1[4].num, ar1[4].den.i, ar1[5].num, ar1[5].den.i, ar1[6].num, ar1[6].den.i, ar1[7].num, ar1[7].den.i, ar1[8].num, ar1[8].den.i, ar1[9].num, ar1[9].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[10].num, ar1[10].den.i, ar1[11].num, ar1[11].den.i, ar1[12].num, ar1[12].den.i, ar1[13].num, ar1[13].den.i, ar1[14].num, ar1[14].den.i, ar1[15].num, ar1[15].den.i, ar1[16].num, ar1[16].den.i, ar1[17].num, ar1[17].den.i, ar1[18].num, ar1[18].den.i, ar1[19].num, ar1[19].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[20].num, ar1[20].den.i, ar1[21].num, ar1[21].den.i, ar1[22].num, ar1[22].den.i, ar1[23].num, ar1[23].den.i, ar1[24].num, ar1[24].den.i, ar1[25].num, ar1[25].den.i, ar1[26].num, ar1[26].den.i, ar1[27].num, ar1[27].den.i, ar1[28].num, ar1[28].den.i, ar1[29].num, ar1[29].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[30].num, ar1[30].den.i, ar1[31].num, ar1[31].den.i, ar1[32].num, ar1[32].den.i, ar1[33].num, ar1[33].den.i, ar1[34].num, ar1[34].den.i, ar1[35].num, ar1[35].den.i, ar1[36].num, ar1[36].den.i, ar1[37].num, ar1[37].den.i, ar1[38].num, ar1[38].den.i, ar1[39].num, ar1[39].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[40].num, ar1[40].den.i, ar1[41].num, ar1[41].den.i, ar1[42].num, ar1[42].den.i, ar1[43].num, ar1[43].den.i, ar1[44].num, ar1[44].den.i, ar1[45].num, ar1[45].den.i, ar1[46].num, ar1[46].den.i, ar1[47].num, ar1[47].den.i, ar1[48].num, ar1[48].den.i, ar1[49].num, ar1[49].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[50].num, ar1[50].den.i, ar1[51].num, ar1[51].den.i, ar1[52].num, ar1[52].den.i, ar1[53].num, ar1[53].den.i, ar1[54].num, ar1[54].den.i, ar1[55].num, ar1[55].den.i, ar1[56].num, ar1[56].den.i, ar1[57].num, ar1[57].den.i, ar1[58].num, ar1[58].den.i, ar1[59].num, ar1[59].den.i);
        
		//printf("element number 32 : %li/%i\n", ar1[32].num, ar1[32].den.i);
		//printf("element number 56 : %li/%i\n", ar1[56].num, ar1[56].den.i);
		
		//printf("dim = %i\n", dim);
		//printf("nel_ar1 = %li\n", nel_ar1);
		
        //for (j = 0; j < dim+1; j++)
        				//printf("ar6[%i] is %li/%i\n",j,ar6[j].num,ar6[j].den.i);*/

        //gentableau(ar1,1,&rowl_inar,&indx); 
        //if(is_set(Long_arithmetic))
        //{
            //RAT a, b;
            //SET_MP_ready;
            //memset( &a, 0, sizeof(a) );
            //memset( &b, 0, sizeof(b) );
            //arith_overflow_func(0,0,a,b,0);
        //}
        //ineq = (cone == points) ? dim : dim + 1;
        //ineq_out = ineq;  /*not used further */
        //gauss(1, points+dim+1,dim+1,dim,ineq,&ineq_out, &equa, indx);
        ///* make indx point to the system-variable section */
        //for (; (*indx) < 0; indx++);
        //if (is_set(Dim)) 
        //{
            //char str[100];
            //fprintf (prt,"\nDIMENSION OF THE POLYHEDRON : %i\n\n",dim-equa);
            
            ///* 17.01.1994: include logging on file porta.log */
            //porta_log( "\nDIMENSION OF THE POLYHEDRON : %i\n\n", dim-equa );
            
			///*sprintf (str,"echo 'DIMENSION OF THE POLYHEDRON : %i' | cat >> %s",
                     //dim-equa,argv[0]);*/
            //sprintf (str,"echo 'DIMENSION OF THE POLYHEDRON : %i'",
                     //dim-equa);
            //system(str);
            //if (equa) 
            //{
                //fprintf(prt,"equations :\n");

                ///* 17.01.1994: include logging on file porta.log */
                //porta_log( "equations :\n");
                
                //listptoar(ar4,equa,dim+1,0); 
                //if ((MP_realised && return_from_mp()) || !MP_realised) 
                    //sort(no_denom(dim+1,0,equa,1), dim+1,0,equa);
                //start = 1;
                ///* last argument of writesys was lost? 
                   //writesys(prt,0,equa,dim+1,0,0,'=');
                   //*/
                //writesys(prt,0,equa,dim+1,0,0,'=', &start);

                ///* log equation system */
                //start = 1;
                //writesys(logfile,0,equa,dim+1,0,0,'=', &start);
            //}
        //}
        //else  
        //{
            ///* POINTS TO INEQUALITIES */
            ///*sprintf(fname,"%s.ieq",*argv);*/
            //fourier_motzkin(fname,ineq-equa,points+dim+1-ineq,
                            //points-ineq+equa,poi_file,indx,0);
            //if (is_set(Validity_table_out)) 
                //red_test(indx,ar1,&rowl_inar);
            //if ((MP_realised && return_from_mp()) || !MP_realised) 
            //{
                //if (equa) sort(no_denom(dim+1, ineq, ineq+equa,1), 
                               //dim+1, ineq, ineq+equa);
                //sort(no_denom(dim+1-equa, 0, ineq,1), 
                     //dim+1-equa, 0, ineq);
            //}
            ///*write_ieq_file("output",outfp,equa,ineq,
                           //dim+1,0,ineq,0,dim+1-equa,indx);*/
            
            //// Here we write the output in the lrs format
			//printf("\n* Result of the conversion:\n");
			//printf("H-representation\n");
			//if (equa > 0)
			//{
				//printf("linearity %i", equa);
				//for (i = 0; i < equa; i++)
					//printf(" %i", i+1);
				//printf("\n");
			//};
			//printf("begin\n");
			//printf("%i %i rational\n", ineq+equa, dim+1);
			//for (i = 0; i < equa; i++)
			//{
				//printf("%li", (porta_list[ineq+i]->sys+dim)->num);
				//if (((porta_list[ineq+i]->sys+dim)->den.i) != 1)
					//printf("/%i ", (porta_list[ineq+dim]->sys+j)->den.i);
				//else
					//printf(" ");
				//for (j = 0; j < dim; j++)
				//{
					//printf("%li", (porta_list[ineq+i]->sys+j)->num);
					//if (((porta_list[ineq+i]->sys+j)->den.i) != 1)
						//printf("/%i ", (porta_list[ineq+i]->sys+j)->den.i);
					//else
						//printf(" ");
				//}
				//printf("\n");
			//};
			
			//// We prepare the inverse index vector to know for each variable
			//// whether it is part of the inequality section or not.
			//int indices[dim];
			//for (j = 0; j < dim; j++)
				//indices[j] = -1;
			//for (j = 0; j < dim-equa; j++)
				//indices[*(indx+j)] = j;
			
			//for (i = 0; i < ineq; i++)
			//{
				//printf("%li", (porta_list[i]->sys+dim-equa)->num);
				//if (((porta_list[i]->sys+dim-equa)->den.i) != 1)
					//printf("/%i ", (porta_list[i]->sys+dim-equa)->den.i);
				//else
					//printf(" ");
				//for (j = 0; j < dim; j++)
				//{
					//if (indices[j] == -1) // In this case the current variable does not appear in the inequality description
						//printf("0 ");
					//else
					//{
						//printf("%li", (porta_list[i]->sys+indices[j])->num);
						//if (((porta_list[i]->sys+indices[j])->den.i) != 1)
							//printf("/%i ", (porta_list[i]->sys+indices[j])->den.i);
						//else
							//printf(" ");
					//}
				//}
				//printf("\n");
			//};
			//printf("end\n\n");

			
			
///*			printf("first few values of porta_list = ...\n");
			//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", (porta_list[0]->sys+0)->num, (porta_list[0]->sys+0)->den.i, (porta_list[0]->sys+1)->num, (porta_list[0]->sys+1)->den.i, (porta_list[0]->sys+2)->num, (porta_list[0]->sys+2)->den.i, (porta_list[0]->sys+3)->num, (porta_list[0]->sys+3)->den.i, (porta_list[0]->sys+4)->num, (porta_list[0]->sys+4)->den.i, (porta_list[0]->sys+5)->num, (porta_list[0]->sys+5)->den.i, (porta_list[0]->sys+6)->num, (porta_list[0]->sys+6)->den.i, (porta_list[0]->sys+7)->num, (porta_list[0]->sys+7)->den.i, (porta_list[0]->sys+8)->num, (porta_list[0]->sys+8)->den.i);
			//printf("Inequalities :\n");
			////int i;
			//for (i = 0; i < ineq; i++)
			//{
				//for (j = 0; j < dim+1-equa; j++)
					//printf("%li, ", (porta_list[i]->sys+j)->num);
				//printf("\n");
			//};
			//printf("Equalities :\n");
			//for (i = 0; i < equa; i++)
			//{
				//for (j = 0; j < dim+1; j++)
					//printf("%li, ", (porta_list[ineq+i]->sys+j)->num);
				//printf("\n");
			//};
			//printf("equa, ineq = %i,%i\n", equa, ineq);
			//printf("dim+1 = %i\n", dim+1);
			//printf("ineq, 0 = %i,%i\n", ineq,0);
			//printf("dim+1-equa = %i\n", dim+1-equa);
			//for (j = 0; j < dim-equa; j++)
				//printf("*(indx+%i) = %i\n", j, *(indx+j));*/
        //}
    //}
    //else if (is_set(Traf) && ieq_file) 
    //{
        ///* INEQUALITIES TO POINTS */
        //RAT *inner,*iep;
        //char *cp1;
        //cp1=strdup("VALID");
        //inner = 0;
        //if(!cp1)
            //msg( "allocation of new space failed", "", 0 );   
            
        ///*points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
                                 //"\0",(int **)&i,
                                 //"\0",(int **)&i,cp1,&inner);*/
        //points = nbPoints;
        ///* Still need to elaborate here... */
        
        //free(cp1);
        //ar6 = inner; if (inner) nel_ar6 = dim;
        //sort_eqie_cvce(ar1,points,dim+2,&equa_in,&ineq_in);
        //iep = ar1+equa_in*(dim+2);
        ///* first equations then inequalities */
        //points = ineq_in;
        //polarformat(iep,&equa_in,ineq_in,inner);
        //gentableau(iep,0,&rowl_inar,&indx);
        //if(is_set(Long_arithmetic)) 
        //{
            //RAT a, b;
            //SET_MP_ready;
            //memset( &a, 0, sizeof(a) );
            //memset( &b, 0, sizeof(b) );
            //arith_overflow_func(0,0,a,b,0);
        //}
        //ineq = (cone == points) ? dim : dim + 1;
        //ineq_out = ineq; /* not used further */
        //gauss(1, points+dim+1,dim+1,dim,ineq,&ineq_out, &equa, indx);
        ///* make indx point to the x-variable section */
        //for (; (*indx) < 0; indx++);
        //fourier_motzkin(0,ineq-equa,points+dim+1-ineq,
                        //points-ineq+equa,poi_file,indx,0);
        //if (is_set(Validity_table_out)) 
            //red_test(indx,iep,&rowl_inar);
        //if (cone >= dim-equa)
            //origin_add(rowl_inar,iep); 
        //resubst(inner,equa_in,indx);
        //if ((MP_realised && return_from_mp()) || !MP_realised) 
        //{
            //if (equa)
                //sort(no_denom(dim+1, ineq, ineq+equa,1), 
                     //dim+1,ineq,ineq+equa);
            //sort(0, dim+1, 0, ineq);
        //}
        //for (cone = 0; !(porta_list[cone]->sys[dim].num); cone++);
        //conv = ineq - cone;
        //if (!MP_realised) no_denom(dim+1, 0, cone,1);
        ///*write_poi_file(*argv,outfp,dim,equa,ineq,cone,0,conv,cone);*/
        ///* Still need to return the computed points here ... */
    //}
    //else 
        //msg( "invalid format of command line", "", 0 );



    ///* 17.01.1994: include logging on file porta.log */
    //fclose( logfile );
    //exit(0);
//}

