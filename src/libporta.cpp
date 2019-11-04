/*******************************************************************************

Copyright (C) 1997-2013 Thomas Christof, Andreas Loebel and Jean-Daniel Bancal
 
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
 

FILENAME: libporta.cpp

AUTHOR: Jean-Daniel bancal

This file provides a c++ library interface to the Porta polytope software of
Thomas Christof and Andreas Loebel. Most of the code here was directly copied
from the porta.c, porta.h and xporta.c source files.

*******************************************************************************/

/* NOTE : It is not clear yet whether porta manages memory in a reversible
 * way or not. Especially with respect to the variable porta_list. So for
 * now the library could be using more memory than really needed... */

#include <iostream>
#include <string>
#include <vector>
#include "libporta.h"
#include "time.h"
extern "C" {
#include "porta.h"
#include "arith.h"
#include "common.h"
#include "log.h"
#include "inout.h"
#include "mp.h"
#include "four_mot.h"
#include "portsort.h"
}

using namespace std;


/* Here are some internal variables */
FILE *logfile;
bool libportaInitialized(false);


// Comparison between two fractions is done without trying to simplify the fraction here.
bool fraction::operator==(const fraction& f) const {return ((num == f.num) && (den == f.den));};

// Basic fraction arithmetic
fraction fraction::operator+ (const fraction& f) const
{
	fraction f2;
	f2.num = num*f.den + den*f.num;
	f2.den = den*f.den;
	f2.simplify();
	return f2;
};

fraction fraction::operator* (const fraction& f) const
{
	fraction f2;
	f2.num = num*f.num;
	f2.den = den*f.den;
	f2.simplify();
	return f2;
};

void fraction::simplify()
{
	long int gcd(longgcd(num, den));
	num = num/gcd;
	den = den/gcd;
};


/* This function initializes the internal variable of porta  and opens
 * the log file. */
void init_libporta()
{
	if (libportaInitialized) fclose( logfile );
	
	logfile = fopen( "porta.log", "a" );
    if( !logfile )
        fprintf( stderr, "can't open logfile porta.log\n" );
    else
    {
        porta_log("\n\n\nlog for libporta function VtoHrep\n\n");
    }

    init_total_time();
    
    initialize();

    prt = stdout;
    setbuf(prt,CP 0);
    
    set_I_functions();

    SET_MP_not_ready;

	// We also initialize some variables for porta.c
	mp_state = 0;
	nel_ar1 = 0;
	nel_ar2 = 0;
	nel_ar3 = 0;
	nel_ar4 = 0;
	nel_ar5 = 0;
	nel_ar6 = 0;
	maxlist = 0;
	total_size = 0;
	dim = 0;
	equa = 0;
	ineq = 0;
	conv = 0;
	cone = 0;
	points = 0;
	blocks = 0;
	option = 0;
	allowed_options = 0;

	// common.c
	comp = 0;
	
	// portsort.c
	comp_ps = 0;
	delay = 0;
	same_vals = 0;
	rowlen = 0;
	
	// four_mot.c
	blocks = 0;
	itr = 0;
	totalineq = 0;

	libportaInitialized = true;
	cout << "Interface with porta initialized.\n"<< endl;
//	long int* junk( logfile);
//	cout << logfile << " : " << (int *) logfile << endl << flush;

}


/* This function closes the log file of porta */
void close_libporta()
{
	if (libportaInitialized) fclose( logfile );
	libportaInitialized = false;
	
	// We also free some pointers which might still not have been freed:
	free(ar1); ar1 = 0;
	free(ar2); ar2 = 0;
	free(ar3); ar3 = 0;
	free(ar4); ar4 = 0;
	free(ar5); ar5 = 0;
	free(ar6); ar6 = 0;
	free(elim_in); elim_in = 0;
//	free(porta_list);
	porta_list = 0; // For now it seems not clear whether we can clear porta_list really (or whether some data were already lost... so we just lose the address to make sure we start with a new allocation.)
	
	cout << "Interface with porta closed.\n"<< endl;
}



void VtoHrep(const Vrep& vrep, Hrep& hrep)
{
    int ieq_file, start;
    char outfname[20];
    char fname[20];
    int   poi_file;
    int   rowl_inar, ierl;
    int  *indx = (int *)0;
    int  *indx0 = (int *)0; // We remember how indx was initializes to allow freeing the memory later...
    int equa_in,ineq_in, ineq_out;
    FILE *outfp;

    int i,j; // Added by J-D B, used as summation indices later.

    init_libporta();
	
    printf("\nPORTA - a POlyhedron Representation Transformation Algorithm\n");
    printf(  "Version %s\n\n", VERSION );

    printf( "Written by Thomas Christof (Uni Heidelberg)\n" );
    printf( "Revised by Andreas Loebel (ZIB Berlin)\n\n" );

    printf( "PORTA is free software and comes with ABSOLUTELY NO WARRENTY! You are welcome\n" );
    printf( "to use, modify, and redistribute it under the GNU General Public Lincese.\n\n" ); 
    
    printf( "This is the program XPORTA from the PORTA package.\n\n" );

    /*Here we will do a traf, so we manually set the option "-T" :*/
	option |= Traf;
	allowed_options = Traf|
		Chernikov_rule_off|Validity_table_out|
		Redundance_check|Statistic_of_coefficients|
		Protocol_to_file|Opt_elim|Long_arithmetic;
    
//    ieq_file = 0; /*!strcmp(*argv+strlen(*argv)-4,".ieq");*/
//    poi_file = 1; /*!strcmp(*argv+strlen(*argv)-4,".poi");*/
        
//    outfp = 0;

	points = vrep.nbPointsAndRays();
	dim = vrep.dimension();

	/* We allocate the memory needed and initialize the table */
	ar1 = (RAT *) RATallo(ar1, 0, (points+1)*(dim+1)); // Note that for some reason porta seems to allocate memory for one more extremal point... so we also do that here.
	for (i = 0; i < points; i++)
	{
		for (j = 0; j < dim+1; j++)
		{
			if (j != dim) /* We put the constant term at the end of the column. */
			{
				ar1[i*(dim+1)+j].num = vrep.pointsAndRays[i][j+1].num;
				ar1[i*(dim+1)+j].den.i = vrep.pointsAndRays[i][j+1].den;
			}
			else
			{
				ar1[i*(dim+1)+j].num = vrep.pointsAndRays[i][0].num;
				ar1[i*(dim+1)+j].den.i = vrep.pointsAndRays[i][0].den;
			}
		}
	}
	/* We also prepare the valid point */
	ar6 = (RAT *) RATallo(ar6, 0, (points+1)*dim);
	for (j = 0; j < dim; j++)
		ar6[j] = ar1[j];

	gentableau(ar1,1,&rowl_inar,&indx);
	indx0 = indx; // We copy the initial address of indx to later be able to free the memory...
	if(is_set(Long_arithmetic))
	{
		RAT a, b;
		SET_MP_ready;
		memset( &a, 0, sizeof(a) );
		memset( &b, 0, sizeof(b) );
		arith_overflow_func(0,0,a,b,0);
	}
	ineq = (cone == points) ? dim : dim + 1;
	ineq_out = ineq;  /*not used further */
	gauss(1, points+dim+1,dim+1,dim,ineq,&ineq_out, &equa, indx);
	/* make indx point to the system-variable section */
	for (; (*indx) < 0; indx++);

	/* POINTS TO INEQUALITIES */
	/*sprintf(fname,"%s.ieq",*argv);*/
	fourier_motzkin(fname,ineq-equa,points+dim+1-ineq,
					points-ineq+equa,poi_file,indx,0);
	if (is_set(Validity_table_out)) 
		red_test(indx,ar1,&rowl_inar);
	if ((MP_realised && return_from_mp()) || !MP_realised) 
	{
//		cout << "REACH HERE, equa=" << equa << ", dim+1-equa =" << dim+1-equa << ", ineq =" << ineq << endl << flush;
		if (equa) sort(no_denom(dim+1, ineq, ineq+equa,1), 
					   dim+1, ineq, ineq+equa);
		sort(no_denom(dim+1-equa, 0, ineq,1), 
			 dim+1-equa, 0, ineq);
	}
	/*write_ieq_file("output",outfp,equa,ineq,
				   dim+1,0,ineq,0,dim+1-equa,indx);*/
	
//	cout << libportaInitialized << endl << flush;
//	cout << logfile << endl << flush;

	// We feed the output into the hrep variable.
	hrep.fill(porta_list, dim, equa, ineq, indx);
	
	// We free the memory:
	free(indx0); indx0 = 0;
	/*cout << "AAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
	cout << ineq+equa << " " << dim << " " << *indx << " aoihoffs; maxlist = "  << maxlist << endl << flush;
	for (j = 0; j < maxlist; j++)
	{
		cout << "round " << j << " : porta_list[j] = " << porta_list[j] << endl << flush;
		if (porta_list[j] != 0)
		{
			free(porta_list[j]->sys);
			free(porta_list[j]->mark);
			free(porta_list[j]->ptr);
			free(porta_list[j]);
		}
	}
	free(porta_list);
	cout << "BBBBBBBBBBBBBBBBBBBBBBBBBB" << endl;*/
	
//	cout << logfile << endl << flush;

    close_libporta();
//	cout << "OK1" << endl << flush;
}

void HtoVrep(const Hrep& hrep, Vrep& vrep)
{
    int ieq_file, start;
    char outfname[20];
    char fname[20];
    int   poi_file;
    int   rowl_inar, ierl;
    int  *indx = (int *)0;
    int  *indx0 = (int *)0; // We remember how indx was initializes to allow freeing the memory later...
    int equa_in,ineq_in, ineq_out;
    FILE *outfp;

	long int nbEqs, nbIneqs;
    int i,j; // Added by J-D B, used as summation indices later.

    init_libporta();
	
    printf("\nPORTA - a POlyhedron Representation Transformation Algorithm\n");
    printf(  "Version %s\n\n", VERSION );

    printf( "Written by Thomas Christof (Uni Heidelberg)\n" );
    printf( "Revised by Andreas Loebel (ZIB Berlin)\n\n" );

    printf( "PORTA is free software and comes with ABSOLUTELY NO WARRENTY! You are welcome\n" );
    printf( "to use, modify, and redistribute it under the GNU General Public Lincese.\n\n" ); 
    
    printf( "This is the program XPORTA from the PORTA package.\n\n" );

    /*Here we will do a traf, so we manually set the option "-T" :*/
	option |= Traf;
	allowed_options = Traf|
		Chernikov_rule_off|Validity_table_out|
		Redundance_check|Statistic_of_coefficients|
		Protocol_to_file|Opt_elim|Long_arithmetic;
    
//    ieq_file = 0; /*!strcmp(*argv+strlen(*argv)-4,".ieq");*/
//    poi_file = 1; /*!strcmp(*argv+strlen(*argv)-4,".poi");*/
        
//    outfp = 0;

	nbIneqs = hrep.nbInequalities();
	nbEqs = hrep.nbEqualities();
	points = nbIneqs+nbEqs;
	dim = hrep.dimension();

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

	/* INEQUALITIES TO POINTS */

	/* TO BE DONE . */


		
	/*points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
							 "\0",(int **)&i,
							 "\0",(int **)&i,cp1,&inner);*/



	nel_ar1 = (dim+2)*(nbEqs+nbIneqs+1); // This seems to be what porta uses for some reason...

	RAT *inner,*iep;

	// We fill up the data inside the table :
	// First the valid point
	// If it was not provided, we just take the origin.
	ar6 = (RAT *) RATallo(ar6,0,dim);
	for (i = 0; i < dim; i++)
	{
		ar6[i].num = hrep.validPoint[i].num;
		ar6[i].den.i = hrep.validPoint[i].den;
//		ar6[i].num = 0; // To take the origin as a valid point
//		ar6[i].den.i = 1;
	}
	inner = ar6;
	nel_ar6 = dim;
	// Then the equalities and inequalities
	/* We allocate the memory needed and initialize the table */
	ar1 = (RAT *) RATallo(ar1, 0, (nbEqs+nbIneqs+1)*(dim+2));
	for (i = 0; i < nbEqs; i++)
	{
		for (j = 0; j < dim; j++)
		{
			ar1[i*(dim+2)+j].num = hrep.equalities[i][j+1].num;
			ar1[i*(dim+2)+j].den.i = hrep.equalities[i][j+1].den;
		}
		ar1[i*(dim+2)+dim].num = -hrep.equalities[i][0].num;
		ar1[i*(dim+2)+dim].den.i = hrep.equalities[i][0].den;
		ar1[i*(dim+2)+dim+1].num = 0;
		ar1[i*(dim+2)+dim+1].den.i = 0;
	}
	// Now the inequalities
	for (i = 0; i < nbIneqs; i++)
	{
		for (j = 0; j < dim; j++)
		{
			ar1[nbEqs*(dim+2)+i*(dim+2)+j].num = -hrep.inequalities[i][j+1].num;
			ar1[nbEqs*(dim+2)+i*(dim+2)+j].den.i = hrep.inequalities[i][j+1].den;
		}
		ar1[nbEqs*(dim+2)+i*(dim+2)+dim].num = hrep.inequalities[i][0].num;
		ar1[nbEqs*(dim+2)+i*(dim+2)+dim].den.i = hrep.inequalities[i][0].den;
		ar1[nbEqs*(dim+2)+i*(dim+2)+dim+1].num = 1;
		ar1[nbEqs*(dim+2)+i*(dim+2)+dim+1].den.i = 1;
	}
	
	
	
	        //printf("The valid point : ");
        //for (i = 0; i < dim; i++)
			//printf("%li/%i ", ar6[i].num, ar6[i].den.i);
		//printf("\n");
        //printf("The valid point : ");
        //for (i = 0; i < dim; i++)
			//printf("%li/%i ", inner[i].num, inner[i].den.i);
		//printf("\n");
		
		//printf("point = %i\n",points);
		//printf("dim = %i\n",dim);
		//printf("nel_ar1 = %i\n",nel_ar1);
		
       	//printf("address of ar1 = %p\n", ( void * ) ar1);
       	//printf("address of equa_in = %p\n", ( void * ) equa_in);
       	//printf("address of ineq_in = %p\n", ( void * ) ineq_in);
       	//printf("value of equa_in = %i\n", ( void * ) equa_in);
       	//printf("value of ineq_in = %i\n", ( void * ) ineq_in);
		//printf("first few values of ar1 = ...\n");
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[0].num, ar1[0].den.i, ar1[1].num, ar1[1].den.i, ar1[2].num, ar1[2].den.i, ar1[3].num, ar1[3].den.i, ar1[4].num, ar1[4].den.i, ar1[5].num, ar1[5].den.i, ar1[6].num, ar1[6].den.i, ar1[7].num, ar1[7].den.i, ar1[8].num, ar1[8].den.i, ar1[9].num, ar1[9].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[10].num, ar1[10].den.i, ar1[11].num, ar1[11].den.i, ar1[12].num, ar1[12].den.i, ar1[13].num, ar1[13].den.i, ar1[14].num, ar1[14].den.i, ar1[15].num, ar1[15].den.i, ar1[16].num, ar1[16].den.i, ar1[17].num, ar1[17].den.i, ar1[18].num, ar1[18].den.i, ar1[19].num, ar1[19].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[20].num, ar1[20].den.i, ar1[21].num, ar1[21].den.i, ar1[22].num, ar1[22].den.i, ar1[23].num, ar1[23].den.i, ar1[24].num, ar1[24].den.i, ar1[25].num, ar1[25].den.i, ar1[26].num, ar1[26].den.i, ar1[27].num, ar1[27].den.i, ar1[28].num, ar1[28].den.i, ar1[29].num, ar1[29].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[30].num, ar1[30].den.i, ar1[31].num, ar1[31].den.i, ar1[32].num, ar1[32].den.i, ar1[33].num, ar1[33].den.i, ar1[34].num, ar1[34].den.i, ar1[35].num, ar1[35].den.i, ar1[36].num, ar1[36].den.i, ar1[37].num, ar1[37].den.i, ar1[38].num, ar1[38].den.i, ar1[39].num, ar1[39].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[40].num, ar1[40].den.i, ar1[41].num, ar1[41].den.i, ar1[42].num, ar1[42].den.i, ar1[43].num, ar1[43].den.i, ar1[44].num, ar1[44].den.i, ar1[45].num, ar1[45].den.i, ar1[46].num, ar1[46].den.i, ar1[47].num, ar1[47].den.i, ar1[48].num, ar1[48].den.i, ar1[49].num, ar1[49].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[50].num, ar1[50].den.i, ar1[51].num, ar1[51].den.i, ar1[52].num, ar1[52].den.i, ar1[53].num, ar1[53].den.i, ar1[54].num, ar1[54].den.i, ar1[55].num, ar1[55].den.i, ar1[56].num, ar1[56].den.i, ar1[57].num, ar1[57].den.i, ar1[58].num, ar1[58].den.i, ar1[59].num, ar1[59].den.i);


	
	// Ok, now comes porta's original code:
	sort_eqie_cvce(ar1,points,dim+2,&equa_in,&ineq_in);
	iep = ar1+equa_in*(dim+2);
	/* first equations then inequalities */
	points = ineq_in;
	polarformat(iep,&equa_in,ineq_in,inner);
	gentableau(iep,0,&rowl_inar,&indx);
	indx0 = indx; // We copy the initial address of indx to later be able to free the memory...
	if(is_set(Long_arithmetic)) 
	{
		RAT a, b;
		SET_MP_ready;
		memset( &a, 0, sizeof(a) );
		memset( &b, 0, sizeof(b) );
		arith_overflow_func(0,0,a,b,0);
	}
	ineq = (cone == points) ? dim : dim + 1;
	ineq_out = ineq; /* not used further */
	gauss(1, points+dim+1,dim+1,dim,ineq,&ineq_out, &equa, indx);
	/* make indx point to the x-variable section */
	for (; (*indx) < 0; indx++);
	fourier_motzkin(0,ineq-equa,points+dim+1-ineq,
					points-ineq+equa,poi_file,indx,0);
	if (is_set(Validity_table_out)) 
		red_test(indx,iep,&rowl_inar);
	if (cone >= dim-equa)
		origin_add(rowl_inar,iep); 
	resubst(inner,equa_in,indx);
	if ((MP_realised && return_from_mp()) || !MP_realised) 
	{
		if (equa)
			sort(no_denom(dim+1, ineq, ineq+equa,1), 
				 dim+1,ineq,ineq+equa);
		sort(0, dim+1, 0, ineq);
	}
	for (cone = 0; !(porta_list[cone]->sys[dim].num); cone++);
	conv = ineq - cone;
	if (!MP_realised) no_denom(dim+1, 0, cone,1);
	/*write_poi_file(*argv,outfp,dim,equa,ineq,cone,0,conv,cone);*/
	/* Still need to return the computed points here ... */

   
  	// We feed the output into the vrep variable.
	//vrep.fill(...porta_list, dim, equa, ineq, indx);



			//printf("first few values of porta_list = ...\n");
			//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", (porta_list[0]->sys+0)->num, (porta_list[0]->sys+0)->den.i, (porta_list[0]->sys+1)->num, (porta_list[0]->sys+1)->den.i, (porta_list[0]->sys+2)->num, (porta_list[0]->sys+2)->den.i, (porta_list[0]->sys+3)->num, (porta_list[0]->sys+3)->den.i, (porta_list[0]->sys+4)->num, (porta_list[0]->sys+4)->den.i, (porta_list[0]->sys+5)->num, (porta_list[0]->sys+5)->den.i, (porta_list[0]->sys+6)->num, (porta_list[0]->sys+6)->den.i, (porta_list[0]->sys+7)->num, (porta_list[0]->sys+7)->den.i, (porta_list[0]->sys+8)->num, (porta_list[0]->sys+8)->den.i);
			//printf("Inequalities :\n");
			////int i;
			//for (i = 0; i < cone; i++)
			//{
				//for (j = 0; j < dim+1; j++)
					//printf("%li, ", (porta_list[i]->sys+j)->num);
				//printf("\n");
			//};
			//printf("Equalities :\n");
			//for (i = 0; i < conv; i++)
			//{
				//for (j = 0; j < dim+1; j++)
					//printf("%li, ", (porta_list[cone+i]->sys+j)->num);
				//printf("\n");
			//};
			//printf("equa, ineq = %i,%i\n", equa, ineq);
			//printf("cone, conv = %i,%i\n", cone, conv);
			//printf("dim+1 = %i\n", dim+1);
			//printf("ineq, 0 = %i,%i\n", ineq,0);
			//printf("dim+1-equa = %i\n", dim+1-equa);
			//for (j = 0; j < dim-equa; j++)
				//printf("*(indx+%i) = %i\n", j, *(indx+j));


	// We feed the output into the hrep variable.
	vrep.fill(porta_list, dim, cone, conv);
	
	// We free the memory:
	// We free the memory:
	free(indx0); indx0 = 0;
	cout << cone+conv << " " << dim << " asfafw ef; maxlist = " << maxlist << endl;
/*	for (j = 0; j < maxlist; j++)
	{
		cout << "step " << j << " : " << endl << flush;
		free(porta_list[j]->sys);
		free(porta_list[j]->mark);
		free(porta_list[j]->ptr);
		free(porta_list[j]);
	}
	free(porta_list);*/
/*	for (j = 0; j < cone+conv; j++)
	{
		free(porta_list[j]->mark);
		free(porta_list[j]);
	}*/
//	free(porta_list);
	
	close_libporta();
}


// The == operator
bool Vrep::operator==(const Vrep& vrep) const
{
	return (pointsAndRays == vrep.pointsAndRays);
};

ostream& operator<<(ostream& os, const Vrep& vrep)
{
	// Here we write the equalities and inequalities in the lrs format
	os << "V-representation\n";
	os << "begin\n";
	os << vrep.pointsAndRays.size() << " " << vrep.pointsAndRays[0].size() << " rational\n";
	for (long int i(0); i < vrep.pointsAndRays.size(); i++)
	{
		for (long int j(0); j < vrep.pointsAndRays[0].size(); j++)
		{
			os << vrep.pointsAndRays[i][j].num;
			if (vrep.pointsAndRays[i][j].den != 1)
				os << "/" << vrep.pointsAndRays[i][j].den << " ";
			else
				os << " ";
		};
		os << "\n";
	};
	os << "end\n\n";

	return os;
}

// The == operator
bool Hrep::operator==(const Hrep& hrep) const
{
	return ((equalities == hrep.equalities) && (inequalities == hrep.inequalities));
};

ostream& operator<<(ostream& os, const Hrep& hrep)
{
	// Here we write the equalities and inequalities in the lrs format
	os << "H-representation\n";
	if (hrep.equalities.size() > 0)
	{
		os << "linearity " << hrep.equalities.size();
		for (long int i(0); i < hrep.equalities.size(); i++)
			os << " " << i+1;
		os << "\n";
	};
	os << "begin\n";
	os << hrep.equalities.size()+hrep.inequalities.size() << " " << hrep.inequalities[0].size() << " rational\n";
	for (long int i(0); i < hrep.equalities.size(); i++)
	{
		for (long int j(0); j < hrep.equalities[0].size(); j++)
		{
			os << hrep.equalities[i][j].num;
			if (hrep.equalities[i][j].den != 1)
				os << "/" << hrep.equalities[i][j].den << " ";
			else
				os << " ";
		};
		os << "\n";
	};
	for (long int i(0); i < hrep.inequalities.size(); i++)
	{
		for (long int j(0); j < hrep.inequalities[i].size(); j++)
		{
			os << hrep.inequalities[i][j].num;
			if (hrep.inequalities[i][j].den != 1)
				os << "/" << hrep.inequalities[i][j].den << " ";
			else
				os << " ";
		};
		os << "\n";
	};
	os << "end\n\n";	
	
	return os;
}



// Functions from the Vrep class :


void Vrep::fill(const Tnum *data, const long int& nbPointsOrRays, const long int& dimension)
{
	//First we empty the table and reinitialise it.
	pointsAndRays.clear();
	pointsAndRays.resize(nbPointsOrRays, vector < fraction > (dimension+1, fraction()));
	
	for (long int i(0); i < nbPointsOrRays; i++)
		for (long int j(0); j < dimension+1; j++)
			pointsAndRays[i][j].num = data[i*(dimension+1)+j];
};


void Vrep::fill(const Tnum *data_num, const Tden *data_den, const long int& nbPointsOrRays, const long int& dimension)
{
	//First we empty the table and reinitialise it.
	pointsAndRays.clear();
	pointsAndRays.resize(nbPointsOrRays, vector < fraction > (dimension+1, fraction()));

	for (long int i(0); i < nbPointsOrRays; i++)
		for (long int j(0); j < dimension+1; j++)
		{
			pointsAndRays[i][j].num = data_num[i*(dimension+1)+j];
			pointsAndRays[i][j].den = data_den[i*(dimension+1)+j];
		}
};


void Vrep::fill(const string& filename)
{
	// We read the dimension first
	// Then allocate the memory
	// And fill up the table
};


void Vrep::fill(const listp* porta_list, const long int& dimension, const long int& nbCone, const long int& nbConv)
{
	//First we empty the tables and reinitialise them.
	pointsAndRays.clear();
	pointsAndRays.resize(nbCone+nbConv, vector < fraction > (dimension+1, fraction()));

	// First, we copy the extremal rays
	for (int i(0); i < nbCone; i++)
	{
		pointsAndRays[i][0].num = (porta_list[i]->sys+dimension)->num; // The constant is brought at the front of each line // CHECK THAT THIS NUMBER IS ZERO AS EXPECTED!!!
		pointsAndRays[i][0].den = (porta_list[i]->sys+dimension)->den.i;
//		cout << pointsAndRays[i][0].num;
		for (int j(0); j < dimension; j++)
		{
			pointsAndRays[i][1+j].num = (porta_list[i]->sys+j)->num;
			pointsAndRays[i][1+j].den = (porta_list[i]->sys+j)->den.i;
//			cout << pointsAndRays[i][1+j].num;
		};
//		cout << "fin" << endl;
	};
	
	// And copy also the extremal vertices
	for (int i(0); i < nbConv; i++)
	{
		pointsAndRays[nbCone+i][0].num = (porta_list[nbCone+i]->sys+dimension)->num;
		pointsAndRays[nbCone+i][0].den = (porta_list[nbCone+i]->sys+dimension)->den.i;
		for (int j(0); j < dimension; j++)
		{
			pointsAndRays[nbCone+i][1+j].num = (porta_list[nbCone+i]->sys+j)->num;
			pointsAndRays[nbCone+i][1+j].den = (porta_list[nbCone+i]->sys+j)->den.i;
		};
	};
};

// Provides some useful info about the polytope
long int Vrep::nbPointsAndRays() const
{
	return pointsAndRays.size();
};
long int Vrep::dimension() const
{
	if (nbPointsAndRays() > 0)
		return pointsAndRays[0].size()-1;
	else
		return 0;
};



// Functions from the Hrep class:


void Hrep::fill(const Tnum *dataIneq, const long int& nbIneq, const long int& dimension)
{
	//First we empty the tables and reinitialise them.
	equalities.clear();
	inequalities.clear();
	inequalities.resize(nbIneq, vector < fraction > (dimension+1, fraction()));

	for (long int i(0); i < nbIneq; i++)
		for (long int j(0); j < dimension+1; j++)
			inequalities[i][j].num = dataIneq[i*(dimension+1)+j];
	
	// We use the null point as a valid point since it was not provided
	validPoint.clear();
	validPoint.resize(dimension, fraction());
	
	checkValidPoint();
};


void Hrep::fill(const Tnum *dataIneq_num, const Tden *dataIneq_den, const long int& nbIneq, const long int& dimension)
{
	//First we empty the tables and reinitialise them.
	equalities.clear();
	inequalities.clear();
	inequalities.resize(nbIneq, vector < fraction > (dimension+1, fraction()));

	for (long int i(0); i < nbIneq; i++)
		for (long int j(0); j < dimension+1; j++)
		{
			inequalities[i][j].num = dataIneq_num[i*(dimension+1)+j];
			inequalities[i][j].den = dataIneq_den[i*(dimension+1)+j];
		}

	// We use the null point as a valid point since it was not provided
	validPoint.clear();
	validPoint.resize(dimension, fraction());
	
	checkValidPoint();
};


void Hrep::fill(const Tnum *dataEq, const long int& nbEq, const Tnum *dataIneq, const long int& nbIneq, const long int& dimension)
{
	//First we empty the tables and reinitialise them.
	equalities.clear();
	equalities.resize(nbEq, vector < fraction > (dimension+1, fraction()));
	inequalities.clear();
	inequalities.resize(nbIneq, vector < fraction > (dimension+1, fraction()));

	for (long int i(0); i < nbEq; i++)
		for (long int j(0); j < dimension+1; j++)
			equalities[i][j].num = dataEq[i*(dimension+1)+j];
	for (long int i(0); i < nbIneq; i++)
		for (long int j(0); j < dimension+1; j++)
			inequalities[i][j].num = dataIneq[i*(dimension+1)+j];
	
	// We use the null point as a valid point since it was not provided
	validPoint.clear();
	validPoint.resize(dimension, fraction());
	
	checkValidPoint();
};


void Hrep::fill(const Tnum *dataEq_num, const Tden *dataEq_den, const long int& nbEq, const Tnum *dataIneq_num, const Tden *dataIneq_den, const long int& nbIneq, const long int& dimension)
{
	//First we empty the tables and reinitialise them.
	equalities.clear();
	equalities.resize(nbEq, vector < fraction > (dimension+1, fraction()));
	inequalities.clear();
	inequalities.resize(nbIneq, vector < fraction > (dimension+1, fraction()));

	for (long int i(0); i < nbEq; i++)
		for (long int j(0); j < dimension+1; j++)
		{
			equalities[i][j].num = dataEq_num[i*(dimension+1)+j];
			equalities[i][j].den = dataEq_den[i*(dimension+1)+j];
		}
	for (long int i(0); i < nbIneq; i++)
		for (long int j(0); j < dimension+1; j++)
		{
			inequalities[i][j].num = dataIneq_num[i*(dimension+1)+j];
			inequalities[i][j].den = dataIneq_den[i*(dimension+1)+j];
		}
	
	// We use the null point as a valid point since it was not provided
	validPoint.clear();
	validPoint.resize(dimension, fraction());
	
	checkValidPoint();
};


void Hrep::fill(const string& filename)
{
	// We read the dimension first
	// Then allocate the memory
	// And fill up the table
	
	// We use the null point as a valid point since it was not provided
	validPoint.clear();
	validPoint.resize(dimension(), fraction());
	
	checkValidPoint();
};


void Hrep::fill(const listp* porta_list, const long int& dimension, const long int& nbEq, const long int& nbIneq, const int* indx)
{
	//First we empty the tables and reinitialise them.
	equalities.clear();
	equalities.resize(nbEq, vector < fraction > (dimension+1, fraction()));
	inequalities.clear();
	inequalities.resize(nbIneq, vector < fraction > (dimension+1, fraction()));

cout << endl << endl << "Hrep::fill, with dimension = " << dimension << ", nbEq = " << nbEq << ", nbIneq = " << nbIneq << endl << endl << endl << flush;

	// First, we copy the equality constraints
	for (int i(0); i < nbEq; i++)
	{
		equalities[i][0].num = (porta_list[nbIneq+i]->sys+dimension)->num; // The constant is brought at the front of each line
		equalities[i][0].den = (porta_list[nbIneq+i]->sys+dimension)->den.i;
		for (int j(0); j < dimension; j++)
		{
			equalities[i][1+j].num = (porta_list[nbIneq+i]->sys+j)->num;
			equalities[i][1+j].den = (porta_list[nbIneq+i]->sys+j)->den.i;
		};
	};
	
	// We prepare the inverse index vector to know for each variable
	// whether it is part of the inequality section or not.
	int indices[dimension];
	for (int j(0); j < dimension; j++)
		indices[j] = -1;
	for (int j(0); j < dimension-nbEq; j++)
		indices[*(indx+j)] = j;

	// And copy also the inequalities (with a zero coefficient for variables which are already fixed by the equality constraints)
	for (int i(0); i < nbIneq; i++)
	{
		inequalities[i][0].num = (porta_list[i]->sys+dimension-nbEq)->num;
		inequalities[i][0].den = (porta_list[i]->sys+dimension-nbEq)->den.i;
		for (int j(0); j < dimension; j++)
		{
			if (indices[j] != -1) // Otherwise the current variable does not appear in the inequality description
			{
				inequalities[i][1+j].num = -(porta_list[i]->sys+indices[j])->num; // These terms are on the other side of the inequality sign
				inequalities[i][1+j].den = (porta_list[i]->sys+indices[j])->den.i;
			}
		};
	};
	
	// We use the null point as a valid point since it was not provided
//	cout << "0" << libportaInitialized << endl;
	validPoint.clear();
	validPoint.resize(dimension, fraction());

//	cout << "1" << libportaInitialized << endl;

	checkValidPoint();
};


// Provides some useful info about the polytope
long int Hrep::nbEqualities() const
{
	return equalities.size();
};
long int Hrep::nbInequalities() const
{
	return inequalities.size();
};
long int Hrep::dimension() const
{
	if (nbEqualities() > 0)
		return equalities[0].size()-1;
	else if (nbInequalities() > 0)
		return inequalities[0].size()-1;
	else
		return 0;
};

void Hrep::checkValidPoint() const
{
//	cout << "2" << libportaInitialized << endl;
	for (long int i(0); i < nbEqualities(); i++)
	{
		fraction value(equalities[i][0]);
		for (long int j(1); j < dimension(); j++)
			value = value + equalities[i][j]*validPoint[j-1];
			
		if (value.num != 0)
		{
			cout << "Error, the point does not satisfy the equality constraints!" << endl;
			return;
		}
	}
	
	for (long int i(0); i < nbInequalities(); i++)
	{
		fraction value(inequalities[i][0]);
		for (long int j(1); j < dimension(); j++)
			value = value + inequalities[i][j]*validPoint[j-1];
			
		if (value.num < 0)
		{
			cout << "Error, the point does not satisfy the inequality constraints!" << endl;
			return;
		}
	}
};





/* The following added by J-D Bancal on 22.2.2012 to call porta as a c library.
	The data variable is expected to be of the standard *.ext format, i.e.
	
	[1 -1
	1  2]   ,  which is given as data=[1 -1 1 2] (i.e. line by line)
	
	for a 1-dimensional polytope with extreme points x=-1, and x=+2 (dimension = 1, nbPoints = 2).

*/
//void quick_and_dirty_poi_conv_call(long int *data, int dimension, int nbPoints)
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

    //if (~libportaInitialized)
		//init_libporta();

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
//}
