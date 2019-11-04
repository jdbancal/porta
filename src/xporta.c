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
 

FILENAME: xporta.c

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/
/*  LAST EDIT: Fri Sep 20 13:36:44 2002 by Andreas Loebel (opt0.zib.de)  */
/* $Id: porta.c,v 1.5 2009/09/21 07:46:48 bzfloebe Exp $ */


/* This file was introduced as a split of porta.c by J-D B to separate 
 * the definitions necessary to run porta from the code of porta itself */

#define _CRT_SECURE_NO_WARNINGS 1


#include "porta.h"
#include "arith.h"
#include "common.h"
#include "log.h"
#include "inout.h"
#include "mp.h"
#include "four_mot.h"
#include "portsort.h"


FILE *logfile;


#define FILENAME_SIZE 2000 // Extended this to accept long strings as input for long filenames (which could also include paths...) J-D B

int main( int argc, char *argv[] )
{
    int i, ieq_file, start;
    char outfname[FILENAME_SIZE];
    char fname[FILENAME_SIZE];
    int   poi_file;
    int   rowl_inar, ierl;
    int  *indx = (int *)0;      
    int equa_in,ineq_in, ineq_out;
    FILE *outfp;

    

    printf("\nPORTA - a POlyhedron Representation Transformation Algorithm\n");
    printf(  "Version %s\n\n", VERSION );

    printf( "Written by Thomas Christof (Uni Heidelberg)\n" );
    printf( "Revised by Andreas Loebel (ZIB Berlin)\n\n" );

    printf( "PORTA is free software and comes with ABSOLUTELY NO WARRENTY! You are welcome\n" );
    printf( "to use, modify, and redistribute it under the GNU General Public Lincese.\n\n" ); 
    
    printf( "This is the program XPORTA from the PORTA package.\n\n" );

    if( argc <= 2 )
    {
        printf( "For more information read the manpages about porta.\n\n" );
        exit(-1);
    }
            
    /* 17.01.1994: include logging on file porta.log */
    logfile = fopen( "porta.log", "a" );
    if( !logfile )
        fprintf( stderr, "can't open logfile porta.log\n" );
    else
    {
        porta_log( "\n\n\nlog for " );
        for( i = 0; i < argc; i++ )
            porta_log( "%s ", argv[i] );
        porta_log( "\n\n" );
    }
            
    init_total_time();
    
    initialize();

    prt = stdout;
    get_options(&argc,&argv);
    
    if (option & Protocol_to_file) 
    {
        strcat(*argv,".prt");
        prt = fopen(*argv,"w");
        (*argv)[strlen(*argv)-4] = '\0';
    }
    setbuf(prt,CP 0);
    
    set_I_functions();
    SET_MP_not_ready;
    ieq_file = !strcmp(*argv+strlen(*argv)-4,".ieq");
    poi_file = !strcmp(*argv+strlen(*argv)-4,".poi");
    printf("%i, %i",ieq_file, poi_file);
    if (!poi_file && !ieq_file)
        msg( "invalid format of command line", "", 0 );
    
    /*
     * change by M.S. 5.6.92:
     * read_input_file writes to the output file, if is_set(Sort).
     */
    outfp = 0;
    strcpy(outfname,*argv);
    if (is_set(Sort) && poi_file) 
    {
        strcat(outfname,".poi");
        outfp = wfopen(outfname);
    }
    if (is_set(Sort) && ieq_file) 
    {
        strcat(outfname,".ieq");
        fprintf(prt,"outfname = %s\n",outfname);
        fflush(stdout);

        /* 17.01.1994: include logging on file porta.log */
        porta_log( "outfname = %s\n",outfname );
        fflush(logfile);

        outfp = wfopen(outfname);
    }

    if (is_set(Fmel) && ieq_file) 
    {  
        /* ONLY FM-ELIMINATON */
        
        int *elim_ord,nel;
            char *cp1;
        elim_ord = 0;
        cp1=strdup("ELIMINATION_ORDER");
        if(!cp1)
            msg( "allocation of new space failed", "", 0 );   
        ineq = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,cp1,
                               &elim_ord,"\0",(int **)&i,"\0",(RAT **)&i);
        free(cp1);
        sort_eqie_cvce(ar1,ineq,dim+2,&equa_in,&ineq_in);
        ineq = ineq_in+equa_in;
        /*     elim_ord = check_and_reorder_elim_ord(elim_ord,&nel); */
        reorder_var(ineq,ar1,&ar2,(int *)&nel_ar2,&nel,&elim_ord,&indx);
        /*     indx = elim_ord+nel; */
        /* 
         * Transform all inequalities into "<= 1" or "<=-1" inequalities,
         * if possible.
         * If the right-hand side is 0, divide the numerators by their gcd
         * and the denominators by their gcd.
         * (This is not really necessary).
         */
        /*
           for (i = 0; i < ineq; i++) 
           (* RAT_row_prim)(ar2+(dim+2)*i,ar2+(dim+1)*i,ar2+(dim+2)*i+dim,dim+1);
        */
        if(is_set(Long_arithmetic))
        {
            RAT a, b;
            SET_MP_ready;
            memset( &a, 0, sizeof(a) );
            memset( &b, 0, sizeof(b) );
            arith_overflow_func(0,0,a,b,0);
        }
        for (i = 0; i < ineq; i++) 
            (* RAT_row_prim)(ar2+(dim+1)*i,ar2+(dim+1)*i,
                             ar2+(dim+1)*i+dim,dim+1);
        /*     nel_ar2 = ineq*(dim+1);  
         */
        equa = 0;
        ineq_out = ineq;
        gauss(0, dim+1, dim+1, dim-nel, equa_in, &ineq_out, &equa, indx);
        for (; (*indx) < 0; indx++);
        ierl = dim-nel-equa +1;
        /* row-length of inequalities after fourier_motzkin
         * =  number of noneliminated variables+1 */
        nel = nel - (ineq - ineq_out);  /* number of variables to be elim. */
        ineq = ineq_out;
        fourier_motzkin(0,ineq-equa,dim+1-equa_in,nel,poi_file,indx,elim_ord);
        if ((MP_realised && return_from_mp()) || !MP_realised) 
            sort(no_denom(ierl, 0, ineq,1), ierl, 0, ineq);
        write_ieq_file(*argv,outfp,equa,ineq,dim+1,0,
                       ineq,0,ierl,indx);
        
    }
    else if (is_set(Sort)) 
    {
        points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,"\0",(int **)&i,
                                 "\0",(int **)&i,"\0",(RAT **)&i);
        if (ieq_file)
            sort_eqie_cvce(ar1,points,dim+2,&equa,&ineq);
        listptoar(ar1,points,ieq_file?dim+2:dim+1,0); 
        if (ieq_file) 
        { 
            if (equa) sort(1,dim+1,0,equa);
            if (ineq) sort(1,dim+1,equa,points); 
            write_ieq_file(*argv,outfp,equa,0,dim+1,0,ineq,equa,dim+1,0);
        }
        else 
        {
            sort(1,dim+1,0,points);
            for (cone = 0; !(porta_list[cone]->sys[dim].num); cone++);
            write_poi_file(*argv,outfp,dim,0,0,cone,0,points-cone,cone);
        }
    }  
    else if ((is_set(Traf) || is_set(Dim)) && poi_file) 
    {
        points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
                                 "\0",(int **)&i,"\0",(int **)&i,"\0",
                                 (RAT **)&i);
        gentableau(ar1,1,&rowl_inar,&indx); 
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
        if (is_set(Dim)) 
        {
			fprintf(prt, "\n\n way1\n\n");
            char str[100];
            fprintf (prt,"\nDIMENSION OF THE POLYHEDRON : %i\n\n",dim-equa);
            
            /* 17.01.1994: include logging on file porta.log */
            porta_log( "\nDIMENSION OF THE POLYHEDRON : %i\n\n", dim-equa );
            
            sprintf (str,"echo 'DIMENSION OF THE POLYHEDRON : %i' | cat >> %s",
                     dim-equa,argv[0]);
            system(str);
            if (equa) 
            {
                fprintf(prt,"equations :\n");

                /* 17.01.1994: include logging on file porta.log */
                porta_log( "equations :\n");
                
                listptoar(ar4,equa,dim+1,0); 
                if ((MP_realised && return_from_mp()) || !MP_realised) 
                    sort(no_denom(dim+1,0,equa,1), dim+1,0,equa);
                start = 1;
                /* last argument of writesys was lost? 
                   writesys(prt,0,equa,dim+1,0,0,'=');
                   */
                writesys(prt,0,equa,dim+1,0,0,'=', &start);

                /* log equation system */
                start = 1;
                writesys(logfile,0,equa,dim+1,0,0,'=', &start);
            }
        }
        else  
        {
			fprintf(prt, "\n\n way2\n\n");
            /* POINTS TO INEQUALITIES */
            sprintf(fname,"%s.ieq",*argv);
            fourier_motzkin(fname,ineq-equa,points+dim+1-ineq,
                            points-ineq+equa,poi_file,indx,0);
			fprintf(prt, "is_set(Validity_table_out) = %i\n", is_set(Validity_table_out));
            if (is_set(Validity_table_out)) 
                red_test(indx,ar1,&rowl_inar);
			if MP_realised
				fprintf(prt, "MP_realised = %i, return_from_mp() = %i\n", MP_realised, return_from_mp());
			else
				fprintf(prt, "MP_realised = %i\n", MP_realised);
            if ((MP_realised && return_from_mp()) || !MP_realised)// || 1) // last condition added by J-D B for test purposes... for real, we should find a better way to exit the MP section...
            {
				fprintf(prt, "equa = %i\n", equa);
                if (equa) sort(no_denom(dim+1, ineq, ineq+equa,1), 
                               dim+1, ineq, ineq+equa);
                sort(no_denom(dim+1-equa, 0, ineq,1), 
                     dim+1-equa, 0, ineq);
            }
            write_ieq_file(*argv,outfp,equa,ineq,
                           dim+1,0,ineq,0,dim+1-equa,indx);
        }
    }
    else if (is_set(Traf) && ieq_file) 
    {
        /* INEQUALITIES TO POINTS */
        RAT *inner,*iep;
        char *cp1;
        cp1=strdup("VALID");
        inner = 0;
        if(!cp1)
            msg( "allocation of new space failed", "", 0 );   
            
        points = read_input_file(*argv,outfp,&dim,&ar1,(int *)&nel_ar1,
                                 "\0",(int **)&i,
                                 "\0",(int **)&i,cp1,&inner);
        free(cp1);
        ar6 = inner; if (inner) nel_ar6 = dim;
        printf("The valid point : ");
        for (i = 0; i < dim; i++)
			printf("%li/%i ", ar6[i].num, ar6[i].den.i);
		printf("\n");
        printf("The valid point : ");
        for (i = 0; i < dim; i++)
			printf("%li/%i ", inner[i].num, inner[i].den.i);
		printf("\n");
		
		printf("point = %i\n",points);
		printf("dim = %i\n",dim);
		printf("nel_ar1 = %li\n",nel_ar1);
		
       	printf("address of ar1 = %p\n", ( void * ) ar1);
       	printf("address of equa_in = %p\n", ( void * ) equa_in);
       	printf("address of ineq_in = %p\n", ( void * ) ineq_in);
       	printf("value of equa_in = %i\n", ( void * ) equa_in);
       	printf("value of ineq_in = %i\n", ( void * ) ineq_in);
		printf("first few values of ar1 = ...\n");
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[0].num, ar1[0].den.i, ar1[1].num, ar1[1].den.i, ar1[2].num, ar1[2].den.i, ar1[3].num, ar1[3].den.i, ar1[4].num, ar1[4].den.i, ar1[5].num, ar1[5].den.i, ar1[6].num, ar1[6].den.i, ar1[7].num, ar1[7].den.i, ar1[8].num, ar1[8].den.i, ar1[9].num, ar1[9].den.i);
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[10].num, ar1[10].den.i, ar1[11].num, ar1[11].den.i, ar1[12].num, ar1[12].den.i, ar1[13].num, ar1[13].den.i, ar1[14].num, ar1[14].den.i, ar1[15].num, ar1[15].den.i, ar1[16].num, ar1[16].den.i, ar1[17].num, ar1[17].den.i, ar1[18].num, ar1[18].den.i, ar1[19].num, ar1[19].den.i);
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[20].num, ar1[20].den.i, ar1[21].num, ar1[21].den.i, ar1[22].num, ar1[22].den.i, ar1[23].num, ar1[23].den.i, ar1[24].num, ar1[24].den.i, ar1[25].num, ar1[25].den.i, ar1[26].num, ar1[26].den.i, ar1[27].num, ar1[27].den.i, ar1[28].num, ar1[28].den.i, ar1[29].num, ar1[29].den.i);
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[30].num, ar1[30].den.i, ar1[31].num, ar1[31].den.i, ar1[32].num, ar1[32].den.i, ar1[33].num, ar1[33].den.i, ar1[34].num, ar1[34].den.i, ar1[35].num, ar1[35].den.i, ar1[36].num, ar1[36].den.i, ar1[37].num, ar1[37].den.i, ar1[38].num, ar1[38].den.i, ar1[39].num, ar1[39].den.i);
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[40].num, ar1[40].den.i, ar1[41].num, ar1[41].den.i, ar1[42].num, ar1[42].den.i, ar1[43].num, ar1[43].den.i, ar1[44].num, ar1[44].den.i, ar1[45].num, ar1[45].den.i, ar1[46].num, ar1[46].den.i, ar1[47].num, ar1[47].den.i, ar1[48].num, ar1[48].den.i, ar1[49].num, ar1[49].den.i);
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[50].num, ar1[50].den.i, ar1[51].num, ar1[51].den.i, ar1[52].num, ar1[52].den.i, ar1[53].num, ar1[53].den.i, ar1[54].num, ar1[54].den.i, ar1[55].num, ar1[55].den.i, ar1[56].num, ar1[56].den.i, ar1[57].num, ar1[57].den.i, ar1[58].num, ar1[58].den.i, ar1[59].num, ar1[59].den.i);
        
        
        sort_eqie_cvce(ar1,points,dim+2,&equa_in,&ineq_in);
        iep = ar1+equa_in*(dim+2);
        /* first equations then inequalities */
        points = ineq_in;
        polarformat(iep,&equa_in,ineq_in,inner);
        
        
       	printf("address of iep = %p\n", ( void * ) iep);
		printf("first few values of iep = ...\n");
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", iep[0].num, iep[0].den.i, iep[1].num, iep[1].den.i, iep[2].num, iep[2].den.i, iep[3].num, iep[3].den.i, iep[4].num, iep[4].den.i, iep[5].num, iep[5].den.i, iep[6].num, iep[6].den.i, iep[7].num, iep[7].den.i, iep[8].num, iep[8].den.i, iep[9].num, iep[9].den.i);
		printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", iep[10].num, iep[10].den.i, iep[11].num, iep[11].den.i, iep[12].num, iep[12].den.i, iep[13].num, iep[13].den.i, iep[14].num, iep[14].den.i, iep[15].num, iep[15].den.i, iep[16].num, iep[16].den.i, iep[17].num, iep[17].den.i, iep[18].num, iep[18].den.i, iep[19].num, iep[19].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[20].num, ar1[20].den.i, ar1[21].num, ar1[21].den.i, ar1[22].num, ar1[22].den.i, ar1[23].num, ar1[23].den.i, ar1[24].num, ar1[24].den.i, ar1[25].num, ar1[25].den.i, ar1[26].num, ar1[26].den.i, ar1[27].num, ar1[27].den.i, ar1[28].num, ar1[28].den.i, ar1[29].num, ar1[29].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[30].num, ar1[30].den.i, ar1[31].num, ar1[31].den.i, ar1[32].num, ar1[32].den.i, ar1[33].num, ar1[33].den.i, ar1[34].num, ar1[34].den.i, ar1[35].num, ar1[35].den.i, ar1[36].num, ar1[36].den.i, ar1[37].num, ar1[37].den.i, ar1[38].num, ar1[38].den.i, ar1[39].num, ar1[39].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[40].num, ar1[40].den.i, ar1[41].num, ar1[41].den.i, ar1[42].num, ar1[42].den.i, ar1[43].num, ar1[43].den.i, ar1[44].num, ar1[44].den.i, ar1[45].num, ar1[45].den.i, ar1[46].num, ar1[46].den.i, ar1[47].num, ar1[47].den.i, ar1[48].num, ar1[48].den.i, ar1[49].num, ar1[49].den.i);
		//printf("%li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i, %li/%i\n", ar1[50].num, ar1[50].den.i, ar1[51].num, ar1[51].den.i, ar1[52].num, ar1[52].den.i, ar1[53].num, ar1[53].den.i, ar1[54].num, ar1[54].den.i, ar1[55].num, ar1[55].den.i, ar1[56].num, ar1[56].den.i, ar1[57].num, ar1[57].den.i, ar1[58].num, ar1[58].den.i, ar1[59].num, ar1[59].den.i);
        
		//printf("element number 32 : %li/%i\n", ar1[32].num, ar1[32].den.i);
		//printf("element number 56 : %li/%i\n", ar1[56].num, ar1[56].den.i);
		
		//printf("dim = %i\n", dim);
		//printf("nel_ar1 = %li\n", nel_ar1);
 
        
        gentableau(iep,0,&rowl_inar,&indx);
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
        write_poi_file(*argv,outfp,dim,equa,ineq,cone,0,conv,cone);
    }
    else 
        msg( "invalid format of command line", "", 0 );



    /* 17.01.1994: include logging on file porta.log */
    fclose( logfile );
    exit(0);
}

