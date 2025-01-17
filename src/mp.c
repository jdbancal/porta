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
 

FILENAME: mp.c

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/
/*  LAST EDIT: Tue Aug 13 16:10:51 2002 by Andreas Loebel (opt0.zib.de)  */
/* $Id: mp.c,v 1.2 2009/09/21 07:05:11 bzfloebe Exp $ */


#include "common.h"
#include "mp.h"
#include "largecalc.h"


// Note by J-D B: The following conditions maybe could be revised in view of long int and long long int data types...
#define m_31        (1 << 31)
#define get_numlen(x)   ((abs((x)->num)>>20) & m0_9)
#define get_denlen(x)   ((abs((x)->num)>>10) & m0_9)
#define get_len(x)  (abs((x)->num) & m0_9)
#define set_len(x,nl,dl,l)  ((x)->num=(nl<<20)|(dl<<10)|l)
#define m0_9    1023
#define m10_19  1047552
#define m0_3    15
#define m0_7    255
#define m0_11   4095
#define lord(la,lb) (lorder(la.val,lb.val,la.len,lb.len))







void arith_overflow_func( int cf, void (*call_func)(), RAT a, RAT b, RAT *c )
{ 
    if (MP_not_ready)  /* not ready for multiple precision artihmetic */
        msg("Arithmetic overflow - try again with long integer arithmetic (`-l' option)", "", 0 );


    if (cf && RAT_add == call_func)
        call_func = L_RAT_add;
    else if (cf && RAT_sub == call_func)
        call_func = L_RAT_sub;
    else if (cf && RAT_mul == call_func)
        call_func = L_RAT_mul;
    
    RAT_add = L_RAT_add;
    RAT_sub = L_RAT_sub;
    RAT_mul = L_RAT_mul;
    RAT_assign = L_RAT_assign;
    RAT_row_prim = L_RAT_row_prim;
    writeline = (void (*)())L_RAT_writeline;
    
    RAT_to_L_RAT(ar1,nel_ar1);
    RAT_to_L_RAT(ar2,nel_ar2);
    RAT_to_L_RAT(ar3,nel_ar3);
    RAT_to_L_RAT(ar4,nel_ar4);
    RAT_to_L_RAT(ar5,nel_ar5);
    RAT_to_L_RAT(ar6,nel_ar6);
    RAT_to_L_RAT(RAT_const,2);
    RAT_to_L_RAT(var,4);
    
    if (cf && call_func) 
    {
        RAT_to_L_RAT(&a,1);
        RAT_to_L_RAT(&b,1);
        (*call_func)(a,b,c);
        L_RAT_to_RAT(&a,1);
        L_RAT_to_RAT(&b,1);
    }
    
    SET_MP_realised; 
    
	fprintf(prt, "\n\n\n --- SET_MP_realised = %i --- \n\n\n", SET_MP_realised);
    
}









void L_RAT_assign( RAT *a, RAT *b )
{ 
    int len;
    unsigned *pa,*pb; 
    
    a->den.p = UP allo(CP a->den.p,U get_len(a)*szU,U ((len = get_len(b))*szU));
    for (pa = a->den.p,pb = b->den.p; len > 0; len--,pa++,pb++)
        *pa = *pb;
    a->num = b->num;
    
}





void RAT_to_L_RAT( RAT *x, int n )
/* transforms the n int-RAT elements of array x to mp-RAT values */
{ 
    int num,den;
    
    for (; n; n--,x++) 
        if (x->num) 
        {
            num = x->num;
            den = x->den.i;
            set_len(x,1,1,2);
            if (num < 0) x->num *= -1;
            x->den.p = UP allo(CP x->den.p,0,2*szU);
            x->den.p[0] = abs(num);
            x->den.p[1] = den;
        }
    
}








int return_from_mp()
{

	fprintf(prt, "vals_lt_MAXINT(ar1,nel_ar1) = %i\n", vals_lt_MAXINT(ar1,nel_ar1));
	fprintf(prt, "vals_lt_MAXINT(ar2,nel_ar2) = %i\n", vals_lt_MAXINT(ar2,nel_ar2));
	fprintf(prt, "vals_lt_MAXINT(ar3,nel_ar3) = %i\n", vals_lt_MAXINT(ar3,nel_ar3));
	fprintf(prt, "vals_lt_MAXINT(ar4,nel_ar4) = %i\n", vals_lt_MAXINT(ar4,nel_ar4));
	fprintf(prt, "vals_lt_MAXINT(ar5,nel_ar5) = %i\n", vals_lt_MAXINT(ar5,nel_ar5));
	fprintf(prt, "vals_lt_MAXINT(ar6,nel_ar6) = %i\n", vals_lt_MAXINT(ar6,nel_ar6));

    if (!vals_lt_MAXINT(ar1,nel_ar1) || !vals_lt_MAXINT(ar2,nel_ar2) ||
        !vals_lt_MAXINT(ar3,nel_ar3) || !vals_lt_MAXINT(ar4,nel_ar4) ||
        !vals_lt_MAXINT(ar5,nel_ar5) || !vals_lt_MAXINT(ar6,nel_ar6))
        return(0); 
    
    L_RAT_to_RAT(ar1,nel_ar1);
    L_RAT_to_RAT(ar2,nel_ar2);
    L_RAT_to_RAT(ar3,nel_ar3);
    L_RAT_to_RAT(ar4,nel_ar4);
    L_RAT_to_RAT(ar5,nel_ar5);
    L_RAT_to_RAT(ar6,nel_ar6);
    L_RAT_to_RAT(RAT_const,2);
    L_RAT_to_RAT(var,4);
    
    set_I_functions(); 
    SET_MP_ready;
    
    return(1);
    
}








void L_RAT_to_RAT( RAT *x, int n )
/* transforms the n mp-RAT elements of array x to int-RAT values */
{ 
    int num,den;
    
    for (; n > 0; n--,x++) 
        if (!x->num)
            x->den.i = 1;
        else {
        den = x->den.p[1];
        num = x->den.p[0]*((x->num<0)?-1:1);
        allo(CP x->den.p,U get_len(x)*szU,0);
        x->num = num;
        x->den.i = den;
    }
    
}







int vals_lt_MAXINT( RAT *x, int n )
/* n mp-RAT values of array x less than max int ?
   yes: returns 1
   no : returns 0

   global: mx = MAXINT, my */
{ 
    
    for (; n > 0; n--, x++) 
        if (x->num) 
        {
            if (get_numlen(x) > 1 || get_denlen(x) > 1)
                return(0);
            if ((x->den.p[0] & m_31) || (x->den.p[1] & m_31))
                return(0);
        }
    
    
    return(1); 
    
}









void expand( unsigned *in, int inlen, loint *out )
{
    int i,r;
    unsigned *pout,x;
    
    for (i = 0, pout = out->val; i < inlen; i++) 
    {
        x = in[i];
        if ((r=i%3) == 1) 
        {
            *pout++ |= ((x & m0_3) << 8);
            x >>= 4;
        }
        else if (r == 2) 
        {
            *pout++ |= ((x & m0_7) << 4);
            x >>= 8;
        }
        *pout++ = x & m0_11; x >>= 12;
        *pout++ = x & m0_11; x >>= 12;
        if (r == 1) 
        {
            *pout = x; 
        }
        else if (!r)
            *pout = x; 
    }
    if (inlen%3 == 0) --pout;
    for ( ; (pout >= out->val) && !(*pout); pout--);
    out->len = pout - out->val+1;
    out->neg = 0;
    
}










void L_RAT_to_lorat( RAT a, lorat *la )
{
    int nl; 

    expand(a.den.p,(nl = get_numlen(&a)),&la->num);
    la->num.neg = (a.num < 0);
    la->den.neg = 0;
    if (nl)
        expand(a.den.p+nl,get_denlen(&a),&la->den);
    else
        la->den.val[0] = la->den.len = 1;
    
}










void comprim( loint in, unsigned *out )
{ 
    int i,r;
    unsigned *pout,x;
    
    for (i = 0,pout = out; i < in.len; i++) 
    {
        x = in.val[i];
        r = i%8;
        switch(r) 
        {
        case 0 : *pout = x; break;
        case 1 : *pout |= (x << 12); break;
        case 2 : *pout |= ((x & m0_7) << 24); *++pout = (x >> 8); break;
        case 3 : *pout |= (x << 4); break;
        case 4 : *pout |= (x << 16); break;
        case 5 : *pout |= ((x & m0_3) << 28); *++pout = (x >> 4); break;
        case 6 : *pout |= (x << 8); break;
        case 7 : *pout |= (x << 20);pout++; break;
        }
    }
    
}






#define storelen(x) ((3*x)/8 + (((3*x)%8)?1:0))


void lorat_to_L_RAT( lorat la, RAT *a )
{ 
    int oldsz,cnstlen,cdstlen,stlen,i;
    
    oldsz = get_len(a)*szU;
    cnstlen = storelen(la.num.len);
    cdstlen = ((cnstlen) ? storelen(la.den.len) : 0);
    set_len(a,cnstlen,cdstlen,(stlen = cnstlen+cdstlen)); 
    a->den.p = UP allo(CP a->den.p,U oldsz,U stlen*szU);
    if (stlen) 
    {
        comprim(la.num,a->den.p);
        comprim(la.den,a->den.p+cnstlen);
        if (!a->den.p[cnstlen-1]) 
        { 
            cnstlen--;
            set_len(a,cnstlen,cdstlen, (stlen=cnstlen+cdstlen)); 
            for (i=0; i < cdstlen; i++)
                a->den.p[cnstlen+i] = a->den.p[cnstlen+i+1]; 
            a->den.p = UP allo(CP a->den.p,U (stlen+1)*szU,U stlen*szU);
        }
        if (!a->den.p[stlen-1]) 
        { 
        cdstlen--;
        set_len(a,cnstlen,cdstlen, (stlen=cnstlen+cdstlen)); 
        a->den.p = UP allo(CP a->den.p,U (stlen+1)*szU,U stlen*szU);
      }
    }

    if (la.num.neg)
      a->num *= -1;

}







void L_RAT_add( RAT a, RAT b, RAT *c )
{ 
  loint lx,ly;
  lorat la,lb,lc;

  if (!a.num)
    L_RAT_assign(c,&b);

  else if (!b.num)
    L_RAT_assign(c,&a);

  else {
    L_RAT_to_lorat(a,&la);
    L_RAT_to_lorat(b,&lb);

    if (!lord(la.den,lb.den)) {
      lc.den = la.den;
      ladd(la.num,lb.num,&lc.num);
    }  
    else {
      lmul(la.den,lb.num,&lx);
      lmul(la.num,lb.den,&ly);
      lmul(la.den,lb.den,&lc.den);
      ladd(lx,ly,&lc.num);
    }  

    L_RAT_kue(&lc.num,&lc.den);
    lorat_to_L_RAT(lc,c);
  }
   
}





void L_RAT_sub( RAT a, RAT b, RAT *c )
{ 
  loint lx,ly;
  lorat la,lb,lc;

  if (!a.num) {
    L_RAT_assign(c,&b);
    c->num *= -1;
  }

  else if (!b.num)
    L_RAT_assign(c,&a);

  else {
    L_RAT_to_lorat(a,&la);
    L_RAT_to_lorat(b,&lb);

    if (!lord(la.den,lb.den)) {
      lc.den = la.den;
      lsub(la.num,lb.num,&lc.num);
    }  
    else {
      lmul(la.num,lb.den,&lx);
      lmul(lb.num,la.den,&ly);
      lmul(la.den,lb.den,&lc.den);
      lsub(lx,ly,&lc.num);
    }  

    L_RAT_kue(&lc.num,&lc.den);
    lorat_to_L_RAT(lc,c);
  }
   
}








void L_RAT_mul( RAT a, RAT b, RAT *c )
{ 
  lorat la,lb,lc;

  if (!a.num || !b.num) 
    L_RAT_assign(c,RAT_const);

  else {
    L_RAT_to_lorat(a,&la);
    L_RAT_to_lorat(b,&lb);
  
    lmul(la.num,lb.num,&lc.num);
    lmul(la.den,lb.den,&lc.den);

    L_RAT_kue(&lc.num,&lc.den);
    lorat_to_L_RAT(lc,c);
  }

}







void L_RAT_kue( loint *x, loint *y )
{
  loint z,r;

  lgcd(*x,*y,&z);
  
  if (z.len > 1 || (z.len == 1 && z.val[0] != 1)) {
    porta_ldiv(*x,z,x,&r); 
    porta_ldiv(*y,z,y,&r); 
  }

}
   
   







loint lgcdrow( loint *x, int m )
{

  int im;
  loint r,q;

  while (m > 0) {

    porta_ldiv(x[0],x[1],&q,&r);
    if (r.len > 1 || (r.len == 1 && r.val[0] > 1)) {
      x[0] = x[1];
      for(im = 1; im < m && lord(x[im+1],r) > 0; im++)
        x[im] = x[im+1];
      if (m > 1 && im < m && !lord(x[im+1],r)) {
        for (;im < m; im++) 
          x[im] = x[im+1];
        m--;
      }
      else
        x[im] = r;
    }
    else if (r.len == 0) {
      for(im = 0; im < m; im++)
        x[im] = x[im+1];
      m--;
    }
    else return(r);
  }

  return(x[0]);

}








void L_RAT_row_prim( RAT *old, RAT *new, RAT *p, int n )
{
  RAT rgcd, *pold;
  loint *x,r;
  int m=0,ii,i,j;
  lorat lgcd;

  x = (loint *) allo(CP NULL,0,U n*sizeof(loint));
  rgcd.num = 0;

  if (p->num != 0) {
    L_RAT_to_lorat(*p,&lgcd);
    r = lgcd.num; 
    x[0] = lgcd.den; lgcd.num = x[0];
    lgcd.den = r;
    lgcd.den.neg = 0;
    lorat_to_L_RAT(lgcd,&rgcd);
    for (i = n; i; i--)
      L_RAT_mul(*old++,rgcd,new++);
  }

  else 
  {

      lgcd.num.neg = lgcd.den.neg = lgcd.num.len = lgcd.den.len = 0;
      
      /* GCD OF NUMERATOR */
      
      for (j = 0,pold = old; j < n; j++,pold++) 
      {
          expand(pold->den.p,get_numlen(pold),&r);
          if (r.len == 1 && r.val[0] == 1) 
          { 
              lgcd.den.len = lgcd.den.val[0] = 1;
              break;
          }
          else if (r.len ) 
          {
              for (i = 0; i != m && lord(r,x[i]) < 0; i++); 
              if (lord(r,x[i]) || i == m) 
              {
                  for (ii = m; ii > i; ii--) 
                      x[ii] = x[ii-1];
                  x[i] = r;
                  m++;
              }
          }
      }
      m--;
      if (!lgcd.den.len ) lgcd.den = lgcdrow(x,m);
      
      /* GCD OF DENUMERATOR */
      
      m = 0;
      for (j = 0,pold = old; j < n; j++,pold++) 
      {
          if (pold->num) 
          {
              expand(pold->den.p+get_numlen(pold),get_denlen(pold),&r);
              if (r.len == 1 && r.val[0] == 1) 
              { 
                  lgcd.num.len = lgcd.num.val[0] = 1;
                  break;
              }
              else 
              {
                  for (i = 0; i != m && lord(r,x[i]) < 0; i++); 
                  if (lord(r,x[i]) || i == m) 
                  {
                      for (ii = m; ii > i; ii--) 
                          x[ii] = x[ii-1];
                      x[i] = r;
                      m++;
                  }
              }
          }
      }
      m--;
      if (!lgcd.num.len) lgcd.num = lgcdrow(x,m);
      
      if (lgcd.num.len > 1 || (lgcd.num.len == 1 && lgcd.num.val[0] > 1)
          || lgcd.den.len > 1 || (lgcd.den.len == 1 && lgcd.den.val[0] > 1)) 
      {
          lorat_to_L_RAT(lgcd,&rgcd);
          for (j = 0; j < n; j++) 
              L_RAT_mul(rgcd,old[j],new+j);
      }
      
      else if (old != new)  
          for (j = 0; j < n; j++) 
              L_RAT_assign(new+j,old+j);
      
  }
  
  x = (loint *) allo(CP x,U n*sizeof(loint),0);
  
}


// This function was added by J-D Bancal on 29.10.2014 to allow
// removal of fractions in presence of large integer fractions
//
// THIS FUNCTION NEEDS TO BE ADAPTED TO L_RAT TYPE AND lorat. THEN
// IT WILL BE POSSIBLE TO CALL IT IN xporta.c AFTER THE FOURIER-
// MOTZKIN ELIMINATION...
//
int L_RAT_no_denom( int sysrow, int first, int last, int outmsg )
/*
 * Make a fractional point integer by multiplying it with a positive number
 * that is as small as possible.
 */
{
	fprintf(prt, "Arrived here with sysrow = %i, first = %i, last = %i, outmsg = %i\n", sysrow, first, last, outmsg);
    int ret=1,i,j,ie;//,*denom;
    int t;
    
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
        
//        qsort(CP denom,sysrow,sizeof(intDenom), // It shouldn't be necessary to order the denominators...
//              (intDenom(*)(const void*,const void*))intcompare);
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



void hexprint( FILE *fp, loint lx )
{
    int i;
    
    if (!lx.len)
        fprintf(fp,"0");
    else  {
    if (lx.neg)
        fprintf(fp,"-");
    if (lx.len == 1 && !(lx.val[0] & m_31))
        fprintf(fp,"%u",lx.val[0]);
    else {
    fprintf(fp,"(hex)");
    for(i = lx.len-1; i >= 0; i--)
        fprintf(fp,"%x%x%x",15&(lx.val[i]>>8),15&(lx.val[i]>>4),15&(lx.val[i]));
}
}
}


void hexprintWithOppositeSign( FILE *fp, loint lx )
{
    int i;
    
    if (!lx.len)
        fprintf(fp,"0");
    else  {
    if (!lx.neg)
        fprintf(fp,"-");
    if (lx.len == 1 && !(lx.val[0] & m_31))
        fprintf(fp,"%u",lx.val[0]);
    else {
    fprintf(fp,"(hex)");
    for(i = lx.len-1; i >= 0; i--)
        fprintf(fp,"%x%x%x",15&(lx.val[i]>>8),15&(lx.val[i]>>4),15&(lx.val[i]));
}
}
}






void size_info( RAT *x, int *nz, int *ld )
{
    int l,nl,dl,b; 
    unsigned y;
#define lb(x) for(l = 0,y = x; y; y >>= 1,l++)
    
    if (x->num) 
    {
        (*nz)++;
        if (MP_realised) 
        {
            nl = get_numlen(x)-1;
            lb((int) x->den.p[nl]);
            if ( (b=l+nl*32) > *ld)
                (*ld) = b;
            dl = get_denlen(x);
            lb((int) x->den.p[nl+dl]);
            if ( (b=l+(dl-1)*32) > *ld)
                (*ld) = b;
        }
        else 
        {
            lb(abs(x->num));
            if (l > *ld)
                (*ld) = l;
            lb(abs(x->den.i));
            if (l > *ld)
                (*ld) = l;
        }
    }
    
}









void L_RAT_writeline( FILE *fp, int rowl, RAT *ptr, int format, RAT* max, char eqie, int *indx )
{
    int j;
    lorat x;
    
	// This function was updated to support a new table format (when format=2), comparable to
	// the one obtained in .ine files.
	
	if (format != 2) // This is just the usual code for format=0 or format=1
		for ( j = 0; j < rowl;ptr++,j++) 
		{
		    if (j == rowl-1 && !format)
		        fprintf(fp," %c= ",eqie);
		    if (!ptr->num && (format || (!format && j == rowl-1)))
		        fprintf(fp,"0");
		    else if (ptr->num && vals_lt_MAXINT(ptr,1)) 
		    {
		        if (ptr->num < 0)
		            fprintf(fp,"-");
		        else if (ptr->num > 0 && ! format)
		            fprintf(fp,"+");
		        fprintf(fp,"%u",ptr->den.p[0]);
		        if ( ptr->den.p[1] > 1) 
		            fprintf(fp,"/%u",ptr->den.p[1]);
		    }       
		    else  if (ptr->num) 
		    {
		        L_RAT_to_lorat(*ptr,&x);
		        if (ptr->num > 0 && ! format)
		            fprintf(fp,"+");
		        hexprint(fp,x.num);
		        if (x.den.len > 1 || x.den.val[0] > 1) 
		        {
		            fprintf(fp,"/");
		            hexprint(fp,x.den);
		        }  
		    }
		    if (format)
		        fprintf(fp," ");
		    else if (j != rowl-1 && ptr->num)
		        fprintf(fp,"x%i",(indx)?(indx[j]+1):j+1);        
		}
	else // Now we define the new table format
	{
		// First we write the constant. For this, we increase the ptr as much as needed
		for ( j = 0; j < rowl-1;ptr++,j++) {};

	    if (!ptr->num)
	        fprintf(fp,"0");
	    else if (ptr->num && vals_lt_MAXINT(ptr,1)) 
	    {
	        if (ptr->num < 0)
	            fprintf(fp,"-");
	        fprintf(fp,"%u",ptr->den.p[0]);
	        if ( ptr->den.p[1] > 1) 
	            fprintf(fp,"/%u",ptr->den.p[1]);
	    }       
	    else  if (ptr->num) 
	    {
	        L_RAT_to_lorat(*ptr,&x);
	        hexprint(fp,x.num);
	        if (x.den.len > 1 || x.den.val[0] > 1) 
	        {
	            fprintf(fp,"/");
	            hexprint(fp,x.den);
	        }  
	    }
        fprintf(fp," ");

		// Now we get back to the beginning of the list to obtain the coefficients...
		for ( j = 0; j < rowl-1;ptr--,j++) {};

		for ( j = 0; j < rowl-1;ptr++,j++) 
		{
		    if (!ptr->num)
		        fprintf(fp,"0");
		    else if (ptr->num && vals_lt_MAXINT(ptr,1)) 
		    {
//				fprintf(prt, "1");
		        if (-(ptr->num) < 0)
		            fprintf(fp,"-");
		        fprintf(fp,"%u",ptr->den.p[0]);
		        if ( ptr->den.p[1] > 1) 
		            fprintf(fp,"/%u",ptr->den.p[1]);
		    }       
		    else  if (ptr->num) 
		    {
//				fprintf(prt, "2");
		        L_RAT_to_lorat(*ptr,&x);
		        hexprintWithOppositeSign(fp,x.num);
		        if (x.den.len > 1 || x.den.val[0] > 1) 
		        {
		            fprintf(fp,"/");
		            hexprint(fp,x.den);
		        }  
		    }
	        fprintf(fp," ");
		}
		ptr++; j++;
	}
}



