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
 

FILENAME: libporta.h

AUTHOR: Jean-Daniel bancal

This file provides a c++ library interface to the Porta polytope software of
Thomas Christof and Andreas Loebel. Most of the code here was directly copied
from the porta.c, porta.h and xporta.c source files.

*******************************************************************************/

#ifndef _LIBPORTA_H
#define _LIBPORTA_H

#include <iostream>
#include <string>
#include <vector>
extern "C" { // These are c libraries.
#include "porta.h"
#include "log.h"
}

using namespace std;


/* The two following classes are meant to be easy ways to describe polytopes 
 * either in terms of half-planes (inequality and equality constraints), or
 * in terms of extremal vertices and rays. */
class Hrep;
class Vrep;


/* Let us define the types of numbers we want to use : the first is the
 * type for numerator coefficients, and the second one for denominator
 * coefficients. */
typedef long int Tnum;
typedef int Tden;


/* Let's define a fraction */
class fraction{
	public:
		Tnum num; // The numerator of the fraction
		Tden den; // The denominator
		
		fraction() : num(0), den(1) {}; // By default, a fraction is initialized to zero.

		bool operator==(const fraction& f) const;
		
		// Basic arithmetic
		fraction operator+ (const fraction& f) const;
		fraction operator* (const fraction& f) const;
		
		void simplify();
};


/* The following defines a V-representation of a polytope. */
class Vrep {
	public:
		/* V-representation constructors :
		 * We accept that V-representations be constructed
		 *  - empty,
		 *  - from a table of integer numbers (a la lrs)
		 *  - from two tables of numbers (for the numerator and denominators)
		 *  - from a *.ext lrs file
		 *  - from the variables produced by porta */
		Vrep() : pointsAndRays(0, vector < fraction > (0, fraction())) {};
		Vrep(const Tnum *data, const long int& nbPointsOrRays, const long int& dimension) {fill(data, nbPointsOrRays, dimension);};
		Vrep(const Tnum *data_num, const Tden *data_den, const long int& nbPointsOrRays, const long int& dimension){fill(data_num, data_den, nbPointsOrRays, dimension);};
		Vrep(const string& filename){fill(filename);};
		
		// Destructor :
		virtual ~Vrep() {};
		
		// Here are the functions that fill the table with different souces of data.
		void fill(const Tnum *data, const long int& nbPointsOrRays, const long int& dimension);
		void fill(const Tnum *data_num, const Tden *data_den, const long int& nbPointsOrRays, const long int& dimension);
		void fill(const string& filename);
		void fill(const listp* porta_list, const long int& dimension, const long int& nbEq, const long int& nbIneq);

		// The == operator
		/* WARNING : This operator checks if the "data structure" is identical for the two representations, so it can
		 * return false for two polytopes that are mathematically identical but listed in different orders, with additional
		 * non-extremal points, or with points defined with a different normalization...*/
		bool operator==(const Vrep& vrep) const;

		// V-representation display (a la lrs)
		friend ostream& operator<<(ostream& os, const Vrep& vrep);
		
		// This function converts to a H-representation
		friend void VtoHrep(const Vrep& vrep, Hrep& hrep);
		
		// Provides some useful info about the polytope
		long int nbPointsAndRays() const;
		long int dimension() const;
		
	protected:
		vector < vector < fraction > > pointsAndRays;
};


/* The following defines an H-representation of a polytope. */
class Hrep {
	public:
		/* H-representation Constructors :
		 * We accept that H-representations be constructed
		 *  - empty,
		 *  - from a table of integers representing inequalities (a la lrs)
		 *  - from two tables of integers representing numerator and denominators of inequalities coefficients
		 *  - from two tables of integer numbers for inequalities and equalities
		 *  - from four tables of numbers for the numerator and denominators of inequalities and equalities coefficients
		 *  - from a *.ine lrs file
		 *  - from the variables produced by porta */
		Hrep() : equalities(0, vector < fraction > (0, fraction())), inequalities(0, vector < fraction > (0, fraction())), validPoint(0, fraction()) {};
		Hrep(const Tnum *dataIneq, const long int& nbIneq, const long int& dimension) {fill(dataIneq, nbIneq, dimension);};
		Hrep(const Tnum *dataIneq_num, const Tden *dataIneq_den, const long int& nbIneq, const long int& dimension) {fill(dataIneq_num, dataIneq_den, nbIneq, dimension);};
		Hrep(const Tnum *dataEq, const long int& nbEq, const Tnum *dataIneq, const long int& nbIneq, const long int& dimension)
			{fill(dataEq, nbEq, dataIneq, nbIneq, dimension);};
		Hrep(const Tnum *dataEq_num, const Tden *dataEq_den, const long int& nbEq, const Tnum *dataIneq_num, const Tden *dataIneq_den, const long int& nbIneq, const long int& dimension)
			{fill(dataEq_num, dataEq_den, nbEq, dataIneq_num, dataIneq_den, nbIneq, dimension);};
		Hrep(const string &filename) {fill(filename);};
		Hrep(const listp* porta_list, const long int& dimension, const long int& nbEq, const long int& nbIneq, const int* indx)
			{fill(porta_list, dimension, nbEq, nbIneq, indx);};
		
		// Destructor
		virtual ~Hrep() {};

		// These functions fill the content of the table with data of given form
		void fill(const Tnum *dataIneq, const long int& nbIneq, const long int& dimension);
		void fill(const Tnum *dataIneq_num, const Tden *dataIneq_den, const long int& nbIneq, const long int& dimension);
		void fill(const Tnum *dataEq, const long int& nbEq, const Tnum *dataIneq, const long int& nbIneq, const long int& dimension);
		void fill(const Tnum *dataEq_num, const Tden *dataEq_den, const long int& nbEq, const Tnum *dataIneq_num, const Tden *dataIneq_den, const long int& nbIneq, const long int& dimension);
		void fill(const string& filename);
		void fill(const listp* porta_list, const long int& dimension, const long int& nbEq, const long int& nbIneq, const int* indx);
		
		// The == operator
		/* WARNING : This operator checks if the "data structure" is identical for the two representations, so it can
		 * return false for two polytopes that are mathematically identical but with inequalities listed in different orders,
		 * with redundant inequalities, inequalities defined with a different normalization...*/
		bool operator==(const Hrep& hrep) const;
		
		// H-representation display (a la lrs)
		friend ostream& operator<<(ostream& os, const Hrep& hrep);
		
		// This function converts to a V-representation
		friend void HtoVrep(const Hrep& hrep, Vrep& vrep);

		// Provides some useful info about the polytope
		long int nbEqualities() const;
		long int nbInequalities() const;
		long int dimension() const;

	protected:
		// This function tests if the validPoint satisfies all current constraints
		void checkValidPoint() const;
	
		vector < vector < fraction > > equalities;
		vector < vector < fraction > > inequalities;
		vector < fraction > validPoint;
};



/* The following added by J-D Bancal on 22.2.2012 to call porta as a c library. */
//extern void quick_and_dirty_poi_conv_call(long int *, int, int);


// Definition of the display functions
extern ostream& operator<<(ostream& os, const Vrep& vrep);
extern ostream& operator<<(ostream& os, const Hrep& hrep);



#endif // _LIBPORTA_H
