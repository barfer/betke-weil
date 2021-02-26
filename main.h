#pragma once
// Betke-Weil Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// std
#include <stddef.h>
#include <iostream>
#include <iomanip>

// CAPD
#include <capd/capdlib.h>

//-----------------------------------------
// Numbers

/** Matrix Type */
typedef capd::IMatrix MatrixType;

/** Vector Type
	Technical necessity for the CAPD package that is capable to handle multi-dimensional maps.
	In our setting the maps are one-dimensional, hence each vector consists of one single component.
*/
typedef capd::IVector VectorType;

/** Scalar Type
	This is the rigorous representation of a value.
	We have the option to work with intervals with double endpoints (capd::DInterval),
	or intervals with higher precisional endpoints.
*/
typedef capd::DInterval ScalarType;

/** Map Type
	Rigorous function evaluation.
*/
typedef capd::IMap MapType;

/** Jet Type
	Rigorous storage and propagation of partial derivatives up to a certain degree.
*/
typedef capd::IJet JetType;

//-----------------------------------------

/** Problem Definition */
typedef enum ProblemDefinition {

	normComparison,
	quadraticFormComparison,
	directInequality

} ProblemDefinition;

/** Problem Configuration
	The configuration captures the precision treshold that is used for localising crossings.
*/
typedef struct ProblemConfiguration {

	/** Printing precision */
	size_t printingPrecision;

	/** Log file */
	std::ofstream logFile;

	/** Debug file */
	std::ofstream debugFile;

	/** Problem definition */
	ProblemDefinition objective;

} ProblemConfiguration;

// Export the problem configuration to be globally accesible
extern ProblemConfiguration* betkeWeil;
