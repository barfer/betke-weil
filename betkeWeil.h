#pragma once
// Betke-Weil Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

#include "main.h"

//-----------------------------------------

namespace BetkeWeil {

	/** Betke-Weil f map: (a2, a3, t1, t2, t3) -> f
		- two parameters: a1, sqrt(3)
	*/
	void f(
		capd::autodiff::Node,
		capd::autodiff::Node in[],     // in[0] = a2, in[1] = a3, in[2] = t1, in[3] = t2, in[4] = t3
		int,
		capd::autodiff::Node out[],
		int,
		capd::autodiff::Node params[], // p[0] = a1, p[1] = sqrt(3)
		int
	) {
		// parameters
		capd::autodiff::Node a1    = params[0];
		capd::autodiff::Node sqrt3 = params[1];

		// variables: (a2, a3, t1, t2, t3) \in W
		capd::autodiff::Node a2 = in[0];
		capd::autodiff::Node a3 = in[1];
		capd::autodiff::Node t1 = in[2];
		capd::autodiff::Node t2 = in[3];
		capd::autodiff::Node t3 = in[4];

		// f_2
		capd::autodiff::Node f2 =
			(  a1 + a2 + a3) *
			(- a1 + a2 + a3) *
			(  a1 - a2 + a3) *
			(  a1 + a2 - a3);

		// f_1
		capd::autodiff::Node f1 =
			sqr(a1 + a2 + a3) + f2 * sqr(t1 / a1 + t2 / a2 + t3 / a3);

		// result
		out[0] = f1 - 3 * sqrt3 * sqrt(f2) * (1 + t1 * t2 + t2 * t3 + t3 * t1);
	}

	/** Betke-Weil f~ map: (y1, y2, y3, y4, y5) -> f~
		- three parameters: a1, sqrt(3), sqrt(2)
		- delegates to the map f after a coordinate transformation
	*/
	void ftilde(
		capd::autodiff::Node,
		capd::autodiff::Node in[],
		int,
		capd::autodiff::Node out[],
		int,
		capd::autodiff::Node params[], // p[0] = a1, p[1] = sqrt(3), p[2] = sqrt(2)
		int
	) {
		// parameters
		capd::autodiff::Node a1    = params[0];
		capd::autodiff::Node sqrt3 = params[1];
		capd::autodiff::Node sqrt2 = params[2];

		// variables: y -> y[i] = (i + 1)th element
		capd::vectalg::Vector<capd::autodiff::Node, 5> y;

		for (int idx = 0; idx < 5; ++idx) {

			y[idx] = in[idx];
		}

		// Orthogonal matrix: S -> S[i][j] = (i + 1)th row, (j + 1)th column
		capd::vectalg::Matrix<capd::autodiff::Node, 5, 5> S;
		
		// first column of S is e1 = (1, 0, 0, 0, 0)
		S[0][0] = capd::autodiff::Node(1);
		S[1][0] = capd::autodiff::Node(0);
		S[2][0] = capd::autodiff::Node(0);
		S[3][0] = capd::autodiff::Node(0);
		S[4][0] = capd::autodiff::Node(0);

		// second column of S is e1 = (0, 1, 0, 0, 0)
		S[0][1] = capd::autodiff::Node(0);
		S[1][1] = capd::autodiff::Node(1);
		S[2][1] = capd::autodiff::Node(0);
		S[3][1] = capd::autodiff::Node(0);
		S[4][1] = capd::autodiff::Node(0);

		// third column of S is e3 = (0, 0, 1/sqrt(2), -1/sqrt(2), 0)
		S[0][2] = capd::autodiff::Node(0);
		S[1][2] = capd::autodiff::Node(0);
		S[2][2] = capd::autodiff::Node( 1) / sqrt2;
		S[3][2] = capd::autodiff::Node(-1) / sqrt2;
		S[4][2] = capd::autodiff::Node(0);

		// fourth column of S is e4 = (0, 0, 1/sqrt(6), 1/sqrt(6), -2/sqrt(6))
		S[0][3] = capd::autodiff::Node(0);
		S[1][3] = capd::autodiff::Node(0);
		S[2][3] = capd::autodiff::Node( 1) / sqrt2 / sqrt3;
		S[3][3] = capd::autodiff::Node( 1) / sqrt2 / sqrt3;
		S[4][3] = capd::autodiff::Node(-2) / sqrt2 / sqrt3;

		// fifth column of S is e5 = (0, 0, 1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
		S[0][4] = capd::autodiff::Node(0);
		S[1][4] = capd::autodiff::Node(0);
		S[2][4] = capd::autodiff::Node(1) / sqrt3;
		S[3][4] = capd::autodiff::Node(1) / sqrt3;
		S[4][4] = capd::autodiff::Node(1) / sqrt3;

		// compute S * y
		capd::vectalg::Vector<capd::autodiff::Node, 5> Sy;

		Sy = S * y;

		// pass w = (a2, a3, t1, t2, t3) = S * y to f
		capd::autodiff::Node w[5];

		for (int idx = 0; idx < 5; ++idx) {
			
			w[idx] = Sy[idx];
		}
		
		// compute f(w)
		capd::autodiff::Node t; // not used, just needed for technical reasons

		BetkeWeil::f(t, w, 5, out, 1, params, 2);
	}

}
