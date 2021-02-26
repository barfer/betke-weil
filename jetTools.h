#pragma once
// Jet Tools
//
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

#include "main.h"

//-----------------------------------------

namespace JetTools {

	/** Transforms <fJet> the jet (at least degree 2) of a map f: R^n -> R into
	the jet (two degrees less than of <fJet>) of (part of) the hessian Hf: R^n -> d^f/dx_i dx_j with i <= j
	- the remaining components of the hessian are obtainable using
	  the symmetry of the matrix i.e. d^f/dx_i dx_j = d^f/dx_j dx_i
*/
	JetType jetToHessianJet(JetType& fJet) {

		if (fJet.imageDimension() != 1) {

			throw std::runtime_error("jetToHessianJet: The dimension of f's image has to be 1.");
		}

		if (fJet.degree() < 2) {

			throw std::runtime_error("jetToHessianJet: The degree of fJet has to be at least 2.");
		}

		JetType hJet(
			(fJet.dimension() * (fJet.dimension() + 1)) / 2,	// n(n+1)/2 possibly unique entries in the hessian
			fJet.dimension(),									// n variables as in f
			fJet.degree() - 2									// degree is two less than of f
		);

		// second order derivative of f
		capd::Multipointer indexOfSecondDerivF = fJet.first(2);

		// index to value of the hessian
		int indexToValueH = 0;

		// process all elements of the hessian
		do {

			// partial derivatives of the value in the hessian up to the proper degree
			for (int order = 0; order <= hJet.degree(); ++order) {

				// the element in the jet of the hessian (a vector of given partial derivatives for all component functions)
				capd::Multipointer indexOfDerivH = hJet.first(order);

				do {

					// the corresponding element in the jet of f
					capd::Multipointer indexOfDerivF =
						capd::vectalg::sumMultipointers(indexOfSecondDerivF, indexOfDerivH);

					// we need to rescale the partial derivative of f
					ScalarType scalingFactor =
						ScalarType(indexOfDerivF.factorial()) /
						ScalarType(indexOfDerivH.factorial());

					//std::cout << indexOfDerivF << " = " << indexOfSecondDerivF << " + " << indexOfDerivH <<  "   " << scalingFactor << " " << indexToValueH << std::endl;

					// f is a scalar, so each partial derivative is a 1D vector
					hJet(indexOfDerivH)[indexToValueH] = scalingFactor * fJet(indexOfDerivF)[0];

					//std::cout << hJet(indexOfDerivH) << std::endl;

				} while (hJet.hasNext(indexOfDerivH)); // next element in the jet of the hessian of the same degree
			}

			++indexToValueH; // go to next element in hessian

		} while (fJet.hasNext(indexOfSecondDerivF)); // go to next element in hessian

		return hJet;
	}

	/** Transforms <jet> representing f: R^n -> R^m with
		inverse information based on the value <x> \in R^n f was evaluated at
		- the result will be a jet with additional n image dimensions in the front,
		  these represent the identity that is, the returned jet is of (Id, f): R^n -> R^(n+m)
	*/
	JetType jetToVariableEnrichedJet(JetType& jet, VectorType x) {

		if (jet.dimension() != x.dimension()) {

			throw std::runtime_error("jetToVariableEnrichedJet: The dimension of the jet's domain and of x has to be the same.");
		}

		JetType enrichedJet(jet.dimension() + jet.imageDimension(), jet.dimension(), jet.degree());

		// process all partial derivatives
		for (int order = 0; order <= enrichedJet.degree(); ++order) {

			capd::Multipointer indexOfDeriv = enrichedJet.first(order);

			do {

				// complete the jet of Identity
				for (int componentInRn = 0; componentInRn < x.dimension(); ++componentInRn) {

					// higher order derivatives are zero
					if (order > 1) {

						enrichedJet(indexOfDeriv)[componentInRn] = 0;

					// first order derivatives
					} else if (order == 1) {

						// derivative is one w.r.t. the itself
						enrichedJet(indexOfDeriv)[componentInRn] = (indexOfDeriv[0] == componentInRn ? 1 : 0);

					// value
					} else {

						enrichedJet(indexOfDeriv)[componentInRn] = x[componentInRn];
					}
				}

				// copy the original jet elements to their new place
				for (int componentInRm = 0; componentInRm < jet.imageDimension(); ++componentInRm) {

					enrichedJet(indexOfDeriv)[jet.dimension() + componentInRm] = jet(indexOfDeriv)[componentInRm];
				}

			} while (enrichedJet.hasNext(indexOfDeriv)); // process next partial derivative of the same degree
		}

		return enrichedJet;
	}

	/** Transforms a <jet> of f:R^n -> R^m with m >= n to
		a jet f~:R^m -> R^m with all new entries being zero i.e.
		partial derivatives w.r.t. the artificially created new variables
	*/
	JetType jetToBloatedSquareJet(JetType& jet) {

		if (jet.dimension() > jet.imageDimension()) {

			throw std::runtime_error("jetToBloatedSquareJet: The dimension of the jet's domain cannot be larger than the dimension of its image.");
		}

		// if the jet is already of square shape, there is nothing to do
		if (jet.dimension() == jet.imageDimension()) {

			return JetType(jet);
		}

		JetType bloatedJet(jet.imageDimension(), jet.imageDimension(), jet.degree());

		// process all partial derivatives in jet (others are initialized to zero)
		for (int order = 0; order <= jet.degree(); ++order) {

			capd::Multipointer indexOfDeriv = jet.first(order);

			// copy all partial derivatives of a given degree
			do {

					bloatedJet(indexOfDeriv) = jet(indexOfDeriv);

			} while (jet.hasNext(indexOfDeriv)); // process next partial derivative of the same degree
		}

		return bloatedJet;
	}

	/** Transforms two jets (<thickJet> and <thinJet>) of the function f:R^n -> R^m into a new jet.
		- <thickJet> is valid for f(<thickX>)
		- <thinJet>  is valid for f(<thinX>) with <thinX> \subseteq <thickX> \in R^n
		The result is an improved bound for Df over <thickX>, where the
		derivative operator is given by <derivative>.
		The procedure sums the corresponding multivariate Taylor expansions with remainders.
	*/
	VectorType taylorSumJetsAtDerivative(
		JetType& thickJet, VectorType thickX,
		JetType& thinJet, VectorType thinX,
		capd::Multipointer indexOfObjectiveDeriv
	) {

		if (thickX.dimension() != thinX.dimension()) {

			throw std::runtime_error("taylorSumJetsAtDerivative: thickX and thinX must be of the same dimension.");

		} else if (thickX.dimension() != thickJet.dimension()) {

			throw std::runtime_error("taylorSumJetsAtDerivative: thickX and thickJet's domain must be of the same dimension.");

		} else if (
			(thickJet.dimension()      != thinJet.dimension()      ) ||
			(thickJet.imageDimension() != thinJet.imageDimension() ) ||
			(thickJet.degree()         != thinJet.degree()         )) {

			throw std::runtime_error("taylorSumJetsAtDerivative: thickJet and thinJet must be of the same shape.");

		} else {

			try {

				thickJet(indexOfObjectiveDeriv);

			} catch (std::runtime_error& e) {

				throw std::runtime_error(
					std::string("taylorSumJetsAtDerivative: The objective derivative has to be present in the jets.\n") +
					std::string(e.what())
				);
			}

			for (int component = 0; component < thickX.dimension(); ++component) {

				if (!(thickX[component].contains(thinX[component]))) {

					throw std::runtime_error("taylorSumJetsAtDerivative: thickX must contain thinX.");
				}
			}
		}

		// we build the multivariate Taylor expansion, this jet is just used to obtain the indices of derivatives
		JetType taylorExpansionJet(
			thickJet.imageDimension(),
			thickJet.dimension(),
			thickJet.degree() - indexOfObjectiveDeriv.dimension()	// the degree of the expansion is limited by the degree of the original jet and the order (dimension) of the sought derivative
		);

		// appropriately sized result
		VectorType result(thickJet.imageDimension());

		// relative domain of the expansion
		VectorType dX = thickX - thinX;

		// sum the Taylor expansion
		for (int order = 0; order <= taylorExpansionJet.degree(); ++order) {

			capd::Multipointer indexOfTaylorElement = taylorExpansionJet.first(order);

			// process elements of a given order
			do {

				// the corresponding element in the other jets
				capd::Multipointer indexOfDeriv = capd::vectalg::sumMultipointers(indexOfObjectiveDeriv, indexOfTaylorElement);

				// we need to rescale the partial derivative
				ScalarType scalingFactor =
					ScalarType(indexOfDeriv.factorial()) /
					ScalarType(indexOfTaylorElement.factorial()) /
					ScalarType(indexOfObjectiveDeriv.factorial());

				result = result + 
					scalingFactor *
					(order == taylorExpansionJet.degree() ? // term in the Taylor expansion(if there are no more terms, this equals to the corresponding element in thickJet)
						thickJet(indexOfDeriv) : 
						thinJet(indexOfDeriv)
					) * 
					power(dX, indexOfTaylorElement);

			} while (taylorExpansionJet.hasNext(indexOfTaylorElement)); // process next element of the same order
		}

		// make sure we do not make things worse
		try {

			result = capd::vectalg::intersection(result, (VectorType)thickJet(indexOfObjectiveDeriv));

		// intersection must be non-empty
		} catch (std::runtime_error& e) {

			throw std::runtime_error(
				std::string("taylorSumJetsAtDerivative: The Taylor sum and the naive evaluation should have non-empty intersection!\n") +
				std::string(e.what())
			);
		}

		return result;
	}

	/** Improves two jets (<thickJet> and <thinJet>) of the function f:R^n -> R^m into a new jet.
		- <thickJet> is valid for f(<thickX>)
		- <thinJet>  is valid for f(<thinX>) with <thinX> \subseteq <thickX> \in R^n
		The result is an improved bound for Df over <thickX> for all permissable D derivatives.
		The procedure sums the corresponding multivariate Taylor expansions with remainders.
	*/
	JetType jetsToImprovedJet(
		JetType& thickJet, VectorType thickX,
		JetType& thinJet, VectorType thinX
	)
	{
		// There's no need to check for errors here as all computation is delegated

		JetType result(thickJet.imageDimension(), thickJet.dimension(), thickJet.degree());

		// for all elements of the jets
		for (int order = 0; order <= result.degree(); ++order) {

			capd::Multipointer indexOfDeriv = result.first(order);

			// for all elements of a given order
			do {

				result(indexOfDeriv) = JetTools::taylorSumJetsAtDerivative(thickJet, thickX, thinJet, thinX, indexOfDeriv);

			} while (result.hasNext(indexOfDeriv)); // process next element of the same order
		}

		return result;
	}

	JetType joinJets(JetType& first, JetType& second) {

		if (first.dimension() != second.dimension()) {
			std::cout << first.dimension() << " " << second.dimension() << std::endl;
			throw std::runtime_error("joinJets: jets of equal dimension required.");
		
		}
		else if (first.degree() != second.degree()) {

			throw std::runtime_error("joinJets: jets of equal degree required.");

		}

		JetType result(
			first.imageDimension() + second.imageDimension(),
			first.dimension(),
			first.degree()
		);

		for (int order = 0; order <= result.degree(); ++order) {

			capd::Multipointer idxToDeriv = result.first(order);

			do {

				for (int index = 0; index < first.imageDimension(); ++index) {

					result(idxToDeriv)[index] = first(idxToDeriv)[index];
				}

				for (int index = 0; index < second.imageDimension(); ++index) {

					result(idxToDeriv)[first.imageDimension() + index] = second(idxToDeriv)[index];
				}

			} while (result.hasNext(idxToDeriv));

		}

		return result;
	}
	
}
