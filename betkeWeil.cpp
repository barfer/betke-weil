// Betke-Weil Rigorous Computation
// 
// Peak Computations
//
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// std
#include <algorithm>
#include <ctime>
#include <iomanip>

// pugi
#include "pugixml.hpp"

// local
#include "main.h"
#include "betkeWeil.h"
#include "tools.h"
#include "jetTools.h"

//-----------------------------------------

/** Global singleton configuration */
ProblemConfiguration* betkeWeil = new ProblemConfiguration;

//-----------------------------------------

/** Parse the config file */
bool parseConfig(ProblemConfiguration* configuration, std::string configFile) {

	// Load configuration file
	pugi::xml_document config;
	if (!config.load_file(configFile.c_str())) return false;

	configuration->printingPrecision =
		config.child("Problem").child("Printing").attribute("Precision").as_int();

	// Log files
	configuration->logFile.open(
		config.child("Problem").child("Log").attribute("Output").as_string(),
		std::ios::out
	);

	configuration->debugFile.open(
		config.child("Problem").child("Log").attribute("Debug").as_string(),
		std::ios::out
	);

	configuration->logFile << std::setprecision(16);
	configuration->debugFile << std::setprecision(16);

	// Problem Definition
	std::string problemDefinition = config.child("Problem").child("ProblemDefinition").attribute("Objective").as_string();

	if (problemDefinition.compare("normComparison") == 0) {

		configuration->objective = ProblemDefinition::normComparison;

	} else if (problemDefinition.compare("quadraticFormComparison") == 0) {

		configuration->objective = ProblemDefinition::quadraticFormComparison;

	} else {

		configuration->objective = ProblemDefinition::directInequality;

	}
	
	return true;
}

//-----------------------------------------

/** Betke-Weil Computations */
int main() {

	// read in technical parameters
	parseConfig(betkeWeil, "betkeWeil.xml");

	// order of expansions
	const int orderOfExpansions = 4;

	// try: catch errors
	try {

		// Directly check f(x) >= ||x - z_t||^2
		if (betkeWeil->objective == ProblemDefinition::directInequality) {

			// -------------
			// f

			// 5 variables: a2, a3, t1, t2, t3
			// 1 dim. output
			// 2 parameters: a1, sqrt(3)
			// up to (orderOfExpansions + 2) order differentiation
			MapType f(BetkeWeil::f, 5, 1, 2, orderOfExpansions + 2);

			// 2 parameters
			ScalarType
				a1(2),
				value_sqrt3(sqrt(ScalarType(3)));

			f.setParameter(0, a1);
			f.setParameter(1, value_sqrt3);

			// result of f
			capd::IJet fJet(          1, 5, orderOfExpansions + 2);
			capd::IJet fJetAtMidpoint(1, 5, orderOfExpansions + 2);

			// -------------
			//  (queue) list of subsets of W to check during computations
			std::list<VectorType> domainsToCheck;

			// (a2, a3, t1, t2, t3) \in W
			ScalarType a2, a3, t1, t2, t3;

			a2 = 2 + (ScalarType(0, 1) / ScalarType(6));                       // 2 + [0, 1]/6     = [2, 2 + 1/6]
			a3 = 2 + (ScalarType(0, 1) / ScalarType(6));                       // 2 + [0, 1]/6     = [2, 2 + 1/6]
			t1 =     (ScalarType(0, 1) / ScalarType(6));                       //     [0, 1]/6     = [0, 1/6]
			t2 =     (ScalarType(0, 1) / ScalarType(6));                       //     [0, 1]/6     = [0, 1/6]
			t3 =     (ScalarType(0, 1) / ScalarType(6));                       //     [0, 1]/6     = [0, 1/6]

			domainsToCheck.push_back(VectorType({ a2, a3, t1, t2, t3 }));

			// -------------
			// computation
		
			// we check the first element of the queue
			while (!domainsToCheck.empty()) {

				// pick working domain: x \subseteq W
				VectorType x = domainsToCheck.front();

				domainsToCheck.pop_front();

				Tools::logValue(betkeWeil->debugFile, "(a2, a3, t1, t2, t3)", x);
			
				// evaluate f
				f(x, fJet);

				// improve the result of the computation using centered form
				f(Tools::mid(x), fJetAtMidpoint);

				JetType fJetImproved =
					JetTools::jetsToImprovedJet(
						fJet,
						x,
						fJetAtMidpoint,
						Tools::mid(x)
					);

				ScalarType fx = fJetImproved(fJetImproved.first(0))[0];

				// compute ||x - z_t||
				VectorType xMinusZt = x;

				xMinusZt[0] -= ScalarType(2);
				
				xMinusZt[1] -= ScalarType(2);
				
				xMinusZt[2] *= ScalarType(2) / ScalarType(3); // t1 - t1 / 3 = t1 * 2/3
				xMinusZt[2] -= (x[3] + x[4]) / ScalarType(3);

				xMinusZt[3] *= ScalarType(2) / ScalarType(3); // t2 - t2 / 3 = t2 * 2/3
				xMinusZt[3] -= (x[2] + x[4]) / ScalarType(3);

				xMinusZt[4] *= ScalarType(2) / ScalarType(3); // t3 - t3 / 3 = t3 * 2/3
				xMinusZt[4] -= (x[2] + x[3]) / ScalarType(3);

				ScalarType normXminusZt = xMinusZt.euclNorm();

				Tools::logValue(betkeWeil->debugFile, "x - z_t",     xMinusZt,     "\\subseteq");
				Tools::logValue(betkeWeil->debugFile, "f(x)",        fx,           "\\subseteq");
				Tools::logValue(betkeWeil->debugFile, "||x - z_t||", normXminusZt, "\\subseteq");

				// the inequality holds, processing is over for x
				if ( ( fx >= sqr(normXminusZt) ) ) {

					Tools::log(betkeWeil->debugFile, "result: success!");
					Tools::log(betkeWeil->debugFile, "----------------------------");

					Tools::logValue(betkeWeil->logFile, "(a2, a3, t1, t2, t3)", x);
					Tools::logValue(betkeWeil->logFile, "f(x)",                 fx,           "\\subseteq");
					Tools::logValue(betkeWeil->logFile, "||x - z_t||",          normXminusZt, "\\subseteq");
					Tools::log(     betkeWeil->logFile, "f(x) >= ||x - z_t||^2 holds.");
					Tools::log(     betkeWeil->logFile, "----------------------------");

				// the inequality could not be verified 
				// we cut x into two parts
				// and add these new parts into the queue for repeated analysis
				} else {
					Tools::log(betkeWeil->debugFile, "result: no success, subdividing...");
					Tools::log(betkeWeil->debugFile, "----------------------------");

					// find the widest component of x
					int idxWidestComponent;

					// if the maximal width is less than 0.001 then we report a failure
					if (Tools::maxDiam(x, idxWidestComponent) < 0.001) {

						Tools::log(betkeWeil->debugFile, "subdivision: widest component is smaller than tolerance, aborting the computation...");

						throw std::runtime_error("Subdivision is too small.");
					}

					ScalarType widestComponent = x[idxWidestComponent];

					// midpoint of the widest component
					double widestComponentMid =
						(
							(widestComponent.right() + widestComponent.left()) / ScalarType(2)
						).rightBound();

					// put "right" half of x into the queue
					x[idxWidestComponent] = ScalarType(widestComponentMid, widestComponent.rightBound());
					domainsToCheck.push_front(x);

					// put "left" half of x into the queue
					x[idxWidestComponent] = ScalarType(widestComponent.leftBound(), widestComponentMid);
					domainsToCheck.push_front(x);
				}

			} // end while (processing the queue)

		// Compute the Hessian and check the inequality either for the norms or for the quadratic form
		} else {

			// -------------
			// f~

			// 5 variables: y1, y2, y3, y4, y5
			// 1 dim. output
			// 3 parameters: a1, sqrt(3), sqrt(2)
			// up to (orderOfExpansions + 2) order differentiation
			MapType ftilde(BetkeWeil::ftilde, 5, 1, 3, orderOfExpansions + 2);

			// 3 parameters
			ScalarType
				a1(2),
				value_sqrt3(sqrt(ScalarType(3))),
				value_sqrt2(sqrt(ScalarType(2)));

			ftilde.setParameter(0, a1);
			ftilde.setParameter(1, value_sqrt3);
			ftilde.setParameter(2, value_sqrt2);

			// result of ftilde
			capd::IJet ftildeJet(          1, 5, orderOfExpansions + 2);
			capd::IJet ftildeJetAtMidpoint(1, 5, orderOfExpansions + 2);

			// -------------
			//  (queue) list of subsets of (W~ \times \Xi) to check during computations
			std::list<VectorType> domainsToCheck;

			// y \in W~
			ScalarType y1, y2, y3, y4, y5;

			y1 = 2 + (ScalarType(0, 1) / ScalarType(6));                       // 2 + [0, 1]/6     = [2, 2 + 1/6]
			y2 = 2 + (ScalarType(0, 1) / ScalarType(6));                       // 2 + [0, 1]/6     = [2, 2 + 1/6]
			y3 =     (ScalarType(14) / ScalarType(100)) * ScalarType(-1, 1);   // 14/100 * [-1, 1] = [-0.14, 0.14]
			y4 =     (ScalarType(14) / ScalarType(100)) * ScalarType(-1, 1);   // 14/100 * [-1, 1] = [-0.14, 0.14]
			y5 =     (ScalarType(0, 3) / ScalarType(10));                      // [0, 3]/10        = [0, 0.3]

			// bar{v} \in \Xi

			// bar{v}(1) = 1, other components are varying within [-1, 1]
			ScalarType
				v1(1),
				v2(-1, 1),
				v3(-1, 1),
				v4(-1, 1);

			domainsToCheck.push_back(
				VectorType({ y1, y2, y3, y4, y5, v1, v2, v3, v4 })
			);

			// bar{v}(2) = 1, other components are varying within [-1, 1]
			v1 = ScalarType(-1, 1);
			v2 = ScalarType(1);

			domainsToCheck.push_back(
				VectorType({ y1, y2, y3, y4, y5, v1, v2, v3, v4 })
			);

			// bar{v}(3) = 1, other components are varying within [-1, 1]
			v2 = ScalarType(-1, 1);
			v3 = ScalarType(1);

			domainsToCheck.push_back(
				VectorType({ y1, y2, y3, y4, y5, v1, v2, v3, v4 })
			);

			// bar{v}(4) = 1, other components are varying within [-1, 1]
			v3 = ScalarType(-1, 1);
			v4 = ScalarType(1);

			domainsToCheck.push_back(
				VectorType({ y1, y2, y3, y4, y5, v1, v2, v3, v4 })
			);

			// -------------
			// computation

			// we check the first element of the queue
			while (!domainsToCheck.empty()) {

				// pick working domain
				VectorType currentJointDomain = domainsToCheck.front();

				domainsToCheck.pop_front();

				// select y \in W~ and bar{v} \in \Xi
				VectorType y(5), v(4);

				// y
				for (int idx = 0; idx < 5; ++idx) {

					y[idx] = currentJointDomain[idx];
				}

				// bar{v}
				for (int idx = 0; idx < 4; ++idx) {

					v[idx] = currentJointDomain[5 + idx];
				}

				Tools::logValue(betkeWeil->debugFile, "(y1, y2, y3, y4, y5)", y);
				Tools::logValue(betkeWeil->debugFile, "(v1, v2, v3, v4)",     v);

				// evaluate f~
				ftilde(y, ftildeJet);

				// improve the result of the computation using centered form
				ftilde(Tools::mid(y), ftildeJetAtMidpoint);

				JetType ftildeJetImproved =
					JetTools::jetsToImprovedJet(
						ftildeJet,
						y,
						ftildeJetAtMidpoint,
						Tools::mid(y)
					);

				// get the Hessian
				JetType hessianJet = JetTools::jetToHessianJet(ftildeJetImproved);

				// Hessian matrix
				capd::IMatrix Hessian(5, 5); // H[i][j] = (i + 1)th row, (j + 1)th column

				int idx = 0;

				for (int row = 0; row < 5; ++row) {
					for (int col = row; col < 5; ++col) {

						Hessian[row][col] = hessianJet(hessianJet.first(0))[idx];
						Hessian[col][row] = hessianJet(hessianJet.first(0))[idx]; // derivatives are stored with symmetry

						++idx;
					}
				}

				// principal minor of the Hessian
				capd::IMatrix H(4, 4);

				for (int row = 0; row < 4; ++row) {
					for (int col = 0; col < 4; ++col) {

						H[row][col] = Hessian[row][col];
					}
				}

				Tools::logValue(betkeWeil->debugFile, "Hessian",                        Hessian, "\\subseteq");
				Tools::logValue(betkeWeil->debugFile, "principal minor of the Hessian", H,       "\\subseteq");

				// compute the norms
				VectorType Hv = H * v;

				ScalarType normHv, normV;

				normHv = Hv.euclNorm();
				normV  =  v.euclNorm();

				Tools::logValue(betkeWeil->debugFile, "||Hv||", normHv, "\\subseteq");
				Tools::logValue(betkeWeil->debugFile, "||v|| ", normV,  "\\subseteq");

				// compute the quadratic form
				ScalarType vHv = v * Hv;

				Tools::logValue(betkeWeil->debugFile, "<v, H v>", vHv, "\\subseteq");

				// the inequality holds, processing is over for this domain of (y \times v) \subset (W~ \times \Xi)
				if (
					((betkeWeil->objective == ProblemDefinition::normComparison)          && (normHv >= ScalarType(2) * normV     )) ||
					((betkeWeil->objective == ProblemDefinition::quadraticFormComparison) && (vHv    >= ScalarType(2) * sqr(normV)))
					) {

					Tools::log(betkeWeil->debugFile, "result: success!");
					Tools::log(betkeWeil->debugFile, "----------------------------");

					Tools::logValue(betkeWeil->logFile, "(y1, y2, y3, y4, y5)", y);
					Tools::logValue(betkeWeil->logFile, "(v1, v2, v3, v4)",     v);
					Tools::logValue(betkeWeil->logFile, "||Hv||",               normHv, "\\subseteq");
					Tools::logValue(betkeWeil->logFile, "||v|| ",               normV,  "\\subseteq");
					Tools::logValue(betkeWeil->logFile, "<v, H v>",             vHv,    "\\subseteq");
					Tools::log(betkeWeil->logFile,
						(betkeWeil->objective == ProblemDefinition::normComparison) ?
						"||Hv||  >= 2 ||v||   holds." :
						"<v, Hv> >= 2 ||v||^2 holds."
					);
					Tools::log(betkeWeil->logFile, "----------------------------");

					// the inequality could not be verified 
					// we cut this domain (y \times v) \subset (W~ \times \Xi) into two parts
					// and add these new parts into the queue for repeated analysis
				}
				else {
					Tools::log(betkeWeil->debugFile, "result: no success, subdividing...");
					Tools::log(betkeWeil->debugFile, "----------------------------");

					// find the widest component of (y \times v)
					int idxWidestComponent;

					// if the maximal width is less than 0.001 then we report a failure
					if (Tools::maxDiam(currentJointDomain, idxWidestComponent) < 0.001) {

						Tools::log(betkeWeil->debugFile, "subdivision: widest component is smaller than tolerance, aborting the computation...");

						throw std::runtime_error("Subdivision is too small.");
					}

					ScalarType widestComponent = currentJointDomain[idxWidestComponent];

					// midpoint of the widest component
					double widestComponentMid =
						(
							(widestComponent.right() + widestComponent.left()) / ScalarType(2)
							).rightBound();

					// put "right" half of (y \times v) into the queue
					currentJointDomain[idxWidestComponent] = ScalarType(widestComponentMid, widestComponent.rightBound());
					domainsToCheck.push_front(currentJointDomain);

					// put "left" half of (y \times v) into the queue
					currentJointDomain[idxWidestComponent] = ScalarType(widestComponent.leftBound(), widestComponentMid);
					domainsToCheck.push_front(currentJointDomain);
				}

			} // end while (processing the queue)

		} // end if (problem definition: which inequality to check)

	// error handling
	} catch (std::exception& e) {

		std::cout << "Exception caught: " << e.what() << std::endl;

		betkeWeil->logFile << std::endl << "Exception caught: " << e.what() << std::endl;

		betkeWeil->logFile.close();
		betkeWeil->debugFile.close();
	}

	// the queue is empty, processing is done
	betkeWeil->logFile.close();
	betkeWeil->debugFile.close();

	// program is done
	return 0;
}
