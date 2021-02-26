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
// Tools

/** Tools */
namespace Tools {

	ScalarType diam(ScalarType interval) {

		return interval.right() - interval.left();
	}

	ScalarType maxDiam(const VectorType& vector, int& index) {

		ScalarType result = 0.0;

		for (int i = 0; i < vector.dimension(); ++i)
			
			if (Tools::diam(vector[i]) > result) {
				
				result = Tools::diam(vector[i]);
				
				index = i;
			}

		return result;
	}

	VectorType mid(VectorType intervalVector) {

		VectorType result(intervalVector);

		for (int i = 0; i < result.dimension(); ++i) {
			
			result[i] = (result[i].left() + result[i].right()) / 2;
		}

		return result;
	}

	std::string stringify(bool logical) {

		std::stringstream buffer;

		logical ? buffer << "True" : buffer << "False";

		return buffer.str();
	}

	std::string stringify(ScalarType interval) {

		std::stringstream buffer;

		buffer << std::setprecision(betkeWeil->printingPrecision) << "[" << interval.leftBound() << ", " << interval.rightBound() << "]";

		return buffer.str();
	}

	template <typename StreamType>
	void log(StreamType& logStream, std::string text) {

		logStream << text << std::endl << std::flush;
	}

	template <typename StreamType>
	void logValue(StreamType& logStream, std::string text, VectorType value, std::string equalityType = "=") {

		logStream << text << " " << equalityType << " " << value << std::endl << std::flush;
	}

	template <typename StreamType>
	void logValue(StreamType& logStream, std::string text, MatrixType value, std::string equalityType = "=") {

		logStream << text << " " << equalityType << " " << value << std::endl << std::flush;
	}

	template <typename StreamType, typename T>
	void logValue(StreamType &logStream, std::string text, T value, std::string equalityType = "=") {

		logStream << text << " " << equalityType << " " << Tools::stringify(value) << std::endl << std::flush;
	}

}
