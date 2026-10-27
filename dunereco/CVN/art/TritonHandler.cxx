////////////////////////////////////////////////////////////////////////
/// \file    TritonHandler.cxx
/// \brief   TritonHandler for CVN
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include  <iostream>
#include  <string>
#include  <iomanip>
#include  <cstdint>
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dunereco/CVN/art/TritonHandler.h"
#include "dunereco/CVN/func/CVNImageUtils.h"


/*
namespace
{
	// Order-sensitive, bit-exact hash over the raw float data (FNV-1a).
	//   // Any difference in content OR ordering will change this value.
	uint64_t ChargeVecHash(const std::vector<float>& v)
	{
		uint64_t hash = 14695981039346656037ULL; // FNV offset basis
		const uint64_t prime = 1099511628211ULL;
		const unsigned char* bytes = reinterpret_cast<const unsigned char*>(v.data());
		size_t nbytes = v.size() * sizeof(float);
		for (size_t i = 0; i < nbytes; ++i)
		{
			hash ^= bytes[i];
			hash *= prime;
		}
		return hash;
	}

	void PrintChargeChecksum(const char* label, const std::vector<float>& v)
	{
		double sum = 0.0;
		for (float x : v) sum += x;

		std::cout << "[ChargeChecksum] " << label
			<< " size=" << v.size()
			<< " sum=" << std::setprecision(17) << sum
			<< " hash=0x" << std::hex << ChargeVecHash(v) << std::dec
			<< std::endl;
	}
}
*/

namespace cvn
{

	TritonHandler::TritonHandler(const fhicl::ParameterSet& pset):
		fUseLogChargeScale(pset.get<bool>("ChargeLogScale")),
		fImageWires(pset.get<unsigned int>("NImageWires")),
		fImageTDCs(pset.get<unsigned int>("NImageTDCs")),
		fNViews(3),
		fReverseViews(pset.get<std::vector<bool>>("ReverseViews")),
		fInputNames(pset.get<std::vector<std::string>>(
					"InputNames", {"view0", "view1", "view2"})),
		fOutputNames(pset.get<std::vector<std::string>>(
					"OutputNames", {"output_is_antineutrino",
					"output_flavour",
					"output_interaction",
					"output_protons",
					"output_pions",
					"output_pizeros",
					"output_neutrons"}))
		{
			mf::LogInfo("TritonHandler") << "Loading Triton client" << std::endl;

			std::cout << "Loading Triton client: ";
			fTritonClient = std::make_unique<lartriton::TritonClient>(pset.get<fhicl::ParameterSet>("TritonConfig"));
			if (!fTritonClient){
				art::Exception(art::errors::Unknown) << "Triton client not created correctly";
			}
		}
	
	// Check the network outputs
	bool check_Triton(const std::vector< std::vector< float > > & outputs)
	{
		if (outputs.size() == 1) return true;
		size_t aux = 0;
		for (size_t o = 0; o < outputs.size(); ++o)
		{
			size_t aux2 = 0;

			for (size_t i = 0; i < outputs[o].size(); ++i)
				if (outputs[o][i] == 0.0 || outputs[o][i] == 1.0)
					aux2++;
			if (aux2 == outputs[o].size()) aux++;
		}
		return aux == outputs.size() ? false : true;
	}

	// Fill outputs with value -3
	void fillEmpty_Triton(std::vector< std::vector< float > > & outputs)
	{
		std::cout << "Inside fillEmpty_Triton: ";
		for (size_t o = 0; o < outputs.size(); ++o)
		{
			for (size_t i = 0; i < outputs[o].size(); ++i)
				outputs[o][i] = -3.0;
		}
		return;
	}
	

	void TritonHandler::SplitViews(const ImageVectorF& fullImage,
			std::vector<std::vector<float>>& viewData) const
	{
		viewData.clear();
		viewData.resize(fNViews);

		const size_t nWires = fullImage.size();
		for (auto& vd : viewData) vd.reserve(nWires * (nWires ? fullImage[0].size() : 0));

		for (size_t w = 0; w < nWires; ++w)
		{
			const size_t nTdcs = fullImage[w].size();
			for (size_t t = 0; t < nTdcs; ++t)
			{
				const auto& channelVec = fullImage[w][t]; // size fNViews (3): {v0, v1, v2}
			for (unsigned v = 0; v < fNViews; ++v)
				viewData[v].push_back(channelVec[v]);
			}
		}
	}


	std::vector< std::vector<float> > TritonHandler::Predict(const PixelMap& pm)
	{
		///====

		std::cout << "Inside Predict: ";
		CVNImageUtils imageUtils(fImageWires,fImageTDCs, fNViews);
		// Configure the image utility
		imageUtils.SetViewReversal(fReverseViews);
		imageUtils.SetImageSize(fImageWires,fImageTDCs,fNViews);
		imageUtils.SetLogScale(fUseLogChargeScale);
		imageUtils.SetPixelMapSize(pm.NWire(), pm.NTdc());


		// --- DEBUG: checksum the raw PixelMap charge vectors ---
		//PrintChargeChecksum("view0 (fPEX)", pm.fPEX);
		//PrintChargeChecksum("view1 (fPEY)", pm.fPEY);
		//PrintChargeChecksum("view2 (fPEZ)", pm.fPEZ);
		// --------------------------------------------------------



		ImageVectorF thisImage;
		imageUtils.ConvertPixelMapToImageVectorF(pm,thisImage);

		// Split the [view][wire][tdc] image into one flat buffer per view
		std::vector<std::vector<float>> viewData;
		SplitViews(thisImage, viewData);

		bool status = false;
		int counter = 0;
		std::vector< std::vector<float> > cvnResults; // shape(#outputs, output_size)

		do{ // do until it gets a correct result
			fTritonClient->reset();

			// model has max_batch_size: 0, fixed shape [1,500,500,1] per view,
			// so no setBatchSize()/setShape() calls are needed here
			for (unsigned v = 0; v < fNViews; ++v)
			{
				//std::cout << "[SplitViewsCheck] view " << v << " size=" << viewData[v].size() << std::endl;
				auto& input = fTritonClient->input().at(fInputNames[v]);
				auto inputData = std::make_shared<lartriton::TritonInput<float>>();
				inputData->push_back(viewData[v]);
				input.toServer(inputData);
			}

			fTritonClient->dispatch();

			cvnResults.clear();
			for (auto const& outName : fOutputNames)
			{
				auto const& output = fTritonClient->output().at(outName);
				lartriton::TritonOutput<float> outputData = output.fromServer<float>();
				cvnResults.emplace_back(outputData[0].begin(), outputData[0].end());
			}

			status = check_Triton(cvnResults);
			counter++;
			if(counter==10){
				std::cout << "Error, CVN never outputing a correct result. Filling result with zeros.";
				std::cout << std::endl;
				fillEmpty_Triton(cvnResults);
				break;
			}
		}while(status == false);

		std::cout << "Classifier summary: ";
		std::cout << "Triton is working: ";
		std::cout << std::endl;
		int output_index = 0;
		for(auto const & output : cvnResults)
		{
			std::cout << "Output " << output_index++ << ": ";
			for(auto const v : output)
				std::cout << v << ", ";
			std::cout << std::endl;
		}
		std::cout << std::endl;

		return cvnResults;
	}

	int TritonHandler::NOutput() const
	{
		return static_cast<int>(fOutputNames.size());
	}

	int TritonHandler::NFeatures() const
	{
		return static_cast<int>(fNViews);
	}

	// The output_flavour head is already a 4-element vector
	// [numu, nue, nutau, NC]-style per config.pbtxt (dims: [1,4])
	std::vector<float> TritonHandler::PredictFlavour(const PixelMap& pm){
		std::cout << "Inside Flavour: ";

		std::vector<std::vector<float>> fullResults = this->Predict(pm);
		return fullResults.at(1); // output_flavour

	}

}

