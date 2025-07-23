///////////////////////////////////////////////////////////////////////
// TorchHandler
//
// A module to evaluate a PyTorch model in ART. 
// Based on TFHandler for tensorflow.
//
// \author $Author: Alejandro Yankelevich
////////////////////////////////////////////////////////////////////////

#include <torch/script.h>
#include "TorchHandler.h"

namespace torch
{
  TorchHandler::TorchHandler(std::string ModelPath):
    fModelPath(ModelPath)
  {
    Initialize();
  }

  TorchHandler::~TorchHandler() {}

  void TorchHandler::Initialize()
  {
    try
    {
      // Deserialize the ScriptModule from a file using torch::jit::load().
      std::cout << "Loading model " << fModelPath << std::endl;
      *fModule = torch::jit::load(fModelPath, c10::kCPU);
    }
    catch (const c10::Error& e)
    {
      std::cout << "Error loading the PyTorch model:" << fModelPath << std::endl;
      return;
    }

    std::cout << "Successfully loaded PyTorch model." << std::endl;
    return;
  }

  torch::jit::IValue TorchHandler::Predict(std::vector<torch::jit::IValue> inputs)
  {
    torch::jit::IValue outputs = fModule->forward(inputs);

    return outputs;
  }

}