///////////////////////////////////////////////////////////////////////
// TorchHandler
//
// A module to evaluate a PyTorch model in ART. 
// Based on TFHandler for tensorflow.
//
// \author $Author: Alejandro Yankelevich
////////////////////////////////////////////////////////////////////////

#ifndef TORCHHANDLER_H
#define TORCHHANDLER_H

//#include <torch/csrc/jit/serialization/import.h>
#include <vector>
#include  <iostream>
#include  <string>

namespace torch {
namespace jit {
  class Module;  // 👈 forward declare the dependency
   } }

namespace c10{
  class IValue;  // 👈 forward declare the dependency
}

namespace torch
{
  /// Wrapper for Tensorflow which handles construction and prediction
  class TorchHandler
  {
  public:
    /// Basic constructor, takes path to model
    TorchHandler(std::string model);
    ~TorchHandler();
    void Initialize();

    c10::IValue Predict(std::vector<c10::IValue> inputs);

  private:
    torch::jit::Module *fModule;
    std::string fModelPath;
  };
}

#endif