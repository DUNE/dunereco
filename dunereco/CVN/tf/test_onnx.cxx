#include "onnx_graph.h"

int main() {

    std::string model_file = "/exp/dune/app/users/calcuttj/tf2onnx_venv/model.onnx";
    std::vector<std::vector<std::vector<std::vector<float>>>> input = {{{{1.}}}};
    auto model = onnx::Model::create(model_file);
    auto result = model->run(input);
    std::cout << result.size() << std::endl;
}
