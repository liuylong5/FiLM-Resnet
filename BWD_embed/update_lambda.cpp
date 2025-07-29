#include <iostream>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <torch/script.h>
#include <../libtorch/include/torch/csrc/autograd/autograd.h>  // autograd library in LibTorch
#include <cstdint>
#include <chrono>
#include <vector>
#include <cstring>  // For std::memcpy
#include <fstream>
#include <sstream>
#include <iomanip>



namespace {

char const* set_model_path() {
	const auto* env = std::getenv("MODEL_PATH");
	if (env) {
		return env;
	} else {
		return "./";
	}
}

struct Auxiliary {
    std::string model_path;
    std::string dnn_conc_filename{"CONC_MODEL.pt"};
    torch::jit::script::Module conc_module;

    Auxiliary()
    : model_path(set_model_path()) {
    }

    void load_model() {
        const std::string conc_file = model_path + "/" + dnn_conc_filename;
        if (!std::filesystem::exists(conc_file)) {
            std::cerr << conc_file << " doesn't exist." << std::endl;
            exit(-1);
        }

        try {
            conc_module = torch::jit::load(conc_file);
            conc_module.eval();  
        } catch (const c10::Error& e) {
            std::cerr << "Error loading the conc_module, due to: " << e.what() << std::endl;
            exit(-1);
        }
    }

    ~Auxiliary() {
    }

    void update_lambda(float* CHKGRID_IN,
                    float* lambda_in,
                    float* lambda_out,
                    int NUMB_MECH_SPC,
                    int NPHOTAB,
                    int NCOUNT,
                    int MYPE,
                    int NG) 
    {
        try {
            torch::TensorOptions options = torch::TensorOptions()
                .dtype(torch::kFloat32)
                .requires_grad(true);
            
            torch::Tensor conc_tensor = torch::from_blob(
                CHKGRID_IN, 
                {NCOUNT, NUMB_MECH_SPC + NPHOTAB - 11}, 
                options
            );

            torch::Tensor x_example = conc_tensor.slice(1, 0, 56);
            torch::Tensor params_example = conc_tensor.slice(1, 56, 82).detach(); 
            torch::Tensor if_light = conc_tensor.slice(1, 82, 83).detach();      // float32
            torch::Tensor time_step = conc_tensor.slice(1, 83, 84).detach();     // float32

            ////////////////////////////////////////////////////////////////////
            // 2. forward pass
            ////////////////////////////////////////////////////////////////////
            
            std::vector<torch::jit::IValue> inputs;
            inputs.push_back(x_example);
            inputs.push_back(params_example);
            inputs.push_back(if_light);
            inputs.push_back(time_step);

            auto output_tuple = conc_module.forward(inputs).toTuple();
            at::Tensor output = output_tuple->elements()[1].toTensor();






            std::vector<int64_t> indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 54, 55, 60, 61, 64};  
            torch::Tensor selected_output = output.index_select(1, torch::tensor(indices));

            torch::Tensor lambda_in_tensor = torch::from_blob(
                lambda_in,
                {NCOUNT, 56},
                torch::TensorOptions().dtype(torch::kFloat32)
            );



            // Jacobian-free embed
            auto gradients = torch::autograd::grad(
                {selected_output},
                {x_example},
                {lambda_in_tensor},
                /*retain_graph=*/false,
                /*create_graph=*/false
            );


            torch::Tensor lambda_out_tensor = gradients[0];


            if (lambda_out_tensor.isnan().any().item<bool>()) {
                std::cerr << "NaN detected in lambda_out_tensor." << std::endl;

                torch::Tensor nan_mask = lambda_out_tensor.isnan();  

            }


            lambda_out_tensor = lambda_out_tensor.contiguous();
            std::memcpy(
                lambda_out,
                lambda_out_tensor.data_ptr<float>(),
                NCOUNT * 56 * sizeof(float)
            );

        } catch (const c10::Error& e) {
            std::cerr << "Torch error: " << e.what() << std::endl;
            throw std::runtime_error("Failed in update_lambda");
        } catch (const std::exception& e) {
            std::cerr << "Standard error: " << e.what() << std::endl;
            throw;
        }
    }


	void initialize() {
		auto load_func = [this](){
			this->load_model();
			return 0;
		};
		static int __ = load_func();
		(void)__;
	}
};


Auxiliary* auxiliary() {
    static Auxiliary aux;
    aux.initialize(); 
    return &aux;
}


} // anonymous namespace


#ifdef __cplusplus
extern "C" {
#endif

void UP_LAMBDA (float* CHKGRID_IN, float* lambda_in, float* lambda_out, int NUMB_MECH_SPC, int NPHOTAB, int NCOUNT, int MYPE, int NG){
	auxiliary()->update_lambda(CHKGRID_IN,lambda_in, lambda_out, NUMB_MECH_SPC, NPHOTAB, NCOUNT, MYPE, NG);
}


#ifdef __cplusplus
}
#endif
