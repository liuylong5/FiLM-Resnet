#include <iostream>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <torch/script.h>
// #include <cuda_runtime.h>
#include <cstdint>
#include <chrono>
#include <vector>
#include <cstring>  // For std::memcpy


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
            // std::cout << "Model loaded successfully from " << conc_file << std::endl;
        } catch (const c10::Error& e) {
            std::cerr << "Error loading the conc_module, due to: " << e.what() << std::endl;
            exit(-1);
        }
    }

    ~Auxiliary() {
    }


    void conc_forward(float* CHKGRID_IN, // C
                    float* CHKGRID_OUT,
                    int NUMB_MECH_SPC,
                    int NPHOTAB,
                    int NCOUNT,
                    int MYPE,
                    int NG) {
        
        auto start = std::chrono::high_resolution_clock::now();

        // Create a tensor from the input array on the CPU with double precision.
        auto options = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
        torch::Tensor conc_tensor = torch::from_blob(CHKGRID_IN, {NCOUNT, NUMB_MECH_SPC + NPHOTAB - 11}, options);

        auto tensor_creation_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const tensor_creation_duration = tensor_creation_end - start;



        torch::Tensor x_example = conc_tensor.slice(1, 0, 56);
        torch::Tensor params_example = conc_tensor.slice(1, 56, 82); 
        torch::Tensor if_light_float = conc_tensor.slice(1, 82, 83);
        torch::Tensor time_step_float = conc_tensor.slice(1, 83, 84);
        torch::Tensor if_light = torch::round(conc_tensor.slice(1, 82, 83));      // float32
        torch::Tensor time_step = torch::round(conc_tensor.slice(1, 83, 84));     // float32


        start = std::chrono::high_resolution_clock::now();


        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(x_example);
        inputs.push_back(params_example);
        inputs.push_back(if_light);
        inputs.push_back(time_step);

        // Assuming `conc_module` is a loaded PyTorch model for CPU use.
        auto output_tuple = conc_module.forward(inputs).toTuple();
        at::Tensor output = output_tuple->elements()[1].toTensor();

        auto forward_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const forward_duration = forward_end - start;
        //std::cout << "Forward pass duration: " << forward_duration.count() << " seconds." << std::endl;
        start = std::chrono::high_resolution_clock::now();

        // Copy the tensor back to the input buffer. Since the tensor was created
        std::memcpy(CHKGRID_OUT, output.data_ptr<float>(), NCOUNT * 66 * sizeof(float));
        auto data_copy_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const data_copy_duration = data_copy_end - start;
        //std::cout << "Data copy duration: " << data_copy_duration.count() << " seconds." << std::endl;
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

void DNN_CONC (float* CHKGRID_IN, float* CHKGRID_OUT, int NUMB_MECH_SPC, int NPHOTAB, int NCOUNT, int MYPE, int NG){
	auxiliary()->conc_forward(CHKGRID_IN, CHKGRID_OUT, NUMB_MECH_SPC, NPHOTAB, NCOUNT, MYPE, NG);
}


#ifdef __cplusplus
}
#endif
