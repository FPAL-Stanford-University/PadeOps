#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <typeinfo>
#include <iostream>
#include <sys/resource.h>


//const int spacialdim = 64;
namespace py = pybind11;
extern "C" float* loadpytorchspace(double* delta_in, double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, int* space1, int* space2, int* space3);
extern "C" void free(void* ptr);


float* loadpytorchspace(double* delta_in, double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, int* space1, int* space2, int* space3) {
	//const rlim_t kStackSize = 512 * 1024 * 1024;   // min stack size = 16 MB
	//struct rlimit rl;
	//int result;

	//result = getrlimit(RLIMIT_STACK, &rl);
	//if (result == 0)
	//{
	//	if (rl.rlim_cur < kStackSize)
	//	{
	//		rl.rlim_cur = kStackSize;
	//		result = setrlimit(RLIMIT_STACK, &rl);
	//		if (result != 0)
	//		{
	//			fprintf(stderr, "setrlimit returned result = %d\n", result);
	//		}
	//	}
	//}
	std::cout << "GOT HERE" << '\n';
	int space1a = *space1;
	int space2a = *space2;
	int space3a = *space3;
	std::cout << "GOT HERE15" << '\n';
	py::scoped_interpreter guard{};
	//py::initialize_interpreter();
	std::cout << "GOT HERE2" << '\n';
	py::module m = py::module::import("interfacetest");
	std::cout << "GOT HERE3" << '\n';
	py::gil_scoped_release release;

	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> deltatensor_mapped(1, std::vector<std::vector<std::vector<std::vector<double>>>>(space1a, std::vector<std::vector<std::vector<double>>>(space2a, std::vector<std::vector<double>>(space3a, std::vector<double>(8)))));
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> train_loadtensor_mapped(1, std::vector<std::vector<std::vector<std::vector<double>>>>(space1a, std::vector<std::vector<std::vector<double>>>(space2a, std::vector<std::vector<double>>(space3a, std::vector<double>(8)))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> strain2tensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> strain2_mtensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> rot2tensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> strainrot_mtensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> rotstrainrot_mtensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> strain2rot_mtensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> rotstrain2_mtensor_mapped(1, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(space1a, std::vector<std::vector<std::vector<std::vector<double>>>>(space2a, std::vector<std::vector<std::vector<double>>>(space3a, std::vector<std::vector<double>>(3, std::vector<double>(3))))));

	//std::array<std::array<std::array<std::array<std::array<double, 8>, space3a>, space2a>, space1a>, 1> deltatensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<double, 7>, space3a>, space2a>, space1a>, 1> train_loadtensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> strain2tensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> strain2_mtensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> rot2tensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> strainrot_mtensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> rotstrainrot_mtensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> strain2rot_mtensor_mapped;
	//std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, space3a>, space2a>, space1a>, 1> rotstrain2_mtensor_mapped;

	// fill test arrays with values
	int count = 0;
	for (int j = 0; j < space1a; j++) {
		for (int k = 0; k < space2a; k++) {
			for (int l = 0; l < space3a; l++) {
				for (int i = 0; i < 8; i++) {
					deltatensor_mapped[0][j][k][l][i] = (float)delta_in[count];
					count = count + 1;
				}
			}
		}
	}

	count = 0;
	for (int j = 0; j < space1a; j++) {
		for (int k = 0; k < space2a; k++) {
			for (int l = 0; l < space3a; l++) {
				for (int h = 0; h < 7; h++) {
					train_loadtensor_mapped[0][j][k][l][h] = (float)train_in[count];
					count = count + 1;
				}
			}
		}
	}

	count = 0;
	for (int j = 0; j < space1a; j++) {
		for (int k = 0; k < space2a; k++) {
			for (int l = 0; l < space3a; l++) {
				for (int i = 0; i < 3; i++) {
					for (int h = 0; h < 3; h++) {
						strain2tensor_mapped[0][j][k][l][i][h] = (float)strain2_in[count];
						strain2_mtensor_mapped[0][j][k][l][i][h] = (float)strain2_m_in[count];
						rot2tensor_mapped[0][j][k][l][i][h] = (float)rot2_in[count];
						strainrot_mtensor_mapped[0][j][k][l][i][h] = (float)strainrot_m_in[count];
						rotstrainrot_mtensor_mapped[0][j][k][l][i][h] = (float)rotstrainrot_m_in[count];
						strain2rot_mtensor_mapped[0][j][k][l][i][h] = (float)strain2rot_m_in[count];
						rotstrain2_mtensor_mapped[0][j][k][l][i][h] = (float)rotstrain2_m_in[count];
						count = count + 1;
					}
				}
			}
		}
	}

	//std::cout << "Strain2 1-1" << '\n';
	//for (int j = 0; j < 4; j++) {
	//	for (int k = 0; k < 4; k++) {
	//		for (int l = 0; l < 4; l++) {
	//			std::cout << strain2tensor_mapped(0, j, k, l, 0, 0) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Strain2 1-1" << '\n';
	//for (int j = 59; j < 64; j++) {
	//	for (int k = 59; k < 64; k++) {
	//		for (int l = 59; l < 64; l++) {
	//			std::cout << strain2tensor_mapped(0, j, k, l, 0, 0) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Strain2 2-1" << '\n';
	//for (int j = 0; j < 4; j++) {
	//	for (int k = 0; k < 4; k++) {
	//		for (int l = 0; l < 4; l++) {
	//			std::cout << strain2tensor_mapped(0, j, k, l, 1, 0) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Strain2 2-1" << '\n';
	//for (int j = 59; j < 64; j++) {
	//	for (int k = 59; k < 64; k++) {
	//		for (int l = 59; l < 64; l++) {
	//			std::cout << strain2tensor_mapped(0, j, k, l, 1, 0) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Invariants 3" << '\n';
	//for (int j = 0; j < 4; j++) {
	//	for (int k = 0; k < 4; k++) {
	//		for (int l = 0; l < 4; l++) {
	//			std::cout << train_loadtensor_mapped(0, j, k, l, 2) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Invariants 3" << '\n';
	//for (int j = 59; j < 64; j++) {
	//	for (int k = 59; k < 64; k++) {
	//		for (int l = 59; l < 64; l++) {
	//			std::cout << train_loadtensor_mapped(0, j, k, l, 2) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Delta" << '\n';
	//for (int j = 0; j < 4; j++) {
	//	for (int k = 0; k < 4; k++) {
	//		for (int l = 0; l < 4; l++) {
	//			std::cout << deltatensor_mapped(0, j, k, l, 2) << '\n';
	//		}
	//	}
	//}

	//std::cout << "Delta" << '\n';
	//for (int j = 59; j < 64; j++) {
	//	for (int k = 59; k < 64; k++) {
	//		for (int l = 59; l < 64; l++) {
	//			std::cout << deltatensor_mapped(0, j, k, l, 2) << '\n';
	//		}
	//	}
	//}
	//std::cout << test3[0, 0, 0, 0, 0, 0] << '\n';
	//std::cout << test3[0, 0, 0, 1, 0, 1] << '\n';
	//std::cout << test3[0, 31, 31, 31, 2, 2] << '\n';
	
	//std::cout << strain2tensor_mapped(0, 31, 31, 31, 2, 2) << '\n';
	//std::cout << "invariants" << '\n';
	//std::cout << train_loadtensor_mapped(0, 0, 0, 0, 0) << '\n';
	//std::cout << train_loadtensor_mapped(0, 0, 0, 1, 1) << '\n';
	//std::cout << train_loadtensor_mapped(0, 31, 31, 31, 6) << '\n';
	//std::cout << "Delta" << '\n';
	//std::cout << deltatensor_mapped(0, 0, 0, 0, 0) << '\n';
	//std::cout << deltatensor_mapped(0, 0, 0, 1, 1) << '\n';
	//std::cout << deltatensor_mapped(0, 31, 31, 31, 7) << '\n';
	//std::cout << "rotstrainrot" << '\n';
	//std::cout << rotstrainrot_mtensor_mapped(0, 0, 0, 0, 0, 0) << '\n';
	//std::cout << rotstrainrot_mtensor_mapped(0, 0, 0, 0, 0, 2) << '\n';
	//std::cout << rotstrainrot_mtensor_mapped(0, 0, 0, 0, 1, 0) << '\n';
	//std::cout << rotstrainrot_mtensor_mapped(0, 0, 0, 1, 0, 1) << '\n';
	//std::cout << rotstrainrot_mtensor_mapped(0, 31, 31, 31, 2, 2) << '\n';

	py::function test_func3 = m.attr("test_func3");
	py::buffer_info buf1 = test_func3(train_loadtensor_mapped, strain2tensor_mapped, strain2_mtensor_mapped, rot2tensor_mapped, strainrot_mtensor_mapped, rotstrainrot_mtensor_mapped, strain2rot_mtensor_mapped, rotstrain2_mtensor_mapped, deltatensor_mapped).cast<py::array_t<float>>().reshape({6*space1a*space2a*space3a}).request();
	//py::buffer_info buf1 = callPythonFunctionAsync(test_func3(train_loadtensor_mapped, strain2tensor_mapped, strain2_mtensor_mapped, rot2tensor_mapped, strainrot_mtensor_mapped, rotstrainrot_mtensor_mapped, strain2rot_mtensor_mapped, rotstrain2_mtensor_mapped, deltatensor_mapped)).cast<py::array_t<float>>().reshape({ 6 * space1a * space2a * space3a }).request();
	float* ptr1 = static_cast<float*>(buf1.ptr);

	//for (int j = 0; j < space1a; j++) {
	//	for (int k = 0; k < space2a; k++) {
	//		std::cout << strain2tensor_mapped(0, j, k, 0, 0, 0) << '\n';
	//	}
	//}

	//std::cout << "new" << '\n';
	//for (int j = 0; j < space1a; j++) {
	//	for (int k = 0; k < space2a; k++) {
	//		std::cout << out[0, 0, j, k, 0] << '\n';
	//	}
	//}

	//std::cout << "new" << '\n';
	//for (int j = 0; j < space1a; j++) {
	//	for (int k = 0; k < space2a; k++) {
	//		std::cout << finalOutput(0, 0, j, k, 0) << '\n';
	//	}
	//}

	//std::cout << out_array[0, 0, 0, 0, 0] << '\n';
	//std::cout << out_array[0, 0, 0, 0, 1] << '\n';
	//std::cout << out_array[0, 0, 0, 0, 2] << '\n';
	//std::cout << out_array[0, 0, 0, 1, 2] << '\n';
	//std::cout << out_array[0, 5, 4, 4, 4] << '\n';
	//std::cout << out_array[0, 5, 31, 31, 31] << '\n';

	float* out_arraypass = (float*)malloc(sizeof(float) * 7 * space1a * space2a * space3a);
	count = space1a * space2a * space3a;
	for (int count2 = 0; count2 < 6*space1a*space2a*space3a; count2++) {
		out_arraypass[count] = ptr1[count2];
		count = count + 1;
	}

	//py::finalize_interpreter();
	//std::cout << "Final" << '\n';
	//for (int j = 0; j < 4; j++) {
	//	for (int k = 0; k < 4; k++) {
	//		for (int l = 0; l < 4; l++) {
	//			std::cout << finalOutput(0, 1, j, k, l) << '\n';
	//			count = count + 1;
	//		}
	//	}
	//}
	//std::cout << "Final" << '\n';
	//for (int j = 59; j < 64; j++) {
	//	for (int k = 59; k < 64; k++) {
	//		for (int l = 59; l < 64; l++) {
	//			std::cout << finalOutput(0, 1, j, k, l) << '\n';
	//			count = count + 1;
	//		}
	//	}
	//}
	//out_arraypass = static_cast<float*>(outtensor_mapped.data());
	//out_arraypass = outtensor_mapped.data();
	//out_arraypass = reinterpret_cast<double*>(outtensor_mapped.data());
	//std::cout << out_arraypass[0 + space1a * space2a * space3a] << '\n';
	//std::cout << out_arraypass[1 + space1a * space2a * space3a] << '\n';
	//std::cout << out_arraypass[2 + space1a * space2a * space3a] << '\n';
	//std::cout << out_arraypass[7*space1a * space2a * space3a-1] << '\n';
	//std::cout << out_arraypass << '\n';
	//std::cout << out_arraypass[4 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[5 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[6 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[7 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[8 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[9 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[10 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[17 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[18 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[500 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[1000 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass[7*8*8*8-1] << '\n';
	return out_arraypass;
}

/*
float* loadtensorflowlstm(float* delta_in, float* train_in, float* strain2_in, float* strain2_m_in, float* rot2_in, float* strainrot_m_in, float* rotstrainrot_m_in, float* strain2rot_m_in, float* rotstrain2_m_in, int* space1, int* space2, int* space3) {
	// Initialize a tensorflow session
	//Session* session;
	//Status status = NewSession(SessionOptions(), &session);
	//if (!status.ok()) {
	//    std::cout << status.ToString() << "\n";
	//    return 1;
	//}
	int space1a = *space1;
	int space2a = *space2;
	int space3a = *space3;
	// We need to use SaveModelBundleLite as a in-memory model object for tensorflow's model bundle.
	const auto savedModelBundle = std::make_unique<tensorflow::SavedModelBundleLite>();

	// Create dummy options.
	tensorflow::SessionOptions sessionOptions;
	tensorflow::RunOptions runOptions;

	// Load the model bundle.
	const auto loadResult = tensorflow::LoadSavedModel(
		sessionOptions,
		runOptions,
		"unet2", //std::string containing path of the model bundle
		{ tensorflow::kSavedModelTagServe },
		savedModelBundle.get());

	// Check if loading was successful
	TF_CHECK_OK(loadResult);

	tensorflow::Tensor train_loadtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 7 }));
	tensorflow::Tensor deltatensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 8 }));
	tensorflow::Tensor strain2tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));
	tensorflow::Tensor strain2_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));
	tensorflow::Tensor rot2tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));
	tensorflow::Tensor strainrot_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));
	tensorflow::Tensor rotstrainrot_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));
	tensorflow::Tensor strain2rot_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));
	tensorflow::Tensor rotstrain2_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));


	auto train_loadtensor_mapped = train_loadtensor.tensor<float, 5>();
	auto deltatensor_mapped = deltatensor.tensor<float, 5>();
	auto strain2tensor_mapped = strain2tensor.tensor<float, 6>();
	auto strain2_mtensor_mapped = strain2_mtensor.tensor<float, 6>();
	auto rot2tensor_mapped = rot2tensor.tensor<float, 6>();
	auto strainrot_mtensor_mapped = strainrot_mtensor.tensor<float, 6>();
	auto rotstrainrot_mtensor_mapped = rotstrainrot_mtensor.tensor<float, 6>();
	auto strain2rot_mtensor_mapped = strain2rot_mtensor.tensor<float, 6>();
	auto rotstrain2_mtensor_mapped = rotstrain2_mtensor.tensor<float, 6>();



	// fill test arrays with values
	int count = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < space1a; j++) {
			for (int k = 0; k < space2a; k++) {
				for (int l = 0; l < space3a; l++) {
					deltatensor_mapped(0, j, k, l, i) = delta_in[count];
					count = count + 1;
				}
			}
		}
	}

	count = 0;
	for (int h = 0; h < 7; h++) {
		for (int j = 0; j < space1a; j++) {
			for (int k = 0; k < space2a; k++) {
				for (int l = 0; l < space3a; l++) {
					train_loadtensor_mapped(0, j, k, l, h) = train_in[count];
					count = count + 1;
				}
			}
		}
	}

	count = 0;
	for (int i = 0; i < 3; i++) {
		for (int h = 0; h < 3; h++) {
			for (int j = 0; j < space1a; j++) {
				for (int k = 0; k < space2a; k++) {
					for (int l = 0; l < space3a; l++) {
						strain2tensor_mapped(0, j, k, l, i, h) = strain2_in[count];
						strain2_mtensor_mapped(0, j, k, l, i, h) = strain2_m_in[count];
						rot2tensor_mapped(0, j, k, l, i, h) = rot2_in[count];
						strainrot_mtensor_mapped(0, j, k, l, i, h) = strainrot_m_in[count];
						rotstrainrot_mtensor_mapped(0, j, k, l, i, h) = rotstrainrot_m_in[count];
						strain2rot_mtensor_mapped(0, j, k, l, i, h) = strain2rot_m_in[count];
						rotstrain2_mtensor_mapped(0, j, k, l, i, h) = rotstrain2_m_in[count];
						count = count + 1;
					}
				}
			}
		}
	}

	std::cout << "Boop" << '\n';
	std::cout << deltatensor_mapped(0, 0, 0, 0, 0) << '\n';
	std::cout << deltatensor_mapped(0, 6, 6, 6, 6) << '\n';
	std::cout << deltatensor_mapped(0, 7, 7, 7, 7) << '\n';
	std::cout << "Boop" << '\n';
	std::cout << train_loadtensor_mapped(0, 0, 0, 0, 0) << '\n';
	std::cout << train_loadtensor_mapped(0, 7, 7, 7, 5) << '\n';
	std::cout << train_loadtensor_mapped(0, 7, 7, 7, 6) << '\n';
	std::cout << "Boop" << '\n';
	std::cout << strain2rot_mtensor_mapped(0, 0, 2, 0, 0, 0) << '\n';
	std::cout << strain2rot_mtensor_mapped(0, 5, 5, 5, 1, 1) << '\n';
	std::cout << strain2rot_mtensor_mapped(0, 7, 7, 7, 2, 2) << '\n';

	std::vector<std::pair<string, tensorflow::Tensor>> inputs = {
	{ "serving_default_input_1", train_loadtensor },
	{ "serving_default_input_2", strain2tensor },
	{ "serving_default_input_3", strain2_mtensor },
	{ "serving_default_input_4", rot2tensor },
	{ "serving_default_input_5", strainrot_mtensor },
	{ "serving_default_input_6", rotstrainrot_mtensor },
	{ "serving_default_input_7", strain2rot_mtensor },
	{ "serving_default_input_8", rotstrain2_mtensor },
	{ "serving_default_input_9", deltatensor },
	};

	std::vector<std::string> fetches = { "StatefulPartitionedCall" };

	auto status = savedModelBundle->GetSession()->Run(inputs, fetches, {}, &outputs);
	TF_CHECK_OK(status);

	for (const auto& record : outputs) {
		LOG(INFO) << record.DebugString();
	}

	auto finalOutput = outputs[0].tensor<float, 5>();
	auto out = finalOutput.data();
	//float* out_array = new float;
	auto out_array = static_cast<float*>(out);
	std::cout << out_array[0, 0, 0, 0, 0] << '\n';
	std::cout << out_array[0, 0, 0, 0, 1] << '\n';
	std::cout << out_array[0, 0, 0, 0, 2] << '\n';
	std::cout << out_array[0, 0, 0, 0, 3] << '\n';
	return out_array;
}
*/

//int main()
//{
//	std::cout << "Running Main Loop";
//}