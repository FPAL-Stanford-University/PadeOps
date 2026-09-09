#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/core/framework/logging.h"
#include <typeinfo>

//const int spacialdim = 64;
std::vector<tensorflow::Tensor> outputs;
using namespace tensorflow;
extern "C" float* loadtensorflowstruc(double* delta_in, double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, int* space1, int* space2, int* space3);
extern "C" float* loadtensorflowmag(double* delta_in, double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, int* space1, int* space2, int* space3);
extern "C" float* loadtensorflowlstm(float* delta_in, float* train_in, float* strain2_in, float* strain2_m_in, float* rot2_in, float* strainrot_m_in, float* rotstrainrot_m_in, float* strain2rot_m_in, float* rotstrain2_m_in, int* space1, int* space2, int* space3);
extern "C" void free(void* ptr);


float* loadtensorflowstruc(double* delta_in, double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, int* space1, int* space2, int* space3) {
	// Initialize a tensorflow session

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
		"unet2normpareto2", //std::string containing path of the model bundle
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
	tensorflow::Tensor outtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({7*space1a*space2a*space3a}));


	auto train_loadtensor_mapped = train_loadtensor.tensor<float, 5>();
	auto deltatensor_mapped = deltatensor.tensor<float, 5>();
	auto strain2tensor_mapped = strain2tensor.tensor<float, 6>();
	auto strain2_mtensor_mapped = strain2_mtensor.tensor<float, 6>();
	auto rot2tensor_mapped = rot2tensor.tensor<float, 6>();
	auto strainrot_mtensor_mapped = strainrot_mtensor.tensor<float, 6>();
	auto rotstrainrot_mtensor_mapped = rotstrainrot_mtensor.tensor<float, 6>();
	auto strain2rot_mtensor_mapped = strain2rot_mtensor.tensor<float, 6>();
	auto rotstrain2_mtensor_mapped = rotstrain2_mtensor.tensor<float, 6>();
	auto outtensor_mapped = outtensor.tensor<float, 1>();



	// fill test arrays with values
	int count = 0;
	for (int j = 0; j < space1a; j++) {
		for (int k = 0; k < space2a; k++) {
			for (int l = 0; l < space3a; l++) {
				for (int i = 0; i < 8; i++) {
					deltatensor_mapped(0, j, k, l, i) = (float)delta_in[count];
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
					train_loadtensor_mapped(0, j, k, l, h) = (float)train_in[count];
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
						strain2tensor_mapped(0, j, k, l, i, h) = (float)strain2_in[count];
						strain2_mtensor_mapped(0, j, k, l, i, h) = (float)strain2_m_in[count];
						rot2tensor_mapped(0, j, k, l, i, h) = (float)rot2_in[count];
						strainrot_mtensor_mapped(0, j, k, l, i, h) = (float)strainrot_m_in[count];
						rotstrainrot_mtensor_mapped(0, j, k, l, i, h) = (float)rotstrainrot_m_in[count];
						strain2rot_mtensor_mapped(0, j, k, l, i, h) = (float)strain2rot_m_in[count];
						rotstrain2_mtensor_mapped(0, j, k, l, i, h) = (float)rotstrain2_m_in[count];
						count = count + 1;
					}
				}
			}
		}
	}

	//for (int j = 126; j < 128; j++) {
	//	for (int k = 30; k < 32; k++) {
	//		for (int l = 94; l < 96; l++) {
	//			std::cout << train_loadtensor_mapped(0, j, k, l, 0) << '\n';
	//		}
	//	}
	//}


	//std::cout << "Strain2_m" << '\n';
	//for (int j = 0; j < 2; j++) {
	//	for (int k = 0; k < 2; k++) {
	//		for (int l = 0; l < 2; l++) {
	//			std::cout << strain2_mtensor_mapped(0, j, k, l, 0, 0) << '\n';
	//		}
	//	}
	//}

	//for (int j = 126; j < 128; j++) {
	//	for (int k = 30; k < 32; k++) {
	//		for (int l = 94; l < 96; l++) {
	//			std::cout << strain2_mtensor_mapped(0, j, k, l, 0, 0) << '\n';
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
	//auto out = finalOutput.data();
	//auto out_array = static_cast<float*>(out);


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
	for (int h = 0; h < 6; h++) {
		for (int j = 0; j < space1a; j++) {
			for (int k = 0; k < space2a; k++) {
				for (int l = 0; l < space3a; l++) {
					out_arraypass[count] = finalOutput(0, h, j, k, l);
					count = count + 1;
				}
			}
		}
	}
	
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

float* loadtensorflowmag(double* delta_in, double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, int* space1, int* space2, int* space3) {

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
		"unet2norm2test", //std::string containing path of the model bundle
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
	tensorflow::Tensor outtensor2(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 2 * space1a * space2a * space3a }));
	tensorflow::Tensor buffer1(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, 6, space1a, space2a, space3a }));
	tensorflow::Tensor buffer2(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, space1a, space2a, space3a, 3, 3 }));


	auto train_loadtensor_mapped = train_loadtensor.tensor<float, 5>();
	auto deltatensor_mapped = deltatensor.tensor<float, 5>();
	auto strain2tensor_mapped = strain2tensor.tensor<float, 6>();
	auto strain2_mtensor_mapped = strain2_mtensor.tensor<float, 6>();
	auto rot2tensor_mapped = rot2tensor.tensor<float, 6>();
	auto strainrot_mtensor_mapped = strainrot_mtensor.tensor<float, 6>();
	auto rotstrainrot_mtensor_mapped = rotstrainrot_mtensor.tensor<float, 6>();
	auto strain2rot_mtensor_mapped = strain2rot_mtensor.tensor<float, 6>();
	auto rotstrain2_mtensor_mapped = rotstrain2_mtensor.tensor<float, 6>();
	auto outtensor_mapped2 = outtensor2.tensor<float, 1>();
	auto buffer1_mapped = buffer1.tensor<float, 5>();
	auto buffer2_mapped = buffer2.tensor<float, 6>();

	int count = 0;
	for (int j = 0; j < space1a; j++) {
		for (int k = 0; k < space2a; k++) {
			for (int l = 0; l < space3a; l++) {
				for (int i = 0; i < 8; i++) {
					deltatensor_mapped(0, j, k, l, i) = (float)delta_in[count];
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
					train_loadtensor_mapped(0, j, k, l, h) = (float)train_in[count];
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
						strain2tensor_mapped(0, j, k, l, i, h) = (float)strain2_in[count];
						strain2_mtensor_mapped(0, j, k, l, i, h) = (float)strain2_m_in[count];
						rot2tensor_mapped(0, j, k, l, i, h) = (float)rot2_in[count];
						strainrot_mtensor_mapped(0, j, k, l, i, h) = (float)strainrot_m_in[count];
						rotstrainrot_mtensor_mapped(0, j, k, l, i, h) = (float)rotstrainrot_m_in[count];
						strain2rot_mtensor_mapped(0, j, k, l, i, h) = (float)strain2rot_m_in[count];
						rotstrain2_mtensor_mapped(0, j, k, l, i, h) = (float)rotstrain2_m_in[count];
						count = count + 1;
					}
				}
			}
		}
	}

	//std::cout << "Strain2" << '\n';
	//std::cout << strain2tensor_mapped(0, 0, 0, 0, 0, 0) << '\n';
	//std::cout << strain2tensor_mapped(0, 0, 0, 1, 0, 1) << '\n';
	//std::cout << strain2tensor_mapped(0, 31, 31, 31, 2, 2) << '\n';

	//std::cout << strain2tensor(0, 0, 0, 0, 0, 0) << '\n';
	//std::cout << strain2tensor(0, 0, 0, 1, 0, 1) << '\n';
	//std::cout << strain2tensor(0, 31, 31, 31, 2, 2) << '\n';

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
	//std::cout << rotstrainrot_mtensor_mapped(0, 0, 0, 1, 0, 1) << '\n';
	//std::cout << rotstrainrot_mtensor_mapped(0, 31, 31, 31, 2, 2) << '\n';

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
	{ "serving_default_input_10", buffer2 },
	{ "serving_default_input_11", buffer1 },
	};

	std::vector<std::string> fetches = { "StatefulPartitionedCall" };

	auto status = savedModelBundle->GetSession()->Run(inputs, fetches, {}, &outputs);
	TF_CHECK_OK(status);

	for (const auto& record : outputs) {
		LOG(INFO) << record.DebugString();
	}

	auto finalOutput = outputs[0].tensor<float, 4>();
	//auto out = finalOutput.data();
	//auto out_array2 = static_cast<float*>(out);
	//std::cout << out_array2[0, 0, 0, 0] << '\n';
	//std::cout << out_array2[0, 0, 0, 1] << '\n';
	//std::cout << out_array2[0, 0, 0, 2] << '\n';
	//std::cout << out_array2[0, 0, 1, 2] << '\n';
	//std::cout << out_array2[0, 4, 4, 4] << '\n';
	//std::cout << out_array2[0, 31, 31, 31] << '\n';

	float* out_arraypass2 = (float*)malloc(sizeof(float) * 2 * space1a * space2a * space3a);
	count = space1a * space2a * space3a;
	for (int j = 0; j < space1a; j++) {
		for (int k = 0; k < space2a; k++) {
			for (int l = 0; l < space3a; l++) {
				out_arraypass2[count] = finalOutput(0, j, k, l);
				count = count + 1;
			}
		}
	}
	//std::cout << "Final" << '\n';
	//for (int j = 0; j < 2; j++) {
	//	for (int k = 0; k < 2; k++) {
	//		for (int l = 0; l < 2; l++) {
	//			std::cout << finalOutput(0, j, k, l) << '\n';
	//			count = count + 1;
	//		}
	//	}
	//}
	//for (int j = 126; j < 128; j++) {
	//	for (int k = 30; k < 32; k++) {
	//		for (int l = 94; l < 96; l++) {
	//			std::cout << finalOutput(0, j, k, l) << '\n';
	//			count = count + 1;
	//		}
	//	}
	//}
	//std::cout << "Final" << '\n';
	//for (int j = 0; j < 4; j++) {
	//	for (int k = 0; k < 4; k++) {
	//		for (int l = 0; l < 4; l++) {
	//			std::cout << finalOutput(0, j, k, l) << '\n';
	//			count = count + 1;
	//		}
	//	}
	//}
	//std::cout << "Final" << '\n';
	//for (int j = 59; j < 64; j++) {
	//	for (int k = 59; k < 64; k++) {
	//		for (int l = 59; l < 64; l++) {
	//			std::cout << finalOutput(0, j, k, l) << '\n';
	//			count = count + 1;
	//		}
	//	}
	//}
	//std::cout << out_arraypass2[0+8*8*8] << '\n';
	//std::cout << out_arraypass2[1 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[2 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[3 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[4 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[5 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[6 + 8 * 8 * 8] << '\n';
    //std::cout << out_arraypass2[7 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[8 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[9 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[10 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[17 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[18 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[500 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[1000 + 8 * 8 * 8] << '\n';
	//std::cout << out_arraypass2[2*8*8*8-1] << '\n';

	return out_arraypass2;
}

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

//int main()
//{
//	std::cout << "Running Main Loop";
//}