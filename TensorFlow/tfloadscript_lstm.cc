#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/core/framework/logging.h"
#include <typeinfo>

const int spacialdim = 8;
std::vector<tensorflow::Tensor> outputs;
using namespace tensorflow;
extern "C" float* loadtensorflow(float* delta_in, float* train_in, float* strain2_in, float* strain2_m_in, float* rot2_in, float* strainrot_m_in, float* rotstrainrot_m_in, float* strain2rot_m_in, float* rotstrain2_m_in);

float* loadtensorflow(float* delta_in, float* train_in, float* strain2_in, float* strain2_m_in, float* rot2_in, float* strainrot_m_in, float* rotstrainrot_m_in, float* strain2rot_m_in, float* rotstrain2_m_in) {
	// Initialize a tensorflow session
	//Session* session;
	//Status status = NewSession(SessionOptions(), &session);
	//if (!status.ok()) {
	//    std::cout << status.ToString() << "\n";
	//    return 1;
	//}

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

	tensorflow::Tensor train_loadtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 7 }));
	tensorflow::Tensor deltatensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 8 }));
	tensorflow::Tensor strain2tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));
	tensorflow::Tensor strain2_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));
	tensorflow::Tensor rot2tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));
	tensorflow::Tensor strainrot_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));
	tensorflow::Tensor rotstrainrot_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));
	tensorflow::Tensor strain2rot_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));
	tensorflow::Tensor rotstrain2_mtensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, spacialdim, spacialdim, spacialdim, 3, 3 }));


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
		for (int j = 0; j < spacialdim; j++) {
			for (int k = 0; k < spacialdim; k++) {
				for (int l = 0; l < spacialdim; l++) {
					deltatensor_mapped(0, j, k, l, i) = delta_in[count];
					count = count + 1;
				}
			}
		}
	}

	count = 0;
	for (int h = 0; h < 7; h++) {
		for (int j = 0; j < spacialdim; j++) {
			for (int k = 0; k < spacialdim; k++) {
				for (int l = 0; l < spacialdim; l++) {
					train_loadtensor_mapped(0, j, k, l, h) = train_in[count];
					count = count + 1;
				}
			}
		}
	}

	count = 0;
	for (int i = 0; i < 3; i++) {
		for (int h = 0; h < 3; h++) {
			for (int j = 0; j < spacialdim; j++) {
				for (int k = 0; k < spacialdim; k++) {
					for (int l = 0; l < spacialdim; l++) {
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