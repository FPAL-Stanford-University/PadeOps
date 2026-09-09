#include <torch/script.h> // One-stop header.
#include <iostream>
#include <memory>
#include <cmath>

extern "C" float* loadpytorchgpu(double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, double* rot_in, int* space1, int* space2, int* space3, double* dx_p, double* dy_p, double* dz_p);
extern "C" void free(void* ptr);

float* loadpytorchgpu(double* train_in, double* strain2_in, double* strain2_m_in, double* rot2_in, double* strainrot_m_in, double* rotstrainrot_m_in, double* strain2rot_m_in, double* rotstrain2_m_in, double* rot_in, int* space1, int* space2, int* space3, double* dx_p, double* dy_p, double* dz_p) {
	// Initialize a tensorflow session
	// std::cout << "Entered Function!\n";

	int space1a = *space1;
	int space2a = *space2;
	int space3a = *space3;

	float dx = (float)*dx_p;
	float dy = (float)*dy_p;
	float dz = (float)*dz_p;

	float myFloat = -1.0f;
	float* rtype = &myFloat;

	torch::NoGradGuard no_grad;

	// create torch module opbject
	torch::jit::script::Module module;
	try {
		// Deserialize the ScriptModule from a file using torch::jit::load().
		module = torch::jit::load("gnnmodel.pt");
	}
	catch (const c10::Error& e) {
		std::cerr << "Error loading the model\n";
		return rtype;
	}

	module.eval();

	// std::cout << "GNN Model is loaded!\n";
	// tensor options
	auto options = torch::TensorOptions().dtype(torch::kFloat32);
	auto options2 = torch::TensorOptions().dtype(torch::kInt64);

	// torch::Tensor train_loadtensor = torch::zeros({ space1a * space2a * space3a, 5 }, options);
	torch::Tensor train_loadtensor = torch::zeros({ space1a * space2a * space3a, 3 }, options);
	torch::Tensor edge_indextensor_long = torch::zeros({ 2, 27 * space1a * space2a * space3a }, options2);
	torch::Tensor edge_featurestensor_long = torch::zeros({ 27 * space1a * space2a * space3a, 1 }, options);
	torch::Tensor strain2tensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor strain2_mtensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor rot2tensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor strainrot_mtensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor rotstrainrot_mtensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor strain2rot_mtensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor rotstrain2_mtensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);
	torch::Tensor rottensor = torch::zeros({ space1a * space2a * space3a, 3, 3 }, options);

	std::vector<int64_t> v(space1a * space2a * space3a);
	std::iota(v.begin(), v.end(), 0);

	torch::Tensor test = torch::from_blob(v.data(), { space1a, space2a, space3a }, options2);
	// std::cout << "Test_idx is completed!\n";

	// auto test = test_v.view({ space1a, space2a, space3a });

	auto train_loadtensor_a = train_loadtensor.accessor<float, 2>();
	auto strain2tensor_a = strain2tensor.accessor<float, 3>();
	auto strain2_mtensor_a = strain2_mtensor.accessor<float, 3>();
	auto rot2tensor_a = rot2tensor.accessor<float, 3>();
	auto strainrot_mtensor_a = strainrot_mtensor.accessor<float, 3>();
	auto rotstrainrot_mtensor_a = rotstrainrot_mtensor.accessor<float, 3>();
	auto strain2rot_mtensor_a = strain2rot_mtensor.accessor<float, 3>();
	auto rotstrain2_mtensor_a = rotstrain2_mtensor.accessor<float, 3>();
	auto rottensor_a = rottensor.accessor<float, 3>();
	auto edge_indextensor_long_a = edge_indextensor_long.accessor<int64_t, 2>();
	auto edge_featurestensor_long_a = edge_featurestensor_long.accessor<float, 2>();
	auto test_a = test.accessor<int64_t, 3>();


	// generate edge indexes and edge features tensor 
	int idx = 0;
	for (int i = 0; i < space1a; i++) {
		for (int j = 0; j < space2a; j++) {
			for (int k = 0; k < space3a; k++) {
				// check x
				if (i + 1 < space1a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j][k];
					edge_featurestensor_long_a[idx][0] = dx;
					idx = idx + 1;
				}

				if (i - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j][k];
					edge_featurestensor_long_a[idx][0] = dx;
					idx = idx + 1;
				}
				// check y
				if (j + 1 < space2a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j + 1][k];
					edge_featurestensor_long_a[idx][0] = dy;
					idx = idx + 1;
				}

				if (j - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j - 1][k];
					edge_featurestensor_long_a[idx][0] = dy;
					idx = idx + 1;
				}
				// check z
				if (k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j][k + 1];
					edge_featurestensor_long_a[idx][0] = dz;
					idx = idx + 1;
				}

				if (k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j][k - 1];
					edge_featurestensor_long_a[idx][0] = dz;
					idx = idx + 1;
				}
				// now do diagonals, check x y
				if (i + 1 < space1a && j + 1 < space2a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j + 1][k];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && j - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j - 1][k];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2));
					idx = idx + 1;
				}

				if (i + 1 < space1a && j - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j - 1][k];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && j + 1 < space2a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j + 1][k];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2));
					idx = idx + 1;
				}

				// check x z
				if (i + 1 < space1a && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i + 1 < space1a && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				// check y z

				if (j + 1 < space2a && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j + 1][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (j - 1 >= 0 && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j - 1][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (j + 1 < space2a && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j + 1][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (j - 1 >= 0 && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i][j - 1][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}
				// now do x yand z
				if (i + 1 < space1a && j + 1 < space2a && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j + 1][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && j + 1 < space2a && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j + 1][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i + 1 < space1a && j - 1 >= 0 && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j - 1][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && j - 1 >= 0 && k + 1 < space3a) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j - 1][k + 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i + 1 < space1a && j + 1 < space2a && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j + 1][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && j + 1 < space2a && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j + 1][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i + 1 < space1a && j - 1 >= 0 && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i + 1][j - 1][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}

				if (i - 1 >= 0 && j - 1 >= 0 && k - 1 >= 0) {
					edge_indextensor_long_a[0][idx] = test_a[i][j][k];
					edge_indextensor_long_a[1][idx] = test_a[i - 1][j - 1][k - 1];
					edge_featurestensor_long_a[idx][0] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
					idx = idx + 1;
				}
			}
		}
	}

	// std::cout << "Edge features all generated!\n";

	torch::Tensor edge_indextensor = edge_indextensor_long.index({ torch::indexing::Slice(0, torch::indexing::None, 1), torch::indexing::Slice(0, idx, 1) });
	torch::Tensor edge_featurestensor = edge_featurestensor_long.index({ torch::indexing::Slice(0, idx, 1), torch::indexing::Slice(0, torch::indexing::None, 1) });

	int count = 0;
	for (int j = 0; j < space1a*space2a*space3a; j++) {
		// changed 5->3
		for (int h = 0; h < 3; h++) {
			train_loadtensor_a[j][h] = (float)train_in[count];
			count = count + 1;
		}
	}

	count = 0;
	for (int j = 0; j < space1a * space2a * space3a; j++) {
		for (int i = 0; i < 3; i++) {
			for (int h = 0; h < 3; h++) {
				strain2tensor_a[j][i][h] = (float)strain2_in[count];
				strain2_mtensor_a[j][i][h] = (float)strain2_m_in[count];
				rot2tensor_a[j][i][h] = (float)rot2_in[count];
				rottensor_a[j][i][h] = (float)rot_in[count];
				strainrot_mtensor_a[j][i][h] = (float)strainrot_m_in[count];
				rotstrainrot_mtensor_a[j][i][h] = (float)rotstrainrot_m_in[count];
				strain2rot_mtensor_a[j][i][h] = (float)strain2rot_m_in[count];
				rotstrain2_mtensor_a[j][i][h] = (float)rotstrain2_m_in[count];
				count = count + 1;
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


	// std::cout << "start eval" << std::endl;

	train_loadtensor = train_loadtensor.to(at::kCUDA);
	strain2tensor = strain2tensor.to(at::kCUDA);
	strain2_mtensor = strain2_mtensor.to(at::kCUDA);
	rot2tensor = rot2tensor.to(at::kCUDA);
	strainrot_mtensor = strainrot_mtensor.to(at::kCUDA);
	rotstrainrot_mtensor = rotstrainrot_mtensor.to(at::kCUDA);
	strain2rot_mtensor = strain2rot_mtensor.to(at::kCUDA);
	rotstrain2_mtensor = rotstrain2_mtensor.to(at::kCUDA);
	rottensor = rottensor.to(at::kCUDA);
	edge_indextensor = edge_indextensor.to(at::kCUDA);
	edge_featurestensor = edge_featurestensor.to(at::kCUDA);

	module.to(at::kCUDA);

	std::vector<torch::jit::IValue> inputs;
	inputs.push_back(train_loadtensor);
	inputs.push_back(strain2tensor);
	inputs.push_back(strain2_mtensor);
	inputs.push_back(rot2tensor);
	inputs.push_back(strainrot_mtensor);
	inputs.push_back(rotstrainrot_mtensor);
	inputs.push_back(strain2rot_mtensor);
	inputs.push_back(rotstrain2_mtensor);
	inputs.push_back(rottensor);
	inputs.push_back(edge_indextensor);
	inputs.push_back(edge_featurestensor);

	at::Tensor output = module.forward(inputs).toTensor();
	at::Tensor output_cpu = output.to(at::kCPU);
	auto output_a = output_cpu.accessor<float, 2>();
	// std::cout << "Model done" << std::endl;

	float* out_arraypass = (float*)malloc(sizeof(float) * 7 * space1a * space2a * space3a);
	count = space1a * space2a * space3a;
	for (int h = 0; h < 6; h++) {
		for (int j = 0; j < space1a * space2a * space3a; j++) {
			out_arraypass[count] = output_a[h][j];
			count = count + 1;
		}
	}
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
	return out_arraypass;
}
