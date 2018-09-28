#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <memory>
#include <algorithm>
#include <vector>
#include <map>
#include <omp.h>

#include "AbstractModel.hpp"

namespace GRANSAC
{
	
// T - AbstractModel
template <class T, int t_num_params>
class RANSAC
{
private:
	std::vector<std::shared_ptr<AbstractParameter>> m_data; // All the data points

	std::vector<std::shared_ptr<T> > m_sampled_models; // Vector of all sampled models
	std::shared_ptr<T> m_best_model; // Pointer to the best model, valid only after Estimate() is called
	std::vector<std::shared_ptr<AbstractParameter>> m_best_inliers;

	int m_max_iterations; // Number of iterations before termination
	float m_threshold; // The threshold for computing model consensus
	float m_best_model_score; // The score of the best model
	int m_best_model_idx;

	std::vector<std::mt19937> m_rand_engines; // Mersenne twister high quality RNG that support *OpenMP* multi-threading

	std::map<std::string, float> m_additional_params;

public:
	RANSAC(void)
	{
		int n_threads = std::max(1, omp_get_max_threads());
		std::cout << "[ INFO ]: Maximum usable threads: " << n_threads << std::endl;
		for (int i = 0; i < n_threads; ++i)
		{
			std::random_device seed_device;
			m_rand_engines.push_back(std::mt19937(seed_device()));
		}

		reset();
	};

	virtual ~RANSAC(void) {};

	void reset(void)
	{
		// Clear sampled models, etc. and prepare for next call. Reset RANSAC estimator state
		m_data.clear();
		m_sampled_models.clear();

		m_best_model_idx = -1;
		m_best_model_score = 0.0;
	};

	void initialize(float threshold, int max_iterations = 1000, 
	                const std::map<std::string, float>& addtional_params = std::map<std::string, float>())
	{
		m_threshold = threshold;
		m_max_iterations = max_iterations;
		m_additional_params = addtional_params;
	};

	std::shared_ptr<T> getBestModel(void) { return m_best_model; };
	const std::vector<std::shared_ptr<AbstractParameter>>& getBestInliers(void) { return m_best_inliers; };

	bool estimate(const std::vector<std::shared_ptr<AbstractParameter>> &data)
	{
		if (data.size() <= t_num_params){
			std::cerr << "[ WARN ]: RANSAC - Number of data points is too less. Not doing anything." << std::endl;
			return false;
		}

		m_data = data;
		int data_size = m_data.size();
		std::uniform_int_distribution<int> uni_dist(0, int(data_size - 1)); // Both inclusive

		std::vector<float> inlier_fraction_accum(m_max_iterations);
		std::vector<std::vector<std::shared_ptr<AbstractParameter>>> inliers_accum(m_max_iterations);
		m_sampled_models.resize(m_max_iterations);

		int sample_num = t_num_params; // deterministic sample
		if ( m_additional_params.count("sample_num") == 1 ) {
			// found
			sample_num = m_additional_params.at("sample_num"); // over-deterministic sample. under-deterministic sample not suppported yet.
		}

		int n_threads = std::max(1, omp_get_max_threads());
		omp_set_dynamic(0); // Explicitly disable dynamic teams
		omp_set_num_threads(n_threads);
#pragma omp parallel for
		for (int i = 0; i < m_max_iterations; ++i)
		{
			// Select t_NumParams random samples
			std::vector<std::shared_ptr<AbstractParameter>> random_samples(sample_num);
			std::vector<std::shared_ptr<AbstractParameter>> remainder_samples = m_data; // Without the chosen random samples

			std::shuffle(remainder_samples.begin(), remainder_samples.end(), m_rand_engines[omp_get_thread_num()]); // To avoid picking the same element more than once
			std::copy(remainder_samples.begin(), remainder_samples.begin() + sample_num, random_samples.begin());

			std::shared_ptr<T> random_model = std::make_shared<T>(random_samples, m_additional_params);

			// Check if the sampled model is the best so far
			std::pair<float, std::vector<std::shared_ptr<AbstractParameter> > > eval_pair = random_model->evaluate(remainder_samples, m_threshold);
			inlier_fraction_accum[i] = eval_pair.first;
			inliers_accum[i] = eval_pair.second;

			// Push back into history. Could be removed later
			m_sampled_models[i] = random_model;
		}

		for (int i = 0; i < m_max_iterations; ++i)
		{
			if (inlier_fraction_accum[i] > m_best_model_score)
			{
				m_best_model_score = inlier_fraction_accum[i];
				m_best_model_idx = m_sampled_models.size() - 1;
				m_best_model = m_sampled_models[i];
				m_best_inliers = inliers_accum[i];
			}
		}

		//std::cerr << "BestInlierFraction: " << m_BestModelScore << std::endl;

		reset();

		return true;
	};
};

} // namespace GRANSAC
