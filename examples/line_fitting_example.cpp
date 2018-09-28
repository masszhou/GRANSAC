#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <random>

#include "GRANSAC.hpp"
#include "line_model.hpp"

float Slope(int x0, int y0, int x1, int y1)
{
	return (float)(y1 - y0) / (x1 - x0);
}

void drawFullLine(cv::Mat& img, cv::Point a, cv::Point b, cv::Scalar color, int line_width)
{
	float slope = Slope(a.x, a.y, b.x, b.y);

	cv::Point p(0, 0), q(img.cols, img.rows);

	p.y = -(a.x - p.x) * slope + a.y;
	q.y = -(b.x - q.x) * slope + b.y;

	cv::line(img, p, q, color, line_width, cv::LINE_AA, 0);
}

int main(int argc, char * argv[])
{
	if (argc != 1 && argc != 3){
		std::cout << "[ USAGE ]: " << argv[0] << " [<Image Size> = 1000] [<nPoints> = 500]" << std::endl;
		return -1;
	}

	int side = 400;
	int n_points = 50;
	if (argc == 3){
		side = std::atoi(argv[1]);
		n_points = std::atoi(argv[2]);
	}

	cv::Mat img_canvas(side, side, CV_8UC3);
	img_canvas.setTo(255);

	// Randomly generate points in a 2D plane roughly aligned in a line for testing
	std::random_device seed_device;
	std::mt19937 RNG = std::mt19937(seed_device());

	std::uniform_int_distribution<int> uni_dist(0, side - 1); // [Incl, Incl]
	int perturb = 25;
	std::normal_distribution<float> perturb_dist(0, perturb);

	std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> cand_points;
	for (int i = 0; i < n_points; ++i)
	{
		int diag = uni_dist(RNG);
		cv::Point pt(floor(diag + perturb_dist(RNG)), floor(diag + perturb_dist(RNG)));
		cv::circle(img_canvas, pt, floor(side / 100) + 3, cv::Scalar(0, 0, 0), 2, cv::LINE_AA);

		std::shared_ptr<GRANSAC::AbstractParameter> cand_pt = std::make_shared<GRANSAC::Point2D>(pt.x, pt.y);
		cand_points.push_back(cand_pt);
	}

	GRANSAC::RANSAC<GRANSAC::Line2DModel, 2> estimator;
	estimator.initialize(20, 100); // Threshold, iterations
    
	float time_average = 0;
    for (int i=0; i<100; i++){
        int start = cv::getTickCount();
        estimator.estimate(cand_points);
        int end = cv::getTickCount();
        float time_of_trial = float(end - start) / float(cv::getTickFrequency()) * 1000.0; // [ms]
        time_average += time_of_trial;
    }
    std::cout << "RANSAC took, average time in 100 trials: " << time_average/100 << " ms." << std::endl;

	auto best_inliers = estimator.getBestInliers();
	if (best_inliers.size() > 0)
	{
		for (auto& inlier : best_inliers)
		{
			auto RPt = std::dynamic_pointer_cast<GRANSAC::Point2D>(inlier);
			cv::Point pt(floor(RPt->m_point2D[0]), floor(RPt->m_point2D[1]));
			cv::circle(img_canvas, pt, floor(side / 100), cv::Scalar(0, 255, 0), -1, cv::LINE_AA);
		}
	}

	auto best_line = estimator.getBestModel();
	if (best_line)
	{
		auto best_line_pt1 = std::dynamic_pointer_cast<GRANSAC::Point2D>(best_line->getModelDefParams()[0]);
		auto best_line_pt2 = std::dynamic_pointer_cast<GRANSAC::Point2D>(best_line->getModelDefParams()[1]);
		if (best_line_pt1 && best_line_pt2)
		{
			cv::Point pt1(best_line_pt1->m_point2D[0], best_line_pt1->m_point2D[1]);
			cv::Point pt2(best_line_pt2->m_point2D[0], best_line_pt2->m_point2D[1]);
			drawFullLine(img_canvas, pt1, pt2, cv::Scalar(0, 0, 255), 2);
		}
	}

	while (true)
	{
		cv::imshow("RANSAC Example", img_canvas);

		char Key = cv::waitKey(1);
		if (Key == 27)
			return 0;
		if (Key == ' ')
			cv::imwrite("LineFitting.png", img_canvas);
	}

	return 0;
}
