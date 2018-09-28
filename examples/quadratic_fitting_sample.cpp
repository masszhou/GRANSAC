#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <random>
#include <limits>

#include "GRANSAC.hpp"
#include "quadratic_model.hpp"

#include <vector>
#include <Eigen/QR>

using namespace std;


int main(int argc, char *argv[])
{
    if (argc != 1 && argc != 3)
    {
        std::cout << "[ USAGE ]: " << argv[0] << " [<Image Size> = 1000] [<nPoints> = 500]" << std::endl;
        return -1;
    }

    int side = 400;
    int n_points = 50;
    if (argc == 3)
    {
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
    while (n_points > 0){
        int diag = uni_dist(RNG);
        double x = diag + perturb_dist(RNG);
        // y = 2100-20x+0.05*x*x
        double y = 2100 - 20 * x + 0.05 * x * x + perturb_dist(RNG);

        if (x<=side && x > 0 && y<=side && y > 0){
            cv::Point pt(floor(x), floor(y));
            cv::circle(img_canvas, pt, floor(side / 100) + 2, cv::Scalar(0, 0, 0), 2, cv::LINE_AA);

            std::shared_ptr<GRANSAC::AbstractParameter> cand_pt = std::make_shared<GRANSAC::Point2D>(pt.x, pt.y);
            cand_points.push_back(cand_pt);

            n_points -= 1;
        }
    }

    //set grid
    std::map<string, float> additional_params = {{"img_width", float(side)}, {"img_height", float(side)}, {"grid_num_x", float(side/10)}, {"grid_num_y", float(side/10)}};
    //draw grid
    cv::Mat img_overlay;
    double alpha = 0.3;
    img_canvas.copyTo(img_overlay);
    for (int i=0; i < side/10; i++){
        cv::Point pt_1(i*10, 0);
        cv::Point pt_2(i*10, side);
        cv::line(img_overlay, pt_1, pt_2, cv::Scalar(0, 100, 0), 1, cv::LINE_AA, 0);
        cv::Point pt_3(0, i*10);
        cv::Point pt_4(side, i*10);
        cv::line(img_overlay, pt_3, pt_4, cv::Scalar(0, 100, 0), 1, cv::LINE_AA, 0);
    }
    cv::addWeighted(img_overlay, alpha, img_canvas, 1 - alpha, 0, img_canvas);
    
    GRANSAC::RANSAC<GRANSAC::QuadraticModel, 3> estimator;
    estimator.initialize(1, 100, additional_params); // Threshold, iterations, the threshold is p-1 distance on grid, 1 means 1 grid cell distance

    float time_average = 0;
    // for (int i=0; i<500; i++){
        int start = cv::getTickCount();
        estimator.estimate(cand_points);
        int end = cv::getTickCount();
        float time_of_trial = float(end - start) / float(cv::getTickFrequency()) * 1000.0; // [ms]
        time_average += time_of_trial;
    // }
    // std::cout << "RANSAC took, average time in 500 trials: " << time_average/500 << " ms." << std::endl;
    std::cout << "RANSAC took, average time in 500 trials: " << time_average << " ms." << std::endl;

    auto best_inliers = estimator.getBestInliers();
    if (best_inliers.size() > 0)
    {
        for (auto &inlier : best_inliers)
        {
            auto RPt = std::dynamic_pointer_cast<GRANSAC::Point2D>(inlier);
            cv::Point pt(floor(RPt->m_point2D[0]), floor(RPt->m_point2D[1]));
            cv::circle(img_canvas, pt, floor(side / 100), cv::Scalar(0, 255, 0), -1, cv::LINE_AA);
        }
    }

    auto best_line = estimator.getBestModel();
    if (best_line)
    {
        std::vector<double> x_values, y_values, coeff;
        for (int i=0; i<3; i++){
            auto best_line_pt = std::dynamic_pointer_cast<GRANSAC::Point2D>(best_line->getModelDefParams()[i]);
            x_values.push_back(best_line_pt->m_point2D[0]); //x
            y_values.push_back(best_line_pt->m_point2D[1]); //x
            cv::Point pt(floor(best_line_pt->m_point2D[0]), floor(best_line_pt->m_point2D[1]));
            cv::circle(img_canvas, pt, floor(side / 100), cv::Scalar(0, 0, 255), -1, cv::LINE_AA);
        }
        GRANSAC::QuadraticModel::polyfit(x_values, y_values, coeff, 2);

        cout << "original  coefficients: a0=2100, a1=-20, a2=0.05" << endl;
        cout << "estimated coefficients: a0="<<coeff[0]<<", a1="<<coeff[1]<<", a2="<<coeff[2]<<endl;

        // draw polynomial line
        float start_point_x = 0;
        float end_point_x = side;
        vector<cv::Point2f> curve_points;
        //Define the curve through equation. In this example, a simple parabola
        for (float x = start_point_x; x <= end_point_x; x+=1){
            float y = coeff[2]*x*x + coeff[1]*x + coeff[0];
            cv::Point2f new_point = cv::Point2f(x, y);                  //resized to better visualize
            curve_points.push_back(new_point);                       //add point to vector/list
        }
        for (int i = 0; i < curve_points.size() - 1; i++){
            cv::line(img_canvas, curve_points[i], curve_points[i + 1], cv::Scalar(0,255,0), 2, CV_AA);
        }
    }

    while (true)
    {
        cv::imshow("RANSAC Example", img_canvas);

        char Key = cv::waitKey(1);
        if (Key == 27)
            return 0;
    }

    return 0;
}