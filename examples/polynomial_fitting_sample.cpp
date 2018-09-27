#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <random>
#include <limits>

#include "GRANSAC.hpp"
#include "LineModel.hpp"
#include "quadratic_model.hpp"
#include "polynomial_model.hpp"

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

    int Side = 400;
    int nPoints = 50;
    if (argc == 3)
    {
        Side = std::atoi(argv[1]);
        nPoints = std::atoi(argv[2]);
    }

    cv::Mat Canvas(Side, Side, CV_8UC3);
    Canvas.setTo(255);

    // Randomly generate points in a 2D plane roughly aligned in a line for testing
    std::random_device SeedDevice;
    std::mt19937 RNG = std::mt19937(SeedDevice());

    std::uniform_int_distribution<int> UniDist(0, Side - 1); // [Incl, Incl]
    int Perturb = 25;
    std::normal_distribution<GRANSAC::VPFloat> PerturbDist(0, Perturb);

    std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> CandPoints;
    while (nPoints > 0){
        int Diag = UniDist(RNG);
        double x = Diag + PerturbDist(RNG);
        // y = 2100-20x+0.05*x*x
        double y = 2100 - 20 * x + 0.05 * x * x + PerturbDist(RNG);

        if (x<=Side && x > 0 && y<=Side && y > 0){
            cv::Point Pt(floor(x), floor(y));
            cv::circle(Canvas, Pt, floor(Side / 100) + 2, cv::Scalar(0, 0, 0), 2, cv::LINE_AA);

            std::shared_ptr<GRANSAC::AbstractParameter> CandPt = std::make_shared<GRANSAC::Point2D>(Pt.x, Pt.y);
            CandPoints.push_back(CandPt);

            nPoints -= 1;
        }
    }

    //set grid
    std::map<string, float> additional_params = {{"img_width", float(Side)}, {"img_height", float(Side)}, {"grid_num_x", float(Side/10)}, {"grid_num_y", float(Side/10)}};
    //draw grid
    cv::Mat overlay;
    double alpha = 0.3;
    Canvas.copyTo(overlay);
    for (int i=0; i < Side/10; i++){
        cv::Point pt_1(i*10, 0);
        cv::Point pt_2(i*10, Side);
        cv::line(overlay, pt_1, pt_2, cv::Scalar(0, 100, 0), 1, cv::LINE_AA, 0);
        cv::Point pt_3(0, i*10);
        cv::Point pt_4(Side, i*10);
        cv::line(overlay, pt_3, pt_4, cv::Scalar(0, 100, 0), 1, cv::LINE_AA, 0);
    }
    cv::addWeighted(overlay, alpha, Canvas, 1 - alpha, 0, Canvas);
    
    // estimate parameters with RANSAC
    int param_num = 4; // 3rd order
    GRANSAC::RANSAC<GRANSAC::PolynomialModel<4>, 4> Estimator;
    Estimator.Initialize(1, 100, additional_params); // Threshold, iterations, the threshold is p-1 distance on grid, 1 means 1 grid cell distance

    float time_average = 0;
    for (int i=0; i<100; i++){
        int start = cv::getTickCount();
        Estimator.Estimate(CandPoints);
        int end = cv::getTickCount();
        float time_of_trial = GRANSAC::VPFloat(end - start) / GRANSAC::VPFloat(cv::getTickFrequency()) * 1000.0; // [ms]
        time_average += time_of_trial;
    }
    std::cout << "RANSAC took, average time in 100 trials: " << time_average/100 << " ms." << std::endl;

    // draw definition points from best model
    auto BestInliers = Estimator.GetBestInliers();
    if (BestInliers.size() > 0)
    {
        for (auto &Inlier : BestInliers)
        {
            auto RPt = std::dynamic_pointer_cast<GRANSAC::Point2D>(Inlier);
            cv::Point Pt(floor(RPt->m_Point2D[0]), floor(RPt->m_Point2D[1]));
            cv::circle(Canvas, Pt, floor(Side / 100), cv::Scalar(0, 255, 0), -1, cv::LINE_AA);
        }
    }

    // draw curve from best model
    auto BestLine = Estimator.GetBestModel();
    if (BestLine)
    {
        auto BestLinePt1 = std::dynamic_pointer_cast<GRANSAC::Point2D>(BestLine->GetModelParams()[0]);
        auto BestLinePt2 = std::dynamic_pointer_cast<GRANSAC::Point2D>(BestLine->GetModelParams()[1]);
        auto BestLinePt3 = std::dynamic_pointer_cast<GRANSAC::Point2D>(BestLine->GetModelParams()[2]);

        std::vector<double> x_values, y_values, coeff;
        for (int i=0; i<BestLine->GetModelParams().size(); i++){
            auto BestLinePt = std::dynamic_pointer_cast<GRANSAC::Point2D>(BestLine->GetModelParams()[i]);
            x_values.push_back(BestLinePt->m_Point2D[0]); //x
            y_values.push_back(BestLinePt->m_Point2D[1]); //x
            cv::Point Pt(floor(BestLinePt->m_Point2D[0]), floor(BestLinePt->m_Point2D[1]));
            cv::circle(Canvas, Pt, floor(Side / 100), cv::Scalar(0, 0, 255), -1, cv::LINE_AA);
        }

        GRANSAC::QuadraticModel::polyfit(x_values, y_values, coeff, param_num-1);
        cout << "original  coefficients: a0=2100, a1=-20, a2=0.05" << endl;
        cout << "estimated coefficients: " ;
        for (int i=0; i < coeff.size(); i++){
            cout << "a" << i << " = " << coeff[i] << ", ";
        }
        cout << endl;

        // draw polynomial line
        float start_point_x = 0;
        float end_point_x = Side;
        vector<cv::Point2f> curvePoints;
        //Define the curve through equation. In this example, a simple parabola
        for (float x = start_point_x; x <= end_point_x; x+=1){
            float y = 0;
            for (int j = 0; j < coeff.size(); j++){
                y += coeff[j]*std::pow(x, float(j));
            } 

            cv::Point2f new_point = cv::Point2f(x, y);                  //resized to better visualize
            curvePoints.push_back(new_point);                       //add point to vector/list
        }
        for (int i = 0; i < curvePoints.size() - 1; i++){
            cv::line(Canvas, curvePoints[i], curvePoints[i + 1], cv::Scalar(0,255,0), 2, CV_AA);
        }
    }

    while (true)
    {
        cv::imshow("RANSAC Example", Canvas);

        char Key = cv::waitKey(1);
        if (Key == 27)
            return 0;
    }

    return 0;
}