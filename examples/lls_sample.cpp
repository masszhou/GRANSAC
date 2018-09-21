#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <random>

#include "GRANSAC.hpp"
#include "LineModel.hpp"

#include <vector>
#include <Eigen/QR>

using namespace std;

void polyfit(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff, int order)
{
    Eigen::MatrixXd A(xv.size(), order+1);
    Eigen::VectorXd yv_mapped = Eigen::VectorXd::Map(&yv.front(), yv.size());
    Eigen::VectorXd result;

    assert(xv.size() == yv.size());
    assert(xv.size() >= order+1);

    // create matrix
    for (size_t i = 0; i < xv.size(); i++)
    for (size_t j = 0; j < order+1; j++)
        A(i, j) = pow(xv.at(i), j);

    // solve for linear least squares fit
    result = A.householderQr().solve(yv_mapped);

    coeff.resize(order+1);
    for (size_t i = 0; i < order+1; i++)
        coeff[i] = result[i];
}

GRANSAC::VPFloat Slope(int x0, int y0, int x1, int y1)
{
    return (GRANSAC::VPFloat)(y1 - y0) / (x1 - x0);
}

void DrawFullLine(cv::Mat &img, cv::Point a, cv::Point b, cv::Scalar color, int LineWidth)
{
    GRANSAC::VPFloat slope = Slope(a.x, a.y, b.x, b.y);

    cv::Point p(0, 0), q(img.cols, img.rows);

    p.y = -(a.x - p.x) * slope + a.y;
    q.y = -(b.x - q.x) * slope + b.y;

    cv::line(img, p, q, color, LineWidth, cv::LINE_AA, 0);
}

int main(int argc, char *argv[])
{
    if (argc != 1 && argc != 3)
    {
        std::cout << "[ USAGE ]: " << argv[0] << " [<Image Size> = 1000] [<nPoints> = 500]" << std::endl;
        return -1;
    }

    int Side = 400;
    int nPoints = 100;
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

    // for polyfit
    std::vector<double> x_values, y_values, coeff;

    std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> CandPoints;
    for (int i = 0; i < nPoints; ++i)
    {
        int Diag = UniDist(RNG);
        double x = Diag + PerturbDist(RNG);
        // y = 2100-20x+0.05*x*x
        double y = 2100 - 20 * x + 0.05 * x * x + PerturbDist(RNG);
        // cv::Point Pt(floor(Diag + PerturbDist(RNG)), floor(Diag + PerturbDist(RNG)));
        cv::Point Pt(floor(x), floor(y));
        cv::circle(Canvas, Pt, floor(Side / 100) + 3, cv::Scalar(0, 0, 0), 2, cv::LINE_AA);

        std::shared_ptr<GRANSAC::AbstractParameter> CandPt = std::make_shared<GRANSAC::Point2D>(Pt.x, Pt.y);
        CandPoints.push_back(CandPt);

        x_values.push_back(x);
        y_values.push_back(y);
    }
    polyfit(x_values, y_values, coeff, 2);
    double m_a0 = coeff[0];
	double m_a1 = coeff[1];
	double m_a2 = coeff[2];
    cout << m_a0 << ", " << m_a1 << ", " << m_a2 << endl;

    // draw polynomial line
    float start_point_x = 0;
    float end_point_x = 400;
    vector<cv::Point2f> curvePoints;
    //Define the curve through equation. In this example, a simple parabola
    for (float x = start_point_x; x <= end_point_x; x+=1){
        float y = m_a2*x*x + m_a1*x + m_a0;
        cv::Point2f new_point = cv::Point2f(x, y);                  //resized to better visualize
        curvePoints.push_back(new_point);                       //add point to vector/list
    }
    for (int i = 0; i < curvePoints.size() - 1; i++){
        cv::line(Canvas, curvePoints[i], curvePoints[i + 1], cv::Scalar(0,255,0), 2, CV_AA);
    }


    //if size=400x400, grid=40x40,
    std::vector<std::vector<int> > m_dist_lookup;
    vector<int> row_init(40, 0);
    for (int i = 0; i < 40; i++)
        m_dist_lookup.push_back(row_init);
    
    float y_0 = m_a0;
    for (int i=1; i<=40; i++){
        int y_0_ind = floor(y_0 / 400) + 1;

        float x = float(i*10); // multiple cell size
        float y_1 = m_a0 + m_a1 * x + m_a2 * x * x;
        int y_1_ind = floor(y_1 / 400);

        if (y_0 <= 400 && y_1 >=0){
            cout << "x = " << x << ", y = " << y_1 << endl;
            cv::line(Canvas, cv::Point(int(x),0), cv::Point(int(x), 400), cv::Scalar(255, 0,0), 1, CV_AA);
            if (y_0 > y_1){
                for (int j=y_1; j<=y_0_ind; j++){
                    m_dist_lookup[j][i] = 1;
                }
            }
        }
        y_0 = y_1;
    }


    for (int i = 0; i < 40; i++){
        for (int j=0; j<40; j++){
            cout << m_dist_lookup[i][j] << ", ";
        }
        cout << endl;
    }

    // cv::Mat Canvas2(Side, Side, CV_8UC3);
    // Canvas2.setTo(255);

    // GRANSAC::RANSAC<GRANSAC::Line2DModel, 2> Estimator;
    // Estimator.Initialize(20, 100); // Threshold, iterations
    // int start = cv::getTickCount();
    // Estimator.Estimate(CandPoints);
    // int end = cv::getTickCount();
    // std::cout << "RANSAC took: " << GRANSAC::VPFloat(end - start) / GRANSAC::VPFloat(cv::getTickFrequency()) * 1000.0 << " ms." << std::endl;

    // auto BestInliers = Estimator.GetBestInliers();
    // if (BestInliers.size() > 0)
    // {
    //     for (auto &Inlier : BestInliers)
    //     {
    //         auto RPt = std::dynamic_pointer_cast<GRANSAC::Point2D>(Inlier);
    //         cv::Point Pt(floor(RPt->m_Point2D[0]), floor(RPt->m_Point2D[1]));
    //         cv::circle(Canvas, Pt, floor(Side / 100), cv::Scalar(0, 255, 0), -1, cv::LINE_AA);
    //     }
    // }

    // auto BestLine = Estimator.GetBestModel();
    // if (BestLine)
    // {
    //     auto BestLinePt1 = std::dynamic_pointer_cast<GRANSAC::Point2D>(BestLine->GetModelParams()[0]);
    //     auto BestLinePt2 = std::dynamic_pointer_cast<GRANSAC::Point2D>(BestLine->GetModelParams()[1]);
    //     if (BestLinePt1 && BestLinePt2)
    //     {
    //         cv::Point Pt1(BestLinePt1->m_Point2D[0], BestLinePt1->m_Point2D[1]);
    //         cv::Point Pt2(BestLinePt2->m_Point2D[0], BestLinePt2->m_Point2D[1]);
    //         DrawFullLine(Canvas, Pt1, Pt2, cv::Scalar(0, 0, 255), 2);
    //     }
    // }

    while (true)
    {
        cv::imshow("RANSAC Example", Canvas);

        char Key = cv::waitKey(1);
        if (Key == 27)
            return 0;
    }

    return 0;
}