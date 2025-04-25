#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdexcept>

// 函数声明
std::vector<double> solveThetaCM(double gamma, double theta_lab_rad);
void printResults(double gamma, double theta_lab_deg, const std::vector<double>& solutions);

int main() {
    double Ap, At, theta_lab_deg, gamma;
    // 用户输入
    std::cout << "请输入Ap, At 的值: ";
    std::cin >> Ap >>At;
    gamma = Ap/At;
    std::cout << "请输入θ_lab (theta_lab) 的值（单位：度）: ";
    std::cin >> theta_lab_deg;

    // 将角度转换为弧度
    double theta_lab_rad = theta_lab_deg * M_PI / 180.0;

    try {
        // 计算θ_c.m.的解
        std::vector<double> solutions = solveThetaCM(gamma, theta_lab_rad);

        // 输出结果
        printResults(gamma, theta_lab_deg, solutions);
    } catch (const std::runtime_error& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }

    return 0;
}

// 函数定义：解θ_c.m.
std::vector<double> solveThetaCM(double gamma, double theta_lab_rad) {
    double cos_theta_lab = cos(theta_lab_rad);
    double cos2_theta_lab = cos_theta_lab * cos_theta_lab;
    double sin2_theta_lab = 1.0 - cos2_theta_lab;

    // 二次方程的系数
    double A = 1.0;
    double B = 2.0 * gamma * sin2_theta_lab;
    double C = gamma * gamma * sin2_theta_lab - cos2_theta_lab;

    // 判别式
    double discriminant = B * B - 4.0 * A * C;

    if (discriminant < 0) {
        throw std::runtime_error("无实数解，输入的参数可能导致判别式为负。");
    }

    // 计算两个解
    double sqrt_discriminant = sqrt(discriminant);
    double x1 = (-B + sqrt_discriminant) / (2.0 * A);
    double x2 = (-B - sqrt_discriminant) / (2.0 * A);

    // 检查解是否在[-1, 1]范围内
    std::vector<double> valid_solutions;
    if (x1 >= -1.0 && x1 <= 1.0) {
        valid_solutions.push_back(acos(x1));
    }
    if (x2 >= -1.0 && x2 <= 1.0) {
        valid_solutions.push_back(acos(x2));
    }

    if (valid_solutions.empty()) {
        throw std::runtime_error("无有效解，可能输入的参数导致解超出[-1, 1]范围。");
    }

    return valid_solutions;
}

// 函数定义：打印结果
void printResults(double gamma, double theta_lab_deg, const std::vector<double>& solutions) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "对于 γ = " << gamma << " 和 θ_lab = " << theta_lab_deg << " 度，可能的 θ_c.m. 值为：" << std::endl;

    for (size_t i = 0; i < solutions.size(); ++i) {
        double theta_cm_deg = solutions[i] * 180.0 / M_PI;
        std::cout << "解 " << i + 1 << ": " << theta_cm_deg << " 度" << std::endl;
    }
}
