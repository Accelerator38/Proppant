// Standalone illustration of the tridiagonal sweep used by the historical Qt model.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<double> solve_tridiagonal(
    const std::vector<double>& lower,
    const std::vector<double>& diagonal,
    const std::vector<double>& upper,
    const std::vector<double>& rhs) {
    const auto n = diagonal.size();
    if (n == 0 || lower.size() != n || upper.size() != n || rhs.size() != n)
        throw std::invalid_argument("inconsistent tridiagonal dimensions");
    std::vector<double> c(n), d(n), x(n);
    if (std::abs(diagonal[0]) < 1e-14) throw std::invalid_argument("zero pivot");
    c[0] = upper[0] / diagonal[0];
    d[0] = rhs[0] / diagonal[0];
    for (std::size_t i = 1; i < n; ++i) {
        const double pivot = diagonal[i] - lower[i] * c[i - 1];
        if (std::abs(pivot) < 1e-14) throw std::invalid_argument("zero pivot");
        c[i] = upper[i] / pivot;
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / pivot;
    }
    x[n - 1] = d[n - 1];
    for (std::size_t i = n - 1; i-- > 0;)
        x[i] = d[i] - c[i] * x[i + 1];
    return x;
}

int main(int argc, char** argv) {
    try {
        const int cells = argc > 1 ? std::stoi(argv[1]) : 21;
        if (cells < 2 || cells > 100000) throw std::invalid_argument("cells must be in [2, 100000]");
        const std::string path = argc > 2 ? argv[2] : "pressure.csv";
        std::vector<double> lower(cells, -1), diagonal(cells, 2), upper(cells, -1), rhs(cells, 0);
        // Steady 1D constant-conductivity pressure: P'' = 0, P(0)=100, P(L)=0.
        diagonal.front() = diagonal.back() = 1;
        lower.front() = upper.front() = lower.back() = upper.back() = 0;
        rhs.front() = 100;
        const auto pressure = solve_tridiagonal(lower, diagonal, upper, rhs);
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot open output file");
        output << "cell,pressure\n";
        for (int i = 0; i < cells; ++i) output << i << ',' << pressure[i] << '\n';
        double max_error = 0;
        for (int i = 0; i < cells; ++i)
            max_error = std::max(max_error, std::abs(pressure[i] - 100.0 * (cells - 1 - i) / (cells - 1)));
        std::cout << "max error against linear profile: " << max_error << '\n';
        return max_error < 1e-8 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
