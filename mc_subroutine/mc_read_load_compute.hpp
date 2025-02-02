//
// Created by adada on 2/2/2025.
//

#ifndef MC_READ_LOAD_COMPUTE_HPP
#define MC_READ_LOAD_COMPUTE_HPP
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#endif //MC_READ_LOAD_COMPUTE_HPP

namespace fs = boost::filesystem;
namespace py = boost::python;
namespace np = boost::python::numpy;
class evolution
{


public:
    double T;// temperature
    double beta;
    double a;
    double J;
    int N;
    double h;// step size
    int sweepToWrite;
    int newFlushNum;
    int flushLastFile;
    std::string TDirRoot;
    std::string U_dipole_dataDir;
    std::ranlux24_base e2;
    std::uniform_real_distribution<> distUnif01;
    int sweep_multiple;
    std::string out_U_path;
    std::string out_Px_path;

    std::string out_Py_path;
    std::string out_Qx_path;
    std::string out_Qy_path;

    double * U_data_ptr;//all U data
    double * Px_ptr;//all Px data
    double * Py_ptr;//all Py data
    double * Qx_ptr;//all Qx data
    double * Qy_ptr;//all Qy data

};