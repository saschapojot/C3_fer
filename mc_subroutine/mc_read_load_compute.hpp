//
// Created by adada on 2/2/2025.
//

#ifndef MC_READ_LOAD_COMPUTE_HPP
#define MC_READ_LOAD_COMPUTE_HPP
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cfenv> // for floating-point exceptions
#include <fstream>
#include <Python.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#endif //MC_READ_LOAD_COMPUTE_HPP

namespace fs = boost::filesystem;
namespace py = boost::python;
namespace np = boost::python::numpy;

class mc_computation
{
public:
    mc_computation(const std::string& cppInParamsFileName): e2(std::random_device{}()), distUnif01(0.0, 1.0)
    {
        std::ifstream file(cppInParamsFileName);
        if (!file.is_open())
        {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }
        std::string line;
        int paramCounter = 0;
        while (std::getline(file, line))
        {
            // Check if the line is empty
            if (line.empty())
            {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);

            //read T
            if (paramCounter == 0)
            {
                iss >> T;
                if (T <= 0)
                {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                } //end if
                std::cout << "T=" << T << std::endl;
                this->beta = 1.0 / T;
                std::cout << "beta=" << beta << std::endl;
                paramCounter++;
                continue;
            } //end T
            // read a
            if (paramCounter == 1)
            {
                iss >> a;
                if (a <= 0)
                {
                    std::cerr << "a must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "a=" << a << std::endl;
                paramCounter++;
                continue;
            } //end a

            //read J
            if (paramCounter == 2)
            {
                iss >> J;
                std::cout << "J=" << J << std::endl;
                paramCounter++;
                continue;
            } //end J

            //read N
            if (paramCounter == 3)
            {
                iss >> N0;
                N1 = N0;
                if (N0 <= 0)
                {
                    std::cerr << "N must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "N0=N1=" << N0 << std::endl;
                paramCounter++;
                continue;
            } //end N

            //read q
            if (paramCounter == 4)
            {
                iss >> q;
                if (q <= 0)
                {
                    std::cerr << "q must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "q=" << q << std::endl;
                paramCounter++;
                continue;
            } //end q
            //read alpha1
            if (paramCounter == 5)
            {
                iss >> alpha1;
                std::cout << "alpha1=" << alpha1 << std::endl;
                paramCounter++;
                continue;
            } //end alpha1

            //read alpha2
            if (paramCounter == 6)
            {
                iss >> alpha2;
                std::cout << "alpha2=" << alpha2 << std::endl;
                paramCounter++;
                continue;
            } //end alpha2

            //read alpha3
            if (paramCounter == 7)
            {
                iss >> alpha3;
                std::cout << "alpha3=" << alpha3 << std::endl;
                paramCounter++;
                continue;
            } //end alpha3

            //read alpha4
            if (paramCounter == 8)
            {
                iss >> alpha4;
                std::cout << "alpha4=" << alpha4 << std::endl;
                paramCounter++;
                continue;
            } //end alpha4

            //read alpha5
            if (paramCounter == 9)
            {
                iss >> alpha5;
                std::cout << "alpha5=" << alpha5 << std::endl;
                paramCounter++;
                continue;
            } //end alpha5

            //read alpha6
            if (paramCounter == 10)
            {
                iss >> alpha6;
                std::cout << "alpha6=" << alpha6 << std::endl;
                paramCounter++;
                continue;
            } //end alpha6

            //read alpha7
            if (paramCounter == 11)
            {
                iss >> alpha7;
                std::cout << "alpha7=" << alpha7 << std::endl;
                paramCounter++;
                continue;
            } //end alpha7

            //read sweepToWrite
            if (paramCounter == 12)
            {
                iss >> sweepToWrite;
                if (sweepToWrite <= 0)
                {
                    std::cerr << "sweepToWrite must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "sweepToWrite=" << sweepToWrite << std::endl;
                paramCounter++;
                continue;
            } //end sweepToWrite

            //read newFlushNum
            if (paramCounter == 13)
            {
                iss >> newFlushNum;
                if (newFlushNum <= 0)
                {
                    std::cerr << "newFlushNum must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "newFlushNum=" << newFlushNum << std::endl;
                paramCounter++;
                continue;
            } //end newFlushNum

            //read flushLastFile
            if (paramCounter == 14)
            {
                iss >> flushLastFile;
                std::cout << "flushLastFile=" << flushLastFile << std::endl;
                paramCounter++;
                continue;
            } //end flushLastFile

            //read TDirRoot
            if (paramCounter == 15)
            {
                iss >> TDirRoot;
                std::cout << "TDirRoot=" << TDirRoot << std::endl;
                paramCounter++;
                continue;
            } //end TDirRoot

            //read U_dipole_dataDir
            if (paramCounter == 16)
            {
                iss >> U_dipole_dataDir;
                std::cout << "U_dipole_dataDir=" << U_dipole_dataDir << std::endl;
                paramCounter++;
                continue;
            } //end U_dipole_dataDir

            //read h
            if (paramCounter == 17)
            {
                iss >> h;
                if (h <= 0)
                {
                    std::cerr << "h must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "h=" << h << std::endl;
                paramCounter++;
                continue;
            } //end h

            //read sweep_multiple
            if (paramCounter == 18)
            {
                iss >> sweep_multiple;
                if (sweep_multiple <= 0)
                {
                    std::cerr << "sweep_multiple must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "sweep_multiple=" << sweep_multiple << std::endl;
                paramCounter++;
                continue;
            } //end sweep_multiple
        } //end while
        //allocate memory for data
        try
        {
            this->U_data_all_ptr = std::shared_ptr<double[]>(new double[sweepToWrite],
                                                             std::default_delete<double[]>());

            this->Px_all_ptr = std::shared_ptr<double[]>(new double[sweepToWrite * N0 * N1],
                                                         std::default_delete<double[]>());

            this->Py_all_ptr = std::shared_ptr<double[]>(new double[sweepToWrite * N0 * N1],
                                                         std::default_delete<double[]>());

            this->Qx_all_ptr = std::shared_ptr<double[]>(new double[sweepToWrite * N0 * N1],
                                                         std::default_delete<double[]>());

            this->Qy_all_ptr = std::shared_ptr<double[]>(new double[sweepToWrite * N0 * N1],
                                                         std::default_delete<double[]>());

            this->Px_init = std::shared_ptr<double[]>(new double[N0 * N1],
                                                      std::default_delete<double[]>());
            this->Py_init = std::shared_ptr<double[]>(new double[N0 * N1],
                                                      std::default_delete<double[]>());

            this->Qx_init = std::shared_ptr<double[]>(new double[N0 * N1],
                                                      std::default_delete<double[]>());
            this->Qy_init = std::shared_ptr<double[]>(new double[N0 * N1],
                                                      std::default_delete<double[]>());
        }
        catch (const std::bad_alloc& e)
        {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception& e)
        {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }
        this->out_U_path = this->U_dipole_dataDir + "/U/";
        if (!fs::is_directory(out_U_path) || !fs::exists(out_U_path))
        {
            fs::create_directories(out_U_path);
        }
        this->out_Px_path = this->U_dipole_dataDir + "/Px/";
        if (!fs::is_directory(out_Px_path) || !fs::exists(out_Px_path))
        {
            fs::create_directories(out_Px_path);
        }
        this->out_Py_path = this->U_dipole_dataDir + "/Py/";
        if (!fs::is_directory(out_Py_path) || !fs::exists(out_Py_path))
        {
            fs::create_directories(out_Py_path);
        }
        this->out_Qx_path = this->U_dipole_dataDir + "/Qx/";
        if (!fs::is_directory(out_Qx_path) || !fs::exists(out_Qx_path))
        {
            fs::create_directories(out_Qx_path);
        }
        this->out_Qy_path = this->U_dipole_dataDir + "/Qy/";
        if (!fs::is_directory(out_Qy_path) || !fs::exists(out_Qy_path))
        {
            fs::create_directories(out_Qy_path);
        }
        this->A = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->B = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->C = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->G = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->R = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->Gamma = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->Theta = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);
        this->Lambda = arma::zeros<arma::dmat>(N0 * N1, N0 * N1);

        this->unif_in_0_N0N1 = std::uniform_int_distribution<int>(0, N0 * N1-1);

        this->dipole_upper_bound = 1.0 / 4.0 * 1.0 / std::sqrt(3.0) * a * q;
        this->dipole_lower_bound = -dipole_upper_bound;
        std::cout << "dipole_upper_bound=" << dipole_upper_bound << std::endl;
        std::cout << "dipole_lower_bound=" << dipole_lower_bound << std::endl;
        this->a_squared = std::pow(a, 2.0);
        this->J_over_a_squared = J / a_squared;
        std::cout << "J_over_a_squared=" << J_over_a_squared << std::endl;
        // int random_number=unif_in_0_N0N1(e2);
        //
        // std::cout << "Random number: " << random_number << std::endl;
    } //end constructor
public:
    void init_and_run();
    void execute_mc(const std::shared_ptr<double[]>& Px_vec,
                    const std::shared_ptr<double[]>& Py_vec,
                    const std::shared_ptr<double[]>& Qx_vec,
                    const std::shared_ptr<double[]>& Qy_vec,
                    const int& flushNum);

    void execute_mc_one_sweep(arma::dvec& Px_arma_vec_curr,
                              arma::dvec& Py_arma_vec_curr,
                              arma::dvec& Qx_arma_vec_curr,
                              arma::dvec& Qy_arma_vec_curr,
                              double& UCurr,
                              arma::dvec& Px_arma_vec_next,
                              arma::dvec& Py_arma_vec_next,
                              arma::dvec& Qx_arma_vec_next,
                              arma::dvec& Qy_arma_vec_next);

    ///
    /// @param n0
    /// @param n1
    /// @param Px_arma_vec_curr
    /// @param Px_arma_vec_next
    /// @param Py_arma_vec_curr
    /// @param Qx_arma_vec_curr
    /// @param Qy_arma_vec_curr
    /// @param UCurr
    /// @param UNext
    void HPx_update(const int& flattened_ind, const arma::dvec& Px_arma_vec_curr,
                    const arma::dvec& Px_arma_vec_next,
                    const arma::dvec& Py_arma_vec_curr,
                    const arma::dvec& Qx_arma_vec_curr,
                    const arma::dvec& Qy_arma_vec_curr,
                    double& UCurr, double& UNext);

    void HPy_update(const int& flattened_ind,
                    const arma::dvec& Py_arma_vec_curr, const arma::dvec& Py_arma_vec_next,
                    const arma::dvec& Px_arma_vec_curr, const arma::dvec& Qx_arma_vec_curr,
                    const arma::dvec& Qy_arma_vec_curr,
                    double& UCurr, double& UNext);

    void HQx_update(const int& flattened_ind,
                    const arma::dvec& Qx_arma_vec_curr, const arma::dvec& Qx_arma_vec_next,
                    const arma::dvec& Px_arma_vec_curr,
                    const arma::dvec& Py_arma_vec_curr,
                    const arma::dvec& Qy_arma_vec_curr,
                    double& UCurr, double& UNext);

    void HQy_update(const int& flattened_ind,
                    const arma::dvec& Qy_arma_vec_curr, const arma::dvec& Qy_arma_vec_next,
                    const arma::dvec& Px_arma_vec_curr,
                    const arma::dvec& Py_arma_vec_curr,
                    const arma::dvec& Qx_arma_vec_curr,
                    double& UCurr, double& UNext);

    void init_Px_Py_Qx_Qy();
    void load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr, std::size_t size);

    void save_array_to_pickle(const std::shared_ptr<double[]>& ptr, int size, const std::string& filename);
    void init_mats();
    ///
    /// @param n0
    /// @param n1
    /// @return flatenned index
    int double_ind_to_flat_ind(const int& n0, const int& n1);

    ///
    /// @param x
    /// @param leftEnd
    /// @param rightEnd
    /// @param eps
    /// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
    double generate_uni_open_interval(const double& x, const double& leftEnd, const double& rightEnd,
                                      const double& eps);

    void proposal_uni(const arma::dvec& arma_vec_curr, arma::dvec& arma_vec_next,
                      const int& flattened_ind);

    ///
    /// @param x proposed value
    /// @param y current value
    /// @param a left end of interval
    /// @param b right end of interval
    /// @param epsilon half length
    /// @return proposal probability S(x|y)
    double S_uni(const double& x, const double& y, const double& a, const double& b, const double& epsilon);

    ///
    /// @param n0 index
    /// @param n1 index
    /// @param Px_arma_vec Px
    /// @param Py_arma_vec Py
    /// @return self energy H1
    double H1(const int& flattened_ind, const arma::dvec& Px_arma_vec, const arma::dvec& Py_arma_vec);

    ///
    /// @param n0 index
    /// @param n1 index
    /// @param Qx_arma_vec Qx
    /// @param Qy_arma_vec Qy
    /// @return self energy H2
    double H2(const int& flattened_ind, const arma::dvec& Qx_arma_vec, const arma::dvec& Qy_arma_vec);

    double acceptanceRatio_uni(const arma::dvec& arma_vec_curr,
                               const arma::dvec& arma_vec_next, const int& flattened_ind,
                               const double& UCurr, const double& UNext);

    template <class T>
    void print_shared_ptr(const std::shared_ptr<T>& ptr, const int& size)
    {
        if (!ptr)
        {
            std::cout << "Pointer is null." << std::endl;
            return;
        }

        for (int i = 0; i < size; i++)
        {
            if (i < size - 1)
            {
                std::cout << ptr[i] << ",";
            }
            else
            {
                std::cout << ptr[i] << std::endl;
            }
        }
    } //end print_shared_ptr
public:
    double T; // temperature
    double beta;
    double a;
    double a_squared;
    double J;
    double J_over_a_squared;

    double h; // step size
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
    //data in 1 flush
    std::shared_ptr<double[]> U_data_all_ptr; //all U data
    std::shared_ptr<double[]> Px_all_ptr; //all Px data
    std::shared_ptr<double[]> Py_all_ptr; //all Py data
    std::shared_ptr<double[]> Qx_all_ptr; //all Qx data
    std::shared_ptr<double[]> Qy_all_ptr; //all Qy data

    //initial value
    std::shared_ptr<double[]> Px_init;
    std::shared_ptr<double[]> Py_init;
    std::shared_ptr<double[]> Qx_init;
    std::shared_ptr<double[]> Qy_init;

    // //for computation
    // arma::dvec Px_arma_vec;
    // arma::dvec Py_arma_vec;
    //
    // arma::dvec Qx_arma_vec;
    // arma::dvec Qy_arma_vec;

    double q;
    double alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7;
    int N0, N1;

    arma::dmat A, B, C, G, R, Gamma, Theta, Lambda;

    std::uniform_int_distribution<int> unif_in_0_N0N1;

    double dipole_lower_bound;
    double dipole_upper_bound;
};
