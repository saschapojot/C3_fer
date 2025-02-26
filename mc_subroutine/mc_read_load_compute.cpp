//
// Created by adada on 2/2/2025.
//
# include "mc_read_load_compute.hpp"

///
/// @param n0
/// @param n1
/// @return flatenned index
int mc_computation::double_ind_to_flat_ind(const int& n0, const int& n1)
{
    return n0 * N1 + n1;
}

void mc_computation::init_mats()
{
    //init A
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    double tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0);
                    if (!(n0 == m0 and n1 == m1))
                    {
                        int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                        int flat_col_num = this->double_ind_to_flat_ind(m0, m1);
                        this->A(flat_row_num, flat_col_num) = 1.0 / tmp;
                    } // end if
                } //end m1
            } //end m0
        }
    } //end n0
    //end A init

    // init B
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    if (!(n0 == m0 and n1 == m1))
                    {
                        double up = std::pow(m0 - n0 - 0.5 * m1 + 0.5 * n1, 2.0);
                        double down_tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0);
                        double down = std::pow(down_tmp, 2.0);
                        int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                        int flat_col_num = this->double_ind_to_flat_ind(m0, m1);
                        this->B(flat_row_num, flat_col_num) = up / down;
                    } //end if
                } //end m1
            } //end m0
        } //end n1
    } //end n0
    //end B init


    //init C
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    if (!(n0 == m0 and n1 == m1))
                    {
                        double up = (m0 - n0 - 0.5 * m1 + 0.5 * n1) * (m1 - n1);
                        double down_tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0);
                        double down = std::pow(down_tmp, 2.0);
                        int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                        int flat_col_num = this->double_ind_to_flat_ind(m0, m1);
                        this->C(flat_row_num, flat_col_num) = up / down;
                    } //end if
                } //end m1
            } //end m0
        } //end n1
    } //end n0
    //end C init

    //init G
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    if (!(n0 == m0 and n1 == m1))
                    {
                        double up = std::pow(m1 - n1, 2.0);
                        double down_tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0);
                        double down = std::pow(down_tmp, 2.0);
                        int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                        int flat_col_num = this->double_ind_to_flat_ind(m0, m1);
                        this->G(flat_row_num, flat_col_num) = up / down;
                    } //end if
                } //end m1
            } //end m0
        } //end n1
    } //end n0

    //end G init

    //init R
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    double down = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0) + m1 - n1 +
                        1.0 / 3.0;
                    int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                    int flat_col_num = this->double_ind_to_flat_ind(m0, m1);
                    this->R(flat_row_num, flat_col_num) = 1.0 / down;
                } //end m1
            } //end m0
        } //end n1
    } //end n0

    //end R init

    //init Gamma
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    double up = std::pow(m0 - n0 - 0.5 * m1 + 0.5 * n1, 2.0);
                    double down_tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0) + m1 - n1
                        + 1.0 / 3.0;
                    double down = std::pow(down_tmp, 2.0);
                    int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                    int flat_col_num = this->double_ind_to_flat_ind(m0, m1);
                    this->Gamma(flat_row_num, flat_col_num) = up / down;
                } //end m1
            } //end m0
        } //end n1
    } //end n0
    //end Gamma init

    //init Theta
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    double up = (m0 - n0 - 0.5 * m1 + 0.5 * n1) * (std::sqrt(3.0) / 2.0 * m1 - std::sqrt(3.0) / 2.0 * n1
                        + std::sqrt(3.0) / 2.0);

                    double down_tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0) + m1 - n1
                        + 1.0 / 3.0;

                    double down = std::pow(down_tmp, 2.0);
                    int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                    int flat_col_num = this->double_ind_to_flat_ind(m0, m1);

                    this->Theta(flat_row_num, flat_col_num) = up / down;
                } //end m1
            } //end m0
        } //end n1
    } //end n0
    //end Theta init

    //init Lambda
    for (int n0 = 0; n0 < N0; n0++)
    {
        for (int n1 = 0; n1 < N1; n1++)
        {
            for (int m0 = 0; m0 < N0; m0++)
            {
                for (int m1 = 0; m1 < N1; m1++)
                {
                    double up = std::pow(std::sqrt(3.0) / 2.0 * m1 - std::sqrt(3.0) / 2.0 * n1 + std::sqrt(3.0) / 3.0,
                                         2.0);

                    double down_tmp = std::pow(m0 - n0, 2.0) - (m0 - n0) * (m1 - n1) + std::pow(m1 - n1, 2.0) + m1 - n1
                        + 1.0 / 3.0;
                    double down = std::pow(down_tmp, 2.0);

                    int flat_row_num = this->double_ind_to_flat_ind(n0, n1);
                    int flat_col_num = this->double_ind_to_flat_ind(m0, m1);

                    this->Lambda(flat_row_num, flat_col_num) = up / down;
                } //end m1
            } //end m0
        } //end n1
    } //end n0
    //end Lambda init
}

void mc_computation::save_array_to_pickle(const std::shared_ptr<double[]>& ptr, int size, const std::string& filename)
{
    using namespace boost::python;
    namespace np = boost::python::numpy;

    // Initialize Python interpreter if not already initialized
    if (!Py_IsInitialized())
    {
        Py_Initialize();
        if (!Py_IsInitialized())
        {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }
        np::initialize(); // Initialize NumPy
    }

    try
    {
        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Convert C++ array to NumPy array using shared_ptr
        np::ndarray numpy_array = np::from_data(
            ptr.get(), // Use shared_ptr's raw pointer
            np::dtype::get_builtin<double>(), // NumPy data type (double)
            boost::python::make_tuple(size), // Shape of the array (1D array)
            boost::python::make_tuple(sizeof(double)), // Strides
            object() // Optional base object
        );

        // Serialize the NumPy array using pickle.dumps
        object serialized_array = pickle_dumps(numpy_array);

        // Extract the serialized data as a string
        std::string serialized_str = extract<std::string>(serialized_array);

        // Write the serialized data to a file
        std::ofstream file(filename, std::ios::binary);
        if (!file)
        {
            throw std::runtime_error("Failed to open file for writing");
        }
        file.write(serialized_str.data(), serialized_str.size());
        file.close();

        // Debug output (optional)
        // std::cout << "Array serialized and written to file successfully." << std::endl;
    }
    catch (const error_already_set&)
    {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}

void mc_computation::load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr,
                                      std::size_t size)
{
    // Initialize Python and NumPy
    Py_Initialize();
    np::initialize();


    try
    {
        // Use Python's 'io' module to open the file directly in binary mode
        py::object io_module = py::import("io");
        py::object file = io_module.attr("open")(filename, "rb"); // Open file in binary mode

        // Import the 'pickle' module
        py::object pickle_module = py::import("pickle");

        // Use pickle.load to deserialize from the Python file object
        py::object loaded_data = pickle_module.attr("load")(file);

        // Close the file
        file.attr("close")();

        // Check if the loaded object is a NumPy array
        if (py::extract<np::ndarray>(loaded_data).check())
        {
            np::ndarray np_array = py::extract<np::ndarray>(loaded_data);

            // Convert the NumPy array to a Python list using tolist()
            py::object py_list = np_array.attr("tolist")();

            // Ensure the list size matches the expected size
            ssize_t list_size = py::len(py_list);
            if (static_cast<std::size_t>(list_size) > size)
            {
                throw std::runtime_error("The provided shared_ptr array size is smaller than the list size.");
            }

            // Copy the data from the Python list to the shared_ptr array
            for (ssize_t i = 0; i < list_size; ++i)
            {
                data_ptr[i] = py::extract<double>(py_list[i]);
            }
        }
        else
        {
            throw std::runtime_error("Loaded data is not a NumPy array.");
        }
    }
    catch (py::error_already_set&)
    {
        PyErr_Print();
        throw std::runtime_error("Python error occurred.");
    }
}

void mc_computation::init_Px_Py_Qx_Qy()
{
    std::string name;

    std::string Px_inFileName, Py_inFileName, Qx_inFileName, Qy_inFileName;
    if (this->flushLastFile == -1)
    {
        name = "init";

        Px_inFileName = out_Px_path + "/Px_" + name + ".pkl";

        Py_inFileName = out_Py_path + "/Py_" + name + ".pkl";

        Qx_inFileName = out_Qx_path + "/Qx_" + name + ".pkl";

        Qy_inFileName = out_Qy_path + "/Qy_" + name + ".pkl";

        this->load_pickle_data(Px_inFileName, Px_init, N0 * N1);
        this->load_pickle_data(Py_inFileName, Py_init, N0 * N1);
        this->load_pickle_data(Qx_inFileName, Qx_init, N0 * N1);
        this->load_pickle_data(Qy_inFileName, Qy_init, N0 * N1);
    } //end flushLastFile==-1
}


///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
double mc_computation::generate_uni_open_interval(const double& x, const double& leftEnd, const double& rightEnd,
                                                  const double& eps)
{
    double xMinusEps = x - eps;
    double xPlusEps = x + eps;

    double unif_left_end = xMinusEps < leftEnd ? leftEnd : xMinusEps;
    double unif_right_end = xPlusEps > rightEnd ? rightEnd : xPlusEps;

    //    std::random_device rd;
    //    std::ranlux24_base e2(rd());

    double unif_left_end_double_on_the_right = std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());


    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right, unif_right_end);
    //(unif_left_end_double_on_the_right, unif_right_end)

    double xNext = distUnif(e2);
    return xNext;
}


///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double& x, const double& y, const double& a, const double& b, const double& epsilon)
{
    if (a < y and y < a + epsilon)
    {
        return 1.0 / (y - a + epsilon);
    }
    else if (a + epsilon <= y and y < b + epsilon)
    {
        return 1.0 / (2.0 * epsilon);
    }
    else if (b - epsilon <= y and y < b)
    {
        return 1.0 / (b - y + epsilon);
    }
    else
    {
        std::cerr << "value out of range." << std::endl;
        std::exit(10);
    }
}

///
/// @param n0 index
/// @param n1 index
/// @param Px_arma_vec Px
/// @param Py_arma_vec Py
/// @return self energy H1
double mc_computation::H1(const int& flattened_ind, const arma::dvec& Px_arma_vec, const arma::dvec& Py_arma_vec)
{
    // int flat_ind = double_ind_to_flat_ind(n0, n1);

    double px_n0n1 = Px_arma_vec(flattened_ind);

    double py_n0n1 = Py_arma_vec(flattened_ind);


    double squared_px_n0n1 = std::pow(px_n0n1, 2.0);
    double squared_py_n0n1 = std::pow(py_n0n1, 2.0);

    double part1 = alpha1 * px_n0n1 * (squared_px_n0n1 - 3.0 * squared_py_n0n1);

    double part2 = alpha2 * py_n0n1 * (3.0 * squared_px_n0n1 - squared_py_n0n1);

    double part3 = alpha3 * (
        squared_px_n0n1 * std::pow(squared_px_n0n1 - 3.0 * squared_py_n0n1, 2.0)
        - squared_py_n0n1 * std::pow(3.0 * squared_px_n0n1 - squared_py_n0n1, 2.0)
    );

    double part4 = alpha4 * px_n0n1 * py_n0n1 * (squared_px_n0n1 - 3.0 * squared_py_n0n1)
        * (3.0 * squared_px_n0n1 - squared_py_n0n1);

    double part5 = alpha5 * (squared_px_n0n1 + squared_py_n0n1);

    double part6 = alpha6 * std::pow(squared_px_n0n1 + squared_py_n0n1, 2.0);

    double part7 = alpha7 * std::pow(squared_px_n0n1 + squared_py_n0n1, 3.0);

    return part1 + part2 + part3 + part4 + part5 + part6 + part7;
}

///
/// @param n0 index
/// @param n1 index
/// @param Qx_arma_vec Qx
/// @param Qy_arma_vec Qy
/// @return self energy H2
double mc_computation::H2(const int& flattened_ind, const arma::dvec& Qx_arma_vec, const arma::dvec& Qy_arma_vec)
{
    // int flat_ind = double_ind_to_flat_ind(n0, n1);

    double qx_n0n1 = Qx_arma_vec(flattened_ind);

    double qy_n0n1 = Qy_arma_vec(flattened_ind);

    double squared_qx_n0n1 = std::pow(qx_n0n1, 2.0);

    double squared_qy_n0n1 = std::pow(qy_n0n1, 2.0);

    double part1 = alpha1 * qx_n0n1 * (squared_qx_n0n1 - 3.0 * squared_qy_n0n1);

    double part2 = alpha2 * qy_n0n1 * (3.0 * squared_qx_n0n1 - squared_qy_n0n1);

    double part3 = alpha3 * (
        squared_qx_n0n1 * std::pow(squared_qx_n0n1 - 3.0 * squared_qy_n0n1, 2.0)
        - squared_qy_n0n1 * std::pow(3.0 * squared_qx_n0n1 - squared_qy_n0n1, 2.0)
    );

    double part4 = alpha4 * qx_n0n1 * qy_n0n1 * (squared_qx_n0n1 - 3.0 * squared_qy_n0n1)
        * (3.0 * squared_qx_n0n1 - squared_qy_n0n1);

    double part5 = alpha5 * (squared_qx_n0n1 + squared_qy_n0n1);

    double part6 = alpha6 * std::pow(squared_qx_n0n1 + squared_qy_n0n1, 2.0);

    double part7 = alpha7 * std::pow(squared_qx_n0n1 + squared_qy_n0n1, 3.0);

    return part1 + part2 + part3 + part4 + part5 + part6 + part7;
}

void mc_computation::init_and_run()
{
    this->init_mats();
    this->init_Px_Py_Qx_Qy();
}


void mc_computation::execute_mc(const std::shared_ptr<double[]>& Px_vec,
                                const std::shared_ptr<double[]>& Py_vec,
                                const std::shared_ptr<double[]>& Qx_vec,
                                const std::shared_ptr<double[]>& Qy_vec,
                                const int& flushNum)
{
    arma::dvec Px_arma_vec_curr(Px_vec.get(), N0 * N1);
    arma::dvec Px_arma_vec_next(N0 * N1, arma::fill::zeros);


    arma::dvec Py_arma_vec_curr(Py_vec.get(), N0 * N1);
    arma::dvec Py_arma_vec_next(N0 * N1, arma::fill::zeros);

    arma::dvec Qx_arma_vec_curr(Qx_vec.get(), N0 * N1);
    arma::dvec Qx_arma_vec_next(N0 * N1, arma::fill::zeros);

    arma::dvec Qy_arma_vec_curr(Qy_vec.get(), N0 * N1);
    arma::dvec Qy_arma_vec_next(N0 * N1, arma::fill::zeros);

    double UCurr=0;
    int flushThisFileStart=this->flushLastFile+1;
    int sweepStart =flushThisFileStart*sweepToWrite*sweep_multiple;
    for (int fls = 0; fls < flushNum; fls++)
    {
        for (int swp = 0; swp < sweepToWrite*sweep_multiple; swp++)
        {
        }//end sweep for
    }//end flush for loop

}

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
void mc_computation::HPx_update(const int& flattened_ind, const arma::dvec& Px_arma_vec_curr,
                                const arma::dvec& Px_arma_vec_next,
                                const arma::dvec& Py_arma_vec_curr,
                                const arma::dvec& Qx_arma_vec_curr,
                                const arma::dvec& Qy_arma_vec_curr,
                                double& UCurr, double& UNext)
{
    // int flattened_ind=this->double_ind_to_flat_ind(n0,n1);

    double H1_self_curr = this->H1(flattened_ind, Px_arma_vec_curr, Py_arma_vec_curr);
    double H1_self_next = this->H1(flattened_ind, Px_arma_vec_next, Py_arma_vec_curr);


    double right_factor1 = J_over_a_squared * arma::dot(A.row(flattened_ind), Px_arma_vec_curr);

    double right_factor2 = -2.0 * J_over_a_squared * arma::dot(B.row(flattened_ind), Px_arma_vec_curr);

    double right_factor3 = -std::sqrt(3.0) * J_over_a_squared * arma::dot(C.row(flattened_ind), Py_arma_vec_curr);

    double right_factor4 = J_over_a_squared * arma::dot(R.row(flattened_ind), Qx_arma_vec_curr);

    double right_factor5 = -2.0 * J_over_a_squared * arma::dot(Gamma.row(flattened_ind), Qx_arma_vec_curr);

    double right_factor6 = -2.0 * J_over_a_squared * arma::dot(Theta.row(flattened_ind), Qy_arma_vec_curr);

    double right_factor = right_factor1 + right_factor2 + right_factor3
        + right_factor4 + right_factor5 + right_factor6;

    double E_int_curr = Px_arma_vec_curr(flattened_ind) * right_factor;

    double E_int_next = Px_arma_vec_next(flattened_ind) * right_factor;

    UCurr = H1_self_curr + E_int_curr;

    UNext = H1_self_next + E_int_next;
}

void mc_computation::HPy_update(const int& flattened_ind,
                                const arma::dvec& Py_arma_vec_curr, const arma::dvec& Py_arma_vec_next,
                                const arma::dvec& Px_arma_vec_curr, const arma::dvec& Qx_arma_vec_curr,
                                const arma::dvec& Qy_arma_vec_curr,
                                double& UCurr, double& UNext)
{
    double H1_self_curr = this->H1(flattened_ind, Px_arma_vec_curr, Py_arma_vec_curr);

    double H1_self_next = this->H1(flattened_ind, Px_arma_vec_curr, Py_arma_vec_next);

    double right_factor1 = J_over_a_squared * arma::dot(A.row(flattened_ind), Py_arma_vec_curr);

    double right_factor2 = -std::sqrt(3.0) * J_over_a_squared * arma::dot(C.row(flattened_ind), Px_arma_vec_curr);
    double right_factor3 = -3.0 / 2.0 * J_over_a_squared * arma::dot(G.row(flattened_ind), Py_arma_vec_curr);

    double right_factor4 = J_over_a_squared * arma::dot(R.row(flattened_ind), Qy_arma_vec_curr);

    double right_factor5 = -2.0 * J_over_a_squared * arma::dot(Theta.row(flattened_ind), Qx_arma_vec_curr);

    double right_factor6 = -2.0 * J_over_a_squared * arma::dot(Lambda.row(flattened_ind), Qy_arma_vec_curr);

    double right_factor = right_factor1 + right_factor2 + right_factor3
        + right_factor4 + right_factor5 + right_factor6;

    double E_int_curr = Py_arma_vec_curr(flattened_ind) * right_factor;

    double E_int_next = Py_arma_vec_next(flattened_ind) * right_factor;

    UCurr = H1_self_curr + E_int_curr;

    UNext = H1_self_next + E_int_next;
}


double mc_computation::acceptanceRatio_uni(const arma::dvec& arma_vec_curr,
                                           const arma::dvec& arma_vec_next, const int& flattened_ind,
                                           const double& UCurr, const double& UNext)
{
    double numerator = -this->beta * UNext;
    double denominator = -this->beta * UCurr;
    double R = std::exp(numerator - denominator);

    double S_curr_next = S_uni(arma_vec_curr(flattened_ind), arma_vec_next(flattened_ind),
                               dipole_lower_bound, dipole_upper_bound, h);

    double S_next_curr = S_uni(arma_vec_next(flattened_ind), arma_vec_curr(flattened_ind),
                               dipole_lower_bound, dipole_upper_bound, h);

    double ratio = S_curr_next / S_next_curr;

    if (std::fetestexcept(FE_DIVBYZERO))
    {
        std::cout << "Division by zero exception caught." << std::endl;
        std::exit(15);
    }
    if (std::isnan(ratio))
    {
        std::cout << "The result is NaN." << std::endl;
        std::exit(15);
    }
    R *= ratio;

    return std::min(1.0, R);
}

void mc_computation::execute_mc_one_sweep(arma::dvec& Px_arma_vec_curr,
                                          arma::dvec& Py_arma_vec_curr,
                                          arma::dvec& Qx_arma_vec_curr,
                                          arma::dvec& Qy_arma_vec_curr,
                                          double& UCurr,
                                          arma::dvec& Px_arma_vec_next,
                                          arma::dvec& Py_arma_vec_next,
                                          arma::dvec& Qx_arma_vec_next,
                                          arma::dvec& Qy_arma_vec_next)
{
    double UNext = 0;
    //update Px
    for (int i = 0; i < N0 * N1; i++)
    {
        int flattened_ind = unif_in_0_N0N1(e2);
        this->proposal_uni(Px_arma_vec_curr, Px_arma_vec_next, flattened_ind);
        this->HPx_update(flattened_ind, Px_arma_vec_curr, Px_arma_vec_next,
                         Py_arma_vec_curr, Qx_arma_vec_curr, Qy_arma_vec_curr, UCurr, UNext);

        double r = this->acceptanceRatio_uni(Px_arma_vec_curr, Px_arma_vec_next,
                                             flattened_ind, UCurr, UNext);
        double u = distUnif01(e2);
        if (u <= r)
        {
            Px_arma_vec_curr = Px_arma_vec_next;
            UCurr = UNext;
        } //end of accept-reject
    } //end updating Px

    //update Py
    for (int i = 0; i < N0 * N1; i++)
    {
        int flattened_ind = unif_in_0_N0N1(e2);
        this->proposal_uni(Py_arma_vec_curr, Py_arma_vec_next, flattened_ind);
        this->HPy_update(flattened_ind, Py_arma_vec_curr, Py_arma_vec_next,
                         Px_arma_vec_curr, Qx_arma_vec_curr, Qy_arma_vec_curr, UCurr, UNext);

        double r = this->acceptanceRatio_uni(Py_arma_vec_curr, Py_arma_vec_next, flattened_ind, UCurr, UNext);

        double u = distUnif01(e2);
        if (u <= r)
        {
            Py_arma_vec_curr = Py_arma_vec_next;
            UCurr = UNext;
        } //end of accept-reject
    } //end updating Py

    //update Qx
    for (int i = 0; i < N0 * N1; i++)
    {
        int flattened_ind = unif_in_0_N0N1(e2);
        this->proposal_uni(Qx_arma_vec_curr, Qx_arma_vec_next, flattened_ind);

        this->HQx_update(flattened_ind, Qx_arma_vec_curr, Qx_arma_vec_next,
                         Px_arma_vec_curr, Py_arma_vec_curr,
                         Qy_arma_vec_curr,
                         UCurr, UNext);

        double r = this->acceptanceRatio_uni(Qx_arma_vec_curr, Qx_arma_vec_next, flattened_ind, UCurr, UNext);
        double u = distUnif01(e2);
        if (u <= r)
        {
            Qx_arma_vec_curr = Qx_arma_vec_next;
            UCurr = UNext;
        } //end of accept-reject
    } //end updating Qx

    //update Qy
    for (int i = 0; i < N0 * N1; i++)
    {
        int flattened_ind = unif_in_0_N0N1(e2);
        this->proposal_uni(Qy_arma_vec_curr, Qy_arma_vec_next, flattened_ind);
        this->HQy_update(flattened_ind, Qy_arma_vec_curr, Qy_arma_vec_next,
                         Px_arma_vec_curr,
                         Py_arma_vec_curr,
                         Qx_arma_vec_curr,
                         UCurr, UNext);
    } //end updating Qy
}

void mc_computation::HQy_update(const int& flattened_ind,
                                const arma::dvec& Qy_arma_vec_curr, const arma::dvec& Qy_arma_vec_next,
                                const arma::dvec& Px_arma_vec_curr,
                                const arma::dvec& Py_arma_vec_curr,
                                const arma::dvec& Qx_arma_vec_curr,
                                double& UCurr, double& UNext)
{
    double H2_self_curr = this->H2(flattened_ind, Qx_arma_vec_curr, Qy_arma_vec_curr);

    double H2_self_next = this->H2(flattened_ind, Qx_arma_vec_curr, Qy_arma_vec_next);

    double factor1 = J_over_a_squared * arma::dot(Py_arma_vec_curr, R.col(flattened_ind));

    double factor2 = -2.0 * J_over_a_squared * arma::dot(Px_arma_vec_curr, Theta.col(flattened_ind));

    double factor3 = -2.0 * J_over_a_squared * arma::dot(Py_arma_vec_curr, Lambda.col(flattened_ind));

    double factor4 = J_over_a_squared * arma::dot(A.row(flattened_ind), Qy_arma_vec_curr);

    double factor5 = -std::sqrt(3.0) * J_over_a_squared * arma::dot(C.row(flattened_ind), Qx_arma_vec_curr);

    double factor6 = -3.0 / 2.0 * J_over_a_squared * arma::dot(G.row(flattened_ind), Qy_arma_vec_curr);

    double factor = factor1 + factor2 + factor3
        + factor4 + factor5 + factor6;

    double E_int_curr = Qy_arma_vec_curr(flattened_ind) * factor;

    double E_int_next = Qy_arma_vec_next(flattened_ind) * factor;

    UCurr = H2_self_curr + E_int_curr;

    UNext = H2_self_next + E_int_next;
}

void mc_computation::HQx_update(const int& flattened_ind,
                                const arma::dvec& Qx_arma_vec_curr, const arma::dvec& Qx_arma_vec_next,
                                const arma::dvec& Px_arma_vec_curr,
                                const arma::dvec& Py_arma_vec_curr,
                                const arma::dvec& Qy_arma_vec_curr,
                                double& UCurr, double& UNext)
{
    double H2_self_curr = this->H2(flattened_ind, Qx_arma_vec_curr, Qy_arma_vec_curr);
    double H2_self_next = this->H2(flattened_ind, Qx_arma_vec_next, Qy_arma_vec_curr);

    double factor1 = J_over_a_squared * arma::dot(Px_arma_vec_curr, R.col(flattened_ind));

    double factor2 = -2.0 * J_over_a_squared * arma::dot(Px_arma_vec_curr, Gamma.col(flattened_ind));

    double factor3 = -2.0 * J_over_a_squared * arma::dot(Py_arma_vec_curr, Theta.col(flattened_ind));

    double factor4 = J_over_a_squared * arma::dot(A.row(flattened_ind), Qx_arma_vec_curr);

    double factor5 = -2.0 * J_over_a_squared * arma::dot(B.row(flattened_ind), Qx_arma_vec_curr);

    double factor6 = -std::sqrt(3.0) * J_over_a_squared * arma::dot(C.row(flattened_ind), Qy_arma_vec_curr);


    double factor = factor1 + factor2 + factor3
        + factor4 + factor5 + factor6;

    double E_int_curr = Qx_arma_vec_curr(flattened_ind) * factor;

    double E_int_next = Qx_arma_vec_next(flattened_ind) * factor;

    UCurr = H2_self_curr + E_int_curr;

    UNext = H2_self_next + E_int_next;
}

void mc_computation::proposal_uni(const arma::dvec& arma_vec_curr, arma::dvec& arma_vec_next,
                                  const int& flattened_ind)
{
    double dp_val_new = this->generate_uni_open_interval(arma_vec_curr(flattened_ind), dipole_lower_bound,
                                                         dipole_upper_bound, h);
    arma_vec_next = arma_vec_curr;
    arma_vec_next(flattened_ind) = dp_val_new;
}
