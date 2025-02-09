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
        return n0*N1+n1;
}
void mc_computation::init_mats()
{
//init A
    for (int n0=0;n0<N0;n0++)
    {
        for (int n1=0;n1<N1;n1++)
        {
            for (int m0=0;m0<N0;m0++)
            {
                for (int m1=0;m1<N1;m1++)
                {
                    double tmp=std::pow(m0-n0,2.0)-(m0-n0)*(m1-n1)+std::pow(m1-n1,2.0);
                }//end m1
            }//end m0
        }
    }//end n0

}