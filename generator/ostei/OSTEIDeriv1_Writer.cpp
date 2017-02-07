#include <algorithm>

#include "generator/Types.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"
#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/ostei/OSTEI_Writer.hpp"


/////////////////////////////
// Basic OSTEI Writer
/////////////////////////////
bool OSTEIDeriv1_Writer::IsSpecialPermutation_(QAM am) const
{
    // is this a special permutation? Handle it that way if it is
    // NOTE: We never switch bra/ket, since that would affect the
    //       vectorization (since we vectorize on the ket)
    if( ( am[0] == 0 && am[1] > 0 ) && ( am[2] == 0 || am[3] == 0 ) )
        return true;
    if( ( am[2] == 0 && am[3] > 0 ) && ( am[0] == 0 || am[1] == 0 ) )
        return true;

    return false;
}


void OSTEIDeriv1_Writer::DeclareContwork(void) const
{
    if(info_.ContMemoryReq() == 0)
        return;

    os_ << indent1 << "// partition workspace into shells\n";
    size_t ptidx = 0;

    for(const auto & it : info_.GetContQ())
    {
        if(!info_.IsFinalAM(it))
        {
            os_ << indent1 << "double * const " << ArrVarName(it) << " = contwork + (SIMINT_NSHELL_SIMD * " << ptidx << ");\n";
            ptidx += NCART(it);
        }
    }

    os_ << "\n";
}

void OSTEIDeriv1_Writer::WriteShellOffsets(void) const
{
    os_ << indent5 << "// calculate the shell offsets\n";
    os_ << indent5 << "// these are the offset from the shell pointed to by cd\n";
    os_ << indent5 << "// for each element\n";
    os_ << indent5 << "int shelloffsets[SIMINT_SIMD_LEN] = {0};\n";
    os_ << indent5 << "int lastoffset = 0;\n";
    os_ << indent5 << "const int nlane = ( ((j + SIMINT_SIMD_LEN) < jend) ? SIMINT_SIMD_LEN : (jend - j));\n";
    os_ << "\n";
    os_ << indent5 << "if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)\n";
    os_ << indent5 << "{\n";

    os_ << indent6 << "// Handle if the first element of the vector is a new shell\n";
    os_ << indent6 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os_ << indent6 << "{\n";
    os_ << indent7 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto it : info_.GetContQ())
        os_ << indent7 << PrimPtrName(it) << " += " << NCART(it) << ";\n";

    os_ << indent6 << "}\n";
    os_ << indent6 << "iprimcd++;\n";

    os_ << indent6 << "for(n = 1; n < SIMINT_SIMD_LEN; ++n)\n";
    os_ << indent6 << "{\n";
    os_ << indent7 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os_ << indent7 << "{\n";
    os_ << indent8 << "shelloffsets[n] = shelloffsets[n-1] + 1;\n";
    os_ << indent8 << "lastoffset++;\n";
    os_ << indent8 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    os_ << indent7 << "}\n";
    os_ << indent7 << "else\n";
    os_ << indent8 << "shelloffsets[n] = shelloffsets[n-1];\n";
    os_ << indent7 << "iprimcd++;\n";
    os_ << indent6 << "}\n";
    os_ << indent5 << "}\n";
    os_ << indent5 << "else\n";
    os_ << indent6 << "iprimcd += SIMINT_SIMD_LEN;\n\n";
}


void OSTEIDeriv1_Writer::WriteAccumulation(void) const
{
    os_ << "\n\n";
    os_ << indent5 << "////////////////////////////////////\n";
    os_ << indent5 << "// Accumulate contracted integrals\n";
    os_ << indent5 << "////////////////////////////////////\n";

    os_ << indent5 << "if(lastoffset == 0)\n";
    os_ << indent5 << "{\n";

    for(const auto it : info_.GetContQ())
    {
        int ncart = NCART(it);
        if(it.tag.size())
            os_ << indent6 << "contract_all_fac(" << ncart << ", cfac_" << it.tag << ", "
                           << PrimVarName(it.notag()) << ", " << PrimPtrName(it) << ");\n";
        else
            os_ << indent6 << "contract_all    (" << ncart << ", " << PrimVarName(it)
                           << ", " << PrimPtrName(it) << ");\n";
    }
    os_ << indent5 << "}\n";
    os_ << indent5 << "else\n";
    os_ << indent5 << "{\n";

    for(const auto it : info_.GetContQ())
    {
        int ncart = NCART(it);
        if(it.tag.size())
            os_ << indent6 << "contract_fac(" << ncart << ", cfac_" << it.tag << ", shelloffsets, "
                           << PrimVarName(it.notag()) << ", " << PrimPtrName(it) << ");\n";
        else
            os_ << indent6 << "contract    (" << ncart << ", shelloffsets, " << PrimVarName(it)
                           << ", " << PrimPtrName(it) << ");\n";
    }

    for(const auto it : info_.GetContQ())
        os_ << indent6 << PrimPtrName(it) << " += lastoffset*" << NCART(it) << ";\n";

    os_ << indent5 << "}\n";
}

std::string OSTEIDeriv1_Writer::FunctionName_(QAM am) const
{
    return StringBuilder("ostei_deriv1_",
                         amchar[am[0]], "_",
                         amchar[am[1]], "_" ,
                         amchar[am[2]], "_",
                         amchar[am[3]]);
}

std::string OSTEIDeriv1_Writer::FunctionPrototype_(QAM am) const
{
    std::string fname = FunctionName_(am);
    std::string indent(fname.length()+1+4, ' '); // +4 for return type

    std::stringstream ss;
    ss << "int " << fname << "(";
    ss << "struct simint_multi_shellpair const P,\n";
    ss << indent << "struct simint_multi_shellpair const Q,\n";
    ss << indent << "double screen_tol,\n";
    ss << indent << "double * const restrict contwork,\n";
    ss << indent << "double * const restrict " << ArrVarName(am) << ")";
    return ss.str();
}

void OSTEIDeriv1_Writer::Write_Permute_(QAM am, bool swap12, bool swap34) const
{
    QAM permuted = am;
    if(swap12)
        std::swap(permuted[0], permuted[1]);
    if(swap34)
        std::swap(permuted[2], permuted[3]);

    // is this permutation unique?
    if(swap34 && !swap12 && permuted[2] == permuted[3])
        return;
    if(swap12 && !swap34 && permuted[0] == permuted[1])
        return;
    if(swap12 && swap34 && (permuted[0] == permuted[1] || permuted[2] == permuted[3]))
        return;

    // output of the function starts here
    os_ << FunctionPrototype_(permuted) << "\n";
    os_ << "{\n";

    const char * P_var = "P";
    const char * Q_var = "Q";

    if(swap12)
    {
        P_var = "P_tmp";
        os_ << indent1 << "double P_AB[3*P.nshell12];\n";
    	os_ << indent1 << "struct simint_multi_shellpair P_tmp = P;\n";
        os_ << indent1 << "P_tmp.PA_x = P.PB_x;  P_tmp.PA_y = P.PB_y;  P_tmp.PA_z = P.PB_z;\n";
        os_ << indent1 << "P_tmp.PB_x = P.PA_x;  P_tmp.PB_y = P.PA_y;  P_tmp.PB_z = P.PA_z;\n";
        os_ << indent1 << "P_tmp.AB_x = P_AB;\n";
        os_ << indent1 << "P_tmp.AB_y = P_AB + P.nshell12;\n";
        os_ << indent1 << "P_tmp.AB_z = P_AB + 2*P.nshell12;\n";
        os_ << "\n";
        os_ << indent1 << "for(int i = 0; i < P.nshell12; i++)\n";
        os_ << indent1 << "{\n";
        os_ << indent2 << "P_tmp.AB_x[i] = -P.AB_x[i];\n";
        os_ << indent2 << "P_tmp.AB_y[i] = -P.AB_y[i];\n";
        os_ << indent2 << "P_tmp.AB_z[i] = -P.AB_z[i];\n";
        os_ << indent1 << "}\n\n";
    }

    if(swap34)
    {
	    Q_var = "Q_tmp";
        os_ << indent1 << "double Q_AB[3*Q.nshell12];\n";
    	os_ << indent1 << "struct simint_multi_shellpair Q_tmp = Q;\n";
		os_ << indent1 << "Q_tmp.PA_x = Q.PB_x;  Q_tmp.PA_y = Q.PB_y;  Q_tmp.PA_z = Q.PB_z;\n";
        os_ << indent1 << "Q_tmp.PB_x = Q.PA_x;  Q_tmp.PB_y = Q.PA_y;  Q_tmp.PB_z = Q.PA_z;\n";
        os_ << indent1 << "Q_tmp.AB_x = Q_AB;\n";
        os_ << indent1 << "Q_tmp.AB_y = Q_AB + Q.nshell12;\n";
        os_ << indent1 << "Q_tmp.AB_z = Q_AB + 2*Q.nshell12;\n";
        os_ << "\n";
        os_ << indent1 << "for(int i = 0; i < Q.nshell12; i++)\n";
        os_ << indent1 << "{\n";
        os_ << indent2 << "Q_tmp.AB_x[i] = -Q.AB_x[i];\n";
        os_ << indent2 << "Q_tmp.AB_y[i] = -Q.AB_y[i];\n";
        os_ << indent2 << "Q_tmp.AB_z[i] = -Q.AB_z[i];\n";
        os_ << indent1 << "}\n\n";
    }

    std::string fname = FunctionName_(am); 
    os_ << indent1 << "int ret = " << fname
        << "(" << P_var << ", " << Q_var << ", screen_tol, "
        << "contwork, " << ArrVarName(permuted) << ");\n"; 


    if(!IsSpecialPermutation_(permuted))
    {
        size_t ncart_abcd = NCART(am);
        //size_t ncart_a  = NCART(am[0]);
        size_t ncart_b  = NCART(am[1]);
        size_t ncart_c  = NCART(am[2]);
        size_t ncart_d  = NCART(am[3]);
        size_t ncart_a2 = NCART(permuted[0]);
        size_t ncart_b2 = NCART(permuted[1]);
        size_t ncart_c2 = NCART(permuted[2]);
        size_t ncart_d2 = NCART(permuted[3]);

        size_t ncart_bcd = ncart_b * ncart_c * ncart_d;
        size_t ncart_cd = ncart_c * ncart_d; 

        char va = 'a';
        char vb = 'b';
        char vc = 'c';
        char vd = 'd';

        if(swap12)
            std::swap(va, vb);
        if(swap34)
            std::swap(vc, vd);

        std::string idx = StringBuilder("q*12", ncart_abcd,
                                        "+", va, "*", ncart_bcd,
                                        "+", vb, "*", ncart_cd,
                                        "+", vc, "*", ncart_d, "+", vd);

        os_ << indent1 << "double buffer[" << ncart_abcd << "] SIMINT_ALIGN_ARRAY_DBL;\n\n"; 


        os_ << indent1 << "for(int q = 0; q < ret; q++)\n";
        os_ << indent1 << "{\n";
        os_ << indent2 << "int idx = 0;\n";
        os_ << indent2 << "for(int a = 0; a < " << ncart_a2 << "; ++a)\n";
        os_ << indent2 << "for(int b = 0; b < " << ncart_b2 << "; ++b)\n";
        os_ << indent2 << "for(int c = 0; c < " << ncart_c2 << "; ++c)\n";
        os_ << indent2 << "for(int d = 0; d < " << ncart_d2 << "; ++d)\n";
        os_ << indent2 << "{\n";
        os_ << indent3 << "memcpy(buffer, &" << ArrVarName(permuted) << "[" << idx << "], 12*sizeof(double));\n";
        os_ << indent3 << "idx += 12;\n";
        os_ << indent2 << "}\n";
        os_ << "\n";
        os_ << indent2 << "memcpy(" << ArrVarName(permuted) << "+q*" << ncart_abcd << "*12"
                       << ", buffer, " << ncart_abcd << "*12*sizeof(double));\n";
        os_ << indent1 << "}\n";
    }

    os_ << "\n";
    os_ << indent1 << "return ret;\n";
    os_ << "}\n";
    os_ << "\n";

    osh_ << FunctionPrototype_(permuted) << ";\n\n";
}

void OSTEIDeriv1_Writer::WriteFormDeriv(void) const
{
    QAM am = info_.FinalAM();
    const int ncart1 = NCART(am[0]);
    const int ncart2 = NCART(am[1]);
    const int ncart3 = NCART(am[2]);
    const int ncart4 = NCART(am[3]);

    const int ncart2p = NCART(am[1]+1);
    const int ncart3p = NCART(am[2]+1);
    const int ncart2m = NCART(am[1]-1);
    const int ncart3m = NCART(am[2]-1);

    std::string outvar = HRRVarName(am);
    std::string invar_1p = HRRVarName(QAM({am[0]+1, am[1]  , am[2]  , am[3]  }, "2a"));
    std::string invar_2p = HRRVarName(QAM({am[0],   am[1]+1, am[2]  , am[3]  }, "2b"));
    std::string invar_3p = HRRVarName(QAM({am[0],   am[1]  , am[2]+1, am[3]  }, "2c"));
    std::string invar_1m = HRRVarName({am[0]-1, am[1]  , am[2]  , am[3]  });
    std::string invar_2m = HRRVarName({am[0],   am[1]-1, am[2]  , am[3]  });
    std::string invar_3m = HRRVarName({am[0],   am[1]  , am[2]-1, am[3]  });

    os_ << indent4 << "struct RecurInfo const * const aminfo_1 = &recurinfo_array[am_recur_map[" << am[0] << "]];\n";
    os_ << indent4 << "struct RecurInfo const * const aminfo_2 = &recurinfo_array[am_recur_map[" << am[1] << "]];\n";
    os_ << indent4 << "struct RecurInfo const * const aminfo_3 = &recurinfo_array[am_recur_map[" << am[2] << "]];\n";
    os_ << indent4 << "struct RecurInfo const * const aminfo_4 = &recurinfo_array[am_recur_map[" << am[3] << "]];\n";
    os_ << "\n";


    os_ << indent4 << "double * restrict " << outvar << " = " << ArrVarName(am)
                   << " + real_abcd * " << ncart1*ncart2*ncart3*ncart4 << " * 12;\n\n";
    os_ << indent4 << "int idx_x, idx_y, idx_z;\n";
    os_ << indent4 << "int full_idx_x, full_idx_y, full_idx_z;\n";
    os_ << indent4 << "int startidx = 0;\n";
    os_ << indent4 << "for(int n1 = 0; n1 < " << ncart1 << "; n1++)\n";
    os_ << indent4 << "for(int n2 = 0; n2 < " << ncart2 << "; n2++)\n";
    os_ << indent4 << "for(int n3 = 0; n3 < " << ncart3 << "; n3++)\n";
    os_ << indent4 << "for(int n4 = 0; n4 < " << ncart4 << "; n4++)\n";
    os_ << indent4 << "{\n";
    os_ << indent5 << "struct RecurInfo const * aminfo_n1 = aminfo_1 + n1;\n";
    os_ << indent5 << "struct RecurInfo const * aminfo_n2 = aminfo_2 + n2;\n";
    os_ << indent5 << "struct RecurInfo const * aminfo_n3 = aminfo_3 + n3;\n";
    os_ << indent5 << "struct RecurInfo const * aminfo_n4 = aminfo_4 + n4;\n";
    os_ << "\n\n";
    os_ << indent5 << "// First Center\n";
    os_ << indent5 << "idx_x = aminfo_n1->idx[0][2];\n";
    os_ << indent5 << "idx_y = aminfo_n1->idx[1][2];\n";
    os_ << indent5 << "idx_z = aminfo_n1->idx[2][2];\n";
    os_ << indent5 << "full_idx_x = idx_x * " << ncart2*ncart3*ncart4
                   << " + n2 * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
    os_ << indent5 << "full_idx_y = idx_y * " << ncart2*ncart3*ncart4
                   << " + n2 * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
    os_ << indent5 << "full_idx_z = idx_z * " << ncart2*ncart3*ncart4
                   << " + n2 * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             

    os_ << indent5 << outvar << "[startidx + 0] = " << invar_1p << "[full_idx_x];\n";
    os_ << indent5 << outvar << "[startidx + 1] = " << invar_1p << "[full_idx_y];\n";
    os_ << indent5 << outvar << "[startidx + 2] = " << invar_1p << "[full_idx_z];\n";

    if(am[0] > 0)
    {
        os_ << indent5 << "idx_x = aminfo_n1->idx[0][0];\n";
        os_ << indent5 << "idx_y = aminfo_n1->idx[1][0];\n";
        os_ << indent5 << "idx_z = aminfo_n1->idx[2][0];\n";
        os_ << indent5 << "full_idx_x = idx_x * " << ncart2*ncart3*ncart4
                       << " + n2 * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "full_idx_y = idx_y * " << ncart2*ncart3*ncart4
                       << " + n2 * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "full_idx_z = idx_z * " << ncart2*ncart3*ncart4
                       << " + n2 * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "if(idx_x >= 0) " << outvar << "[startidx + 0] -= " << invar_1m << "[full_idx_x];\n";
        os_ << indent5 << "if(idx_y >= 0) " << outvar << "[startidx + 1] -= " << invar_1m << "[full_idx_y];\n";
        os_ << indent5 << "if(idx_z >= 0) " << outvar << "[startidx + 2] -= " << invar_1m << "[full_idx_z];\n";
    }

    os_ << "\n\n";
    os_ << indent5 << "// Second Center\n";
    os_ << indent5 << "idx_x = aminfo_n2->idx[0][2];\n";
    os_ << indent5 << "idx_y = aminfo_n2->idx[1][2];\n";
    os_ << indent5 << "idx_z = aminfo_n2->idx[2][2];\n";
    os_ << indent5 << "full_idx_x = n1 * " << ncart2p*ncart3*ncart4
                   << " + idx_x * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
    os_ << indent5 << "full_idx_y = n1 * " << ncart2p*ncart3*ncart4
                   << " + idx_y * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
    os_ << indent5 << "full_idx_z = n1 * " << ncart2p*ncart3*ncart4
                   << " + idx_z * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             

    os_ << indent5 << outvar << "[startidx + 3] = " << invar_2p << "[full_idx_x];\n";
    os_ << indent5 << outvar << "[startidx + 4] = " << invar_2p << "[full_idx_y];\n";
    os_ << indent5 << outvar << "[startidx + 5] = " << invar_2p << "[full_idx_z];\n";

    if(am[1] > 0)
    {
        os_ << indent5 << "idx_x = aminfo_n2->idx[0][0];\n";
        os_ << indent5 << "idx_y = aminfo_n2->idx[1][0];\n";
        os_ << indent5 << "idx_z = aminfo_n2->idx[2][0];\n";
        os_ << indent5 << "full_idx_x = n1 * " << ncart2m*ncart3*ncart4
                       << " + idx_x * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "full_idx_y = n1 * " << ncart2m*ncart3*ncart4
                       << " + idx_y * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "full_idx_z = n1 * " << ncart2m*ncart3*ncart4
                       << " + idx_z * " << ncart3*ncart4 << " + n3 * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "if(idx_x >= 0) " << outvar << "[startidx + 3] -= " << invar_2m << "[full_idx_x];\n";
        os_ << indent5 << "if(idx_y >= 0) " << outvar << "[startidx + 4] -= " << invar_2m << "[full_idx_y];\n";
        os_ << indent5 << "if(idx_z >= 0) " << outvar << "[startidx + 5] -= " << invar_2m << "[full_idx_z];\n";
    }

    os_ << "\n\n";
    os_ << indent5 << "// Third Center\n";
    os_ << indent5 << "idx_x = aminfo_n3->idx[0][2];\n";
    os_ << indent5 << "idx_y = aminfo_n3->idx[1][2];\n";
    os_ << indent5 << "idx_z = aminfo_n3->idx[2][2];\n";
    os_ << indent5 << "full_idx_x = n1 * " << ncart2*ncart3p*ncart4
                   << " + n2 * " << ncart3p*ncart4 << " + idx_x * " << ncart4 << " + n4;\n";             
    os_ << indent5 << "full_idx_y = n1 * " << ncart2*ncart3p*ncart4
                   << " + n2 * " << ncart3p*ncart4 << " + idx_y * " << ncart4 << " + n4;\n";             
    os_ << indent5 << "full_idx_z = n1 * " << ncart2*ncart3p*ncart4
                   << " + n2 * " << ncart3p*ncart4 << " + idx_z * " << ncart4 << " + n4;\n";             

    os_ << indent5 << outvar << "[startidx + 6] = " << invar_3p << "[full_idx_x];\n";
    os_ << indent5 << outvar << "[startidx + 7] = " << invar_3p << "[full_idx_y];\n";
    os_ << indent5 << outvar << "[startidx + 8] = " << invar_3p << "[full_idx_z];\n";

    if(am[2] > 0)
    {
        os_ << indent5 << "idx_x = aminfo_n3->idx[0][0];\n";
        os_ << indent5 << "idx_y = aminfo_n3->idx[1][0];\n";
        os_ << indent5 << "idx_z = aminfo_n3->idx[2][0];\n";
        os_ << indent5 << "full_idx_x = n1 * " << ncart2*ncart3m*ncart4
                       << " + n2 * " << ncart3m*ncart4 << " + idx_x * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "full_idx_y = n1 * " << ncart2*ncart3m*ncart4
                       << " + n2 * " << ncart3m*ncart4 << " + idx_y * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "full_idx_z = n1 * " << ncart2*ncart3m*ncart4
                       << " + n2 * " << ncart3m*ncart4 << " + idx_z * " << ncart4 << " + n4;\n";             
        os_ << indent5 << "if(idx_x >= 0) " << outvar << "[startidx + 6] -= " << invar_3m << "[full_idx_x];\n";
        os_ << indent5 << "if(idx_y >= 0) " << outvar << "[startidx + 7] -= " << invar_3m << "[full_idx_y];\n";
        os_ << indent5 << "if(idx_z >= 0) " << outvar << "[startidx + 8] -= " << invar_3m << "[full_idx_z];\n";
    }

    os_ << "\n\n";
    os_ << indent5 << "// Fourth Center\n";
    os_ << indent5 << outvar << "[startidx +  9] = -(" << outvar << "[startidx + 0] + " << outvar << "[startidx + 3] + " << outvar << "[startidx + 6]);\n";
    os_ << indent5 << outvar << "[startidx + 10] = -(" << outvar << "[startidx + 1] + " << outvar << "[startidx + 4] + " << outvar << "[startidx + 7]);\n";
    os_ << indent5 << outvar << "[startidx + 11] = -(" << outvar << "[startidx + 2] + " << outvar << "[startidx + 5] + " << outvar << "[startidx + 8]);\n";

    os_ << indent5 << "startidx += 12;\n";
    os_ << indent4 << "} // close loop for forming derivatives\n\n";
    
}

void OSTEIDeriv1_Writer::Write_Full_(void) const
{
    const QAM am = info_.FinalAM();
    const int ncart = NCART(am);

    // some helper bools
    const bool hashrr = hrr_writer_.HasHRR();
    const bool hasbrahrr = hrr_writer_.HasBraHRR();
    const bool haskethrr = hrr_writer_.HasKetHRR();

    const bool hasbravrr = vrr_writer_.HasBraVRR();
    const bool hasketvrr = vrr_writer_.HasKetVRR();
    //const bool hasvrr = (hasbravrr || hasketvrr);

    //const bool hasoneoverp = hasbravrr;
    //const bool hasoneoverq = hasketvrr;
    //const bool hasoneover2p = (hasbravrr && (am[0]+am[1]) > 1); 
    //const bool hasoneover2q = (hasketvrr && (am[2]+am[3]) > 1); 
    //const bool hasoneover2pq = (hasketvrr && (am[0]+am[1]) > 0);
    const bool hasoneoverp = true;
    const bool hasoneoverq = true;
    const bool hasoneover2p = true;
    const bool hasoneover2q = true;
    const bool hasoneover2pq = true;


    // add includes
    IncludeSet includes{"<string.h>",
                        "<math.h>",
                        "\"simint/ostei/gen/ostei_deriv1_generated.h\"",
                        "\"simint/vectorization/vectorization.h\"",
                        "\"simint/ostei/recur_lookup.h\"",
                        "\"simint/boys/boys.h\""};

    // Constants
    ConstantMap cm;
    cm.emplace("const_1", "1");  // for 1/x

    // Note: vrr_writer handles the prim arrays for s_s_s_s
    //       so we always want to run this
    for(const auto & it : vrr_writer_.GetConstants())
        cm.insert(it);

    if(hrr_writer_.HasHRR())
        for(const auto & it : hrr_writer_.GetConstants())
            cm.insert(it);


    // need these factors sometimes
    if(hasoneover2p || hasoneover2q || hasoneover2pq)
        cm.emplace("one_half", "0.5");



    ///////////////////////////////////////
    // Beginning of file writing
    ///////////////////////////////////////

    // Write out all the includes
    for(const auto & it : includes)
        os_ << "#include " << it << "\n";
    os_ << "\n\n";

    //////////////////////////////
    // Function name & signature
    //////////////////////////////
    os_ << FunctionPrototype_(am) << "\n";
    os_ << "{\n";
    os_ << "\n";

    os_ << indent1 << "SIMINT_ASSUME_ALIGN_DBL(contwork);\n";
    os_ << indent1 << "SIMINT_ASSUME_ALIGN_DBL(" << ArrVarName(am) << ");\n";

    ///////////////////////////////////
    // NOW IN THE ACTUAL OSTEI FUNCTION
    ///////////////////////////////////

    // If there is no HRR, integrals are accumulated from inside the primitive loop
    // directly into the final integral array that was passed into this function, so it must be zeroed first
    if(!hashrr)
        os_ << indent1 << "memset(" << ArrVarName(am)
                       << ", 0, P.nshell12_clip * Q.nshell12_clip * "
                       << ncart << " * sizeof(double));\n\n";


    // abcd = index within simd loop, 
    os_ << indent1 << "int ab, cd, abcd;\n";
    os_ << indent1 << "int istart, jstart;\n";
    os_ << indent1 << "int iprimcd, nprim_icd, icd;\n";
    os_ << indent1 << "const int check_screen = (screen_tol > 0.0);\n";
    os_ << indent1 << "int i, j;\n";
    os_ << indent1 << "int n;\n";
    os_ << indent1 << "int not_screened;\n";
    

    // real_abcd is the absolute actual abcd in terms of all the shells that we are doing
    // (only needed if we do HRR)
    if(hashrr)
        os_ << indent1 << "int real_abcd;\n";


    if(hasbrahrr)
        os_ << indent1 << "int iket;\n";
    if(haskethrr)
        os_ << indent1 << "int ibra;\n";
    
    os_ << "\n";


    // Declare the temporary space 
    // Only needed if we are doing HRR
    if(hashrr)
        DeclareContwork();


    os_ << indent1 << "// Create constants\n";
    for(const auto & it : cm)
        os_ << indent1 << "const SIMINT_DBLTYPE " << it.first << " = SIMINT_DBLSET1(" << it.second << ");\n";

    os_ << "\n\n";
    os_ << indent1 << "////////////////////////////////////////\n";
    os_ << indent1 << "// Loop over shells and primitives\n";
    os_ << indent1 << "////////////////////////////////////////\n";
    os_ << "\n";

    if(hashrr)
        os_ << indent1 << "real_abcd = 0;\n";
    else
        os_ << indent1 << "abcd = 0;\n";

    os_ << indent1 << "istart = 0;\n";
    os_ << indent1 << "for(ab = 0; ab < P.nshell12_clip; ++ab)\n";
    os_ << indent1 << "{\n";

    os_ << indent2 << "const int iend = istart + P.nprim12[ab];\n";
    os_ << "\n";

    os_ << indent2 << "cd = 0;\n";
    os_ << indent2 << "jstart = 0;\n";
    os_ << "\n";

    os_ << indent2 << "for(cd = 0; cd < Q.nshell12_clip; cd += SIMINT_NSHELL_SIMD)\n";
    os_ << indent2 << "{\n";
    os_ << indent3 << "const int nshellbatch = ((cd + SIMINT_NSHELL_SIMD) > Q.nshell12_clip) ? Q.nshell12_clip - cd : SIMINT_NSHELL_SIMD;\n";

    os_ << indent3 << "int jend = jstart;\n";
    os_ << indent3 << "for(i = 0; i < nshellbatch; i++)\n";
    os_ << indent4 << "jend += Q.nprim12[cd+i];\n";
    os_ << "\n";


    if(hashrr)
    {
        size_t contmem = info_.ContMemoryReq();
        if(contmem > 0)
            os_ << indent3 << "memset(contwork, 0, SIMINT_NSHELL_SIMD * " << contmem << ");\n";

        os_ << indent3 << "abcd = 0;\n";
        os_ << "\n";
    }

    os_ << "\n";
    os_ << indent3 << "for(i = istart; i < iend; ++i)\n";
    os_ << indent3 << "{\n";
    os_ << "\n";

    os_ << indent4 << "// Skip this whole thing if always insignificant\n";
    os_ << indent4 << "if(check_screen && (P.screen[i] * Q.screen_max) < screen_tol)\n";
    os_ << indent5 << "continue;\n\n";

    os_ << indent4 << "icd = 0;\n";
    os_ << indent4 << "iprimcd = 0;\n";
    os_ << indent4 << "nprim_icd = Q.nprim12[cd];\n";

    DeclarePrimPointers();
    os_ << "\n";

    os_ << indent4 << "// Load these one per loop over i\n";
    os_ << indent4 << "const SIMINT_DBLTYPE P_alpha = SIMINT_DBLSET1(P.alpha[i]);\n";
    os_ << indent4 << "const SIMINT_DBLTYPE P_prefac = SIMINT_DBLSET1(P.prefac[i]);\n";
    os_ << indent4 << "const SIMINT_DBLTYPE Pxyz[3] = { SIMINT_DBLSET1(P.x[i]), SIMINT_DBLSET1(P.y[i]), SIMINT_DBLSET1(P.z[i]) };\n";
                             
    os_ << indent4 << "const SIMINT_DBLTYPE bra_screen_max = SIMINT_DBLSET1(P.screen[i]);\n";
    os_ << "\n";
    os_ << indent4 << "// Contraction factors (needed for derivatives)\n";
    os_ << indent4 << "const SIMINT_DBLTYPE cfac_2a = SIMINT_DBLSET1(P.alpha2[i]);\n";
    os_ << indent4 << "const SIMINT_DBLTYPE cfac_2b = SIMINT_DBLSET1(P.beta2[i]);\n\n";

    if(hasbravrr)
    {
        if(vrr_writer_.HasVRR_I())
            os_ << indent4 << "const SIMINT_DBLTYPE P_PA[3] = { SIMINT_DBLSET1(P.PA_x[i]), SIMINT_DBLSET1(P.PA_y[i]), SIMINT_DBLSET1(P.PA_z[i]) };\n";
        else
            os_ << indent4 << "const SIMINT_DBLTYPE P_PB[3] = { SIMINT_DBLSET1(P.PB_x[i]), SIMINT_DBLSET1(P.PB_y[i]), SIMINT_DBLSET1(P.PB_z[i]) };\n";
    }

    os_ << "\n";


    os_ << indent4 << "for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)\n";
    os_ << indent4 << "{\n";

    WriteShellOffsets();


    os_ << indent5 << "// Do we have to compute this vector (or has it been screened out)?\n";
    os_ << indent5 << "// (not_screened != 0 means we have to do this vector)\n";
    os_ << indent5 << "if(check_screen)\n";
    os_ << indent5 << "{\n";
    os_ << indent6 << "const double vmax = vector_max(SIMINT_MUL(bra_screen_max, SIMINT_DBLLOAD(Q.screen, j)));\n";
    os_ << indent6 << "if(vmax > screen_tol)\n";
    os_ << indent6 << "{\n";
    for(const auto it : info_.GetContQ())
        os_ << indent7 << PrimPtrName(it) << " += lastoffset*" << NCART(it) << ";\n";
    os_ << indent7 << "continue;\n";
    os_ << indent6 << "}\n";
    os_ << indent5 << "}\n\n";


    // Note: vrr_writer handles the prim arrays for s_s_s_s
    //       so we always want to run this
    vrr_writer_.DeclarePrimArrays(os_);

    os_ << indent5 << "const SIMINT_DBLTYPE Q_alpha = SIMINT_DBLLOAD(Q.alpha, j);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE PQalpha_mul = SIMINT_MUL(P_alpha, Q_alpha);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE PQalpha_sum = SIMINT_ADD(P_alpha, Q_alpha);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE one_over_PQalpha_sum = SIMINT_DIV(const_1, PQalpha_sum);\n";
    os_ << "\n";
    os_ << indent5 << "// Contraction factors (needed for derivatives)\n";
    os_ << indent5 << "const SIMINT_DBLTYPE cfac_2c = SIMINT_DBLLOAD(Q.alpha2, j);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE cfac_2d = SIMINT_DBLLOAD(Q.beta2, j);\n\n";
    os_ << "\n";
    os_ << indent5 << "/* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os_ << indent5 << "SIMINT_DBLTYPE PQ[3];\n";
    os_ << indent5 << "PQ[0] = SIMINT_SUB(Pxyz[0], SIMINT_DBLLOAD(Q.x, j));\n";
    os_ << indent5 << "PQ[1] = SIMINT_SUB(Pxyz[1], SIMINT_DBLLOAD(Q.y, j));\n";
    os_ << indent5 << "PQ[2] = SIMINT_SUB(Pxyz[2], SIMINT_DBLLOAD(Q.z, j));\n";


    os_ << indent5 << "SIMINT_DBLTYPE R2 = SIMINT_MUL(PQ[0], PQ[0]);\n";
    os_ << indent5 << "R2 = SIMINT_FMADD(PQ[1], PQ[1], R2);\n";
    os_ << indent5 << "R2 = SIMINT_FMADD(PQ[2], PQ[2], R2);\n";
    os_ << "\n";
    os_ << indent5 << "const SIMINT_DBLTYPE alpha = SIMINT_MUL(PQalpha_mul, one_over_PQalpha_sum); // alpha from MEST\n";

    if(hasoneoverp)
        os_ << indent5 << "const SIMINT_DBLTYPE one_over_p = SIMINT_DIV(const_1, P_alpha);\n";

    if(hasoneoverq)
        os_ << indent5 << "const SIMINT_DBLTYPE one_over_q = SIMINT_DIV(const_1, Q_alpha);\n";

    if(hasoneover2p)    
        os_ << indent5 << "const SIMINT_DBLTYPE one_over_2p = SIMINT_MUL(one_half, one_over_p);\n";

    if(hasoneover2q)    
        os_ << indent5 << "const SIMINT_DBLTYPE one_over_2q = SIMINT_MUL(one_half, one_over_q);\n";

    if(hasoneover2pq)
        os_ << indent5 << "const SIMINT_DBLTYPE one_over_2pq = SIMINT_MUL(one_half, one_over_PQalpha_sum);\n";

    if(hasketvrr)
    {
        if(vrr_writer_.HasVRR_K())
            os_ << indent5 << "const SIMINT_DBLTYPE Q_PA[3] = { SIMINT_DBLLOAD(Q.PA_x, j), SIMINT_DBLLOAD(Q.PA_y, j), SIMINT_DBLLOAD(Q.PA_z, j) };\n"; 
        else
            os_ << indent5 << "const SIMINT_DBLTYPE Q_PB[3] = { SIMINT_DBLLOAD(Q.PB_x, j), SIMINT_DBLLOAD(Q.PB_y, j), SIMINT_DBLLOAD(Q.PB_z, j) };\n"; 
    }

    if(hasbravrr)
    {
        os_ << "\n";
        os_ << indent5 << "// NOTE: Minus sign!\n";
        os_ << indent5 << "const SIMINT_DBLTYPE a_over_p = SIMINT_MUL(SIMINT_NEG(alpha), one_over_p);\n";
        os_ << indent5 << "SIMINT_DBLTYPE aop_PQ[3];\n";
        os_ << indent5 << "aop_PQ[0] = SIMINT_MUL(a_over_p, PQ[0]);\n";
        os_ << indent5 << "aop_PQ[1] = SIMINT_MUL(a_over_p, PQ[1]);\n";
        os_ << indent5 << "aop_PQ[2] = SIMINT_MUL(a_over_p, PQ[2]);\n";
    }

    if(hasketvrr)
    {
        os_ << "\n";
        os_ << indent5 << "SIMINT_DBLTYPE a_over_q = SIMINT_MUL(alpha, one_over_q);\n";
        os_ << indent5 << "SIMINT_DBLTYPE aoq_PQ[3];\n";
        os_ << indent5 << "aoq_PQ[0] = SIMINT_MUL(a_over_q, PQ[0]);\n";
        os_ << indent5 << "aoq_PQ[1] = SIMINT_MUL(a_over_q, PQ[1]);\n";
        os_ << indent5 << "aoq_PQ[2] = SIMINT_MUL(a_over_q, PQ[2]);\n";

        os_ << indent5 << "// Put a minus sign here so we don't have to in RR routines\n";
        os_ << indent5 << "a_over_q = SIMINT_NEG(a_over_q);\n";
    }

    os_ << "\n";
    os_ << "\n";
    os_ << indent5 << "//////////////////////////////////////////////\n";
    os_ << indent5 << "// Fjt function section\n";
    os_ << indent5 << "// Maximum v value: " << (info_.L()+1) << "\n";
    os_ << indent5 << "//////////////////////////////////////////////\n";
    os_ << indent5 << "// The parameter to the Fjt function\n";
    os_ << indent5 << "const SIMINT_DBLTYPE F_x = SIMINT_MUL(R2, alpha);\n";
    os_ << "\n";
    os_ << "\n";

    // we need to zero out any that are beyond the end of the batch (that's been clipped)
    os_ << indent5 << "const SIMINT_DBLTYPE Q_prefac = mask_load(nlane, Q.prefac + j);\n";
    os_ << "\n\n"; 
    os_ << indent5 << "boys_F_split(" << PrimVarName({0,0,0,0})
                   << ", F_x, " << (info_.L()+1) << ");\n";


    // prefac = sqrt(1/PQalpha_sum) * P_prefac * Q_prefac
    os_ << indent5 << "SIMINT_DBLTYPE prefac = SIMINT_SQRT(one_over_PQalpha_sum);\n";
    os_ << indent5 << "prefac = SIMINT_MUL(SIMINT_MUL(P_prefac, Q_prefac), prefac);\n";

    const std::string name0000 = PrimVarName({0,0,0,0});
    const std::string name0000n = name0000 + "[n]";

    os_ << indent5 << "for(n = 0; n <= " << (info_.L()+1) << "; n++)\n"
        << indent6 << name0000n << " = SIMINT_MUL(" << name0000n << ", prefac);\n";


    if(vrr_writer_.HasVRR())
        vrr_writer_.WriteVRR(os_);

    WriteAccumulation();

        
    os_ << "\n";
    os_ << indent4 << "}  // close loop over j\n";
    os_ << indent3 << "}  // close loop over i\n";

    os_ << indent3 << "\n";
    os_ << indent3 << "//Advance to the next batch\n";
    os_ << indent3 << "jstart = SIMINT_SIMD_ROUND(jend);\n";
    if(!hashrr)
        os_ << indent3 << "abcd += nshellbatch;\n";
    os_ << indent3 << "\n";

    hrr_writer_.WriteHRR(os_);
    WriteFormDeriv();

    os_ << "\n";
    os_ << indent3 << "}  // close HRR loop\n";
    os_ << "\n\n";

    os_ << indent2 << "}   // close loop cdbatch\n";

    os_ << "\n";
    os_ << indent2 << "istart = iend;\n";

    os_ << indent1 << "}  // close loop over ab\n";
    os_ << "\n";
    os_ << indent1 << "return P.nshell12_clip * Q.nshell12_clip;\n";
    os_ << "}\n";
    os_ << "\n";
}


void OSTEIDeriv1_Writer::Write_Permutations_(void) const
{
    ///////////////////////////////////////////////////////////
    // Note that we never permute bra, ket. That would
    // affect vectorization since we only vectorize on the ket
    ///////////////////////////////////////////////////////////
    QAM am = info_.FinalAM();

    // permute 1,2
    Write_Permute_(am, true, false);

    // permute 3,4
    Write_Permute_(am, false, true);

    // permute 1,2 and 3,4
    Write_Permute_(am, true, true);
}


void OSTEIDeriv1_Writer::WriteFile(void) const
{
    const QAM am = info_.FinalAM();

    // is this a special permutation? Handle it if so.
    Write_Full_();

    // Add to the header
    osh_ << FunctionPrototype_(am) << ";\n\n";

    // Write out the code for permuting final integrals, if necessary
    if(info_.FinalPermute())
        Write_Permutations_();
}


