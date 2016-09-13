#include "generator/Types.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"
#include "generator/BoysGenerator.hpp"
#include "generator/eri/ERIGeneratorInfo.hpp"
#include "generator/eri/ERI_VRR_Writer.hpp"
#include "generator/eri/ERI_ET_Writer.hpp"
#include "generator/eri/ERI_HRR_Writer.hpp"
#include "generator/eri/ERI_Writer.hpp"


///////////////////////////
// ERI_Writer Base Class //
///////////////////////////
ERI_Writer::ERI_Writer(std::ostream & os,
                       std::ostream & osh,
                       const ERIGeneratorInfo & info,
                       const BoysGenerator & bg,
                       const ERI_VRR_Writer & vrr_writer,
                       const ERI_ET_Writer & et_writer,
                       const ERI_HRR_Writer & hrr_writer)
   : os_(os), osh_(osh), info_(info), vinfo_(info_.GetVectorInfo()), bg_(bg),
     vrr_writer_(vrr_writer), et_writer_(et_writer), hrr_writer_(hrr_writer)
{ }





/////////////////////////////
// Basic ERI Writer
/////////////////////////////
void ERI_Writer_Basic::DeclareContwork(void) const
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


void ERI_Writer_Basic::ZeroContwork(void) const
{
    size_t contmem = info_.ContMemoryReq();
    if(contmem > 0)
        os_ << indent3 << "memset(contwork, 0, SIMINT_NSHELL_SIMD * " << contmem << ");\n";
    
}


void ERI_Writer_Basic::FreeContwork(void) const
{
    size_t contmem = info_.ContMemoryReq();

    if(contmem > 0 && info_.UseHeap())
    {
        os_ << indent1 << "// Free contracted workspace\n";
        os_ << indent1 << "FREE(contwork);\n\n";
    }       
}

void ERI_Writer_Basic::WriteShellOffsets(void) const
{
    os_ << indent5 << "// calculate the shell offsets\n";
    os_ << indent5 << "// these are the offset from the shell pointed to by cd\n";
    os_ << indent5 << "// for each element\n";
    os_ << indent5 << "int shelloffsets[SIMINT_SIMD_LEN] = {0};\n";
    os_ << indent5 << "int lastoffset = 0;\n";
    os_ << indent5 << "if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)\n";
    os_ << indent5 << "{\n";

    os_ << indent6 << "// Handle if the first element of the vector is a new shell\n";
    os_ << indent6 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os_ << indent6 << "{\n";
    os_ << indent7 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto qam : info_.GetContQ())
        os_ << indent7 << PrimPtrName(qam) << " += " << NCART(qam) << ";\n";

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


void ERI_Writer_Basic::WriteShellOffsets_Scalar(void) const
{
    os_ << indent5 << "// Move pointers if this is the end of a shell\n";
    os_ << indent5 << "// Handle if the first element of the vector is a new shell\n";
    os_ << indent5 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os_ << indent5 << "{\n";
    os_ << indent6 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto qam : info_.GetContQ())
        os_ << indent7 << PrimPtrName(qam) << " += " << NCART(qam) << ";\n";

    os_ << indent5 << "}\n";
    os_ << indent5 << "iprimcd++;\n";
}


void ERI_Writer_Basic::WriteAccumulation(void) const
{
    os_ << "\n\n";
    os_ << indent5 << "////////////////////////////////////\n";
    os_ << indent5 << "// Accumulate contracted integrals\n";
    os_ << indent5 << "////////////////////////////////////\n";
    os_ << indent5 << "if(lastoffset == 0)\n";
    os_ << indent5 << "{\n";

    for(const auto qam : info_.GetContQ())
    {
        int ncart = NCART(qam);
        os_ << indent6 << "contract_all(SIMINT_SIMD_LEN, " << ncart << ", " << PrimVarName(qam) << ", " << PrimPtrName(qam) << ");\n";
    }
    os_ << indent5 << "}\n";
    os_ << indent5 << "else\n";
    os_ << indent5 << "{\n";

    for(const auto qam : info_.GetContQ())
    {
        int ncart = NCART(qam);
        os_ << indent6 << "contract(SIMINT_SIMD_LEN, " << ncart << ", shelloffsets, " << PrimVarName(qam) << ", " << PrimPtrName(qam) << ");\n";
    }

    for(const auto qam : info_.GetContQ())
        os_ << indent6 << PrimPtrName(qam) << " += lastoffset*" << NCART(qam) << ";\n";

    os_ << indent5 << "}\n";
}


void ERI_Writer_Basic::WriteFile_Permute_(void) const
{
    QAM am = info_.FinalAM();
    QAM tocall = am;

    bool swap_ab = false;
    bool swap_cd = false;

    // TODO - more thoroughly check?
    if(am[0] == 0)
    {
        swap_ab = true;
        std::swap(tocall[0], tocall[1]);
    }
    if(am[2] == 0)
    {
        swap_cd = true;
        std::swap(tocall[2], tocall[3]);
    }

    std::string funcline = StringBuilder("int eri_sharedwork_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string indent(funcline.length(), ' ');

    // start output to the file
    // we only need this one include
    os_ << "#include \"simint/eri/eri.h\"\n";
    os_ << "\n";

    os_ << "\n\n";
    os_ << funcline;
    os_ << "struct simint_multi_shellpair const P,\n";
    os_ << indent << "struct simint_multi_shellpair const Q,\n";
    os_ << indent << "double * const restrict contwork,\n";
    os_ << indent << "double * const restrict " << ArrVarName(am) << ")\n";
    os_ << "{\n";
    os_ << indent1 << "// Can be accomplished by swapping some variables\n";
    os_ << indent1 << "// and calling another function\n";
    os_ << indent1 << "// Note that the struct was passed by copy\n";
    os_ << "\n";


    const char * P_var = "P";
    const char * Q_var = "Q";

    if(swap_ab)
    {
    	os_ << indent1 << "struct simint_multi_shellpair P_tmp = P;\n";
        os_ << indent1 << "P_tmp.PA_x = P.PB_x;  P_tmp.PA_y = P.PB_y;  P_tmp.PA_z = P.PB_z;\n";
        os_ << indent1 << "P_tmp.PB_x = P.PA_x;  P_tmp.PB_y = P.PA_y;  P_tmp.PB_z = P.PA_z;\n";
        P_var = "P_tmp";
    }

    if(swap_cd)
    {
    	os_ << indent1 << "struct simint_multi_shellpair Q_tmp = Q;\n";
		os_ << indent1 << "Q_tmp.PA_x = Q.PB_x;  Q_tmp.PA_y = Q.PB_y;  Q_tmp.PA_z = Q.PB_z;\n";
        os_ << indent1 << "Q_tmp.PB_x = Q.PA_x;  Q_tmp.PB_y = Q.PA_y;  Q_tmp.PB_z = Q.PA_z;\n";
	    Q_var = "Q_tmp";
    }

 
    os_ << indent1 << "return eri_sharedwork_" << amchar[tocall[0]] << "_" 
                                               << amchar[tocall[1]] << "_" 
                                               << amchar[tocall[2]] << "_" 
                                               << amchar[tocall[3]] << "(" << P_var << ", " << Q_var 
                                                                    << ", contwork, " << ArrVarName(am) << ");\n"; 

    os_ << "}\n";
    os_ << "\n";

}


void ERI_Writer_Basic::WriteFile_NoPermute_(void) const
{
    const QAM am = info_.FinalAM();
    const int ncart = NCART(am);

    // some helper bools
    bool hashrr = hrr_writer_.HasHRR();
    bool hasbrahrr = hrr_writer_.HasBraHRR();
    bool haskethrr = hrr_writer_.HasKetHRR();
    bool inline_hrr = (hashrr && hrr_writer_.IsInline());

    bool hasbravrr = vrr_writer_.HasBraVRR();
    bool hasketvrr = vrr_writer_.HasKetVRR();

    //bool haset = et_writer_.HasET();
    bool hasbraet = et_writer_.HasBraET(); 
    bool hasketet = et_writer_.HasKetET(); 

    bool hasoneoverp = (hasbravrr || hasbraet);
    bool hasoneoverq = (hasketvrr || hasketet);
    bool hasoneover2p = (hasbraet || (hasbravrr && (am[0]+am[1]) > 1)); 
    bool hasoneover2q = (hasketet || (hasketvrr && (am[2]+am[3]) > 1)); 
    bool hasoneover2pq = (hasketvrr && (am[0]+am[1]) > 0);


    std::string dbltype = vinfo_.DoubleType();
    std::string cdbltype = vinfo_.ConstDoubleType();


    // add includes
    IncludeSet includes{"<string.h>", "<math.h>", "\"simint/eri/eri.h\"", "\"simint/eri/eri_contract.h\""};
    for(const auto & it : bg_.GetIncludes())
        includes.insert(it);

    // Constants
    ConstantMap cm;
    cm.emplace("const_1", "1");  // for 1/x
    for(const auto & it : bg_.GetConstants())
        cm.insert(it);

    // Note: vrr_writer handles the prim arrays for s_s_s_s
    //       so we always want to run this
    for(const auto & it : vrr_writer_.GetConstants())
        cm.insert(it);

    if(et_writer_.HasET())
        for(const auto & it : et_writer_.GetConstants())
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


    //////////////////////////////
    // Function name & signature
    //////////////////////////////
    std::string funcline = StringBuilder("int eri_sharedwork_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string indent(funcline.length(), ' ');


    os_ << "\n\n";
    os_ << funcline;
    os_ << "struct simint_multi_shellpair const P,\n";
    os_ << indent << "struct simint_multi_shellpair const Q,\n";
    os_ << indent << "double * const restrict contwork,\n";
    os_ << indent << "double * const restrict " << ArrVarName(am) << ")\n";
    os_ << "{\n";
    os_ << "\n";


    ///////////////////////////////////
    // NOW IN THE ACTUAL ERI FUNCTION
    ///////////////////////////////////

    // If there is no HRR, integrals are accumulated from inside the primitive loop
    // directly into the final integral array that was passed into this function, so it must be zeroed first
    if(!hashrr)
        os_ << indent1 << "memset(" << ArrVarName(am) << ", 0, P.nshell12 * Q.nshell12 * " << ncart << " * sizeof(double));\n";
    os_ << "\n";


    // abcd = index within simd loop, 
    os_ << indent1 << "int ab, cd, cdbatch, abcd;\n";
    os_ << indent1 << "int istart, jstart;\n";
    os_ << indent1 << "int iprimcd, nprim_icd, icd;\n";
    os_ << indent1 << "int i, j, n;\n";

    // real_abcd is the absolute actual abcd in terms of all the shells that we are doing
    // (only needed if we do HRR)
    if(hashrr)
        os_ << indent1 << "int real_abcd;\n";


    // Needed only if we are doing inline HRR
    if(inline_hrr)
    {
        if(hasbrahrr)
            os_ << indent1 << "int iket;\n";
        if(haskethrr)
            os_ << indent1 << "int ibra;\n";
    }

    os_ << "\n";


    // Declare the temporary space 
    // Only needed if we are doing HRR
    if(hashrr)
        DeclareContwork();



    // Write out all the constants 
    // This is only needed if this is NOT scalar code
    if(info_.Vectorized())
    {
        os_ << indent1 << "// Create constants\n";
        for(const auto & it : cm)
            os_ << indent1 << vinfo_.NewConstDoubleSet1(it.first, it.second) << ";\n";
    }

    // some constants used here
    std::string const1_str = vinfo_.IntConstant(1);
    std::string consthalf_str = info_.Vectorized() ? "one_half" : "0.5";


    os_ << "\n\n";
    os_ << indent1 << "////////////////////////////////////////\n";
    os_ << indent1 << "// Loop over shells and primitives\n";
    os_ << indent1 << "////////////////////////////////////////\n";
    os_ << "\n";

    os_ << indent1 << "// Number of batches of Q\n";
    os_ << indent1 << "int nbatch_Q = Q.nshell12_slice / SIMINT_NSHELL_SIMD;\n";
    os_ << indent1 << "if(Q.nshell12_slice % SIMINT_NSHELL_SIMD)\n";
    os_ << indent2 << "    nbatch_Q++;\n\n";


    if(hashrr)
        os_ << indent1 << "real_abcd = 0;\n";
    else
        os_ << indent1 << "abcd = 0;\n";

    os_ << indent1 << "istart = 0;\n";
    os_ << indent1 << "for(ab = 0; ab < P.nshell12; ++ab)\n";
    os_ << indent1 << "{\n";

    os_ << indent2 << "const int iend = istart + P.nprim12[ab];\n";
    os_ << "\n";

    os_ << indent2 << "cd = 0;\n";
    os_ << indent2 << "jstart = 0;\n";
    os_ << "\n";

    os_ << indent2 << "for(cdbatch = 0; cdbatch < nbatch_Q; ++cdbatch)\n";
    os_ << indent2 << "{\n";
    os_ << indent3 << "const int nshellbatch = ((cd + SIMINT_NSHELL_SIMD) > Q.nshell12) ? Q.nshell12 - cd : SIMINT_NSHELL_SIMD;\n";

    os_ << indent3 << "const int jend = jstart + Q.nbatchprim[cdbatch];\n";


    if(hashrr)
    {
        ZeroContwork();
        os_ << indent3 << "abcd = 0;\n";
        os_ << "\n";
    }

    os_ << indent3 << "for(i = istart; i < iend; ++i)\n";
    os_ << indent3 << "{\n";
    os_ << "\n";

    os_ << indent4 << "icd = 0;\n";
    os_ << indent4 << "iprimcd = 0;\n";
    os_ << indent4 << "nprim_icd = Q.nprim12[cd];\n";

    // Note: vrr_writer handles the prim pointers for s_s_s_s
    //       so we always want to run this
    vrr_writer_.DeclarePrimPointers(os_);
    if(et_writer_.HasET())
        et_writer_.DeclarePrimPointers(os_);
    os_ << "\n";

    os_ << indent4 << "// Load these one per loop over i\n";
    os_ << indent4 << vinfo_.NewConstDoubleSet1("P_alpha", "P.alpha[i]") << ";\n";
    os_ << indent4 << vinfo_.NewConstDoubleSet1("P_prefac", "P.prefac[i]") << ";\n";
    os_ << indent4 << vinfo_.NewConstDoubleSet1("P_x", "P.x[i]") << ";\n";
    os_ << indent4 << vinfo_.NewConstDoubleSet1("P_y", "P.y[i]") << ";\n";
    os_ << indent4 << vinfo_.NewConstDoubleSet1("P_z", "P.z[i]") << ";\n";

    if(hasbravrr)
    {
        if(vrr_writer_.HasVRR_I())
        {
            os_ << indent4 << vinfo_.NewConstDoubleSet1("P_PA_x", "P.PA_x[i]") << ";\n";
            os_ << indent4 << vinfo_.NewConstDoubleSet1("P_PA_y", "P.PA_y[i]") << ";\n";
            os_ << indent4 << vinfo_.NewConstDoubleSet1("P_PA_z", "P.PA_z[i]") << ";\n";
        }
        else
        {
            os_ << indent4 << vinfo_.NewConstDoubleSet1("P_PB_x", "P.PB_x[i]") << ";\n";
            os_ << indent4 << vinfo_.NewConstDoubleSet1("P_PB_y", "P.PB_y[i]") << ";\n";
            os_ << indent4 << vinfo_.NewConstDoubleSet1("P_PB_z", "P.PB_z[i]") << ";\n";
        }
    }

    if(hasketet)
    {
        os_ << indent4 << vinfo_.NewConstDoubleSet1("P_bAB_x", "P.bAB_x[i]") << ";\n";
        os_ << indent4 << vinfo_.NewConstDoubleSet1("P_bAB_y", "P.bAB_y[i]") << ";\n";
        os_ << indent4 << vinfo_.NewConstDoubleSet1("P_bAB_z", "P.bAB_z[i]") << ";\n";
    }

    os_ << "\n";


    os_ << indent4 << "for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)\n";
    os_ << indent4 << "{\n";

    if(info_.Vectorized())
        WriteShellOffsets();
    else
        WriteShellOffsets_Scalar();

    os_ << "\n";

    // Note: vrr_writer handles the prim arrays for s_s_s_s
    //       so we always want to run this
    vrr_writer_.DeclarePrimArrays(os_);
    if(et_writer_.HasET())
        et_writer_.DeclarePrimArrays(os_);

    os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_alpha", "Q.alpha", "j") << ";\n";
    os_ << indent5 << cdbltype << " PQalpha_mul = P_alpha * Q_alpha;\n";
    os_ << indent5 << cdbltype << " PQalpha_sum = P_alpha + Q_alpha;\n";
    os_ << indent5 << cdbltype << " one_over_PQalpha_sum = " << const1_str << " / PQalpha_sum;\n";
    os_ << "\n";
    os_ << "\n";
    os_ << indent5 << "/* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os_ << indent5 << cdbltype << " PQ_x = P_x - " << vinfo_.DoubleLoad("Q.x", "j") << ";\n";
    os_ << indent5 << cdbltype << " PQ_y = P_y - " << vinfo_.DoubleLoad("Q.y", "j") << ";\n";
    os_ << indent5 << cdbltype << " PQ_z = P_z - " << vinfo_.DoubleLoad("Q.z", "j") << ";\n";


    os_ << indent5 << cdbltype << " R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os_ << "\n";
    os_ << indent5 << cdbltype << " alpha = PQalpha_mul * one_over_PQalpha_sum;   // alpha from MEST\n";

    if(hasoneoverp)
        os_ << indent5 << cdbltype << " one_over_p = " << const1_str << " / P_alpha;\n";

    if(hasoneoverq)
        os_ << indent5 << cdbltype << " one_over_q = " << const1_str << " / Q_alpha;\n";

    if(hasoneover2p)    
        os_ << indent5 << cdbltype << " one_over_2p = " << consthalf_str << " * one_over_p;  // gets multiplied in VRR\n";

    if(hasoneover2q)    
        os_ << indent5 << cdbltype << " one_over_2q = " << consthalf_str << " * one_over_q;  // gets multiplied in VRR\n";

    if(hasoneover2pq)
        os_ << indent5 << cdbltype << " one_over_2pq = " << consthalf_str << " * one_over_PQalpha_sum;\n";

    if(hasketvrr)
    {
        if(vrr_writer_.HasVRR_K())
        {
            os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_PA_x", "Q.PA_x", "j") << ";\n";
            os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_PA_y", "Q.PA_y", "j") << ";\n";
            os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_PA_z", "Q.PA_z", "j") << ";\n";
        }
        else
        {
            os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_PB_x", "Q.PB_x", "j") << ";\n";
            os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_PB_y", "Q.PB_y", "j") << ";\n";
            os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_PB_z", "Q.PB_z", "j") << ";\n";
        }
    }

    if(hasbravrr)
    {
        os_ << "\n";
        os_ << indent5 << "// NOTE: Minus sign!\n";
        os_ << indent5 << cdbltype << " a_over_p =  -alpha * one_over_p;     // a/p from MEST\n";
        os_ << indent5 << cdbltype << " aop_PQ_x = a_over_p * PQ_x;\n"; 
        os_ << indent5 << cdbltype << " aop_PQ_y = a_over_p * PQ_y;\n"; 
        os_ << indent5 << cdbltype << " aop_PQ_z = a_over_p * PQ_z;\n"; 
    }

    if(hasketvrr)
    {
        os_ << "\n";
        os_ << indent5 << "// NOTE: Minus sign\n";
        os_ << indent5 << "// NOTE2: Plus sign taken care of on aoq_PQ_x!\n";
        os_ << indent5 << cdbltype << " a_over_q =  -alpha * one_over_q;     // a/q from MEST\n";
        os_ << indent5 << cdbltype << " aoq_PQ_x = -a_over_q * PQ_x;\n"; 
        os_ << indent5 << cdbltype << " aoq_PQ_y = -a_over_q * PQ_y;\n"; 
        os_ << indent5 << cdbltype << " aoq_PQ_z = -a_over_q * PQ_z;\n"; 

    }

    if(hasketet || hasbraet)
    {
        os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_bAB_x", "Q.bAB_x", "j") << ";\n";
        os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_bAB_y", "Q.bAB_y", "j") << ";\n";
        os_ << indent5 << vinfo_.NewConstDoubleLoad("Q_bAB_z", "Q.bAB_z", "j") << ";\n";
    }

    if(hasketet)
    {
        os_ << "\n";
        os_ << indent5 << cdbltype << " p_over_q = P_alpha * one_over_q;\n";
        os_ << indent5 << cdbltype << " etfac_k[3] = {\n";
        os_ << indent6 << "-(P_bAB_x + Q_bAB_x) * one_over_q,\n";
        os_ << indent6 << "-(P_bAB_y + Q_bAB_y) * one_over_q,\n";
        os_ << indent6 << "-(P_bAB_z + Q_bAB_z) * one_over_q,\n";
        os_ << indent6 << "};\n";
        os_ << "\n";
    }

    if(hasbraet)
    {
        os_ << "\n";
        os_ << indent5 << cdbltype << " q_over_p = Q_alpha * one_over_p;\n";
        os_ << indent5 << cdbltype << " etfac_b[3] = {\n";
        os_ << indent6 << "-(P_bAB_x + Q_bAB_x) * one_over_p,\n";
        os_ << indent6 << "-(P_bAB_y + Q_bAB_y) * one_over_p,\n";
        os_ << indent6 << "-(P_bAB_z + Q_bAB_z) * one_over_p,\n";
        os_ << indent6 << "};\n";
        os_ << "\n";
    }

    os_ << "\n";
    os_ << "\n";
    os_ << indent5 << "//////////////////////////////////////////////\n";
    os_ << indent5 << "// Boys function section\n";
    os_ << indent5 << "// Maximum v value: " << info_.L() << "\n";
    os_ << indent5 << "//////////////////////////////////////////////\n";
    os_ << indent5 << "// The parameter to the boys function\n";
    os_ << indent5 << cdbltype << " F_x = R2 * alpha;\n";
    os_ << "\n";
    os_ << "\n";

    bg_.WriteBoys(os_);

    if(vrr_writer_.HasVRR())
        vrr_writer_.WriteVRR(os_);

    if(et_writer_.HasET())
        et_writer_.WriteET(os_);

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

    if(hrr_writer_.HasHRR())
        hrr_writer_.WriteHRR(os_);

    os_ << "\n";

    os_ << indent3 << "cd += nshellbatch;\n";

    os_ << indent2 << "}   // close loop cdbatch\n";

    os_ << "\n";
    os_ << indent2 << "istart = iend;\n";

    os_ << indent2 << "// if this is the end of a batch in the bra part, skip the padding\n";
    os_ << indent2 << "if( ((ab+1) % SIMINT_NSHELL_SIMD) == 0)\n";
    os_ << indent3 << "istart = SIMINT_SIMD_ROUND(istart);\n";
    os_ << "\n";

    os_ << indent1 << "}  // close loop over ab\n";
    os_ << "\n";
    os_ << "\n";

    os_ << "\n";


    os_ << indent1 << "return P.nshell12 * Q.nshell12;\n";
    os_ << "}\n";
    os_ << "\n";
}



void ERI_Writer_Basic::WriteFile(void) const
{
    const QAM am = info_.FinalAM();

    // is this a special permutation? Handle it if so.
    if( ( (am[0] == 0 && am[1] > 0)  && ( am[2] == 0 || am[3] == 0 ) ) ||
        ( (am[2] == 0 && am[3] > 0)  && ( am[0] == 0 || am[1] == 0 ) ) )
    {
        WriteFile_Permute_();
    }
    else
        WriteFile_NoPermute_();


    // for header and for non-shared-work version
    // note - no return type
    std::string funcline = StringBuilder("eri_sharedwork_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string funcline2 = StringBuilder("eri_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string funcindent = std::string(funcline.length()+4, ' ');  // +4 for return type
    std::string funcindent2 = std::string(funcline2.length()+4, ' ');
    std::stringstream sscwork, ssig, ssig2;

    // comment out contwork if its not needed
    ssig  << "struct simint_multi_shellpair const P,\n"
          << funcindent << "struct simint_multi_shellpair const Q,\n"
          << funcindent << "double * const restrict contwork,\n"
          << funcindent << "double * const restrict " << ArrVarName(am) << ")";
    ssig2 << "struct simint_multi_shellpair const P,\n"
          << funcindent2 << "struct simint_multi_shellpair const Q,\n"
          << funcindent2 << "double * const restrict " << ArrVarName(am) << ")";


    // create the version that allocates contwork for the user
    os_ << "\n\n";
    os_ << "int " << funcline2 << ssig2.str() << "\n";
    os_ << "{\n";
    size_t contmem = info_.ContMemoryReq();
    if(contmem == 0)
        os_ << indent1 << "int ret = " << StringBuilder(funcline, "P, Q, NULL, ", ArrVarName(am), ");");
    else
    {
        size_t contnel = info_.ContNElements();

        os_ << indent1 << "// Workspace for contracted integrals\n";
        if(info_.UseHeap())
            os_ << indent1 << "double * const contwork = ALLOC(SIMINT_NSHELL_SIMD * " << contmem << ");\n\n";
        else
            os_ << indent1 << "double contwork[SIMINT_NSHELL_SIMD * " << contnel << "] SIMINT_ALIGN_ARRAY_DBL;\n\n";

        os_ << indent1 << "int ret = " << StringBuilder(funcline, "P, Q, contwork, ", ArrVarName(am), ");");

    }

    os_ << "\n\n"; 

    FreeContwork();
    os_ << indent1 << "return ret;\n";
    os_ << "}\n";


    // Add both versions to the header
    osh_ << "int " << funcline << ssig.str() << ";\n\n";
    osh_ << "int " << funcline2 << ssig2.str() << ";\n\n";
}


