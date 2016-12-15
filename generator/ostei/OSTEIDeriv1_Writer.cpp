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
void OSTEIDeriv1_Writer::DeclareContwork(void) const
{
    if(info_.ContMemoryReq() == 0)
        return;

    os_ << indent1 << "// partition workspace into shells\n";
    size_t ptidx = 0;

    for(const auto & it : info_.GetContQ())
    {
        if(!info_.IsFinalAM(it.qam))
        {
            os_ << indent1 << "double * const " << ArrVarName(it) << " = contwork + (SIMINT_NSHELL_SIMD * " << ptidx << ");\n";
            ptidx += NCART(it.qam);
        }
    }

    os_ << "\n";
}


void OSTEIDeriv1_Writer::ZeroContwork(void) const
{
    size_t contmem = info_.ContMemoryReq();
    if(contmem > 0)
        os_ << indent3 << "memset(contwork, 0, SIMINT_NSHELL_SIMD * " << contmem << ");\n";
    
}


void OSTEIDeriv1_Writer::FreeContwork(void) const
{
    size_t contmem = info_.ContMemoryReq();

    if(contmem > 0 && info_.UseHeap())
    {
        os_ << indent1 << "// Free contracted workspace\n";
        os_ << indent1 << "SIMINT_FREE(contwork);\n\n";
    }       
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
        os_ << indent7 << PrimPtrName(it) << " += " << NCART(it.qam) << ";\n";

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


void OSTEIDeriv1_Writer::WriteShellOffsets_Scalar(void) const
{
    os_ << indent5 << "// Move pointers if this is the end of a shell\n";
    os_ << indent5 << "// Handle if the first element of the vector is a new shell\n";
    os_ << indent5 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os_ << indent5 << "{\n";
    os_ << indent6 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto it : info_.GetContQ())
        os_ << indent6 << PrimPtrName(it) << " += " << NCART(it.qam) << ";\n";

    os_ << indent5 << "}\n";
    os_ << indent5 << "iprimcd++;\n\n";
}


void OSTEIDeriv1_Writer::WriteAccumulation(void) const
{
    os_ << "\n\n";
    os_ << indent5 << "////////////////////////////////////\n";
    os_ << indent5 << "// Accumulate contracted integrals\n";
    os_ << indent5 << "////////////////////////////////////\n";

    if(info_.Vectorized())
    {
        os_ << indent5 << "if(lastoffset == 0)\n";
        os_ << indent5 << "{\n";

        for(const auto it : info_.GetContQ())
        {
            int ncart = NCART(it.qam);
            os_ << indent6 << "contract_all(" << ncart << ", " << PrimVarName(it.qam) << ", " << PrimPtrName(it.qam) << ");\n";
        }
        os_ << indent5 << "}\n";
        os_ << indent5 << "else\n";
        os_ << indent5 << "{\n";

        for(const auto it : info_.GetContQ())
        {
            int ncart = NCART(it.qam);
            os_ << indent6 << "contract(" << ncart << ", shelloffsets, " << PrimVarName(it) << ", " << PrimPtrName(it) << ");\n";
        }

        for(const auto it : info_.GetContQ())
            os_ << indent6 << PrimPtrName(it) << " += lastoffset*" << NCART(it.qam) << ";\n";

        os_ << indent5 << "}\n";
    }
    else
    {
        for(const auto it : info_.GetContQ())
        {
            int ncart = NCART(it.qam);
            os_ << indent6 << "contract(" << ncart << ", " << PrimVarName(it) << ", " << PrimPtrName(it.qam) << ");\n";
        }
    }
}


void OSTEIDeriv1_Writer::WriteFile_Permute_(void) const
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

    std::string funcline = StringBuilder("int ostei_sharedwork_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string indent(funcline.length(), ' ');

    // start output to the file
    // we only need this one include
    os_ << "#include \"simint/ostei/gen/ostei_generated.h\"\n";
    os_ << "\n";

    os_ << "\n\n";
    os_ << funcline;
    os_ << "struct simint_multi_shellpair const P,\n";
    os_ << indent << "struct simint_multi_shellpair const Q,\n";
    os_ << indent << "double screen_tol,\n";
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

 
    os_ << indent1 << "return ostei_sharedwork_" << amchar[tocall[0]] << "_" 
                                               << amchar[tocall[1]] << "_" 
                                               << amchar[tocall[2]] << "_" 
                                               << amchar[tocall[3]] << "(" << P_var << ", " << Q_var << ", screen_tol, "
                                                                    << "contwork, " << ArrVarName(am) << ");\n"; 

    os_ << "}\n";
    os_ << "\n";

}


void OSTEIDeriv1_Writer::WriteFile_NoPermute_(void) const
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
                        "\"simint/ostei/gen/ostei_generated.h\"",
                        "\"simint/vectorization/vectorization.h\"",
                        "\"simint/boys/boys.h\"",
                        "\"simint/ostei/ostei_contract.h\""};

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


    //////////////////////////////
    // Function name & signature
    //////////////////////////////
    std::string funcline = StringBuilder("int ostei_sharedwork_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string indent(funcline.length(), ' ');


    os_ << "\n\n";
    os_ << funcline;
    os_ << "struct simint_multi_shellpair const P,\n";
    os_ << indent << "struct simint_multi_shellpair const Q,\n";
    os_ << indent << "double screen_tol,\n";
    os_ << indent << "double * const restrict contwork,\n";
    os_ << indent << "double * const restrict " << ArrVarName(am) << ")\n";
    os_ << "{\n";
    os_ << "\n";


    ///////////////////////////////////
    // NOW IN THE ACTUAL OSTEI FUNCTION
    ///////////////////////////////////

    // If there is no HRR, integrals are accumulated from inside the primitive loop
    // directly into the final integral array that was passed into this function, so it must be zeroed first
    if(!hashrr)
        os_ << indent1 << "memset(" << ArrVarName(am) << ", 0, P.nshell12_clip * Q.nshell12_clip * " << ncart << " * sizeof(double));\n";
    os_ << "\n";


    // abcd = index within simd loop, 
    os_ << indent1 << "int ab, cd, abcd;\n";
    os_ << indent1 << "int istart, jstart;\n";
    os_ << indent1 << "int iprimcd, nprim_icd, icd;\n";
    os_ << indent1 << "const int check_screen = (screen_tol > 0.0);\n";
    os_ << indent1 << "int i, j;\n";
    os_ << indent1 << "int n;\n";

    if(info_.Vectorized())
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



    // Write out all the constants 
    // This is only needed if this is NOT scalar code
    if(info_.Vectorized())
    {
        os_ << indent1 << "// Create constants\n";
        for(const auto & it : cm)
            os_ << indent1 << "const SIMINT_DBLTYPE " << it.first << " = SIMINT_DBLSET1(" << it.second << ");\n";
    }


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
        ZeroContwork();
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

    // Note: vrr_writer handles the prim pointers for s_s_s_s
    //       so we always want to run this
    vrr_writer_.DeclarePrimPointers(os_);
    os_ << "\n";

    os_ << indent4 << "// Load these one per loop over i\n";
    os_ << indent4 << "const SIMINT_DBLTYPE P_alpha = SIMINT_DBLSET1(P.alpha[i]);\n";
    os_ << indent4 << "const SIMINT_DBLTYPE P_prefac = SIMINT_DBLSET1(P.prefac[i]);\n";
    os_ << indent4 << "const SIMINT_DBLTYPE Pxyz[3] = { SIMINT_DBLSET1(P.x[i]), SIMINT_DBLSET1(P.y[i]), SIMINT_DBLSET1(P.z[i]) };\n";
                             
    os_ << indent4 << "const SIMINT_DBLTYPE bra_screen_max = SIMINT_DBLSET1(P.screen[i]);\n";
    os_ << "\n";

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

    if(info_.Vectorized())
        WriteShellOffsets();
    else
        WriteShellOffsets_Scalar();


    if(info_.Vectorized())
    {
        os_ << indent5 << "// Do we have to compute this vector (or has it been screened out)?\n";
        os_ << indent5 << "// (not_screened != 0 means we have to do this vector)\n";
        os_ << indent5 << "if(check_screen)\n";
        os_ << indent5 << "{\n";
        os_ << indent6 << "const SIMINT_DBLTYPE screen_max = SIMINT_MUL(bra_screen_max, SIMINT_DBLLOAD(Q.screen, j));\n";
        os_ << indent6 << "SIMINT_UNIONTYPE screen_max_v = { screen_max };\n";
        os_ << indent6 << "not_screened = 0;\n";
        os_ << indent6 << "for(n = 0; n < SIMINT_SIMD_LEN; n++)\n";
        os_ << indent7 << "not_screened = ( screen_max_v.d[n] >= screen_tol ? 1 : not_screened );\n";

        os_ << indent6 << "if(not_screened == 0)\n";
        os_ << indent6 << "{\n";
        for(const auto it : info_.GetContQ())
            os_ << indent7 << PrimPtrName(it) << " += lastoffset*" << NCART(it.qam) << ";\n";
        os_ << indent7 << "continue;\n";
        os_ << indent6 << "}\n";
        os_ << indent5 << "}\n\n";
    }
    else
    {
        os_ << indent5 << "// Skip if screened out\n";
        os_ << indent5 << "if(check_screen && ( (bra_screen_max * Q.screen[j]) < screen_tol ) )\n";
        os_ << indent6 << "continue;\n\n";
    }

    // Note: vrr_writer handles the prim arrays for s_s_s_s
    //       so we always want to run this
    vrr_writer_.DeclarePrimArrays(os_);

    os_ << indent5 << "const SIMINT_DBLTYPE Q_alpha = SIMINT_DBLLOAD(Q.alpha, j);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE PQalpha_mul = SIMINT_MUL(P_alpha, Q_alpha);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE PQalpha_sum = SIMINT_ADD(P_alpha, Q_alpha);\n";
    os_ << indent5 << "const SIMINT_DBLTYPE one_over_PQalpha_sum = SIMINT_DIV(const_1, PQalpha_sum);\n";
    os_ << "\n";
    os_ << "\n";
    os_ << indent5 << "/* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os_ << indent5 << "SIMINT_DBLTYPE PQ[3];\n";
    os_ << indent5 << "PQ[0] = SIMINT_SUB(Pxyz[0], SIMINT_DBLLOAD(Q.x, j));\n";
    os_ << indent5 << "PQ[1] = SIMINT_SUB(Pxyz[1], SIMINT_DBLLOAD(Q.y, j));\n";
    os_ << indent5 << "PQ[2] = SIMINT_SUB(Pxyz[2], SIMINT_DBLLOAD(Q.z, j));\n";


    os_ << indent5 << "SIMINT_DBLTYPE R2 = SIMINT_MUL(PQ[0], PQ[0]);\n";
    os_ << indent5 << " R2 = SIMINT_FMADD(PQ[1], PQ[1], R2);\n";
    os_ << indent5 << " R2 = SIMINT_FMADD(PQ[2], PQ[2], R2);\n";
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
    os_ << indent5 << "// Maximum v value: " << info_.L() << "\n";
    os_ << indent5 << "//////////////////////////////////////////////\n";
    os_ << indent5 << "// The parameter to the Fjt function\n";
    os_ << indent5 << "const SIMINT_DBLTYPE F_x = SIMINT_MUL(R2, alpha);\n";
    os_ << "\n";
    os_ << "\n";

    // we need to zero out any that are beyond the end of the batch (that's been clipped)
    if(info_.Vectorized())
    {
        os_ << indent5 << "SIMINT_UNIONTYPE Q_prefac_u = { SIMINT_DBLLOAD(Q.prefac, j) };\n";
        os_ << indent5 << "for(n = nlane; n < SIMINT_SIMD_LEN; n++)\n";
        os_ << indent6 << "Q_prefac_u.d[n] = 0.0;\n";
        os_ << indent5 << "const SIMINT_DBLTYPE Q_prefac = SIMINT_UNIONMEMBER(Q_prefac_u);\n";
    }
    else
        os_ << indent5 << "const double Q_prefac = Q.prefac[j];\n";
    os_ << "\n\n"; 


    os_ << indent5 << "boys_F_split(" << PrimVarName({0,0,0,0})
                   << ", &F_x, " << info_.L() << ");\n";


    // prefac = sqrt(1/PQalpha_sum) * P_prefac * Q_prefac
    os_ << indent5 << "SIMINT_DBLTYPE prefac = SIMINT_SQRT(one_over_PQalpha_sum);\n";
    os_ << indent5 << "prefac = SIMINT_MUL(SIMINT_MUL(P_prefac, Q_prefac), prefac);\n";

    const std::string name0000 = PrimVarName({0,0,0,0});
    const std::string name0000n = name0000 + "[n]";

    os_ << indent5 << "for(n = 0; n <= " << info_.L() << "; n++)\n"
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

    if(hrr_writer_.HasHRR())
        hrr_writer_.WriteHRR(os_);

    os_ << indent2 << "}   // close loop cdbatch\n";

    os_ << "\n";
    os_ << indent2 << "istart = iend;\n";

    os_ << indent1 << "}  // close loop over ab\n";
    os_ << "\n";
    os_ << "\n";

    os_ << "\n";


    os_ << indent1 << "return P.nshell12_clip * Q.nshell12_clip;\n";
    os_ << "}\n";
    os_ << "\n";
}



void OSTEIDeriv1_Writer::WriteFile(void) const
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
    std::string funcline = StringBuilder("ostei_sharedwork_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string funcline2 = StringBuilder("ostei_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string funcindent = std::string(funcline.length()+4, ' ');  // +4 for return type
    std::string funcindent2 = std::string(funcline2.length()+4, ' ');
    std::stringstream sscwork, ssig, ssig2;

    // comment out contwork if its not needed
    ssig  << "struct simint_multi_shellpair const P,\n"
          << funcindent << "struct simint_multi_shellpair const Q,\n"
          << funcindent << "double screen_tol,\n"
          << funcindent << "double * const restrict contwork,\n"
          << funcindent << "double * const restrict " << ArrVarName(am) << ")";
    ssig2 << "struct simint_multi_shellpair const P,\n"
          << funcindent2 << "struct simint_multi_shellpair const Q,\n"
          << funcindent2 << "double screen_tol,\n"
          << funcindent2 << "double * const restrict " << ArrVarName(am) << ")";


    // create the version that allocates contwork for the user
    os_ << "\n\n";
    os_ << "int " << funcline2 << ssig2.str() << "\n";
    os_ << "{\n";
    size_t contmem = info_.ContMemoryReq();
    if(contmem == 0)
        os_ << indent1 << "int ret = " << StringBuilder(funcline, "P, Q, screen_tol, NULL, ", ArrVarName(am), ");");
    else
    {
        size_t contnel = info_.ContNElements();

        os_ << indent1 << "// Workspace for contracted integrals\n";
        if(info_.UseHeap())
            os_ << indent1 << "double * const contwork = SIMINT_ALLOC(SIMINT_NSHELL_SIMD * " << contmem << ");\n\n";
        else
            os_ << indent1 << "double contwork[SIMINT_NSHELL_SIMD * " << contnel << "] SIMINT_ALIGN_ARRAY_DBL;\n\n";

        os_ << indent1 << "int ret = " << StringBuilder(funcline, "P, Q, screen_tol, contwork, ", ArrVarName(am), ");");

    }

    os_ << "\n\n"; 

    FreeContwork();
    os_ << indent1 << "return ret;\n";
    os_ << "}\n";


    // Add both versions to the header
    osh_ << "int " << funcline << ssig.str() << ";\n\n";
    osh_ << "int " << funcline2 << ssig2.str() << ";\n\n";
}


