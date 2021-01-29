#include "../include/WienerSVD.h"

// -----------------------------------------------------------------------------
void WienerSVD::DoUnfolding(Int_t C_type, Float_t Norm_type, TH1D *sig, TH1D *mes, TH2D* res, TH2D* cov)
{
    
    std::cout << "Derivative matrix type: "<< C_type <<std::endl;
    std::cout << "Signal normalization type (power exponent): "<< Norm_type <<std::endl;
    
    Int_t n = sig->GetNbinsX();
    Double_t Nuedges[n+1];
    
    for(int i=0; i<n+1; i++){
        Nuedges[i] = sig->GetBinLowEdge(i+1);
    }

    Int_t m = mes->GetNbinsX();
 
    // Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input
    TVectorD signal(n);
    TVectorD measure(m);
    TMatrixD response(m, n);
    TMatrixD covariance(m, m);
   
    // Convert input into mathematical formats, easy and clean to be processed. 
    // Converted defined/implemented in source files, see include/Util.h
    H2V(sig, signal);
    H2V(mes, measure);
    H2M(res, response, kTRUE);
    H2M(cov, covariance, kTRUE);

    // Construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    TMatrixD AddSmear(n,n);
    TVectorD WF(n);
    TMatrixD UnfoldCov(n,n);
    
    TH2D* smear     = new TH2D("smear",     "Additional Smearing Matirx",n,Nuedges,n,Nuedges);
    TH1D* wiener    = new TH1D("wiener",    "Wiener Filter Vector",n,0,n);
    TH2D* unfcov    = new TH2D("unfcov",    "Unfolded spectrum covariance", n, Nuedges, n, Nuedges);
    TH1D* unf       = new TH1D("unf",       "unfolded spectrum",n,Nuedges);
    TH1D* diff      = new TH1D("diff",      "Fractional difference of unf and signal model",n, Nuedges);
    TH1D* bias      = new TH1D("bias",      "intrinsic bias w.r.t. model",n, Nuedges);
    TH1D* bias2     = new TH1D("bias2",     "intrinsic bias2 w.r.t. unfolded result",n, Nuedges);
    TH1D* fracError = new TH1D("fracError", "Fractional uncertainty", n, Nuedges);
    TH1D* absError  = new TH1D("absError",  "Absolute uncertainty", n, Nuedges);
    TH1D* MSE       = new TH1D("MSE",       "Mean Square Error: variance+bias^2", n, Nuedges);
    TH1D* MSE2      = new TH1D("MSE2",      "Mean Square Error: variance", n, Nuedges);


    // Core implementation of Wiener-SVD
    // Input as names read. AddSmear and WF to record the core information in the unfolding.
    TVectorD unfold = WienerSVDUnfold(response, signal, measure, covariance, C_type, Norm_type, AddSmear, WF, UnfoldCov);

    // Output and comparison between true/expected signal spectrum with unfolded one
    V2H(unfold, unf);

    for (Int_t i=1; i<=n; i++) {
        Double_t s2s = 0;
        Double_t u = unf->GetBinContent(i);
        Double_t t = signal(i-1);
        if(t!=0) s2s = u/t - 1;
        else s2s = 1.; // attention to this t=0
        diff->SetBinContent(i, s2s); // in percentage 
    }
   
    // Intrinsic bias (Ac-I)*s_bar formula
    TMatrixD unit(n,n);
    unit.UnitMatrix();
    TVectorD intrinsicbias = (AddSmear - unit)*signal;
    TVectorD intrinsicbias2 = (AddSmear - unit)*unfold;
    
    for (int i=0; i<n; i++) {
        
        if(signal(i)!=0)
            intrinsicbias(i) = intrinsicbias(i)/signal(i);
        else
            intrinsicbias(i)=0.;
        
        if (unfold(i)!=0)
            intrinsicbias2(i) = intrinsicbias2(i)/unfold(i);
        else
            intrinsicbias2(i)=0.;
    }

    V2H(intrinsicbias, bias);
    V2H(intrinsicbias2, bias2);
    
    // Diagonal uncertainty
    for (int i=1; i<=n; i++) {
        fracError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1))/unfold(i-1));
        absError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1)));
    }
    
    // MSE
    for (int i=0; i<n; i++) {
        MSE ->SetBinContent(i+1, TMath::Power(intrinsicbias2(i)*unfold(i),2)+UnfoldCov(i,i));
        MSE2->SetBinContent(i+1, UnfoldCov(i,i));
    }

    // Convert matrix/vector to histogram and save
    M2H(AddSmear, smear);
    V2H(WF, wiener);
    M2H(UnfoldCov, unfcov);

    return;
}

// -----------------------------------------------------------------------------
TMatrixD WienerSVD::Matrix(Int_t row, Int_t column) {
    TMatrixD MMM(row, column);
    for(Int_t i=0; i<row; i++)
    {
        for(Int_t j=0; j<column; j++)
        {
            Double_t eee;
            std::cout<<"Row:"<<i+1<<" Column:"<<j+1<<std::endl;
            std::cin>>eee;
            MMM(i, j) = eee;
        }
    }
    MMM.Print();
    return MMM;
}
// -----------------------------------------------------------------------------
void WienerSVD::MatrixMatrix(TMatrixD M1, TMatrixD M2)
{
    std::cout<<"Matrix 1 ========="<<std::endl;
    M1.Print();
    std::cout<<"Matrix 2 ========="<<std::endl;
    M2.Print();
    TMatrixD M = M1*M2;
    std::cout<<"Matrix 1*2 ========="<<std::endl;
    M.Print();
}
// -----------------------------------------------------------------------------
TVectorD WienerSVD::Vector(Int_t row)
{
    TVectorD V(row);
    for(Int_t i=0; i<row; i++)
    {
        Double_t eee;
        std::cout<<"Row "<<i+1<<std::endl;
        std::cin>>eee;
        V(i) = eee;
    }
    V.Print();
    return V;
}
// -----------------------------------------------------------------------------
void WienerSVD::MatrixVector(TMatrixD M, TVectorD V)
{
    std::cout<<"Matrix 1 ========="<<std::endl;
    M.Print();
    std::cout<<"Vector 1 ========="<<std::endl;
    V.Print();
    std::cout<<"Matrix * Vector =========="<<std::endl;
    TVectorD MM = M*V;
    MM.Print();
}
// -----------------------------------------------------------------------------
void WienerSVD::SVD(TMatrixD M)
{
    Int_t nrow = M.GetNrows();
    Int_t ncol = M.GetNcols();

    TDecompSVD svd(M);
    TMatrixD U = svd.GetU();
    TMatrixD V (TMatrixD::kTransposed, svd.GetV());
    TVectorD S = svd.GetSig();

    std::cout<<"U ========="<<std::endl;
    U.Print();
    TMatrixD uu(U, TMatrixD::kMultTranspose, U);
    uu.Print();
    std::cout<<"V ========="<<std::endl;
    V.Print();
    TMatrixD vv(V, TMatrixD::kMultTranspose, V);
    vv.Print();
    std::cout<<"D ========="<<std::endl;
    TMatrixD D(nrow, ncol);
    for(Int_t i=0; i<nrow; i++)
    {
        for(Int_t j=0; j<ncol; j++)
        {
            D(i, j) = 0;
            if(i == j)
            {
                D(i, j) = S(i);
            }
        }
    }
    D.Print();

    std::cout<<"Orignal ========="<<std::endl;
    M.Print();
    std::cout<<"Reassembly ========"<<std::endl;
    TMatrixD MM = U*D*V;
    MM.Print();
    MM.Draw("colz");
}
// -----------------------------------------------------------------------------
void WienerSVD::H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn)
{
    // Fill 2D histogram into matrix
    // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        for(Int_t j=0; j<histo->GetNbinsY(); j++)
        {
            if(rowcolumn) mat(i, j) = histo->GetBinContent(i+1, j+1);
            else mat(j, i) = histo->GetBinContent(i+1, j+1);
        }
    }
}
// -----------------------------------------------------------------------------
void WienerSVD::H2V(const TH1D* histo, TVectorD& vec)
{
    // Fill 1D histogram into matrix
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        vec(i) = histo->GetBinContent(i+1);
    }
}
// -----------------------------------------------------------------------------
void WienerSVD::M2H(const TMatrixD mat, TH2D* histo)
{
    // Fill matrix to histogram
    for(Int_t i=0; i<mat.GetNrows(); i++)
    {
        for(Int_t j=0; j<mat.GetNcols(); j++)
        {
            histo->SetBinContent(i+1, j+1, mat(i, j));
        }
    }
}
// -----------------------------------------------------------------------------
void WienerSVD::V2H(const TVectorD vec, TH1D* histo)
{
    // Fill vector to histogram,
    for(Int_t i=0; i<vec.GetNrows(); i++)
    {
        histo->SetBinContent(i+1, vec(i));
    }
}
// -----------------------------------------------------------------------------
TMatrixD WienerSVD::Matrix_C(Int_t n, Int_t type)
{
    //Type 0: Unit matrix
    //Type I: First derivative matrix
    //Type II: Second derivative matrix
    Double_t epsilon = 1e-6; // needed for 2nd derivative matrix inversion
    Double_t epsilon2 = 1e-2; // needed for 3rd derivative matrix inversion
    TMatrixD C(n, n);
    for(Int_t i=0; i<n; i++)
    {
        for(Int_t j=0; j<n; j++)
        {
            C(i, j) = 0;
            if(type == 0)
            {
                if(i == j) C(i, j) = 1;
            }

            else if(type == 1)
            {
                if( j-i == 1 ) C(i, j) = 1;
                if(i==j) 
                {
                    C(i, j) = -1;
                }
            }

            else if(type == 2)
            {
                if(TMath::Abs(i-j) == 1) C(i, j) = 1;
                if(i==j) 
                {
                    C(i, j) = -2+epsilon;
                    if(i==0 || i==n-1)
                    {
                        C(i, j) = -1+epsilon;
                    }
                }
            }

            else // third derivative matrix
            {
                if(i-j == 2) 
                { 
                    C(i, j) = -1;
                }
                if(i-j == 1) 
                {
                    C(i, j) = 2;
                    if(i==1)
                    {
                        C(i, j) = 0; 
                    }
                }
                if(i-j == -1) 
                {
                    C(i, j) = -2;
                    if(i==n-2)
                    {
                        C(i, j) = 0;
                    }
                }
                if(i-j == -2) 
                {
                    C(i, j) = 1;
                }
                if(i==j) 
                {
                    C(i, j) = 0 + epsilon2;
                    if(i==0 || i==1)
                    {
                        C(i, j) = 1 + epsilon2;
                    }
                    if(i==n-1 || i==n-2)
                    {
                        C(i, j) = -1 + epsilon2;
                    }
                }
            }
        }
    }

    return C;
}
// -----------------------------------------------------------------------------
TVectorD WienerSVD::WienerSVDUnfold(TMatrixD Response, TVectorD Signal, TVectorD Measure, TMatrixD Covariance, Int_t C_type, Float_t Norm_type, TMatrixD& AddSmear, TVectorD& WF, TMatrixD& UnfoldCov)
{
    Int_t m = Response.GetNrows(); // measure, M
    Int_t n = Response.GetNcols(); // signal, S
    // Decomposition of Covariance Matrix to get orthogonal Q to rotate the current frame, 
    // then make the uncertainty for each bin equal to 1
    TDecompSVD decV(Covariance);
    TMatrixD Q0 (TMatrixD::kTransposed, decV.GetV());
    TVectorD err0 = decV.GetSig();

    TMatrixD err(m, m);
    for(Int_t i=0; i<m; i++)
    {
        for(Int_t j=0; j<m; j++)
        {
            err(i, j) = 0;
            if(i == j)
            {
                if(err0(i)) err(i, j) = 1./TMath::Sqrt( err0(i) );
                else err(i, j) = 0;
            }
        }
    }

    TMatrixD Q = err*Q0;
    // transform Measure and Response
    TVectorD M = Q*Measure;
    TMatrixD R = Q*Response;

    // For addtion of smoothness matrix, e.g. 2nd derivative C
    TMatrixD C0(n, n);
    C0 = Matrix_C(n, C_type);
    TMatrixD normsig(n, n);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            normsig(i, j) = 0;
            if(i==j) normsig(i, j) = 1./TMath::Power(Signal(i), Norm_type);
        }
    }
    C0 = C0*normsig;

    TMatrixD C = C0;
    C0.Invert();
    TMatrixD C_inv = C0;
    Signal = C*Signal;
    R = R*C_inv;
  
    // SVD decomposition of R 
    TDecompSVD udv(R);
    TMatrixD U = udv.GetU();
    TMatrixD U_t (TMatrixD::kTransposed, U);
    TMatrixD V = udv.GetV();
    TMatrixD V_t (TMatrixD::kTransposed, V);
    TVectorD D = udv.GetSig();
    // for matrix D transverse
    TMatrixD D_t(n, m);
    for(Int_t i=0; i<n; i++)
    {
        for(Int_t j=0; j<m; j++)
        {
            D_t(i,j) = 0;
            if(i==j)
            {
                D_t(i, j) = D(i);
            }
        }
    }

    TVectorD S = V_t*Signal;
    // Wiener Filter 
    TMatrixD W(n, n);
    TMatrixD W0(n, n);
    for(Int_t i=0; i<n; i++)
    {
        for(Int_t j=0; j<n; j++)
        {
            W(i, j) = 0;
            W0(i, j) = 0;
            if(i == j)
            {
                //W(i, j) = 1./(D(i)*D(i)+2e-7); //S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                //WF(i) = D(i)*D(i)*W(i, j);//S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                W(i, j) = S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                WF(i) = D(i)*D(i)*W(i, j);
                W0(i, j) = WF(i); 
            }
        }
    }

    TVectorD unfold = C_inv*V*W*D_t*U_t*M;
    AddSmear = C_inv*V*W0*V_t*C;

    // covariance matrix of the unfolded spectrum
    TMatrixD covRotation = C_inv*V*W*D_t*U_t*Q;
    TMatrixD covRotation_t (TMatrixD::kTransposed, covRotation); 
    UnfoldCov = covRotation*Covariance*covRotation_t;  

    return unfold;
}
// -----------------------------------------------------------------------------