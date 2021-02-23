#ifndef WIENERSVD_H
#define WIENERSVD_H

#include "Utility.h"

// Class to carry out the WienerSVD unfolding
// This is a reimplementation of the code in https://github.com/BNLIF/Wiener-SVD-Unfolding


class WienerSVD{

    public:
    // Default constructor
    WienerSVD(){};

    TH2D* smear;
    TH1D* wiener;
    TH2D* unfcov;
    TH1D* unf;
    TH1D* diff;
    TH1D* bias;
    TH1D* bias2;
    TH1D* fracError;
    TH1D* absError;
    TH1D* MSE;
    TH1D* MSE2;

    Utility _util;

    // -------------------------------------------------------------------------
    // Initialise the Class
    void Initialise(Utility _utility);
    // -------------------------------------------------------------------------
    /*
    Main function to run the Wiener Unfolding
    
    C_type: Options of  matrix for smoothness, 0 (unitary matrix), 1 (1st derivative matrix), 2 (2nd derivative matrix), else (3rd derivative matrix)
    sig: true signal spectrum or expected signal spectrum from model, used to construct wiener filter 
    mes: measured spectrum
    res: response matrix in histogram format
    cov: covariance matrix of measured spectrum in histrogram format

    */
    void DoUnfolding(Int_t C_type, Float_t Norm_type, TH1D *sig, TH1D *mes, TH2D* res, TH2D* cov);
    // -------------------------------------------------------------------------
    // Make a comparison of data to MC
    void CompareModel(TH1D *sig);
    // -------------------------------------------------------------------------
    // converters for histrogram <---> matrix/vector
    // former is input, latter is output

    // rowcolumn is TRUE, histo(i+1, j+1) = matrix(i, j)
    // rowcolumn is FALSE, histo(i+1, j+1) = matrix(j, i)
    // -------------------------------------------------------------------------
    // 2D histo to matrix
    void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn); 
    // -------------------------------------------------------------------------
    // 1D histo to vector
    void H2V(const TH1D* histo, TVectorD& vec);
    // -------------------------------------------------------------------------
    // Matrix to 2D histogram
    void M2H(const TMatrixD mat, TH2D* histo);
    // -------------------------------------------------------------------------
    // Vector to 1D histo
    void V2H(const TVectorD vec, TH1D* histo);
    // -------------------------------------------------------------------------
    // print out two matrices M1, M2, and their product M1*M2
    void MatrixMatrix(TMatrixD M1, TMatrixD M2);
    // -------------------------------------------------------------------------
    // print out matrix M and vector V, and their product M*V
    void MatrixVector(TMatrixD M, TVectorD V);
    // -------------------------------------------------------------------------
    // test/validation on SVD decomposition implemented by ROOT
    void SVD(TMatrixD M);
    // -------------------------------------------------------------------------
    // interactive initilization of a matrix/vector
    TMatrixD Matrix(Int_t row, Int_t column);
    TVectorD Vector(Int_t row);
    // -------------------------------------------------------------------------
    // construction of additional matrix for smoothness/improvement
    TMatrixD Matrix_C(Int_t n, Int_t type);
    // -------------------------------------------------------------------------
    // wiener filter unfolding
    // first five parameters are inputs as names read
    // AddSmear is additional smearing matrix after unfolding
    // WF is the elements of Wiener Filter
    // Covariance matrix of the unfolded spectrum
    TVectorD WienerSVDUnfold(TMatrixD Response, TVectorD Signal, TVectorD Measure, TMatrixD Covariance, Int_t C_type, Float_t Norm_type, TMatrixD& AddSmear, TVectorD& WF, TMatrixD& UnfoldCov);
    // -------------------------------------------------------------------------


}; // End Class WienerSVD

#endif