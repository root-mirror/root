/// \file
/// \ingroup tutorial_matrix
/// \notebook -nodraw
/// This tutorial shows how to decompose a matrix A in an orthogonal matrix Q and an upper
/// triangular matrix R using QR Householder decomposition with the TDecompQRH class.
/// The matrix same matrix as in this example is used: https://en.wikipedia.org/wiki/QR_decomposition#Example_2

#include <iostream>
#include "TMath.h"
#include "TDecompQRH.h"

void decomposeQR() {

   const int n = 3;


   double a[] = {12, -51, 4, 6, 167, -68, -4, 24, -41};

   TMatrixT<double> A(3, 3, a);

   std::cout << "initial matrox A " << std::endl;

   A.Print();

   TDecompQRH decomp(A);

   bool ret = decomp.Decompose();

   std::cout << "Orthogonal Q matrix " << std::endl;

   // note that decomp.GetQ()  returns an intenrnal matrix which is not Q defined as A = QR
   auto Q = decomp.GetOrthogonalMatrix();
   Q.Print();

   std::cout << "Upper Triangular R matrix " << std::endl;
   auto R = decomp.GetR();

   R.Print();

   // check that we have a correct Q-R decomposition

   TMatrixT<double> compA = Q * R;

   std::cout << "Computed A matrix from Q * R " << std::endl;
   compA.Print();

   for (int i = 0; i < A.GetNrows(); ++i) {
      for (int j = 0; j < A.GetNcols(); ++j) {
         if (!TMath::AreEqualAbs( compA(i,j), A(i,j), 1.E-6) )
            Error("decomposeQR","Reconstrate matrix is not equal to the original : %f different than %f",compA(i,j),A(i,j));
      }
   }

   // chech also that Q is orthogonal (Q^T * Q = I)
   auto QT = Q;
   QT.Transpose(Q);

   auto qtq = QT * Q;
   for (int i = 0; i < Q.GetNrows(); ++i) {
      for (int j = 0; j < Q.GetNcols(); ++j) {
         if ((i == j && !TMath::AreEqualAbs(qtq(i, i), 1., 1.E-6)) ||
             (i != j && !TMath::AreEqualAbs(qtq(i, j), 0., 1.E-6))) {
            Error("decomposeQR", "Q matrix is not orthogonal ");
            qtq.Print();
            break;
         }
      }
   }
}
