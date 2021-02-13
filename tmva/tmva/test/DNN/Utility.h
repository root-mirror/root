#ifndef TMVA_TEST_DNN_UTILITY
#define TMVA_TEST_DNN_UTILITY

#include <cassert>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <cstdlib>
#include <vector>
#include <string>

#include "TRandom.h"
//#include "TMVA/DNN/Architectures/Reference.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/Net.h"
#include "TMVA/DNN/DeepNet.h"
#include "TMVA/DNN/CNN/ConvLayer.h"
#include "TMVA/DNN/CNN/MaxPoolLayer.h"
#include "TMVA/DNN/DenseLayer.h"

namespace TMVA {
namespace DNN {

/** Construct a convolutional neural network with one convolutional layer,
 *  one pooling layer and two fully connected layers. The dimensions are
 *  predetermined. The activation functions are chosen randomly. */
//______________________________________________________________________________
template <typename AArchitecture>
void constructConvNet(TDeepNet<AArchitecture> &net)
{
   /* For random selection */
   std::vector<EActivationFunction> ActivationFunctions = {EActivationFunction::kIdentity, EActivationFunction::kRelu,
                                                           EActivationFunction::kSigmoid, EActivationFunction::kTanh};

   size_t depth1 = 10;
   size_t filterHeightConv1 = 3;
   size_t filterWidthConv1 = 3;
   size_t strideRowsConv1 = 1;
   size_t strideColsConv1 = 1;
   size_t zeroPaddingHeight1 = 1;
   size_t zeroPaddingWidth1 = 1;

   //EActivationFunction fConv1 = EActivationFunction::kIdentity;

   //TRandom &  r = (AArchitecture::GetRandomGenerator());
   //r.SetSeed(0);

   EActivationFunction fConv1 =  EActivationFunction::kRelu;
   //EActivationFunction fConv1 = ActivationFunctions[r.Uniform(ActivationFunctions.size())];

   net.AddConvLayer(depth1, filterHeightConv1, filterWidthConv1, strideRowsConv1, strideColsConv1, zeroPaddingHeight1,
                    zeroPaddingWidth1, fConv1);

   std::cout << "added Conv layer " <<  net.GetLayerAt(net.GetDepth() - 1)->GetDepth() << " x " <<  net.GetLayerAt(net.GetDepth() - 1)->GetHeight()
             << " x " << net.GetLayerAt(net.GetDepth() - 1)->GetWidth() << std::endl;


   size_t depth2 = 10;
   size_t filterHeightConv2 = 3;
   size_t filterWidthConv2 = 3;
   size_t strideRowsConv2 = 1;
   size_t strideColsConv2 = 1;
   size_t zeroPaddingHeight2 = 1;
   size_t zeroPaddingWidth2 = 1;

   EActivationFunction fConv2 = EActivationFunction::kSigmoid;
   //EActivationFunction fConv2 = ActivationFunctions[r.Uniform(ActivationFunctions.size())];

   net.AddConvLayer(depth2, filterHeightConv2, filterWidthConv2, strideRowsConv2, strideColsConv2, zeroPaddingHeight2,
                    zeroPaddingWidth2, fConv2);

   std::cout << "added Conv layer " <<  net.GetLayerAt(net.GetDepth() - 1)->GetDepth() << " x " <<  net.GetLayerAt(net.GetDepth() - 1)->GetHeight()
             << " x " << net.GetLayerAt(net.GetDepth() - 1)->GetWidth() << std::endl;


   size_t filterHeightPool = 2;
   size_t filterWidthPool = 2;
   size_t strideRowsPool = 1;
   size_t strideColsPool = 1;


   net.AddMaxPoolLayer(filterHeightPool, filterWidthPool, strideRowsPool, strideColsPool);

   std::cout << "added MaxPool layer " <<  net.GetLayerAt(net.GetDepth() - 1)->GetDepth() << " x " <<  net.GetLayerAt(net.GetDepth() - 1)->GetHeight()
             << " x " << net.GetLayerAt(net.GetDepth() - 1)->GetWidth() << std::endl;


   size_t depthReshape = 1;
   size_t heightReshape = 1;
   size_t widthReshape = net.GetLayerAt(net.GetDepth() - 1)->GetDepth() *
                         net.GetLayerAt(net.GetDepth() - 1)->GetHeight() *
                         net.GetLayerAt(net.GetDepth() - 1)->GetWidth();

   net.AddReshapeLayer(depthReshape, heightReshape, widthReshape, true);

   size_t widthFC1 = 20;

   //EActivationFunction fFC1 = EActivationFunction::kIdentity;
   EActivationFunction fFC1 = EActivationFunction::kSigmoid;

   //EActivationFunction fFC1 = ActivationFunctions[r.Uniform(ActivationFunctions.size())];
   net.AddDenseLayer(widthFC1, fFC1);

   size_t widthFC2 = 2;
   EActivationFunction fFC2 = EActivationFunction::kIdentity;
   net.AddDenseLayer(widthFC2, fFC2);
}

/** Construct a linear convolutional neural network with one convolutional layer,
 *  one pooling layer and two fully connected layers. The dimensions are
 *  predetermined. The activation functions are all linear.  */
//______________________________________________________________________________
template <typename AArchitecture>
void constructLinearConvNet(TDeepNet<AArchitecture> &net)
{

   size_t depth1 = 2;
   size_t filterHeightConv1 = 3;
   size_t filterWidthConv1 = 3;
   size_t strideRowsConv1 = 1;
   size_t strideColsConv1 = 1;
   size_t zeroPaddingHeight1 = 1;
   size_t zeroPaddingWidth1 = 1;

   EActivationFunction fConv1 = EActivationFunction::kIdentity;

   net.AddConvLayer(depth1, filterHeightConv1, filterWidthConv1, strideRowsConv1, strideColsConv1, zeroPaddingHeight1,
                    zeroPaddingWidth1, fConv1);

   std::cout << "added Conv layer " <<  net.GetLayerAt(net.GetDepth() - 1)->GetDepth() << " x " <<  net.GetLayerAt(net.GetDepth() - 1)->GetHeight()
             << " x " << net.GetLayerAt(net.GetDepth() - 1)->GetWidth() << std::endl;



   size_t depth2 = 2;
   size_t filterHeightConv2 = 2;
   size_t filterWidthConv2 = 2;
   size_t strideRowsConv2 = 1;
   size_t strideColsConv2 = 1;
   size_t zeroPaddingHeight2 = 0;
   size_t zeroPaddingWidth2 = 0;

   EActivationFunction fConv2 = EActivationFunction::kIdentity;

   net.AddConvLayer(depth2, filterHeightConv2, filterWidthConv2, strideRowsConv2, strideColsConv2, zeroPaddingHeight2,
                    zeroPaddingWidth2, fConv2);

   std::cout << "added Conv layer " <<  net.GetLayerAt(net.GetDepth() - 1)->GetDepth() << " x " <<  net.GetLayerAt(net.GetDepth() - 1)->GetHeight()
             << " x " << net.GetLayerAt(net.GetDepth() - 1)->GetWidth() << std::endl;

   // size_t depth3 = 12;
   // size_t filterHeightConv3 = 3;
   // size_t filterWidthConv3 = 3;
   // size_t strideRowsConv3 = 1;
   // size_t strideColsConv3 = 1;
   // size_t zeroPaddingHeight3 = 1;
   // size_t zeroPaddingWidth3 = 1;

   // EActivationFunction fConv3 = EActivationFunction::kIdentity;

   // net.AddConvLayer(depth3, filterHeightConv3, filterWidthConv3, strideRowsConv3, strideColsConv3, zeroPaddingHeight3,
   //                 zeroPaddingWidth3, fConv3);


   size_t filterHeightPool = 2;
   size_t filterWidthPool = 2;
   size_t strideRowsPool = 1;
   size_t strideColsPool = 1;


   net.AddMaxPoolLayer(filterHeightPool, filterWidthPool, strideRowsPool, strideColsPool);

   std::cout << "added MaxPool layer " <<  net.GetLayerAt(net.GetDepth() - 1)->GetDepth() << " x " <<  net.GetLayerAt(net.GetDepth() - 1)->GetHeight()
             << " x " << net.GetLayerAt(net.GetDepth() - 1)->GetWidth() << std::endl;


   size_t depthReshape = 1;
   size_t heightReshape = 1;
   size_t widthReshape = net.GetLayerAt(net.GetDepth() - 1)->GetDepth() *
                         net.GetLayerAt(net.GetDepth() - 1)->GetHeight() *
                         net.GetLayerAt(net.GetDepth() - 1)->GetWidth();

   net.AddReshapeLayer(depthReshape, heightReshape, widthReshape, true);

   size_t widthFC1 = 3;
   EActivationFunction fFC1 = EActivationFunction::kIdentity;
   net.AddDenseLayer(widthFC1, fFC1);

   size_t widthFC2 = 1;
   EActivationFunction fFC2 = EActivationFunction::kIdentity;
   net.AddDenseLayer(widthFC2, fFC2);
}

//______________________________________________________________________________
template <typename AArchitecture>
void constructMasterSlaveConvNets(TDeepNet<AArchitecture> &master, std::vector<TDeepNet<AArchitecture>> &nets)
{
   /* For random selection */
   std::vector<EActivationFunction> ActivationFunctions = {EActivationFunction::kIdentity, EActivationFunction::kRelu,
                                                           EActivationFunction::kSigmoid, EActivationFunction::kTanh};

   // Add Convolutional Layer
   size_t depth = 12;
   size_t filterHeightConv = 2;
   size_t filterWidthConv = 2;
   size_t strideRowsConv = 1;
   size_t strideColsConv = 1;
   size_t zeroPaddingHeight = 1;
   size_t zeroPaddingWidth = 1;

   EActivationFunction fConv = ActivationFunctions[rand() % ActivationFunctions.size()];

   TConvLayer<AArchitecture> *convLayer =
      master.AddConvLayer(depth, filterHeightConv, filterWidthConv, strideRowsConv, strideColsConv, zeroPaddingHeight,
                          zeroPaddingWidth, fConv);

   convLayer->Initialize();
   TConvLayer<AArchitecture> *copyConvLayer = new TConvLayer<AArchitecture>(*convLayer);

   // add the copy to all slave nets
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].AddConvLayer(copyConvLayer);
   }

   // Add Max Pooling Layer
   size_t filterHeightPool = 6;
   size_t filterWidthPool = 6;
   size_t strideRowsPool = 1;
   size_t strideColsPool = 1;

   // Add the Max pooling layer
   TMaxPoolLayer<AArchitecture> *maxPoolLayer =
      master.AddMaxPoolLayer(filterHeightPool, filterWidthPool, strideRowsPool, strideColsPool);
   TMaxPoolLayer<AArchitecture> *copyMaxPoolLayer = new TMaxPoolLayer<AArchitecture>(*maxPoolLayer);

   // Add the copy to all slave nets
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].AddMaxPoolLayer(copyMaxPoolLayer);
   }

   // Add the reshape layer
   size_t depthReshape = 1;
   size_t heightReshape = 1;
   size_t widthReshape = master.GetLayerAt(master.GetDepth() - 1)->GetDepth() *
                         master.GetLayerAt(master.GetDepth() - 1)->GetHeight() *
                         master.GetLayerAt(master.GetDepth() - 1)->GetWidth();

   TReshapeLayer<AArchitecture> *reshapeLayer = master.AddReshapeLayer(depthReshape, heightReshape, widthReshape, true);
   TReshapeLayer<AArchitecture> *copyReshapeLayer = new TReshapeLayer<AArchitecture>(*reshapeLayer);

   // Add the copy to all slave nets
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].AddReshapeLayer(copyReshapeLayer);
   }

   // Add Dense Layer
   size_t widthFC1 = 20;
   EActivationFunction fFC1 = ActivationFunctions[rand() % ActivationFunctions.size()];
   TDenseLayer<AArchitecture> *denseLayer = master.AddDenseLayer(widthFC1, fFC1);
   denseLayer->Initialize();
   TDenseLayer<AArchitecture> *copyDenseLayer = new TDenseLayer<AArchitecture>(*denseLayer);

   // add the copy to all slave nets
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].AddDenseLayer(copyDenseLayer);
   }

   // Add the final Dense Layer
   size_t widthFC2 = 5;
   EActivationFunction fFC2 = EActivationFunction::kIdentity;
   TDenseLayer<AArchitecture> *finalDenseLayer = master.AddDenseLayer(widthFC2, fFC2);
   finalDenseLayer->Initialize();
   TDenseLayer<AArchitecture> *copyFinalDenseLayer = new TDenseLayer<AArchitecture>(*finalDenseLayer);

   // add the copy to all slave nets
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].AddDenseLayer(copyFinalDenseLayer);
   }
}

/** Construct a random linear neural network with up to five layers.*/
//______________________________________________________________________________
template <typename AArchitecture>
void constructRandomLinearNet(TNet<AArchitecture> &net)
{
   int nlayers = rand() % 5 + 1;

   std::vector<EActivationFunction> ActivationFunctions = {EActivationFunction::kIdentity};

   for (int i = 0; i < nlayers; i++) {
      int width = rand() % 20 + 1;
      EActivationFunction f = ActivationFunctions[rand() % ActivationFunctions.size()];
      net.AddLayer(width, f);
   }
}

/*! Set matrix to the identity matrix */
//______________________________________________________________________________
template <typename AMatrix>
void identityMatrix(AMatrix &X)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         X(i, j) = 0.0;
      }
      if (i < n) {
         X(i, i) = 1.0;
      }
   }
}

/*! Fill matrix with given value.*/
//______________________________________________________________________________
template <typename AMatrix, typename AReal>
void fillMatrix(AMatrix &X, AReal x)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         X(i, j) = x;
      }
   }
}

/*! Fill matrix with random, Gaussian-distributed values. */
//______________________________________________________________________________
template <typename AMatrix>
void randomMatrix(AMatrix &X, double mean = 0.0, double sigma = 1.0,
                  TRandom & rand = *gRandom)
{
   size_t m = X.GetNrows();
   size_t n = X.GetNcols();

   for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < n; ++j)
         X(i,j) = rand.Gaus(mean, sigma);
}

/*! Fill matrix with random, uniform-distributed values in [-1, 1] */
//______________________________________________________________________________
template <typename AMatrix>
void uniformMatrix(AMatrix &X)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   TRandom & rand = *gRandom;

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         X(i, j) = rand.Uniform(-1, 1);
      }
   }
}

/*! Generate a random batch as input for a neural net. */
//______________________________________________________________________________
template <typename ATensor>
void randomBatch(ATensor &X)
{
   // for (size_t i = 0; i < X.GetFirstSize(); ++i) {
   //    auto rX = X.At(i);
   //    if (X.GetShape().size() == 3 || rX.GetFirstSize() == 1) {
   //       auto mX = rX.GetMatrix();
   //       randomMatrix(mX);
   //    }
   //    else {
   //       for (size_t j = 0; j < rX.GetFirstSize(); ++j ) {
   //          auto mX = rX.At(j).GetMatrix();
   //          randomMatrix(mX);
   //       }
   //    }
   // }
   using Scalar_t = typename ATensor::Scalar_t;
   TMatrixT<Scalar_t>  m(1,X.GetSize() ); 
   randomMatrix(m);
   X = ATensor(m.GetMatrixArray(), X.GetShape(), X.GetLayout() );
}

// one should use Architecture::Copy function
#if 0
/*! Generate a random batch as input for a neural net. */
//______________________________________________________________________________
template <typename AMatrix>
void copyMatrix(AMatrix &X, const AMatrix &Y)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         X(i, j) = Y(i, j);
      }
   }
}

template <typename ATensor>
void copyTensor(ATensor &X, const ATensor &Y)
{
   size_t n = Y.GetSize();
   assert (n == X.GetSize());
   std::copy(Y.GetData(), Y.GetData()+n, X.GetData());
}
#endif

/*! Apply functional to each element in the matrix. */
//______________________________________________________________________________
template <typename AMatrix, typename F>
void applyMatrix(AMatrix &X, F f)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         X(i, j) = f(X(i, j));
      }
   }
}

template <typename ATensor, typename F>
void applyTensor(ATensor &X, F f)
{
   size_t m = X.GetFirstSize();

   for (size_t i = 0; i < m; i++) {
      auto mX = X.At(i).GetMatrix();
      applyMatrix(mX,f);
   }
}

/*! Combine elements of two given matrices into a single matrix using
 *  the given function f. */
//______________________________________________________________________________
template <typename AMatrix, typename F>
void zipWithMatrix(AMatrix &Z, F f, const AMatrix &X, const AMatrix &Y)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         Z(i, j) = f(X(i, j), Y(i, j));
      }
   }
}

/** Generate a random batch as input for a neural net. */
//______________________________________________________________________________
template <typename AMatrix, typename AFloat, typename F>
AFloat reduce(F f, AFloat start, const AMatrix &X)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   AFloat result = start;

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         result = f(result, X(i, j));
      }
   }
   return result;
}

/** Apply function to matrix element-wise and compute the mean of the resulting
 *  element values */
//______________________________________________________________________________
template <typename AMatrix, typename AFloat, typename F>
AFloat reduceMean(F f, AFloat start, const AMatrix &X)
{
   size_t m, n;
   m = X.GetNrows();
   n = X.GetNcols();

   AFloat result = start;

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         result = f(result, X(i, j));
      }
   }
   return result / (AFloat)(m * n);
}

/** Compute the relative error of x and y. If their difference is too small,
 *  compute the absolute error instead. */
//______________________________________________________________________________
template <typename T>
inline T relativeError(const T &x, const T &y)
{
   using std::abs;

   if (x == y) return T(0.0);

   T diff = abs(x - y);

   if (x * y == T(0.0) || diff < std::numeric_limits<T>::epsilon()) return diff;

   return diff / (abs(x) + abs(y));
}

/*! Compute the maximum, element-wise relative error of the matrices
*  X and Y normalized by the element of Y. Protected against division
*  by zero. */
//______________________________________________________________________________
template <typename Matrix1, typename Matrix2>
auto maximumRelativeError(const Matrix1 &X, const Matrix2 &Y) -> Double_t
{
   Double_t curError, maxError = 0.0;

   Int_t m = X.GetNrows();
   Int_t n = X.GetNcols();

   assert(m == (Int_t) Y.GetNrows());
   assert(n == (Int_t) Y.GetNcols());

   for (Int_t i = 0; i < m; i++) {
      for (Int_t j = 0; j < n; j++) {
         curError = relativeError<Double_t>(X(i, j), Y(i, j));
         maxError = std::max(curError, maxError);
      }
   }

   return maxError;
}

template <typename Tensor1, typename Tensor2>
auto maximumRelativeErrorTensor(const Tensor1 &X, const Tensor2 &Y) -> Double_t
{

   size_t fsize = X.GetFirstSize();
   assert(fsize == Y.GetFirstSize());

   Double_t curError, maxError = 0.0;

   for (size_t i = 0; i < fsize; i++) {
      curError = maximumRelativeError( X.At(i).GetMatrix(), Y.At(i).GetMatrix() );
      maxError = std::max(curError, maxError);
   }


   return maxError;
}

/** Compute the average element-wise absolute error of the matrices
 *  X and Y.
 */

//______________________________________________________________________________
template <typename Matrix1, typename Matrix2>
auto meanAbsoluteError(const Matrix1 &X, const Matrix2 &Y) -> Double_t
{
   Double_t avgError = 0;

   Int_t m = X.GetNrows();
   Int_t n = X.GetNcols();

   assert(m == Y.GetNrows());
   assert(n == Y.GetNcols());

   for (Int_t i = 0; i < m; i++) {
      for (Int_t j = 0; j < n; j++) {
         avgError += std::abs(X(i, j) - Y(i, j));
      }
   }

   avgError /= (n * m);
   return avgError;
}

/*! Numerically compute the derivative of the functional f using finite
*  differences. */
//______________________________________________________________________________
template <typename F, typename AFloat>
inline AFloat finiteDifference(F f, AFloat dx)
{
   return f(dx) - f(-dx);
}

/*! Color code error. */
//______________________________________________________________________________
template <typename AFloat>
std::string print_error(AFloat &e)
{
   std::ostringstream out{};

   out << ("\033[");

   if (e > 1e-5)
      out << "31m";
   else if (e > 1e-9)
      out << "33m";
   else
      out << "32m";

   out << e;
   out << "\033[39m";

   return out.str();
}
}
}

#endif
