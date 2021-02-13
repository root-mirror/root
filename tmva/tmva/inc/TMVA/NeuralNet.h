/**
 * @file TMVA/NeuralNet.h
 * @author  Peter Speckmayer
 * @version 1.0
 *
 * @section LICENSE
 *
 *
 * @section Neural net implementation
 *
 * An implementation of a neural net for TMVA. This neural net uses multithreading
 * 
 */


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// NeuralNet                                                            //
//                                                                      //
// A neural net implementation                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef TMVA_NEURAL_NET
#define TMVA_NEURAL_NET
#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <functional>
#include <tuple>
#include <cmath>
#include <cassert>
#include <random>
#include <thread>
#include <future>
#include <type_traits>
#include <string>
#include <utility>

#include "Pattern.h"
#include "Monitoring.h"

#include "TApplication.h"
#include "Timer.h"

#include "TH1F.h"
#include "TH2F.h"

#include <fenv.h> // turn on or off exceptions for NaN and other numeric exceptions


namespace TMVA
{

  class IPythonInteractive;

   namespace DNN
   {

      //    double gaussDoubl (edouble mean, double sigma);



      double gaussDouble (double mean, double sigma);
      double uniformDouble (double minValue, double maxValue);
      int randomInt (int maxValue);




      class MeanVariance
      {
      public:
      MeanVariance() 
          : m_n(0)
              , m_sumWeights(0)
              , m_mean(0)
              , m_squared(0)
          {}

          inline void clear() 
          { 
              m_n = 0; 
              m_sumWeights = 0;
              m_mean = 0;
              m_squared = 0;
          }

          template <typename T>
              inline void add(T value, double weight = 1.0)
          {
              ++m_n; // a value has been added

              if (m_n == 1) // initialization
              {
                  m_mean = value;
                  m_squared = 0.0;
                  m_sumWeights = weight;
                  return;
              }

              double tmpWeight = m_sumWeights+weight;
              double Q      = value - m_mean;

              double R = Q*weight/tmpWeight;
              m_mean    += R;
              m_squared += m_sumWeights*R*Q;

              m_sumWeights = tmpWeight;
          }

          template <typename ITERATOR>
              inline void add (ITERATOR itBegin, ITERATOR itEnd)
          {
              for (ITERATOR it = itBegin; it != itEnd; ++it)
                  add (*it);
          }



          inline int    count()      const { return m_n; }
          inline double weights()    const { if(m_n==0) return 0; return m_sumWeights; }
          inline double mean()       const { if(m_n==0) return 0; return m_mean; }
          inline double var() const
          {
              if(m_n==0)
                  return 0;
              if (m_squared <= 0)
                  return 0;
              return (m_squared/m_sumWeights);
          }
    
          inline double var_corr ()   const
          {
              if (m_n <= 1)
                  return var ();
        
              return (var()*m_n/(m_n-1));    // unbiased for small sample sizes
          } 
    
          inline double stdDev_corr () const { return sqrt( var_corr() ); }
          inline double stdDev ()   const { return sqrt( var() ); } // unbiased for small sample sizes

      private:
          size_t m_n;
          double m_sumWeights;
          double m_mean;
          double m_squared;
      };



      enum class EnumFunction
      {
         ZERO = '0',
            LINEAR = 'L',
            TANH = 'T',
            RELU = 'R',
            SYMMRELU = 'r',
            TANHSHIFT = 't',
            SIGMOID = 's',
            SOFTSIGN = 'S',
            GAUSS = 'G',
            GAUSSCOMPLEMENT = 'C'
            };



      enum class EnumRegularization
      {
         NONE, L1, L2, L1MAX
            };


      enum class ModeOutputValues : int
         {
            DIRECT = 0x01,
               SIGMOID = 0x02,
               SOFTMAX = 0x04,
               BATCHNORMALIZATION = 0x08
               };



      inline ModeOutputValues operator| (ModeOutputValues lhs, ModeOutputValues rhs)
      {
         return (ModeOutputValues)(static_cast<std::underlying_type<ModeOutputValues>::type>(lhs) | static_cast<std::underlying_type<ModeOutputValues>::type>(rhs));
      }

      inline ModeOutputValues operator|= (ModeOutputValues& lhs, ModeOutputValues rhs)
      {
         lhs = (ModeOutputValues)(static_cast<std::underlying_type<ModeOutputValues>::type>(lhs) | static_cast<std::underlying_type<ModeOutputValues>::type>(rhs));
         return lhs;
      }

      inline ModeOutputValues operator& (ModeOutputValues lhs, ModeOutputValues rhs)
      {
         return (ModeOutputValues)(static_cast<std::underlying_type<ModeOutputValues>::type>(lhs) & static_cast<std::underlying_type<ModeOutputValues>::type>(rhs));
      }

      inline ModeOutputValues operator&= (ModeOutputValues& lhs, ModeOutputValues rhs)
      {
         lhs = (ModeOutputValues)(static_cast<std::underlying_type<ModeOutputValues>::type>(lhs) & static_cast<std::underlying_type<ModeOutputValues>::type>(rhs));
         return lhs;
      }


      template <typename T>
         bool isFlagSet (T flag, T value)
         {
            return (int)(value & flag) != 0;
         }



      class Net;







      typedef std::vector<char> DropContainer;


      /*! \brief The Batch class encapsulates one mini-batch
       *
       *  Holds a const_iterator to the beginning and the end of one batch in a vector of Pattern
       */
      class Batch 
      {
      public:
         typedef typename std::vector<Pattern>::const_iterator const_iterator;

         Batch (typename std::vector<Pattern>::const_iterator itBegin, typename std::vector<Pattern>::const_iterator itEnd)
            : m_itBegin (itBegin)
            , m_itEnd (itEnd)
         {}

         const_iterator begin () const { return m_itBegin; }
         const_iterator end   () const { return m_itEnd; }

         size_t size () const { return std::distance (begin (), end ()); }
            
      private:
         const_iterator m_itBegin; ///< iterator denoting the beginning of the batch
         const_iterator m_itEnd;   ///< iterator denoting the end of the batch
      };






      template <typename ItSource, typename ItWeight, typename ItTarget>
         void applyWeights (ItSource itSourceBegin, ItSource itSourceEnd, ItWeight itWeight, ItTarget itTargetBegin, ItTarget itTargetEnd);



      template <typename ItSource, typename ItWeight, typename ItPrev>
         void applyWeightsBackwards (ItSource itCurrBegin, ItSource itCurrEnd, ItWeight itWeight, ItPrev itPrevBegin, ItPrev itPrevEnd);





      template <typename ItValue, typename ItFunction>
         void applyFunctions (ItValue itValue, ItValue itValueEnd, ItFunction itFunction);


      template <typename ItValue, typename ItFunction, typename ItInverseFunction, typename ItGradient>
         void applyFunctions (ItValue itValue, ItValue itValueEnd, ItFunction itFunction, ItInverseFunction itInverseFunction, ItGradient itGradient);



      template <typename ItSource, typename ItDelta, typename ItTargetGradient, typename ItGradient>
         void update (ItSource itSource, ItSource itSourceEnd, 
                      ItDelta itTargetDeltaBegin, ItDelta itTargetDeltaEnd, 
                      ItTargetGradient itTargetGradientBegin, 
                      ItGradient itGradient);



      template <EnumRegularization Regularization, typename ItSource, typename ItDelta, typename ItTargetGradient, typename ItGradient, typename ItWeight>
         void update (ItSource itSource, ItSource itSourceEnd, 
                      ItDelta itTargetDeltaBegin, ItDelta itTargetDeltaEnd, 
                      ItTargetGradient itTargetGradientBegin, 
                      ItGradient itGradient, 
                      ItWeight itWeight, double weightDecay);



      // ----- signature of a minimizer -------------
      // class Minimizer
      // {
      // public:

      //     template <typename Function, typename Variables, typename PassThrough>
      //     double operator() (Function& fnc, Variables& vars, PassThrough& passThrough) 
      //     {
      //         // auto itVars = begin (vars);
      //         // auto itVarsEnd = end (vars);

      //         std::vector<double> myweights;
      //         std::vector<double> gradients;

      //         double value = fnc (passThrough, myweights);
      //         value = fnc (passThrough, myweights, gradients);
      //         return value;
      //     } 
      // };



      ///< list all the minimizer types
      enum MinimizerType
      {
         fSteepest ///< SGD
      };





      /*! \brief Steepest Gradient Descent algorithm (SGD)
       *
       *  Implements a steepest gradient descent minimization algorithm
       */
      class Steepest
      {
      public:

         size_t m_repetitions;

    
         /*! \brief c'tor
          *
          *  C'tor
          * 
          * \param learningRate denotes the learning rate for the SGD algorithm
          * \param momentum fraction of the velocity which is taken over from the last step
          * \param repetitions re-compute the gradients each "repetitions" steps
          */
         Steepest (double learningRate = 1e-4, 
                   double momentum = 0.5, 
                   size_t repetitions = 10) 
            : m_repetitions (repetitions)
            , m_alpha (learningRate)
            , m_beta (momentum)
         {}

         /*! \brief operator to call the steepest gradient descent algorithm
          *
          *  entry point to start the minimization procedure
          * 
          * \param fitnessFunction (templated) function which has to be provided. This function is minimized
          * \param weights (templated) a reference to a container of weights. The result of the minimization procedure 
          *                is returned via this reference (needs to support std::begin and std::end
          * \param passThrough (templated) object which can hold any data which the fitness function needs. This object 
          *                    is not touched by the minimizer; This object is provided to the fitness function when
          *                    called
          */
         template <typename Function, typename Weights, typename PassThrough>
            double operator() (Function& fitnessFunction, Weights& weights, PassThrough& passThrough);


         double m_alpha; ///< internal parameter (learningRate)
         double m_beta;  ///< internal parameter (momentum)
         std::vector<double> m_prevGradients; ///< vector remembers the gradients of the previous step

         std::vector<double> m_localWeights; ///< local weights for reuse in thread. 
         std::vector<double> m_localGradients; ///< local gradients for reuse in thread. 
      };


















      template <typename ItOutput, typename ItTruth, typename ItDelta, typename ItInvActFnc>
         double sumOfSquares (ItOutput itOutputBegin, ItOutput itOutputEnd, ItTruth itTruthBegin, ItTruth itTruthEnd, ItDelta itDelta, ItDelta itDeltaEnd, ItInvActFnc itInvActFnc, double patternWeight);



      template <typename ItProbability, typename ItTruth, typename ItDelta, typename ItInvActFnc>
         double crossEntropy (ItProbability itProbabilityBegin, ItProbability itProbabilityEnd, ItTruth itTruthBegin, ItTruth itTruthEnd, ItDelta itDelta, ItDelta itDeltaEnd, ItInvActFnc itInvActFnc, double patternWeight);




      template <typename ItOutput, typename ItTruth, typename ItDelta, typename ItInvActFnc>
         double softMaxCrossEntropy (ItOutput itProbabilityBegin, ItOutput itProbabilityEnd, ItTruth itTruthBegin, ItTruth itTruthEnd, ItDelta itDelta, ItDelta itDeltaEnd, ItInvActFnc itInvActFnc, double patternWeight);





      template <typename ItWeight>
         double weightDecay (double error, ItWeight itWeight, ItWeight itWeightEnd, double factorWeightDecay, EnumRegularization eRegularization);














      /*! \brief LayerData holds the data of one layer
       *
       *     LayerData holds the data of one layer, but not its layout 
       *
       *  
       */
      class LayerData
      {
      public:
         typedef std::vector<double> container_type;

         typedef container_type::iterator iterator_type;
         typedef container_type::const_iterator const_iterator_type;

         typedef std::vector<std::function<double(double)> > function_container_type;
         typedef function_container_type::iterator function_iterator_type;
         typedef function_container_type::const_iterator const_function_iterator_type;

         typedef DropContainer::const_iterator const_dropout_iterator;
    
         /*! \brief c'tor of LayerData
          *
          *  C'tor of LayerData for the input layer
          * 
          * \param itInputBegin iterator to the begin of a vector which holds the values of the nodes of the neural net
          * \param itInputEnd iterator to the end of a vector which holdsd the values of the nodes of the neural net
          * \param eModeOutput indicates a potential tranformation of the output values before further computation
          *                    DIRECT does not further transformation; SIGMOID applies a sigmoid transformation to each
          *                    output value (to create a probability); SOFTMAX applies a softmax transformation to all 
          *                    output values (mutually exclusive probability)
          */
         LayerData (const_iterator_type itInputBegin, const_iterator_type itInputEnd, ModeOutputValues eModeOutput = ModeOutputValues::DIRECT);


         /*! \brief c'tor of LayerData
          *
          *  C'tor of LayerData for the input layer
          * 
          * \param inputSize input size of this layer
          */
         LayerData  (size_t inputSize);
         ~LayerData ()    {}


         /*! \brief c'tor of LayerData
          *
          *  C'tor of LayerData for all layers which are not the input layer; Used during the training of the DNN
          * 
          * \param size size of the layer
          * \param itWeightBegin indicates the start of the weights for this layer on the weight vector
          * \param itGradientBegin indicates the start of the gradients for this layer on the gradient vector
          * \param activationFunction indicates activation functions for this layer
          * \param inverseActivationFunction indicates the inverse activation functions for this layer
          * \param eModeOutput indicates a potential tranformation of the output values before further computation
          *                    DIRECT does not further transformation; SIGMOID applies a sigmoid transformation to each
          *                    output value (to create a probability); SOFTMAX applies a softmax transformation to all 
          *                    output values (mutually exclusive probability)
          */
         LayerData (size_t size, 
                    const_iterator_type itWeightBegin, 
                    iterator_type itGradientBegin, 
                    std::shared_ptr<std::function<double(double)>> activationFunction, 
                    std::shared_ptr<std::function<double(double)>> inverseActivationFunction,
                    ModeOutputValues eModeOutput = ModeOutputValues::DIRECT);

         /*! \brief c'tor of LayerData
          *
          *  C'tor of LayerData for all layers which are not the input layer; Used during the application of the DNN
          * 
          * \param size size of the layer
          * \param itWeightBegin indicates the start of the weights for this layer on the weight vector
          * \param activationFunction indicates the activation function for this layer
          * \param eModeOutput indicates a potential tranformation of the output values before further computation
          *                    DIRECT does not further transformation; SIGMOID applies a sigmoid transformation to each
          *                    output value (to create a probability); SOFTMAX applies a softmax transformation to all 
          *                    output values (mutually exclusive probability)
          */
         LayerData (size_t size, const_iterator_type itWeightBegin, 
                    std::shared_ptr<std::function<double(double)>> activationFunction, 
                    ModeOutputValues eModeOutput = ModeOutputValues::DIRECT);

         /*! \brief copy c'tor of LayerData
          *
          * 
          */
         LayerData (const LayerData& other)
            : m_size (other.m_size)
            , m_itInputBegin (other.m_itInputBegin)
            , m_itInputEnd (other.m_itInputEnd)
            , m_deltas (other.m_deltas)
            , m_valueGradients (other.m_valueGradients)
            , m_values (other.m_values)
        , m_itDropOut (other.m_itDropOut)
        , m_hasDropOut (other.m_hasDropOut)
            , m_itConstWeightBegin   (other.m_itConstWeightBegin)
            , m_itGradientBegin (other.m_itGradientBegin)
            , m_activationFunction (other.m_activationFunction)
            , m_inverseActivationFunction (other.m_inverseActivationFunction)
            , m_isInputLayer (other.m_isInputLayer)
            , m_hasWeights (other.m_hasWeights)
            , m_hasGradients (other.m_hasGradients)
            , m_eModeOutput (other.m_eModeOutput) 
            {}

         /*! \brief move c'tor of LayerData
          *
          * 
          */
         LayerData (LayerData&& other)
            : m_size (other.m_size)
            , m_itInputBegin (other.m_itInputBegin)
            , m_itInputEnd (other.m_itInputEnd)
        , m_deltas (std::move(other.m_deltas))
        , m_valueGradients (std::move(other.m_valueGradients))
        , m_values (std::move(other.m_values))
        , m_itDropOut (other.m_itDropOut)
        , m_hasDropOut (other.m_hasDropOut)
            , m_itConstWeightBegin   (other.m_itConstWeightBegin)
            , m_itGradientBegin (other.m_itGradientBegin)
        , m_activationFunction (std::move(other.m_activationFunction))
        , m_inverseActivationFunction (std::move(other.m_inverseActivationFunction))
            , m_isInputLayer (other.m_isInputLayer)
            , m_hasWeights (other.m_hasWeights)
            , m_hasGradients (other.m_hasGradients)
            , m_eModeOutput (other.m_eModeOutput) 
            {}


         /*! \brief change the input iterators
          *
          * 
          * \param itInputBegin indicates the start of the input node vector
          * \param itInputEnd indicates the end of the input node vector
          *
          */
         void setInput (const_iterator_type itInputBegin, const_iterator_type itInputEnd)
         {
            m_isInputLayer = true;
            m_itInputBegin = itInputBegin;
            m_itInputEnd = itInputEnd;
         }

         /*! \brief clear the values and the deltas
          *
          * 
          */
         void clear ()
         {
            m_values.assign (m_values.size (), 0.0);
            m_deltas.assign (m_deltas.size (), 0.0);
         }

         const_iterator_type valuesBegin () const { return m_isInputLayer ? m_itInputBegin : begin (m_values); } ///< returns const iterator to the begin of the (node) values
         const_iterator_type valuesEnd   () const { return m_isInputLayer ? m_itInputEnd   : end (m_values); } ///< returns iterator to the end of the (node) values
    
         iterator_type valuesBegin () { assert (!m_isInputLayer); return begin (m_values); }  ///< returns iterator to the begin of the (node) values
         iterator_type valuesEnd   () { assert (!m_isInputLayer); return end (m_values); } ///< returns iterator to the end of the (node) values

         ModeOutputValues outputMode () const { return m_eModeOutput; } ///< returns the output mode
    container_type probabilities () const { return computeProbabilities (); } ///< computes the probabilities from the current node values and returns them 

         iterator_type deltasBegin () { return begin (m_deltas); } ///< returns iterator to the begin of the deltas (back-propagation)
         iterator_type deltasEnd   () { return end   (m_deltas); } ///< returns iterator to the end of the deltas (back-propagation)

         const_iterator_type deltasBegin () const { return begin (m_deltas); } ///< returns const iterator to the begin of the deltas (back-propagation)
         const_iterator_type deltasEnd   () const { return end   (m_deltas); } ///< returns const iterator to the end of the deltas (back-propagation)

         iterator_type valueGradientsBegin () { return begin (m_valueGradients); } ///< returns iterator to the begin of the gradients of the node values
         iterator_type valueGradientsEnd   () { return end   (m_valueGradients); } ///< returns iterator to the end of the gradients of the node values

         const_iterator_type valueGradientsBegin () const { return begin (m_valueGradients); } ///< returns const iterator to the begin of the gradients
         const_iterator_type valueGradientsEnd   () const { return end   (m_valueGradients); } ///< returns const iterator to the end of the gradients

         iterator_type gradientsBegin () { assert (m_hasGradients); return m_itGradientBegin; } ///< returns iterator to the begin of the gradients
         const_iterator_type gradientsBegin () const { assert (m_hasGradients); return m_itGradientBegin; } ///< returns const iterator to the begin of the gradients
         const_iterator_type weightsBegin   () const { assert (m_hasWeights); return m_itConstWeightBegin; } ///< returns const iterator to the begin of the weights for this layer

         std::shared_ptr<std::function<double(double)>> activationFunction () const { return m_activationFunction; }
         std::shared_ptr<std::function<double(double)>> inverseActivationFunction () const { return m_inverseActivationFunction; }

         /*! \brief set the drop-out info for this layer
          *
          */
         template <typename Iterator>
            void setDropOut (Iterator itDrop) { m_itDropOut = itDrop; m_hasDropOut = true; }

         /*! \brief clear the drop-out-data for this layer
          *
          * 
          */
         void clearDropOut () { m_hasDropOut = false; }
    
         bool hasDropOut () const { return m_hasDropOut; } ///< has this layer drop-out turned on?
    const_dropout_iterator dropOut () const { assert (m_hasDropOut); return m_itDropOut; } ///< return the begin of the drop-out information
    
         size_t size () const { return m_size; } ///< return the size of the layer

      private:

         /*! \brief compute the probabilities from the node values
          *
          * 
          */
    container_type computeProbabilities () const;

      private:
    
         size_t m_size; ////< layer size

         const_iterator_type m_itInputBegin; ///< iterator to the first of the nodes in the input node vector
         const_iterator_type m_itInputEnd;   ///< iterator to the end of the nodes in the input node vector

         std::vector<double> m_deltas; ///< stores the deltas for the DNN training 
         std::vector<double> m_valueGradients; ///< stores the gradients of the values (nodes) 
         std::vector<double> m_values; ///< stores the values of the nodes in this layer
         const_dropout_iterator m_itDropOut; ///< iterator to a container indicating if the corresponding node is to be dropped
         bool m_hasDropOut; ///< dropOut is turned on?

         const_iterator_type m_itConstWeightBegin; ///< const iterator to the first weight of this layer in the weight vector
         iterator_type       m_itGradientBegin;  ///< iterator to the first gradient of this layer in the gradient vector

         std::shared_ptr<std::function<double(double)>> m_activationFunction; ///< activation function for this layer
         std::shared_ptr<std::function<double(double)>> m_inverseActivationFunction; ///< inverse activation function for this layer

         bool m_isInputLayer; ///< is this layer an input layer
         bool m_hasWeights;  ///< does this layer have weights (it does not if it is the input layer)
         bool m_hasGradients; ///< does this layer have gradients (only if in training mode)
 
         ModeOutputValues m_eModeOutput; ///< stores the output mode (DIRECT, SIGMOID, SOFTMAX)

      };





      /*! \brief Layer defines the layout of a layer
       *
       *     Layer defines the layout of a specific layer in the DNN
       *     Objects of this class don't hold the layer data itself (see class "LayerData")
       *  
       */
      class Layer
      {
      public:

         /*! \brief c'tor for defining a Layer
          *
          * 
          */
         Layer (size_t numNodes, EnumFunction activationFunction, ModeOutputValues eModeOutputValues = ModeOutputValues::DIRECT);

         ModeOutputValues modeOutputValues () const { return m_eModeOutputValues; } ///< get the mode-output-value (direct, probabilities)
         void modeOutputValues (ModeOutputValues eModeOutputValues) { m_eModeOutputValues = eModeOutputValues; } ///< set the mode-output-value

         size_t numNodes () const { return m_numNodes; } ///< return the number of nodes of this layer
         size_t numWeights (size_t numInputNodes) const { return numInputNodes * numNodes (); } ///< return the number of weights for this layer (fully connected)

         std::shared_ptr<std::function<double(double)>> activationFunction  () const { return m_activationFunction; } ///< fetch the activation function for this layer
         std::shared_ptr<std::function<double(double)>> inverseActivationFunction  () const { return m_inverseActivationFunction; } ///< fetch the inverse activation function for this layer

         EnumFunction activationFunctionType () const { return m_activationFunctionType; } ///< get the activation function type for this layer

      private:


         std::shared_ptr<std::function<double(double)>> m_activationFunction;  ///< stores the activation function
         std::shared_ptr<std::function<double(double)>> m_inverseActivationFunction;  ///< stores the inverse activation function


         size_t m_numNodes;

         ModeOutputValues m_eModeOutputValues; ///< do the output values of this layer have to be transformed somehow (e.g. to probabilities) or returned as such
         EnumFunction m_activationFunctionType;

         friend class Net;
      };





      template <typename LAYERDATA>
         void forward (const LAYERDATA& prevLayerData, LAYERDATA& currLayerData);


      template <typename LAYERDATA>
         void backward (LAYERDATA& prevLayerData, LAYERDATA& currLayerData);


      template <typename LAYERDATA>
         void update (const LAYERDATA& prevLayerData, LAYERDATA& currLayerData, double weightDecay, EnumRegularization regularization);



      /*! \brief Settings for the training of the neural net
       *
       * 
       */
      class Settings
      {
      public:

         /*! \brief c'tor
          *
          * 
          */
         Settings (TString name,
                   size_t _convergenceSteps = 15, size_t _batchSize = 10, size_t _testRepetitions = 7, 
                   double _factorWeightDecay = 1e-5, TMVA::DNN::EnumRegularization _regularization = TMVA::DNN::EnumRegularization::NONE,
                   MinimizerType _eMinimizerType = MinimizerType::fSteepest, 
                   double _learningRate = 1e-5, double _momentum = 0.3, 
                   int _repetitions = 3,
                   bool _multithreading = true);
    
         /*! \brief d'tor
          *
          * 
          */
         virtual ~Settings ();


         /*! \brief set the drop-out configuration (layer-wise)
          *
          * \param begin begin of an array or vector denoting the drop-out probabilities for each layer
          * \param end end of an array or vector denoting the drop-out probabilities for each layer 
          * \param _dropRepetitions denotes after how many repetitions the drop-out setting (which nodes are dropped out exactly) is changed
          */
         template <typename Iterator>
            void setDropOut (Iterator begin, Iterator end, size_t _dropRepetitions) { m_dropOut.assign (begin, end); m_dropRepetitions = _dropRepetitions; }

         size_t dropRepetitions () const { return m_dropRepetitions; }
         const std::vector<double>& dropFractions () const { return m_dropOut; }

         void setMonitoring (std::shared_ptr<Monitoring> ptrMonitoring) { fMonitoring = ptrMonitoring; } ///< prepared for monitoring

         size_t convergenceSteps () const { return m_convergenceSteps; } ///< how many steps until training is deemed to have converged
         size_t batchSize () const { return m_batchSize; } ///< mini-batch size
         size_t testRepetitions () const { return m_testRepetitions; } ///< how often is the test data tested
         double factorWeightDecay () const { return m_factorWeightDecay; } ///< get the weight-decay factor

         double learningRate () const { return fLearningRate; } ///< get the learning rate
         double momentum () const { return fMomentum; } ///< get the momentum (e.g. for SGD)
         int repetitions () const { return fRepetitions; } ///< how many steps have to be gone until the batch is changed
         MinimizerType minimizerType () const { return fMinimizerType; } ///< which minimizer shall be used (e.g. SGD)






         virtual void testSample (double /*error*/, double /*output*/, double /*target*/, double /*weight*/) {} ///< virtual function to be used for monitoring (callback)
         virtual void startTrainCycle () ///< callback for monitoring and logging
         {
            m_convergenceCount = 0;
            m_maxConvergenceCount= 0;
            m_minError = 1e10;
         }
         virtual void endTrainCycle (double /*error*/) {} ///< callback for monitoring and logging

         virtual void setProgressLimits (double minProgress = 0, double maxProgress = 100) ///< for monitoring and logging (set the current "progress" limits for the display of the progress)
         { 
            m_minProgress = minProgress;
            m_maxProgress = maxProgress; 
         }
         virtual void startTraining () ///< start drawing the progress bar
         {
            m_timer.DrawProgressBar (Int_t(m_minProgress));
         }
         virtual void cycle (double progress, TString text) ///< advance on the progress bar
         {
            m_timer.DrawProgressBar (Int_t(m_minProgress+(m_maxProgress-m_minProgress)*(progress/100.0)), text);
         }

         virtual void startTestCycle () {} ///< callback for monitoring and loggging
         virtual void endTestCycle () {} ///< callback for monitoring and loggging
         virtual void testIteration () {} ///< callback for monitoring and loggging
         virtual void drawSample (const std::vector<double>& /*input*/, const std::vector<double>& /* output */, const std::vector<double>& /* target */, double /* patternWeight */) {} ///< callback for monitoring and loggging

         virtual void computeResult (const Net& /* net */, std::vector<double>& /* weights */) {} ///< callback for monitoring and loggging

         virtual bool hasConverged (double testError); ///< has this training converged already?

         EnumRegularization regularization () const { return m_regularization; } ///< some regularization of the DNN is turned on?

         bool useMultithreading () const { return m_useMultithreading; } ///< is multithreading turned on?


         void pads (int numPads) { if (fMonitoring) fMonitoring->pads (numPads); } ///< preparation for monitoring
         void create (std::string histoName, int bins, double min, double max) { if (fMonitoring) fMonitoring->create (histoName, bins, min, max); } ///< for monitoring
         void create (std::string histoName, int bins, double min, double max, int bins2, double min2, double max2) { if (fMonitoring) fMonitoring->create (histoName, bins, min, max, bins2, min2, max2); } ///< for monitoring
         void addPoint (std::string histoName, double x) { if (fMonitoring) fMonitoring->addPoint (histoName, x); } ///< for monitoring
         void addPoint (std::string histoName, double x, double y) {if (fMonitoring) fMonitoring->addPoint (histoName, x, y); } ///< for monitoring
         void plot (std::string histoName, std::string options, int pad, EColor color) { if (fMonitoring) fMonitoring->plot (histoName, options, pad, color); } ///< for monitoring
         void clear (std::string histoName) { if (fMonitoring) fMonitoring->clear (histoName); } ///< for monitoring
         bool exists (std::string histoName) { if (fMonitoring) return fMonitoring->exists (histoName); return false; } ///< for monitoring

         size_t convergenceCount () const { return m_convergenceCount; } ///< returns the current convergence count
         size_t maxConvergenceCount () const { return m_maxConvergenceCount; } ///< returns the max convergence count so far
         size_t minError () const { return m_minError; } ///< returns the smallest error so far
    
      public:
         Timer  m_timer; ///< timer for monitoring
         double m_minProgress; ///< current limits for the progress bar
         double m_maxProgress; ///< current limits for the progress bar


         size_t m_convergenceSteps; ///< number of steps without improvement to consider the DNN to have converged
         size_t m_batchSize; ///< mini-batch size
         size_t m_testRepetitions; 
         double m_factorWeightDecay;

         size_t count_E;
         size_t count_dE;
         size_t count_mb_E;
         size_t count_mb_dE;

         EnumRegularization m_regularization;

         double m_dropRepetitions;
         std::vector<double> m_dropOut;

         double fLearningRate;
         double fMomentum;
         int fRepetitions;
         MinimizerType fMinimizerType;

         size_t m_convergenceCount;
         size_t m_maxConvergenceCount;
         double m_minError;


      protected:
         bool m_useMultithreading;

         std::shared_ptr<Monitoring> fMonitoring;
      };























      /*! \brief Settings for classification
       *
       * contains additional settings if the DNN problem is classification
       */
      class ClassificationSettings : public Settings
      {
      public:
         /*! \brief c'tor
          *
          * 
          */
         ClassificationSettings (TString name,
                                 size_t _convergenceSteps = 15, size_t _batchSize = 10, size_t _testRepetitions = 7, 
                                 double _factorWeightDecay = 1e-5, EnumRegularization _regularization = EnumRegularization::NONE, 
                                 size_t _scaleToNumEvents = 0, MinimizerType _eMinimizerType = MinimizerType::fSteepest, 
                                 double _learningRate = 1e-5, double _momentum = 0.3, int _repetitions = 3,
                                 bool _useMultithreading = true)
            : Settings (name, _convergenceSteps, _batchSize, _testRepetitions, _factorWeightDecay, 
                        _regularization, _eMinimizerType, _learningRate, _momentum, _repetitions, _useMultithreading)
            , m_ams ()
            , m_sumOfSigWeights (0)
            , m_sumOfBkgWeights (0)
            , m_scaleToNumEvents (_scaleToNumEvents)
            , m_cutValue (10.0)
            , m_pResultPatternContainer (NULL)
            , m_fileNameResult ()
            , m_fileNameNetConfig ()
            {
            }

         /*! \brief d'tor
          *
          * 
          */
         virtual ~ClassificationSettings () 
            {
            }

         void startTrainCycle ();
         void endTrainCycle (double /*error*/);
         void testIteration () { if (fMonitoring) fMonitoring->ProcessEvents (); }


         /* void createHistograms () */
         /* { */
         /*     std::cout << "is hist ROC existing?" << std::endl; */
         /*     if (m_histROC) */
         /*     { */
         /*         std::cout << "--> yes" << std::endl; */
         /*         fMonitoring->ProcessEvents (); */
         /*         return; */
         /*     } */

         /*     std::cout << "create histograms" << std::endl; */
         /*     TCanvas* canvas = fMonitoring->GetCanvas (); */
         /*     if (canvas) */
         /*     { */
         /*         std::cout << "canvas divide" << std::endl; */
         /*         canvas->cd (); */
         /*         canvas->Divide (2,2); */
         /*     } */
         /*     if (!m_histROC)  */
         /*     {  */
         /*         m_histROC = new TH2F ("ROC","ROC", 1000, 0, 1.0, 1000, 0, 1.0); m_histROC->SetDirectory (0);  */
         /*         m_histROC->SetLineColor (kBlue); */
         /*     } */
         /*     if (!m_histSignificance)  */
         /*     {  */
         /*         m_histSignificance = new TH2F ("Significance", "Significance", 1000, 0,1.0, 5, 0.0, 2.0);  */
         /*         m_histSignificance->SetDirectory (0);  */
         /*         m_histSignificance->SetBit (TH1::kCanRebin);  */
         /*         m_histROC->SetLineColor (kRed); */
         /*     } */
         /*     if (!m_histError)  */
         /*     {  */
         /*         m_histError = new TH1F ("Error", "Error", 100, 0, 100);  */
         /*         m_histError->SetDirectory (0);  */
         /*         m_histError->SetBit (TH1::kCanRebin); */
         /*         m_histROC->SetLineColor (kGreen); */
         /*     } */
         /*     if (!m_histOutputSignal)  */
         /*     {  */
         /*         m_histOutputSignal = new TH1F ("Signal", "Signal", 100, 0, 1.0);   */
         /*         m_histOutputSignal->SetDirectory (0);  */
         /*         m_histOutputSignal->SetBit (TH1::kCanRebin); */
         /*     } */
         /*     if (!m_histOutputBackground)  */
         /*     {  */
         /*         m_histOutputBackground = new TH1F ("Background", "Background", 100, 0, 1.0);  */
         /*         m_histOutputBackground->SetDirectory (0);  */
         /*         m_histOutputBackground->SetBit (TH1::kCanRebin); */
         /*     } */

         /*     fMonitoring->ProcessEvents (); */
         /* } */

         void testSample (double error, double output, double target, double weight);

         virtual void startTestCycle ();
         virtual void endTestCycle ();


         void setWeightSums (double sumOfSigWeights, double sumOfBkgWeights);
         void setResultComputation (std::string _fileNameNetConfig, std::string _fileNameResult, std::vector<Pattern>* _resultPatternContainer);

         std::vector<double> m_input;
         std::vector<double> m_output;
         std::vector<double> m_targets;
         std::vector<double> m_weights;

         std::vector<double> m_ams;
         std::vector<double> m_significances;


         double m_sumOfSigWeights;
         double m_sumOfBkgWeights;
         size_t m_scaleToNumEvents;

         double m_cutValue;
         std::vector<Pattern>* m_pResultPatternContainer;
         std::string m_fileNameResult;
         std::string m_fileNameNetConfig;


         /* TH2F* m_histROC; */
         /* TH2F* m_histSignificance; */

         /* TH1F* m_histError; */
         /* TH1F* m_histOutputSignal; */
         /* TH1F* m_histOutputBackground; */
      };







      ///< used to distinguish between different function signatures
      enum class ModeOutput
      {
         FETCH
            };

      /*! \brief error functions to be chosen from 
       *
       * 
       */
      enum class ModeErrorFunction
      {
         SUMOFSQUARES = 'S',
            CROSSENTROPY = 'C',
            CROSSENTROPY_MUTUALEXCLUSIVE = 'M'
            };

      /*! \brief weight initialization strategies to be chosen from
       *
       * 
       */
      enum class WeightInitializationStrategy
      {
         XAVIER, TEST, LAYERSIZE, XAVIERUNIFORM
            };



      /*! \brief neural net 
       *
       * holds the structure of all layers and some data for the whole net
       * does not know the layer data though (i.e. values of the nodes and weights)
       */
      class Net
      {
      public:

         typedef std::vector<double> container_type;
         typedef container_type::iterator iterator_type;
         typedef std::pair<iterator_type,iterator_type> begin_end_type;


         /*! \brief c'tor
          *
          * 
          */
         Net () 
            : m_eErrorFunction (ModeErrorFunction::SUMOFSQUARES)
            , m_sizeInput (0)
            , m_layers ()
            {
            }

         /*! \brief d'tor
          *
          * 
          */
         Net (const Net& other)
            : m_eErrorFunction (other.m_eErrorFunction)
            , m_sizeInput (other.m_sizeInput)
            , m_layers (other.m_layers)
            {
            }

         void setInputSize (size_t sizeInput) { m_sizeInput = sizeInput; } ///< set the input size of the DNN
         void setOutputSize (size_t sizeOutput) { m_sizeOutput = sizeOutput; } ///< set the output size of the DNN
         void addLayer (Layer& layer) { m_layers.push_back (layer); } ///< add a layer (layout)
         void addLayer (Layer&& layer) { m_layers.push_back (layer); } 
         void setErrorFunction (ModeErrorFunction eErrorFunction) { m_eErrorFunction = eErrorFunction; } ///< which error function is to be used
    
         size_t inputSize () const { return m_sizeInput; } ///< input size of the DNN
         size_t outputSize () const { return m_sizeOutput; } ///< output size of the DNN

         /*! \brief set the drop out configuration
          *
          * 
          */
         template <typename WeightsType, typename DropProbabilities>
            void dropOutWeightFactor (WeightsType& weights,
                                      const DropProbabilities& drops, 
                                      bool inverse = false);

         /*! \brief start the training
          *
          * \param weights weight vector
          * \param trainPattern training pattern 
          * \param testPattern test pattern
          * \param minimizer use this minimizer for training (e.g. SGD)
          * \param settings settings used for this training run
          */
         template <typename Minimizer>
            double train (std::vector<double>& weights, 
                          std::vector<Pattern>& trainPattern, 
                          const std::vector<Pattern>& testPattern, 
                  Minimizer& minimizer,
                  Settings& settings);

         /*! \brief pre-training for future use
          *
          * 
          */
         template <typename Minimizer>
            void preTrain (std::vector<double>& weights,
                           std::vector<Pattern>& trainPattern,
                           const std::vector<Pattern>& testPattern,
                           Minimizer& minimizer, Settings& settings);

    
         /*! \brief executes one training cycle
          *
          * \param minimizer the minimizer to be used
          * \param weights the weight vector to be used
          * \param itPatternBegin the pattern to be trained with
          * \param itPatternEnd the pattern to be trainied with
          * \param settings the settings for the training
          * \param dropContainer the configuration for DNN drop-out
          */
         template <typename Iterator, typename Minimizer>
            inline double trainCycle (Minimizer& minimizer, std::vector<double>& weights, 
			      Iterator itPatternBegin, Iterator itPatternEnd,
                              Settings& settings,
                              DropContainer& dropContainer);

         size_t numWeights (size_t trainingStartLayer = 0) const; ///< returns the number of weights in this net
    size_t numNodes   (size_t trainingStartLayer = 0) const; ///< returns the number of nodes in this net

         template <typename Weights>
            std::vector<double> compute (const std::vector<double>& input, const Weights& weights) const; ///< compute the net with the given input and the given weights

         template <typename Weights, typename PassThrough>
            double operator() (PassThrough& settingsAndBatch, const Weights& weights) const; ///< execute computation of the DNN for one mini-batch (used by the minimizer); no computation of gradients

         template <typename Weights, typename PassThrough, typename OutContainer>
            double operator() (PassThrough& settingsAndBatch, const Weights& weights, ModeOutput eFetch, OutContainer& outputContainer) const; ///< execute computation of the DNN for one mini-batch; helper function
    
         template <typename Weights, typename Gradients, typename PassThrough>
        double operator() (PassThrough& settingsAndBatch, Weights& weights, Gradients& gradients) const;  ///< execute computation of the DNN for one mini-batch (used by the minimizer); returns gradients as well

         template <typename Weights, typename Gradients, typename PassThrough, typename OutContainer>
        double operator() (PassThrough& settingsAndBatch, Weights& weights, Gradients& gradients, ModeOutput eFetch, OutContainer& outputContainer) const;


    template <typename LayerContainer, typename DropContainer, typename ItWeight, typename ItGradient>
        std::vector<std::vector<LayerData>> prepareLayerData (LayerContainer& layers,
                                                              Batch& batch,
                                                              const DropContainer& dropContainer,
                                                              ItWeight itWeightBegin,
                                                              ItWeight itWeightEnd, 
                                                              ItGradient itGradientBegin,
                                                              ItGradient itGradientEnd,
                                                              size_t& totalNumWeights) const;

    template <typename LayerContainer>
        void forwardPattern (const LayerContainer& _layers,
                             std::vector<LayerData>& layerData) const;


    template <typename LayerContainer, typename LayerPatternContainer>
        void forwardBatch (const LayerContainer& _layers,
                           LayerPatternContainer& layerPatternData,
                           std::vector<double>& valuesMean,
                           std::vector<double>& valuesStdDev,
                           size_t trainFromLayer) const;
    
    template <typename OutputContainer>
        void fetchOutput (const LayerData& lastLayerData, OutputContainer& outputContainer) const;

    template <typename OutputContainer>
        void fetchOutput (const std::vector<LayerData>& layerPatternData, OutputContainer& outputContainer) const;


    template <typename ItWeight>
        std::tuple</*sumError*/double,/*sumWeights*/double> computeError (const Settings& settings,
                                                                          std::vector<LayerData>& lastLayerData,
                                                                          Batch& batch,
                                                                          ItWeight itWeightBegin,
                                                                          ItWeight itWeightEnd) const;

    template <typename Settings>
        void backPropagate (std::vector<std::vector<LayerData>>& layerPatternData,
                            const Settings& settings,
                            size_t trainFromLayer,
                            size_t totalNumWeights) const;

    

    /*! \brief main NN computation function
          *
          * 
          */
         template <typename LayerContainer, typename PassThrough, typename ItWeight, typename ItGradient, typename OutContainer>
            double forward_backward (LayerContainer& layers, PassThrough& settingsAndBatch, 
			     ItWeight itWeightBegin, ItWeight itWeightEnd, 
                                     ItGradient itGradientBegin, ItGradient itGradientEnd, 
                                     size_t trainFromLayer, 
                                     OutContainer& outputContainer, bool fetchOutput) const;


    
         double E ();
         void dE ();


         /*! \brief computes the error of the DNN
          *
          * 
          */
         template <typename Container, typename ItWeight>
            double errorFunction (LayerData& layerData,
                                  Container truth,
                                  ItWeight itWeight,
                                  ItWeight itWeightEnd,
                                  double patternWeight,
                                  double factorWeightDecay,
                                  EnumRegularization eRegularization) const;


         const std::vector<Layer>& layers () const { return m_layers; } ///< returns the layers (structure)
         std::vector<Layer>& layers ()  { return m_layers; } ///< returns the layers (structure)

         void removeLayer () { m_layers.pop_back (); } ///< remove one layer
    

         void clear () ///< clear one layer
         {
            m_layers.clear ();
            m_eErrorFunction = ModeErrorFunction::SUMOFSQUARES;
         }


         template <typename OutIterator>
            void initializeWeights (WeightInitializationStrategy eInitStrategy, 
                                    OutIterator itWeight); ///< initialize the weights with the given strategy

      protected:

         void fillDropContainer (DropContainer& dropContainer, double dropFraction, size_t numNodes) const; ///< prepare the drop-out-container (select the nodes which are to be dropped out)
    
    
      private:

         ModeErrorFunction m_eErrorFunction; ///< denotes the error function
         size_t m_sizeInput; ///< input size of this DNN
         size_t m_sizeOutput; ///< outut size of this DNN
         std::vector<Layer> m_layers; ///< layer-structure-data

      protected:
         // variables for JsMVA (interactive training in jupyter notebook)
         IPythonInteractive *fInteractive = nullptr;
         bool * fExitFromTraining = nullptr;
         UInt_t *fIPyMaxIter = nullptr, *fIPyCurrentIter = nullptr;

      public:

        // setup ipython interactive variables
        void SetIpythonInteractive(IPythonInteractive* fI, bool* fE, UInt_t *M, UInt_t *C){
          fInteractive = fI;
          fExitFromTraining = fE;
          fIPyMaxIter = M;
          fIPyCurrentIter = C;
        }
      };




typedef std::tuple<Settings&, Batch&, DropContainer&> pass_through_type;







   } // namespace DNN
} // namespace TMVA


// include the implementations (in header file, because they are templated)
#include "TMVA/NeuralNet.icc"

#endif

