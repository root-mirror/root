// @(#)root/r:$Id$
// Author: Omar Zapata   30/05/2015


/*************************************************************************
 * Copyright (C) 2013-2015, Omar Andres Zapata Mesa                      *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_R_TRDataFrame
#define ROOT_R_TRDataFrame

#include <RExports.h>

#include <TRObject.h>

#include <TRFunctionImport.h>


namespace ROOT {
   namespace R {

      /**
      \class TRDataFrame

      This is a class to create DataFrames from ROOT to R
      <center><h2>TRDataFrame class</h2></center>

      DataFrame is a very important datatype in R and in ROOTR we have a class to manipulate<br>
      dataframes called TRDataFrame, with a lot of very useful operators overloaded to work with TRDataFrame's objects<br>
      in a similar way that in the R environment but from c++ in ROOT.<br>
      Example:<br>
      <br>
      Lets to create need data to play with dataframe features<br>

      <h2>Creating variables</h2><br>
      \code{.cpp}
      TVectorD v1(3);
      std::vector<Double_t> v2(3);
      std::array<Int_t,3>  v3{ {1,2,3} };
      std::list<std::string> names;
      \endcode

      <h2> Assigning values </h2><br>
      \code{.cpp}
      v1[0]=1;
      v1[1]=2;
      v1[2]=3;

      v2[0]=0.101;
      v2[1]=0.202;
      v2[2]=0.303;

      names.push_back("v1");
      names.push_back("v2");
      names.push_back("v3");

      ROOT::R::TRInterface &r=ROOT::R::TRInterface::Instance();
      \endcode

      In R the dataframe have associate to every column a label,
      in ROOTR you can have the same label using the class ROOT::R::Label to create a TRDataFrame where you data
      have a label associate.
      <h2> Creating dataframe object with its labels</h2> <br>
      \code{.cpp}
      using namespace ROOT::R;
      TRDataFrame  df1(Label["var1"]=v1,Label["var2"]=v2,Label["var3"]=v3,Label["strings"]=names);
      \endcode

      <h2>Passing dataframe to R's environment</h2><br>
      \code{.cpp}
      r["df1"]<<df1;
      r<<"print(df1)";
      \endcode
      Output
      \code
      var1  var2 var3 strings
      1    1 0.101    1      v1
      2    2 0.202    2      v2
      3    3 0.303    3      v3
      \endcode

      Manipulating data between dataframes
      <h2>Adding colunms to dataframe</h2><br>
      \code{.cpp}
      TVectorD v4(3);
      //filling the vector fro R's environment
      r["c(-1,-2,-3)"]>>v4;
      //adding new colunm to df1 with name var4
      df1["var4"]=v4;
      //updating df1 in R's environment
      r["df1"]<<df1;
      //printing df1
      r<<"print(df1)";
      \endcode

      Output
      var1  var2 var3 strings var4
      1    1 0.101    1      v1   -1
      2    2 0.202    2      v2   -2
      3    3 0.303    3      v3   -3

      <h2>Getting dataframe from R's environment</h2><br>
      \code{.cpp}
      ROOT::R::TRDataFrame df2;

      r<<"df2<-data.frame(v1=c(0.1,0.2,0.3),v2=c(3,2,1))";
      r["df2"]>>df2;

      TVectorD v(3);
      df2["v1"]>>v;
      v.Print();

      df2["v2"]>>v;
      v.Print();
      \endcode

      Output
      \code
      Vector (3)  is as follows

           |        1  |
      ------------------
         0 |0.1
         1 |0.2
         2 |0.3

      Vector (3)  is as follows

           |        1  |
      ------------------
         0 |3
         1 |2
         2 |1
      \endcode

      <h2>Working with colunms between dataframes</h2><br>
      \code{.cpp}
      df2["v3"]<<df1["strings"];

      //updating df2 in R's environment
      r["df2"]<<df2;
      r<<"print(df2)";
      \endcode
      Output
      \code
      v1 v2 v3
      1 0.1  3 v1
      2 0.2  2 v2
      3 0.3  1 v3
      \endcode

      <h2>Working with colunms between dataframes</h2><br>
      \code{.cpp}
      //passing values from colunm v3 of df2 to var1 of df1
      df2["v3"]>>df1["var1"];
      //updating df1 in R's environment
      r["df1"]<<df1;
      r<<"print(df1)";
      \endcode
      Output
      \code
      var1  var2 var3 strings var4
      1   v1 0.101    1      v1   -1
      2   v2 0.202    2      v2   -2
      3   v3 0.303    3      v3   -3
      \endcode
      <h2>Users Guide </h2>
      <a href="http://oproject.org/tiki-index.php?page=ROOT+R+Users+Guide"> http://oproject.org/tiki-index.php?page=ROOT+R+Users+Guide</a><br>
      <a href="https://root.cern.ch/drupal/content/how-use-r-root-root-r-interface"> https://root.cern.ch/drupal/content/how-use-r-root-root-r-interface</a>
         @ingroup R
      */


      class TRDataFrame: public TObject {
         friend class TRInterface;
         friend SEXP Rcpp::wrap<TRDataFrame>(const TRDataFrame &f);
      protected:
         Rcpp::DataFrame df; //internal Rcpp::DataFrame
      public:
         //Proxy class to use operators for assignation Ex: df["name"]>>object
         class Binding {
            friend class TRDataFrame;
         public:
            /**
            Construct a Binding nestead class for facilities with operators
            \param _df Rcpp::DataFrame (internal from TDataFrame)
            \param name string to use in assignations
            */
            Binding(Rcpp::DataFrame &_df, TString name): fName(name), fDf(_df) {}
            /**
            Copy constructor for Binding nestead class
            \param obj object with Rcpp::DataFame objecta and string with name
            */
            Binding(const Binding &obj): fName(obj.fName), fDf(obj.fDf) {}
            /**
            template method for operator assignation
            \param var any R wrappable datatype
            */
            template <class T> Binding operator=(T var)
            {
               Int_t size = fDf.size(), i = 0;
               Rcpp::CharacterVector names = fDf.attr("names");
               Bool_t found = false;
               while (i < size) {
                  if (names[i] == fName.Data()) {
                     found = true;
                     break;
                  }
                  i++;
               }
               if (found) fDf[fName.Data()] = var;
               else {
                  if (size == 0) {
                     fDf = Rcpp::DataFrame::create(ROOT::R::Label[fName.Data()] = var);
                  } else {
                     Rcpp::List nDf(size + 1);
                     Rcpp::CharacterVector nnames(size + 1);
                     for (i = 0; i < size; i++) {
                        nDf[i] = fDf[i] ;
                        nnames[i] = names[i];
                     }
                     nDf[size] = var;
                     nnames[size] = fName.Data();
                     nDf.attr("class") = fDf.attr("class") ;
                     nDf.attr("row.names") = fDf.attr("row.names") ;
                     nDf.attr("names") = nnames ;
                     fDf = nDf;
                  }
               }
               return *this;
            }
            /**
            method for operator assignation of Binding class
            \param obj other Binding object
            */
            Binding operator=(Binding obj)
            {
               Int_t size = fDf.size(), i = 0;
               Rcpp::CharacterVector names = fDf.attr("names");
               Bool_t found = false;
               while (i < size) {
                  if (names[i] == fName.Data()) {
                     found = true;
                     break;
                  }
                  i++;
               }
               if (found) fDf[fName.Data()] = obj.fDf[obj.fName.Data()];
               else {
                  Rcpp::List nDf(size + 1);
                  Rcpp::CharacterVector nnames(size + 1);
                  for (i = 0; i < size; i++) {
                     nDf[i] = obj.fDf[i] ;
                     nnames[i] = names[i];
                  }
                  nDf[size] = obj.fDf[obj.fName.Data()];
                  nnames[size] = fName.Data();

                  nDf.attr("class") = obj.fDf.attr("class") ;
                  nDf.attr("row.names") = obj.fDf.attr("row.names") ;
                  nDf.attr("names") = nnames ;
                  fDf = nDf;
               }

               return *this;
            }

            /**
            Template method for operator >> that lets to use dataframes like streams
            example: df["v"]>>vector;
            \param var any datatype that can be assigned from dataframe label
            */
            template <class T> Binding &operator >>(T &var)
            {
               var = Rcpp::as<T>(fDf[fName.Data()]);
               return *this;
            }
            Binding operator >>(Binding var)
            {
               var.fDf[var.fName.Data()] = fDf[fName.Data()];
               return var;
            }

            /**
            Template method for operator << that lets to use dataframes like streams
            example: df["v"]<<vector;
            \param var any datatype that can be assigned to dataframe label
            */
            template <class T> Binding &operator <<(T var)
            {
               Int_t size = fDf.size(), i = 0;
               Rcpp::CharacterVector names = fDf.attr("names");
               Bool_t found = false;
               while (i < size) {
                  if (names[i] == fName.Data()) {
                     found = true;
                     break;
                  }
                  i++;
               }
               if (found) fDf[fName.Data()] = var;
               else {
                  Rcpp::List nDf(size + 1);
                  Rcpp::CharacterVector nnames(size + 1);
                  for (i = 0; i < size; i++) {
                     nDf[i] = fDf[i] ;
                     nnames[i] = names[i];
                  }
                  nDf[size] = var;
                  nnames[size] = fName.Data();

                  nDf.attr("class") = fDf.attr("class") ;
                  nDf.attr("row.names") = fDf.attr("row.names") ;
                  nDf.attr("names") = nnames ;
                  fDf = nDf;
               }
               return *this;
            }
            template <class T> operator T()
            {
               return Rcpp::as<T>(fDf[fName.Data()]);
            }
            template <class T> operator T() const
            {
               return Rcpp::as<T>(fDf[fName.Data()]);
            }

         private:
            TString fName; //name of label
            Rcpp::DataFrame &fDf;//internal dataframe
         };

         /**
         Default TDataFrame constructor
         */
         TRDataFrame();
         /**
         TDataFrame constructor
         \param obj raw R object that can be casted to DataFrame
         */
         TRDataFrame(SEXP obj)
         {
            df = Rcpp::as<Rcpp::DataFrame>(obj);
         }
         /**
         TDataFrame copy constructor
         \param _df other TRDataFrame
         */
         TRDataFrame(const TRDataFrame &_df);
         /**
         TDataFrame constructor for Rcpp::DataFrame
         \param _df raw dataframe from Rcpp
         */
         TRDataFrame(const Rcpp::DataFrame &_df): df(_df) {};

#include <TRDataFrame__ctors.h>

         Binding operator[](const TString &name);

         TRDataFrame &operator=(TRDataFrame &obj)
         {
            df = obj.df;
            return *this;
         }

         TRDataFrame &operator=(TRDataFrame obj)
         {
            df = obj.df;
            return *this;
         }

         TRDataFrame &operator=(SEXP obj)
         {
            df = Rcpp::as<Rcpp::DataFrame>(obj);
            return *this;
         }

         operator SEXP()
         {
            return df;
         }

         operator SEXP() const
         {
            return df;
         }

         /**
         Method to get the number of colunms
         \return number of cols
         */
         Int_t GetNcols() { return df.size(); }
         /**
         Method to get the number of rows
         \return number of rows
         */
         Int_t GetNrows() { return df.nrows(); }
         /**
         Method to get labels of dataframe
         \return colunms names
         */
         TVectorString GetColNames()
         {
            Rcpp::CharacterVector names = df.attr("names");
            TVectorString rnames(GetNcols());
            for (Int_t i = 0; i < GetNcols(); i++) rnames[i] = names[i];
            return rnames;
         }

         /**
         Method to get dataframe as matrix
         \note only work on numerical dataframes if some column if string or other it will fail
         \return TMatrixT with a given tamplate data type
         */
         template<class T> TMatrixT<T> AsMatrix()
         {
            TRFunctionImport asMatrix("as.matrix");
            return Rcpp::as<TMatrixT<T> >(asMatrix(df));
         }

         /**
         Method to print the dataframe in stdout or a column given the label
         \param label nomber of the column to print
         */
         void Print(const Char_t *label = "")
         {
            TRFunctionImport print("print");
            if (label && !label[0]) {
               // label is ""
               print(df);
            } else {
               print(df[label]);
            }
         }
         ClassDef(TRDataFrame, 0) //
      };
   }
}



#endif
