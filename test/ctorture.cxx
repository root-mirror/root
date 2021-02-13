// @(#)root/test:$Id$
// Author: Federico Carminati    22/04/2004

// test program for the class TComplex

#include <iostream>
#include <cstdio>

#include <strlcpy.h>
#include <TRandom.h>
#include <TComplex.h>

Double_t Error(TComplex a, TComplex b)
{return 2*TComplex::Abs(a-b)/(a.Rho()+b.Rho());}

void Verify(const TComplex a, const TComplex b,
            Double_t epsmin, Double_t epsmax,
            const char* where, Int_t & ifail, Double_t &serr)
{
  Double_t err=Error(a,b);
  serr+=err;
  if(epsmin<err) {
    ifail++;
    if(err<epsmax) {
      printf("Fail %s %e\n",where,err);
      std::cout << a << std::endl;
      std::cout << b << std::endl;
    }
  }
}

void Summary(const char* title, Int_t ifail, Double_t serr, Int_t np)
{
  printf("Results for %s\n",title);
  printf("Fail= %5.2f%%, av err= %e\n",100.*ifail/np,serr/np);
}


int main () {
  //
  // Torture for complex numbers
  //

  const Int_t np=10000;
  Int_t i;
  Int_t ifail;
  Double_t serr;
  Double_t x;
  char title[20];
  TComplex a,b,c,d,e;

  // Torture Square roots
  serr=0;
  ifail=0;
  strlcpy(title,"Sqrt",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=a*a;
    c=TComplex::Sqrt(b);
    // Cater for the fact that there are two roots!
    if(a.Re()*c.Re()<0)  c*=-1;
    Verify(a,c,1e-14,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture exp and log
  serr=0;
  ifail=0;
  strlcpy(title,"Exp&Log",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::Log(a);
    c=TComplex::Exp(b);
    Verify(a,c,1e-14,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture sin and asin
  serr=0;
  ifail=0;
  strlcpy(title,"Sin&ASin",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::ASin(a);
    c=TComplex::Sin(b);
    Verify(a,c,1e-13,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture cos and acos
  serr=0;
  ifail=0;
  strlcpy(title,"Cos&ACos",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::ACos(a);
    c=TComplex::Cos(b);
    Verify(a,c,1e-13,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture tan and atan
  serr=0;
  ifail=0;
  strlcpy(title,"Tan&ATan",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::ATan(a);
    c=TComplex::Tan(b);
    Verify(a,c,1e-14,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture SinH and ASinH
  serr=0;
  ifail=0;
  strlcpy(title,"SinH&ASinH",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::ASinH(a);
    c=TComplex::SinH(b);
    Verify(a,c,1e-13,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture CosH and ACosH
  serr=0;
  ifail=0;
  strlcpy(title,"CosH&ACosH",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::ACosH(a);
    c=TComplex::CosH(b);
    Verify(a,c,1e-13,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture TanH and ATanH
  serr=0;
  ifail=0;
  strlcpy(title,"TanH&ATanH",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    b=TComplex::ATanH(a);
    c=TComplex::TanH(b);
    Verify(a,c,1e-14,1e10,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture Power complex - complex
  //
  // Important note on the following tests. The operation of raising a complex
  // number to a power does not yet a single value, but rather an infinite
  // number of values, particularly if the number is non rational.
  // For a real number a, you can define a^(b+ic) by writing a = e^(ln a):
  //
  //       b+ic    (ln a)(b+ic)    (b ln a) + i(c ln a)
  //      a     = e             = e
  //
  //               (b ln a)
  //            = e         ( cos (c ln a) + i sin (c ln a) )
  //
  //               b
  //            = a  ( cos (c ln a) + i sin (c ln a) ).
  //
  // Now, if a is a complex number instead of a real number, there is no
  // single value to "ln a": there are lots of different complex numbers z
  // for which e^z = a, and for any such complex number z, you could define
  // a^(b+ic) to be e^(z(b+ic)) and use the above technique to calculate it.
  //
  // In fact, the same thing is true even when a is a real number. The
  // expression a^(b+ic) has many possible values (infinite except when b
  // and c are both rational numbers), because instead of doing the calculation
  // writing a = e^(ln a), you could also do it by writing a = e^(ln a + 2pi i)
  // or by writing a = e^(ln a + 4 pi i), or a = e^(ln a + 6 pi i), and so on.
  // Each of these equalities is true (in fact e^(2pi n i)=1 for integer n).
  //
  // When a is real it is more "natural" to use the ordinary real-valued
  // logarithm ln a rather than than something like ln a + 2 pi i.
  // Technically, this value is called the principal value. This is what
  // the formula up above gives you. Unfortunately this alone does not
  // guarantees that the inverse operation brings you back where you
  // started from.
  //
  // When a is not real there is no one natural choice of logarithm to prefer
  // over any other, so in those cases we have to say that an expression like
  // a^(b+ic) has many different values.
  //
  // This is because in these tests we exclude from the error output the
  // results where we ended up very far from the initial value, and the
  // difference is more than 50%.
  //
  serr=0;
  ifail=0;
  strlcpy(title,"Power C-C",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    while (TComplex::Abs(
                         e=TComplex(10*(1-2*gRandom->Rndm()),
                                    10*(1-2*gRandom->Rndm()))
                         )<0.1) { }
    b=TComplex::Power(a,1./e);
    c=TComplex::Power(b,e);
    Verify(a,c,2e-14,1.,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture Power complex - real
  serr=0;
  ifail=0;
  strlcpy(title,"Power C-R",20);
  for(i=0; i<np; i++) {
    a=TComplex(10*(1-2*gRandom->Rndm()),10*(1-2*gRandom->Rndm()));
    while (TMath::Abs(x=10*(1-2*gRandom->Rndm()))<0.1) { }
    b=TComplex::Power(a,1./x);
    c=TComplex::Power(b,x);
    Verify(a,c,5e-14,.5,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  // Torture Power real - complex
  serr=0;
  ifail=0;
  strlcpy(title,"Power R-C",20);
  for(i=0; i<np; i++) {
    while (TComplex::Abs(
                         a=TComplex(10*(1-2*gRandom->Rndm()),
                                    10*(1-2*gRandom->Rndm()))
                         )<0.1) { }
    x=10*(1-2*gRandom->Rndm());
    b=TComplex::Power(x,1./a);
    c=TComplex::Power(b,a);
    Verify(c,TComplex(x,0),2e-14,1.5,title,ifail,serr);
  }
  Summary(title,ifail,serr,np);

  return 0;
}
