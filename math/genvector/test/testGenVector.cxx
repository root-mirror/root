
#include "Math/Vector3D.h"
#include "Math/Point3D.h"

#include "Math/Vector2D.h"
#include "Math/Point2D.h"

#include "Math/EulerAngles.h"

#include "Math/Transform3D.h"
#include "Math/Translation3D.h"

#include "Math/Rotation3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Quaternion.h"
#include "Math/AxisAngle.h"
#include "Math/RotationZYX.h"

#include "Math/LorentzRotation.h"

#include "Math/VectorUtil.h"
#ifndef NO_SMATRIX
#include "Math/SMatrix.h"
#endif

#include <vector>

using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;



typedef  DisplacementVector3D<Cartesian3D<double>, GlobalCoordinateSystemTag>  GlobalXYZVector;
typedef  DisplacementVector3D<Cartesian3D<double>, LocalCoordinateSystemTag>   LocalXYZVector;
typedef  DisplacementVector3D<Polar3D<double>, GlobalCoordinateSystemTag>      GlobalPolar3DVector;



typedef  PositionVector3D<Cartesian3D<double>, GlobalCoordinateSystemTag>   GlobalXYZPoint;
typedef  PositionVector3D<Cartesian3D<double>, LocalCoordinateSystemTag>    LocalXYZPoint;
typedef  PositionVector3D<Polar3D<double>, GlobalCoordinateSystemTag>       GlobalPolar3DPoint;
typedef  PositionVector3D<Polar3D<double>, LocalCoordinateSystemTag>       LocalPolar3DPoint;


//#define TEST_COMPILE_ERROR


int compare( double v1, double v2, const std::string & name = "", double scale = 1.0) {
  //  ntest = ntest + 1;

  // numerical double limit for epsilon
  double eps = scale* std::numeric_limits<double>::epsilon();
  int iret = 0;
  double delta = v2 - v1;
  double d = 0;
  if (delta < 0 ) delta = - delta;
  if (v1 == 0 || v2 == 0) {
    if  (delta > eps ) {
      iret = 1;
    }
  }
  // skip case v1 or v2 is infinity
  else {
     d = v1;

    if ( v1 < 0) d = -d;
    // add also case when delta is small by default
    if ( delta/d  > eps && delta > eps )
      iret =  1;
  }

  if (iret == 0)
    std::cout << ".";
  else {
    int pr = std::cout.precision (18);
    std::cout << "\nDiscrepancy in " << name << "() : " << v1 << " != " << v2 << " discr = " << int(delta/d/eps)
              << "   (Allowed discrepancy is " << eps  << ")\n";
    std::cout.precision (pr);
    //nfail = nfail + 1;
  }
  return iret;
}

template<class Transform>
bool IsEqual(const Transform & t1, const Transform & t2, unsigned int size)  {
// size should be an enum of the Transform class
   std::vector<double> x1(size);
   std::vector<double> x2(size);
   t1.GetComponents(x1.begin(), x1.end() );
   t2.GetComponents(x2.begin(), x2.end() );
   bool ret = true;
   unsigned int i = 0;
   while (ret && i < size) {
      // from TMath::AreEqualRel(x1,x2,2*eps)
      bool areEqual = std::abs(x1[i]-x2[i]) < std::numeric_limits<double>::epsilon() *
         ( std::abs(x1[i]) + std::abs(x2[i] ) );
      ret &= areEqual;
      i++;
   }
   return ret;
}

int testVector3D() {

  int iret = 0;

  std::cout << "testing Vector3D   \t:\t";

  // test the vector tags

  GlobalXYZVector vg(1.,2.,3.);
  GlobalXYZVector vg2(vg);

  GlobalPolar3DVector vpg(vg);

  iret |= compare(vpg.R(), vg2.R() );

//   std::cout << vg2 << std::endl;

  double r = vg.Dot(vpg);
  iret |= compare(r, vg.Mag2() );

  GlobalXYZVector vcross = vg.Cross(vpg);
  iret |= compare(vcross.R(), 0.0,"cross",10 );

//   std::cout << vg.Dot(vpg) << std::endl;
//   std::cout << vg.Cross(vpg) << std::endl;





  GlobalXYZVector vg3 = vg + vpg;
  iret |= compare(vg3.R(), 2*vg.R() );

  GlobalXYZVector vg4 = vg - vpg;
  iret |= compare(vg4.R(), 0.0,"diff",10 );




#ifdef TEST_COMPILE_ERROR
   LocalXYZVector vl; vl = vg;
  LocalXYZVector vl2(vg2);
  LocalXYZVector vl3(vpg);
  vg.Dot(vl);
  vg.Cross(vl);
  vg3 = vg + vl;
  vg4 = vg - vl;
#endif


  if (iret == 0) std::cout << "\t\t\t\t\tOK\n";
  else std::cout << "\t\t\t\tFAILED\n";


  return iret;
}



int testPoint3D() {

  int iret = 0;

  std::cout << "testing Point3D    \t:\t";

  // test the vector tags

  GlobalXYZPoint pg(1.,2.,3.);
  GlobalXYZPoint pg2(pg);

  GlobalPolar3DPoint ppg(pg);

  iret |= compare(ppg.R(), pg2.R() );
  //std::cout << pg2 << std::endl;




  GlobalXYZVector vg(pg);

  double r = pg.Dot(vg);
  iret |= compare(r, pg.Mag2() );

  GlobalPolar3DVector vpg(pg);
  GlobalXYZPoint pcross = pg.Cross(vpg);
  iret |= compare(pcross.R(), 0.0,"cross",10 );

  GlobalPolar3DPoint pg3 = ppg + vg;
  iret |= compare(pg3.R(), 2*pg.R() );

  GlobalXYZVector vg4 = pg - ppg;
  iret |= compare(vg4.R(), 0.0,"diff",10 );


#ifdef TEST_COMPILE_ERROR
  LocalXYZPoint pl; pl = pg;
  LocalXYZVector pl2(pg2);
  LocalXYZVector pl3(ppg);
  pl.Dot(vg);
  pl.Cross(vg);
  pg3 = ppg + pg;
  pg3 = ppg + pl;
  vg4 = pg - pl;
#endif

  // operator -
  XYZPoint q1(1.,2.,3.);
  XYZPoint q2 = -1.* q1;
  XYZVector v2 = -XYZVector(q1);
  iret |= compare(XYZVector(q2) == v2,true,"reflection");


  if (iret == 0) std::cout << "\t\t\t\t\tOK\n";
  else std::cout << "\t\t\t\tFAILED\n";

  return iret;
}



typedef  DisplacementVector2D<Cartesian2D<double>, GlobalCoordinateSystemTag>  GlobalXYVector;
typedef  DisplacementVector2D<Cartesian2D<double>, LocalCoordinateSystemTag>   LocalXYVector;
typedef  DisplacementVector2D<Polar2D<double>, GlobalCoordinateSystemTag>      GlobalPolar2DVector;




int testVector2D() {

  int iret = 0;

  std::cout << "testing Vector2D   \t:\t";

  // test the vector tags

  GlobalXYVector vg(1.,2.);
  GlobalXYVector vg2(vg);

  GlobalPolar2DVector vpg(vg);

  iret |= compare(vpg.R(), vg2.R() );

//   std::cout << vg2 << std::endl;

  double r = vg.Dot(vpg);
  iret |= compare(r, vg.Mag2() );

//   std::cout << vg.Dot(vpg) << std::endl;


  GlobalXYVector vg3 = vg + vpg;
  iret |= compare(vg3.R(), 2*vg.R() );

  GlobalXYVector vg4 = vg - vpg;
  iret |= compare(vg4.R(), 0.0,"diff",10 );


  double angle = 1.;
  vg.Rotate(angle);
  iret |= compare(vg.Phi(), vpg.Phi() + angle );
  iret |= compare(vg.R(), vpg.R()  );

  GlobalXYZVector v3d(1,2,0);
  GlobalXYZVector vr3d = RotationZ(angle) * v3d;
  iret |= compare(vg.X(), vr3d.X() );
  iret |= compare(vg.Y(), vr3d.Y()  );

  GlobalXYVector vu = vg3.Unit();
  iret |= compare(vu.R(), 1. );


#ifdef TEST_COMPILE_ERROR
   LocalXYVector vl; vl = vg;
  LocalXYVector vl2(vg2);
  LocalXYVector vl3(vpg);
  vg.Dot(vl);
  vg3 = vg + vl;
  vg4 = vg - vl;
#endif


  if (iret == 0) std::cout << "\t\t\t\tOK\n";
  else std::cout << "\t\t\tFAILED\n";


  return iret;
}


typedef  PositionVector2D<Cartesian2D<double>, GlobalCoordinateSystemTag>   GlobalXYPoint;
typedef  PositionVector2D<Cartesian2D<double>, LocalCoordinateSystemTag>    LocalXYPoint;
typedef  PositionVector2D<Polar2D<double>, GlobalCoordinateSystemTag>       GlobalPolar2DPoint;
typedef  PositionVector2D<Polar2D<double>, LocalCoordinateSystemTag>       LocalPolar2DPoint;



int testPoint2D() {

  int iret = 0;

  std::cout << "testing Point2D    \t:\t";

  // test the vector tags

  GlobalXYPoint pg(1.,2.);
  GlobalXYPoint pg2(pg);

  GlobalPolar2DPoint ppg(pg);

  iret |= compare(ppg.R(), pg2.R() );
  //std::cout << pg2 << std::endl;




  GlobalXYVector vg(pg);

  double r = pg.Dot(vg);
  iret |= compare(r, pg.Mag2() );

  GlobalPolar2DVector vpg(pg);

  GlobalPolar2DPoint pg3 = ppg + vg;
  iret |= compare(pg3.R(), 2*pg.R() );

  GlobalXYVector vg4 = pg - ppg;
  iret |= compare(vg4.R(), 0.0,"diff",10 );


#ifdef TEST_COMPILE_ERROR
  LocalXYPoint pl; pl = pg;
  LocalXYVector pl2(pg2);
  LocalXYVector pl3(ppg);
  pl.Dot(vg);
  pl.Cross(vg);
  pg3 = ppg + pg;
  pg3 = ppg + pl;
  vg4 = pg - pl;
#endif

  // operator -
  XYPoint q1(1.,2.);
  XYPoint q2 = -1.* q1;
  XYVector v2 = -XYVector(q1);
  iret |= compare(XYVector(q2) == v2,true,"reflection");



  double angle = 1.;
  pg.Rotate(angle);
  iret |= compare(pg.Phi(), ppg.Phi() + angle );
  iret |= compare(pg.R(), ppg.R()  );

  GlobalXYZVector v3d(1,2,0);
  GlobalXYZVector vr3d = RotationZ(angle) * v3d;
  iret |= compare(pg.X(), vr3d.X() );
  iret |= compare(pg.Y(), vr3d.Y()  );



  if (iret == 0) std::cout << "\t\t\t\tOK\n";
  else std::cout << "\t\t\tFAILED\n";

  return iret;
}


// missing LV test

int testRotations3D() {

  int iret=0;
  std::cout << "testing 3D Rotations\t:\t";


  Rotation3D rot = RotationZ(1.) * RotationY(2) * RotationX(3);
  GlobalXYZVector vg(1.,2.,3);
  GlobalXYZPoint  pg(1.,2.,3);
  GlobalPolar3DVector vpg(vg);

  //  GlobalXYZVector vg2 = rot.operator()<Cartesian3D,GlobalCoordinateSystemTag, GlobalCoordinateSystemTag> (vg);
  GlobalXYZVector vg2 = rot(vg);
  iret |= compare(vg2.R(), vg.R(),"rot3D" );

  GlobalXYZPoint pg2 = rot(pg);
  iret |= compare(pg2.X(), vg2.X(),"x diff");
  iret |= compare(pg2.Y(), vg2.Y(),"y diff");
  iret |= compare(pg2.Z(), vg2.Z(),"z diff");


  Quaternion qrot(rot);

  pg2 = qrot(pg);
  iret |= compare(pg2.X(), vg2.X(),"x diff",10);
  iret |= compare(pg2.Y(), vg2.Y(),"y diff",10);
  iret |= compare(pg2.Z(), vg2.Z(),"z diff",10);

  GlobalPolar3DVector vpg2 = qrot * vpg;

  iret |= compare(vpg2.X(), vg2.X(),"x diff",10 );
  iret |= compare(vpg2.Y(), vg2.Y(),"y diff",10 );
  iret |= compare(vpg2.Z(), vg2.Z(),"z diff",10 );

  AxisAngle arot(rot);
  pg2 = arot(pg);
  iret |= compare(pg2.X(), vg2.X(),"x diff",10 );
  iret |= compare(pg2.Y(), vg2.Y(),"y diff",10 );
  iret |= compare(pg2.Z(), vg2.Z(),"z diff",10 );

  vpg2 = arot (vpg);
  iret |= compare(vpg2.X(), vg2.X(),"x diff",10 );
  iret |= compare(vpg2.Y(), vg2.Y(),"y diff",10 );
  iret |= compare(vpg2.Z(), vg2.Z(),"z diff",10 );

  EulerAngles erot(rot);

  vpg2 = erot (vpg);
  iret |= compare(vpg2.X(), vg2.X(),"x diff",10 );
  iret |= compare(vpg2.Y(), vg2.Y(),"y diff",10 );
  iret |= compare(vpg2.Z(), vg2.Z(),"z diff",10 );

  GlobalXYZVector vrx = RotationX(3) * vg;
  GlobalXYZVector vry = RotationY(2) * vrx;
  vpg2 = RotationZ(1) * GlobalPolar3DVector (vry);
  iret |= compare(vpg2.X(), vg2.X(),"x diff",10 );
  iret |= compare(vpg2.Y(), vg2.Y(),"y diff",10 );
  iret |= compare(vpg2.Z(), vg2.Z(),"z diff",10 );

  // test Get/SetComponents
  XYZVector v1,v2,v3;
  rot.GetComponents(v1,v2,v3);
  const Rotation3D rot2(v1,v2,v3);
  //rot2.SetComponents(v1,v2,v3);
  double r1[9],r2[9];
  rot.GetComponents(r1,r1+9);
  rot2.GetComponents(r2,r2+9);
  for (int i = 0; i < 9; ++i) {
     iret |= compare(r1[i],r2[i],"Get/SetComponents");
  }
  // operator == fails for numerical precision
  //iret |= compare( (rot2==rot),true,"Get/SetComponens");

  // test get/set with a matrix
#ifndef NO_SMATRIX
  SMatrix<double,3> mat;
  rot2.GetRotationMatrix(mat);
  rot.SetRotationMatrix(mat);
  iret |= compare( (rot2==rot),true,"Get/SetRotMatrix");
#endif

  //test inversion
  Rotation3D rotInv = rot.Inverse();
  rot.Invert(); // invert in place
  bool comp = (rotInv == rot );
  iret |= compare(comp,true,"inversion");

  // rotation and scaling of points
  XYZPoint q1(1.,2,3); double a = 3;
  XYZPoint qr1 =  rot( a * q1);
  XYZPoint qr2 =  a * rot( q1);
  iret |= compare(qr1.X(), qr2.X(),"x diff",10 );
  iret |= compare(qr1.Y(), qr2.Y(),"y diff",10 );
  iret |= compare(qr1.Z(), qr2.Z(),"z diff",10 );


  if (iret == 0) std::cout << "\tOK\n";
  else std::cout << "\t FAILED\n";

  return iret;
}


int testTransform3D() {


  std::cout << "testing 3D Transform\t:\t";
  int iret = 0;

  EulerAngles r(1.,2.,3.);

  GlobalPolar3DVector v(1.,2.,3.);
  GlobalXYZVector w(v);

  Transform3D t1( v );
  GlobalXYZPoint pg;
  t1.Transform( LocalXYZPoint(), pg );
  iret |= compare(pg.X(), v.X(),"x diff",10 );
  iret |= compare(pg.Y(), v.Y(),"y diff",10 );
  iret |= compare(pg.Z(), v.Z(),"z diff",10 );


  Transform3D t2( r, v );

  GlobalPolar3DVector vr = r.Inverse()*v;

//   std::cout << GlobalXYZVector(v) << std::endl;
//   std::cout << GlobalXYZVector(vr) << std::endl;
//   std::cout << GlobalXYZVector (r(v)) << std::endl;
//   std::cout << GlobalXYZVector (r(vr)) << std::endl;
//   std::cout << vr << std::endl;
//   std::cout << r(vr) << std::endl;



//   std::cout << r << std::endl;
//   std::cout << r.Inverse() << std::endl;
//   std::cout << r * r.Inverse() << std::endl;
//   std::cout << Rotation3D(r) * Rotation3D(r.Inverse()) << std::endl;
//   std::cout << Rotation3D(r) * Rotation3D(r).Inverse() << std::endl;


  // test Translation3D

  Translation3D tr1(v);
  Translation3D tr2(v.X(),v.Y(),v.Z());
// skip this test on 32 bits architecture. It might fail due to extended precision
#if !defined(__i386__)
   iret |= compare(tr1 ==tr2, 1,"eq transl",1 );
#else
  // add a dummy test to have the same outputfile for roottest
  // otherwise it will complain that the output is different !
  iret |= compare(0, 0,"dummy test",1 );
#endif

  Translation3D tr3 = tr1 * tr1.Inverse();
  GlobalPolar3DVector vp2 = tr3 * v;
  iret |= compare(vp2.X(), v.X(),"x diff",10 );
  iret |= compare(vp2.Y(), v.Y(),"y diff",10 );
  iret |= compare(vp2.Z(), v.Z(),"z diff",10 );


  Transform3D t2b = tr1 * Rotation3D(r);
  // this above fails on Windows - use a comparison with tolerance
  // 12 is size of Transform3D internal vector
  iret |= compare( IsEqual(t2,t2b,12), true,"eq1 transf",1 );
  //iret |= compare(t2 ==t2b, 1,"eq1 transf",1 );
  Transform3D t2c( r, tr1);
  iret |= compare( IsEqual(t2,t2c,12), true,"eq2 transf",1 );
  //iret |= compare(t2 ==t2c, 1,"eq2 transf",1 );


  Transform3D t3 =  Rotation3D(r) * Translation3D(vr);

  Rotation3D rrr;
  XYZVector vvv;
  t2b.GetDecomposition(rrr,vvv);
  iret |= compare(Rotation3D(r) ==rrr, 1,"eq transf rot",1 );
  iret |= compare( tr1.Vect() == vvv, 1,"eq transf vec",1 );
//   if (iret) std::cout << vvv << std::endl;
//   if (iret) std::cout << Translation3D(vr) << std::endl;

  Translation3D ttt;
  t2b.GetDecomposition(rrr,ttt);
  iret |= compare( tr1 == ttt, 1,"eq transf trans",1 );
//   if (iret) std::cout << ttt << std::endl;

  EulerAngles err2;
  GlobalPolar3DVector vvv2;
  t2b.GetDecomposition(err2,vvv2);
  iret |= compare( r.Phi(), err2.Phi(),"transf rot phi",4 );
  iret |= compare( r.Theta(), err2.Theta(),"transf rot theta",1 );
  iret |= compare( r.Psi(), err2.Psi(),"transf rot psi",1 );

//iret |= compare( v == vvv2, 1,"eq transf g vec",1 );
  iret |= compare( v.X(), vvv2.X(),"eq transf g vec",4 );
  iret |= compare( v.Y(), vvv2.Y(),"eq transf g vec",1 );
  iret |= compare( v.Z(), vvv2.Z(),"eq transf g vec",1 );

  // create from other rotations
  RotationZYX rzyx(r);
  Transform3D t4( rzyx);
  iret |= compare( t4.Rotation() == Rotation3D(rzyx), 1,"eq trans rzyx",1 );

  Transform3D trf2 = tr1 * r;
  iret |= compare( trf2 == t2b, 1,"trasl * e rot",1 );
  Transform3D trf3 = r * Translation3D(vr);
  //iret |= compare( trf3 == t3, 1,"e rot * transl",1 );
  // this above fails on i686-slc5-gcc43-opt - use a comparison with tolerance
  iret |= compare( IsEqual(trf3,t3,12), true,"e rot * transl",1 );

  Transform3D t5(rzyx, v);
  Transform3D trf5 = Translation3D(v) * rzyx;
  //iret |= compare( trf5 == t5, 1,"trasl * rzyx",1 );
  iret |= compare( IsEqual(trf5,t5,12), true,"trasl * rzyx",1 );

  Transform3D t6(rzyx, rzyx * Translation3D(v).Vect() );
  Transform3D trf6 = rzyx * Translation3D(v);
  iret |= compare( trf6 == t6, 1,"rzyx * transl",1 );
  if (iret) std::cout << t6 << "\n---\n" << trf6 << std::endl;



  Transform3D trf7 = t4 * Translation3D(v);
  //iret |= compare( trf7 == trf6, 1,"tranf * transl",1 );
  iret |= compare( IsEqual(trf7,trf6,12), true,"tranf * transl",1 );
  Transform3D trf8 = Translation3D(v) * t4;
  iret |= compare( trf8 == trf5, 1,"trans * transf",1 );

  Transform3D trf9 = Transform3D(v) * rzyx;
  iret |= compare( trf9 == trf5, 1,"tranf * rzyx",1 );
  Transform3D trf10 = rzyx * Transform3D(v);
  iret |= compare( trf10 == trf6, 1,"rzyx * transf",1 );
  Transform3D trf11 = Rotation3D(rzyx) * Transform3D(v);
  iret |= compare( trf11 == trf10, 1,"r3d * transf",1 );

  RotationZYX rrr2 = trf10.Rotation<RotationZYX>();
  //iret |= compare( rzyx == rrr2, 1,"gen Rotation()",1 );
  iret |= compare( rzyx.Phi() , rrr2.Phi(),"gen Rotation() Phi",1 );
  iret |= compare( rzyx.Theta(), rrr2.Theta(),"gen Rotation() Theta",10 );
  iret |= compare( rzyx.Psi(), rrr2.Psi(),"gen Rotation() Psi",1 );
  if (iret) std::cout << rzyx << "\n---\n" << rrr2 << std::endl;


  //std::cout << t2 << std::endl;
  //std::cout << t3 << std::endl;

  XYZPoint p1(-1.,2.,-3);

  XYZPoint p2 = t2 (p1);
  Polar3DPoint p3 = t3 ( Polar3DPoint(p1) );
  iret |= compare(p3.X(), p2.X(),"x diff",10 );
  iret |= compare(p3.Y(), p2.Y(),"y diff",10 );
  iret |= compare(p3.Z(), p2.Z(),"z diff",10 );

  GlobalXYZVector v1(1.,2.,3.);
  LocalXYZVector v2;  t2.Transform (v1, v2);
  GlobalPolar3DVector v3;  t3.Transform ( GlobalPolar3DVector(v1), v3 );

  iret |= compare(v3.X(), v2.X(),"x diff",10 );
  iret |= compare(v3.Y(), v2.Y(),"y diff",10 );
  iret |= compare(v3.Z(), v2.Z(),"z diff",10 );

  XYZPoint q1(1,2,3);
  XYZPoint q2(-1,-2,-3);
  XYZPoint q3 = q1 +  XYZVector(q2);
  //std::cout << q3 << std::endl;
  XYZPoint qt3 = t3(q3);
  //std::cout << qt3 << std::endl;
  XYZPoint qt1 = t3(q1);
  XYZVector vt2 = t3( XYZVector(q2) );
  XYZPoint qt4 = qt1 + vt2;
  iret |= compare(qt3.X(), qt4.X(),"x diff",10 );
  iret |= compare(qt3.Y(), qt4.Y(),"y diff",10 );
  iret |= compare(qt3.Z(),  qt4.Z(),"z diff",10 );
     //std::cout << qt4 << std::endl;

  // this fails
//  double a = 3;
  //XYZPoint q4 = a*q1;
//   std::cout << t3( a * q1) << std::endl;
//   std::cout << a * t3(q1) << std::endl;

  // test get/set with a matrix
#ifndef NO_SMATRIX
  SMatrix<double,3,4> mat;
  t3.GetTransformMatrix(mat);
  Transform3D t3b;  t3b.SetTransformMatrix(mat);
  iret |= compare( (t3==t3b),true,"Get/SetTransformMatrix");

  // test LR
  Boost b(0.2,0.4,0.8);
  LorentzRotation lr(b);
  SMatrix<double,4> mat4;
  lr.GetRotationMatrix(mat4);
  LorentzRotation lr2;  lr2.SetRotationMatrix(mat4);
  iret |= compare( (lr==lr2),true,"Get/SetLRotMatrix");
#endif

  {
     // testing ApplyInverse on Point
     XYZPoint point(1., -2., 3.);
     Transform3D tr(EulerAngles(10, -10, 10), XYZVector(10, -10, 0));

     // test that applying transformation + Inverse is identity
     auto r0 = tr.ApplyInverse(tr(point));
     auto r0_2 = tr.Inverse()(tr(point));

     iret |= compare(r0.X(), point.X(), "ApplyInverse/PointX", 100);
     iret |= compare(r0_2.X(), point.X(), "ApplyInverse/PointX", 100);
     iret |= compare(r0.Y(), point.Y(), "ApplyInverse/PointY", 10);
     iret |= compare(r0_2.Y(), point.Y(), "ApplyInverse/PointY", 10);
     iret |= compare(r0.Z(), point.Z(), "ApplyInverse/PointZ", 10);
     iret |= compare(r0_2.Z(), point.Z(), "ApplyInverse/PointZ", 10);

     // compare ApplyInverse with Inverse()
     auto r1 = tr.ApplyInverse(point);
     auto r2 = tr.Inverse()(point);
     iret |= compare(r1.X(), r2.X(), "ApplyInverse/Point", 10);
     iret |= compare(r1.Y(), r2.Y(), "ApplyInverse/Point", 10);
     iret |= compare(r1.Z(), r2.Z(), "ApplyInverse/Point", 10);
  }

  {
     // testing ApplyInverse on Vector
     XYZVector vector(1, -2., 3);
     Transform3D tr(EulerAngles(10, -10, 10), XYZVector(10, -10, 0));
     auto r1 = tr.ApplyInverse(vector);
     auto r2 = tr.Inverse()(vector);
     iret |= compare(r1.X(), r2.X(), "ApplyInverse/Vector", 10);
     iret |= compare(r1.Y(), r2.Y(), "ApplyInverse/Vector", 10);
     iret |= compare(r1.Z(), r2.Z(), "ApplyInverse/Vector", 10);
  }

  if (iret == 0) std::cout << "OK\n";
  else std::cout << "\t\t\tFAILED\n";

  return iret;
}


int testVectorUtil() {

  std::cout << "testing VectorUtil  \t:\t";
   int iret = 0;

   // test new perp functions
   XYZVector v(1.,2.,3.);

   XYZVector vx = ProjVector(v,XYZVector(3,0,0) );
   iret |= compare(vx.X(), v.X(),"x",1 );
   iret |= compare(vx.Y(), 0,"y",1 );
   iret |= compare(vx.Z(), 0,"z",1 );

   XYZVector vpx = PerpVector(v,XYZVector(2,0,0) );
   iret |= compare(vpx.X(), 0,"x",1 );
   iret |= compare(vpx.Y(), v.Y(),"y",1 );
   iret |= compare(vpx.Z(), v.Z(), "z",1 );

   double perpy = Perp(v, XYZVector(0,2,0) );
   iret |= compare(perpy, std::sqrt( v.Mag2() - v.y()*v.y()),"perpy" );

   XYZPoint  u(1,1,1);
   XYZPoint  un = u/u.R();


   XYZVector vl = ProjVector(v,u);
   XYZVector vl2 = XYZVector(un) * ( v.Dot(un ) );

   iret |= compare(vl.X(), vl2.X(),"x",1 );
   iret |= compare(vl.Y(), vl2.Y(),"y",1 );
   iret |= compare(vl.Z(), vl2.Z(),"z",1 );

   XYZVector vp = PerpVector(v,u);
   XYZVector vp2 = v - XYZVector ( un * ( v.Dot(un ) ) );
   iret |= compare(vp.X(), vp2.X(),"x",10 );
   iret |= compare(vp.Y(), vp2.Y(),"y",10 );
   iret |= compare(vp.Z(), vp2.Z(),"z",10 );

   double perp = Perp(v,u);
   iret |= compare(perp, vp.R(),"perp",1 );
   double perp2 = Perp2(v,u);
   iret |= compare(perp2, vp.Mag2(),"perp2",1 );

   // test rotations
   double angle = 1;
   XYZVector vr1 = RotateX(v,angle);
   XYZVector vr2 = RotationX(angle) * v;
   iret |= compare(vr1.Y(), vr2.Y(),"y",1 );
   iret |= compare(vr1.Z(), vr2.Z(),"z",1 );

   vr1 = RotateY(v,angle);
   vr2 = RotationY(angle) * v;
   iret |= compare(vr1.X(), vr2.X(),"x",1 );
   iret |= compare(vr1.Z(), vr2.Z(),"z",1 );

   vr1 = RotateZ(v,angle);
   vr2 = RotationZ(angle) * v;
   iret |= compare(vr1.X(), vr2.X(),"x",1 );
   iret |= compare(vr1.Y(), vr2.Y(),"y",1 );


  if (iret == 0) std::cout << "\t\t\tOK\n";
  else std::cout << "\t\t\t\t\t\tFAILED\n";
  return iret;

}

int testGenVector() {

  int iret = 0;
  iret |= testVector3D();
  iret |= testPoint3D();

  iret |= testVector2D();
  iret |= testPoint2D();

  iret |= testRotations3D();

  iret |= testTransform3D();

  iret |= testVectorUtil();


  if (iret !=0) std::cout << "\nTest GenVector FAILED!!!!!!!!!\n";
  return iret;

}

int main() {
   int ret = testGenVector();
   if (ret)  std::cerr << "test FAILED !!! " << std::endl;
   else   std::cout << "test OK " << std::endl;
   return ret;
}
