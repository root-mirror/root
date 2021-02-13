
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"
#include "Math/EulerAngles.h"

#include "Math/Transform3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Boost.h"

#include "Math/Rotation3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Quaternion.h"
#include "Math/AxisAngle.h"

#include "Math/VectorUtil.h"

#include <iostream>
#include <iterator>
#include <vector>
#include <list>
#include <set>
#include <algorithm>

using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;

template <class Rot>
void printRot(const Rot & rot) {
   std::cout << "rot:  ( " ;
   std::ostream_iterator<double> oi(std::cout,"  ");
   rot.GetComponents(oi);
   std::cout << ") " << std::endl;
}
template<class V>
void printVec(const V & v ) {
   std::cout << "vec : ( " ;
   std::ostream_iterator<double> oi(std::cout,"  ");
   v.GetCoordinates(oi);
   std::cout << ") " << std::endl;
}
template<class L>
void printList(const L & l ) {
   std::cout << "list : ( " ;
   std::ostream_iterator<double> oi(std::cout,"  ");
   std::copy(l.begin(),l.end(),oi);
   std::cout << ") " << std::endl;
}

void testOstreamIter() {

   XYZVector v(1.,2.,3);
   printVec(v);
   XYZPoint p(v); printVec(p);
   XYZTVector q(1.,2,3,4); printVec(q);

   AxisAngle ar(v,4.); printRot(ar);
   EulerAngles er(ar); printRot(er);
   Quaternion qr(er); printRot(qr);
   Rotation3D rr(qr) ; printRot(rr);

   Transform3D t(rr,v); printRot(t);


   Boost b(0.3,0.4,0.8); printRot(b);
   LorentzRotation lr(rr); printRot(lr);
   LorentzRotation lr2(b); printRot(lr2);

}

void testListIter() {

   // test with lists
   double d[10] = {1,2,3,4,5,6,7,8,9,10};
   std::list<double>  inputData(d,d+3);

   XYZVector v; v.SetCoordinates(inputData.begin(), inputData.end() );
   std::list<double> data(3);
   v.GetCoordinates(data.begin(),data.end());
   printList(data);

   inputData = std::list<double>(d+3,d+6);
   XYZPoint p; p.SetCoordinates(inputData.begin(),inputData.end() );
   data.clear();
   data = std::list<double>(3);
   p.GetCoordinates(data.begin(), data.end() );
   printList(data);

   inputData = std::list<double>(d+6,d+10);
   XYZTVector q; q.SetCoordinates(inputData.begin(),inputData.end() );
   data.clear();
   data = std::list<double>(4);
   q.GetCoordinates(data.begin(), data.end() );
   printList(data);

   // test on rotations
   inputData = std::list<double>(d,d+3);
   EulerAngles re(inputData.begin(), inputData.end() );
   data = std::list<double>(3);
   re.GetComponents(data.begin(), data.end() );
   printList(data);

   inputData = std::list<double>(d,d+4);
   AxisAngle ra(inputData.begin(), inputData.end() );
   data = std::list<double>(4);
   ra.GetComponents(data.begin(), data.end() );
   printList(data);

   inputData = std::list<double>(d,d+4);
   Quaternion rq(inputData.begin(), inputData.end() );
   data = std::list<double>(4);
   rq.GetComponents(data.begin(), data.end() );
   printList(data);

   double b[3] = {0.3,0.4,0.8};
   inputData = std::list<double>(b,b+3);
   Boost bst(inputData.begin(), inputData.end() );
   data = std::list<double>(3);
   bst.GetComponents(data.begin(), data.end() );
   printList(data);

   Rotation3D tmp(ra);
   inputData = std::list<double>(9);
   tmp.GetComponents(inputData.begin());  printList(inputData);
   Rotation3D r(inputData.begin(),inputData.end());
   data = std::list<double>(9);
   r.GetComponents(data.begin(), data.end() );
   printList(data);


   Transform3D ttmp(r,XYZVector(1,2,3));
   inputData = std::list<double>(12);
   ttmp.GetComponents(inputData.begin());  printList(inputData);
   Transform3D t(inputData.begin(),inputData.end());
   data = std::list<double>(12);
   t.GetComponents(data.begin(), data.end() );
   printList(data);


   LorentzRotation ltmp(bst);
   inputData = std::list<double>(16);
   ltmp.GetComponents(inputData.begin());  printList(inputData);
   LorentzRotation lr(inputData.begin(),inputData.end());
   data = std::list<double>(16);
   lr.GetComponents(data.begin(), data.end() );
   printList(data);


}

int main() {

   testOstreamIter();
   testListIter();


}
