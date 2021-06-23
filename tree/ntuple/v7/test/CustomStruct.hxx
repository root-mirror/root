#include <string>
#include <vector>

/**
 * Used to test serialization and deserialization of classes in RNTuple with TClass
 */
struct CustomStruct {
   float a = 0.0;
   std::vector<float> v1;
   std::vector<std::vector<float>> v2;
   typedef std::string typedef_string;
   typedef_string s;
};
