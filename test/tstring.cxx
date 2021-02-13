// @(#)root/test:$Id$
// Author: Fons Rademakers   19/08/96

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "TString.h"
#include "TRegexp.h"
#include "TPRegexp.h"


void Ok(int i, int b)
{
   std::cout << "Test " << i;
   if (b)
      std::cout << "   ok" << std::endl;
   else
      std::cout << "   NOT ok" << std::endl;
}



int main()
{
   // create base string
   TString s = "aap noot mies";
   std::cout << s << std::endl;

   // use different constructors and excercise +, += and == operators
   TString s1 = TString("aa") + "p ";
   s1 += TString("noot ") + "mi" + "es";
   std::cout << s1 << std::endl;
   Ok(1, s1==s);

   // use single char constructor and excercise again += and == operators
   TString s2 = 'a';
   s2 += "ap ";
   s2 += "noot ";
   s2 += "mies";
   std::cout << s2 << std::endl;
   Ok(2, s2==s);

   // get and set a subrange (via TSubString)
   TString s3 = s(4,9);
   s3(0,4) = "mama papa";
   std::cout << s3 << std::endl;
   Ok(3, s3=="mama papa mies");

   // create a regular expression, make search in string and replace
   // matched substring
   TRegexp re(" n.*t ");
   TString s4 = s;
   s4(re) = " pipo ";
   std::cout << s4 << std::endl;
   Ok(4, s4=="aap pipo mies");

   // use "const char*" operator and the Data() member function to access
   // the string as a const char*
   const char *a = (const char*)s;
   std::cout << a << std::endl;
   Ok(5, !strcmp(a, s.Data()));

   // return 13th character and replace it by 't', getting 14th character
   // should result in an error message since operator[] uses range checking
   TString s6 = s;
   s6[12] = 't';
   std::cout << s6 << std::endl;
   std::cout << "*** Error message ok, accessing intentionaly out-of-bounds" << std::endl;
   char b = s6[13];
   std::cout << b << std::endl;
   Ok(6, s6=="aap noot miet");

   // return 13th character and replace it by 'p', getting 14th character
   // should NOT result in an error message since operator() does not use
   // range checking
   TString s7 = s;
   s7(12) = 'p';
   std::cout << s7 << std::endl;
   char c = s7(13);
   std::cout << c << std::endl;
   Ok(7, s7=="aap noot miep");

   // use Append, Remove, Prepend, Replace and Insert
   TString s8 = s;
   s8.Append(" tante ");
   s8.Append(s7(9,4));
   std::cout << s8 << std::endl;

   s8.Remove(0,14);
   std::cout << s8 << std::endl;
   s8.Prepend("oom jan");
   std::cout << s8 << std::endl;
   s8.Insert(7," en ");
   std::cout << s8 << std::endl;
   s8.Replace(4,3,"jaap");
   std::cout << s8 << std::endl;
   Ok(8, s8=="oom jaap en tante miep");

   // use CompareTo to compare char * and TString
   TString s9 = s;
   Ok(9,  !s9.CompareTo(s));
   Ok(10, !s9.CompareTo("AAP NOOT MIES", TString::kIgnoreCase));

   // use Contains to find if "string" is contained
   Ok(11, s9.Contains("NooT", TString::kIgnoreCase));
   Ok(12, !s9.Contains("pipo"));

   // return the index to the first and last character 'c'
   Ok(13, s.First('o')==5);
   Ok(14, s.First('z')==kNPOS);
   Ok(15, s.Last('a')==1);
   Ok(16, s.Last('z')==kNPOS);

   // returns the index of the start of match
   Ok(17, s.Index("mies")==9);
   Ok(18, s.Index("aa",1)==kNPOS);

   // returns length of string
   Ok(19, s.Length()==13);

   // test IsNull
   TString s10;

   Ok(20, s10.IsNull());
   Ok(21, !s.IsNull());

   // test IsAscii
   s10 = "\xb9";

   Ok(22, !s10.IsAscii());
   Ok(23, s.IsAscii());

   // some excercises with the Perl Compatible Regular Expressions
   TString s11("Food is on the foo table.");
   TPRegexp("\\b(foo)\\s+(\\w+)").Substitute(s11, "round $2");
   Ok(24, s11=="Food is on the round table.");

   TString s12("pepernotenkoek");
   TPRegexp("peper(.*)koek").Substitute(s12, "wal$1boom");
   Ok(25, s12=="walnotenboom");

   TString s13("hihi haha");
   TPRegexp("^([^ ]*) *([^ ]*)").Substitute(s13, "$2 $1");
   Ok(26, s13=="haha hihi");

   Ok(27, TPRegexp("^(\\w+) *(\\w+)").Match(s13) == 3);

   // test Resize and Strip
   s9.Prepend("   ");
   std::cout << s9 << std::endl;

   s9.Resize(50);
   std::cout << s9 << "<<ends here" << std::endl;

   std::cout << s9.Strip(TString::kBoth) << "<<ends here" << std::endl;

   Printf("Using Print: %s (%d)\n", (const char*) s9, s9.Length());

   // test comparisons
   TString s20 = "abc";
   TString s21 = "abcd";
   TString s22 = "bcde";
   TString s23 = "Bcde";
   TString s24 = "";

   Ok(28, s20 < s21);
   Ok(29, s21 < s22);
   Ok(30, s23 < s22);
   Ok(31, s22.CompareTo(s23, TString::kIgnoreCase) == 0);
   Ok(32, (s23 < s24) != (s24 < s23));

   // test file access
   std::ifstream f("tstring.cxx");
   f.seekg(0, std::ios::end);
   Ssiz_t size = f.tellg();
   f.seekg(0, std::ios::beg);

   TString fs;
   fs.ReadFile(f);
   Ok(33, size == fs.Length());

   std::ifstream f2("tstring.cxx");
   fs.ReadLine(f2);  // read '// @(#)root/test:$Id: tstring.cxx 38977..."
   Ok(34, fs.Contains("root/test"));

   fs.ReadToken(f2);   // read '//'
   fs.ReadToken(f2);   // read 'Author:'
   Ok(35, fs == "Author:");
   
   // test tokenizer
   Ssiz_t from = 0;  
   TString act;
   TString option="-FF 1  -FF  2"; // with 2 spaces
   Int_t results[2] = {1,2};
   Int_t i = 0;
 
   while( option.Tokenize(act, from, " ")) {
      if (act == "-FF") {
         option.Tokenize(act, from, " ");
         Ok(36, act.Atoi() == results[i]);
         i++; 
      } 
   }    

   return 0;
}
