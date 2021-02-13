#include <cctype>
#include <cstdio>
#include <cstdlib>

#include "strlcpy.h"

#include "NdbMTReacDesc.h"

ClassImp(NdbMTReacDesc);

/* ================ NdbMTReacDesc ================== */
NdbMTReacDesc::NdbMTReacDesc(const char *filename)
{
   Init(filename);
} // NdbMTReacDesc

/* ------- ~NdbMTReacDesc ------- */
NdbMTReacDesc::~NdbMTReacDesc()
{
   if (shrt) {
      for (int i=0; i<mt.GetSize(); i++) {
         if (shrt[i]) free(shrt[i]);
         if (desc[i]) free(desc[i]);
         if (comment[i]) free(comment[i]);
      }
      free(shrt);
      free(desc);
      free(comment);
   }
} // ~NdbMTReacDesc

/* ------ Init ------ */
void
NdbMTReacDesc::Init(const char *filename)
{
   shrt = NULL;
   desc = NULL;
   comment = NULL;

   FILE   *f;

   if ((f=fopen(filename,"r"))==NULL) {
      fprintf(stderr,"ERROR: NdbMTReacDesc::NdbMTReacDesc(%s) cannot open file.\n",filename);
      return;
   }

   /* ----- First read total number of MT's ------ */
   Int_t   N;

   fscanf(f,"%d",&N);
   if (N<=0) {
      fprintf(stderr,"ERROR: NdbMTReacDesc::NdbMTReacDesc(%s) error reading from file.\n",filename);
      fclose(f);
   }

   /* ----- Allocate memory ----- */
   mt.Set(N);
   shrt = (char **)malloc(N * sizeof(char*));
   desc = (char **)malloc(N * sizeof(char*));
   comment = (char **)malloc(N * sizeof(char*));

   for (int i=0; i<N; i++) {
      /* --- read mt number and short description --- */
      Int_t   mt_num, ch;
      char   str[512], str2[256];

      fscanf(f,"%d",&mt_num);

      mt.AddAt(mt_num,i);

      /* --- skip blanks --- */
      do {} while (isspace(ch=fgetc(f)));
      ungetc(ch,f);

      /* get rest of line */
      fgets(str,sizeof(str),f);

      /* --- strip last newline char --- */
      str[strlen(str)-1] = 0;

      shrt[i] = strdup(str);

      /* --- Read the description line until the empty line --- */
      str[0] = 0;
      while (1) {
         fgets(str2,sizeof(str2),f);
         if (str2[0] == '\n') break;
         strlcat(str,str2,512);
      }
      str[strlen(str)-1] = 0;
      desc[i] = strdup(str);

      /* --- Read the description line until the empty line --- */
      str[0] = 0;
      while (1) {
         fgets(str2,sizeof(str2),f);
         if (str2[0] == '\n') break;
         strlcat(str,str2,512);
      }
      if (str[0])
         str[strlen(str)-1] = 0;
      comment[i] = strdup(str);
   }
   fclose(f);
} // NdbMTReacDesc

/* ------- FindMT -------- */
Int_t
NdbMTReacDesc::FindMT( Int_t MT )
{
   /* Make a linear search */
   if (shrt){
      for (int i=0; i<mt.GetSize(); i++){
         if (mt[i] == MT){
            return i;
         } else {
            if (mt[i] > MT){
               break;
            }
         }
      }
   }
   return -1;

} // FindMT

/* ------- GetShort -------- */
char *
NdbMTReacDesc::GetShort(Int_t MT)
{
   Int_t   idx = FindMT(MT);
   if (idx<0)
      return NULL;

   return shrt[idx];
} // GetShort

/* ------- GetDescription -------- */
char *
NdbMTReacDesc::GetDescription(Int_t MT)
{
   Int_t   idx = FindMT(MT);
   if (idx<0)
      return NULL;

   return desc[idx];
} // GetDescription

/* ------- GetComment -------- */
char *
NdbMTReacDesc::GetComment(Int_t MT)
{
   Int_t   idx = FindMT(MT);
   if (idx<0)
      return NULL;

   return comment[idx];
} // GetComment
