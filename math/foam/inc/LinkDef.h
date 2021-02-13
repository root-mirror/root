#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TFoamIntegrand+;
#pragma link C++ class TFoamMaxwt+;
#pragma link C++ class TFoamVect+;
#pragma link C++ class TFoamCell+;
#pragma link C++ class TFoam+;
#pragma link C++ class TFoamSampler+;
#pragma read sourceClass="TFoam" targetClass="TFoam" version="[1]" \
  source="Int_t fNCells; TFoamCell **fCells; TRefArray *fCellsAct" target="fNCells,fCells,fCellsAct"\
  include="TRefArray.h" \
  code="{fNCells = onfile.fNCells; \
         fCells = onfile.fCells; \
         onfile.fCells = nullptr; \
         fCellsAct.clear(); \
         for (Int_t i=0; i < onfile.fCellsAct->GetEntries(); ++i) { \
            const TObject* cellp = onfile.fCellsAct->At(i); \
            for (Int_t j=0; j < fNCells; ++j) { \
               if (cellp == fCells[j]) { \
                 fCellsAct.push_back(j); \
                 break; \
               } \
            } \
         }}";
#endif
