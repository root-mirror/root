/* @(#)root/gpad:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ROOT::Experimental::RStyle+;
#pragma link C++ class ROOT::Experimental::RStyle::Block_t+;

#pragma link C++ class ROOT::Experimental::RDrawable+;
#pragma link C++ class ROOT::Experimental::RDisplayItem+;
#pragma link C++ class ROOT::Experimental::RDrawableDisplayItem+;
#pragma link C++ class ROOT::Experimental::RIndirectDisplayItem+;
#pragma link C++ class ROOT::Experimental::TObjectDisplayItem+;
#pragma link C++ class ROOT::Experimental::RDrawableReply+;
#pragma link C++ class ROOT::Experimental::RDrawableRequest+;
#pragma link C++ class ROOT::Experimental::RDrawableExecRequest+;

#pragma link C++ class ROOT::Experimental::Internal::RIOSharedBase+;
#pragma link C++ class ROOT::Experimental::Internal::RIOShared<TObject>+;
#pragma link C++ class ROOT::Experimental::Internal::RIOShared<ROOT::Experimental::RDrawable>+;

#pragma link C++ class ROOT::Experimental::Detail::RMenuItem+;
#pragma link C++ class std::vector<ROOT::Experimental::Detail::RMenuItem*>+;
#pragma link C++ class ROOT::Experimental::Detail::RCheckedMenuItem+;
#pragma link C++ class ROOT::Experimental::Detail::RMenuArgument+;
#pragma link C++ class std::vector<ROOT::Experimental::Detail::RMenuArgument>+;
#pragma link C++ class ROOT::Experimental::Detail::RArgsMenuItem+;
#pragma link C++ class ROOT::Experimental::RMenuItems+;
#pragma link C++ class ROOT::Experimental::RDrawableMenuRequest+;
#pragma link C++ class ROOT::Experimental::TObjectDrawable+;
#pragma link C++ class ROOT::Experimental::RPadUserAxisBase+;
#pragma link C++ class ROOT::Experimental::RPadCartesianUserAxis+;
#pragma link C++ class ROOT::Experimental::RPadExtent+;
#pragma link C++ class ROOT::Experimental::RPadPos+;
#pragma link C++ class ROOT::Experimental::RChangeAttrRequest+;
#pragma link C++ class ROOT::Experimental::RPadBase+;
#pragma link C++ class ROOT::Experimental::RPad+;
#pragma link C++ class ROOT::Experimental::RCanvas+;
#pragma link C++ class ROOT::Experimental::RFrame+;
#pragma link C++ class ROOT::Experimental::RFrame::RUserRanges+;
#pragma link C++ class ROOT::Experimental::RFrame::RZoomRequest+;
#pragma link C++ class ROOT::Experimental::RPave+;
#pragma link C++ class ROOT::Experimental::RPadLength+;
#pragma link C++ class ROOT::Experimental::RPadLength::Pixel+;
#pragma link C++ class ROOT::Experimental::RPadLength::Normal+;
#pragma link C++ class ROOT::Experimental::RPadLength::User+;
#pragma link C++ class ROOT::Experimental::RPadLength::CoordSysBase<ROOT::Experimental::RPadLength::Pixel>+;
#pragma link C++ class ROOT::Experimental::RPadLength::CoordSysBase<ROOT::Experimental::RPadLength::Normal>+;
#pragma link C++ class ROOT::Experimental::RPadLength::CoordSysBase<ROOT::Experimental::RPadLength::User>+;
#pragma link C++ class ROOT::Experimental::RPalette+;
#pragma link C++ struct ROOT::Experimental::RPalette::OrdinalAndColor+;
#pragma link C++ class ROOT::Experimental::RPaletteDrawable+;
#pragma link C++ class ROOT::Experimental::RPadBaseDisplayItem+;
#pragma link C++ class ROOT::Experimental::RPadDisplayItem+;
#pragma link C++ class ROOT::Experimental::RCanvasDisplayItem+;

#pragma read sourceClass="ROOT::Experimental::RCanvas" targetClass="ROOT::Experimental::RCanvas" source="" target="" code="{ newObj->ResolveSharedPtrs() ; }"

#pragma link C++ class ROOT::Experimental::RAttrMap+;
#pragma link C++ class ROOT::Experimental::RAttrMap::Value_t+;
#pragma link C++ class ROOT::Experimental::RAttrMap::NoValue_t+;
#pragma link C++ class ROOT::Experimental::RAttrMap::BoolValue_t+;
#pragma link C++ class ROOT::Experimental::RAttrMap::IntValue_t+;
#pragma link C++ class ROOT::Experimental::RAttrMap::DoubleValue_t+;
#pragma link C++ class ROOT::Experimental::RAttrMap::StringValue_t+;
#pragma link C++ class ROOT::Experimental::RAttrBase+;
#pragma link C++ class ROOT::Experimental::RAttrValue<bool>+;
#pragma link C++ class ROOT::Experimental::RAttrValue<int>+;
#pragma link C++ class ROOT::Experimental::RAttrValue<double>+;
#pragma link C++ class ROOT::Experimental::RAttrValue<std::string>+;
#pragma link C++ class ROOT::Experimental::RAttrValue<ROOT::Experimental::RPadLength>+;

#pragma link C++ class ROOT::Experimental::RColor+;
#pragma link C++ class ROOT::Experimental::RAttrColor+;
#pragma link C++ class ROOT::Experimental::RAttrFill+;
#pragma link C++ class ROOT::Experimental::RAttrLine+;
#pragma link C++ class ROOT::Experimental::RAttrBox+;
#pragma link C++ class ROOT::Experimental::RAttrMarker+;
#pragma link C++ class ROOT::Experimental::RAttrText+;
#pragma link C++ class ROOT::Experimental::RAttrAxis+;
#pragma link C++ class ROOT::Experimental::RAttrMargins+;

#pragma link C++ class ROOT::Experimental::RAxisDrawableBase+;
#pragma link C++ class ROOT::Experimental::RAxisDrawable+;
#pragma link C++ class ROOT::Experimental::RAxisLabelsDrawable+;

#endif
