# Pretty printers for gdb
# \author: Stephan Hageboeck, CERN
# These pretty printers will make ROOT objects more readable when printed in gdb.
# If the pretty-printed output is not sufficient, one can always use "print /r <object>"
# for raw printing.
#
# When a debug build is used, they will be installed next to the ROOT libraries.
# gdb will load them automatically if the auto-load-safe-path is set to ROOT's library directory.
# For this, one has to add `add-auto-load-safe-path <ROOT lib dir>` to .gdbinit
#
# If loaded successfully, typing `info pretty-printer` at the gdb prompt should list the
# printers registered at the end of this file.

import gdb


class RooCollectionPrinter(object):
   "Print a RooAbsCollection"

   def __init__(self, val):
      self.val = val
      self.viz = gdb.default_visualizer(self.val['_list'])
      
   def to_string(self):
      ret = "{" + str(self.val.dynamic_type) + " " + str(self.val['_name']) +": "
      try:
         for name,val in self.viz.children():
            itemName = val.referenced_value()['fName']
            ret += str(itemName) + ","
      except:
         ret += "<exception " + str(sys.exc_info()[0]) + ">,"
            
      ret += "}"
      return ret
      
   def children(self):
      for name,val in self.viz.children():
         try:
            itemName = val.referenced_value()['fName']
            key = name + " " + str(val.address) + " (" + str(val.dynamic_type) +") " + str(itemName)
            yield key, val.referenced_value()
         except:
#            print("<exception " + str(sys.exc_info()[0]) + ">,")
            raise

   def display_hint(self):
      return 'RooAbsCollection printer'


        
class RooSpanPrinter(object):
   "Print a RooSpan"

   def __init__(self, val):
      self.val = val
      
   def to_string(self):
      return "span of length " + str(self.val['_span']['length_'])
      
   def children(self):
      length = self.val['_span']['length_']
      values = ""
      for i in range(0, min(length, 10)):
         values += ' ' + str((self.val['_span']['data_']+i).dereference())
      yield 'Values', values + '...'

   def display_hint(self):
      return 'RooSpan printer'


        
class RooAbsArgPrinter(object):
   "Print a RooAbsArg"

   def __init__(self, val):
      self.val = val
      
   def children(self):
      for name,item in self.val.fields():
         yield name, item

   def to_string(self):
      ret += str(self.val.address) + " " + str(self.val.dynamic_type)
      itemName = self.val['fName']
      ret += " = { <fName> = {" + str(itemName) + "} }"
      return ret

   def display_hint(self):
      return 'RooAbsArg printer'


      
class RooSTLRefCountListPrinter(object):
   "Print ref count lists"
   
   def __init__(self, val):
      self.val = val
      
   def to_string(self):
      ret = "{"
      viz = gdb.default_visualizer(self.val['_storage'])
      vizRC = gdb.default_visualizer(self.val['_refCount'])
      for (name,val), (name2,val2) in zip(viz.children(), vizRC.children()):
         ret += str(val['fName']) + ": " + str(val2) + ", "
      return ret + "}"
      


class NonePrinter(object):
   "Prevent printing an object"

   def __init__(self, val):
      self.val = val

   def to_string(self):
      return ""

   def display_hint(self):
      return 'Disables printing'




def build_pretty_printer():
   pp = gdb.printing.RegexpCollectionPrettyPrinter("libRooFitCore")
   pp.add_printer('RooSpan', '^RooSpan.*$', RooSpanPrinter)
   pp.add_printer('Collections', '^Roo(AbsCollection|ArgList|ArgSet)$', RooCollectionPrinter)
   pp.add_printer('RooSTLRefCountList', '^RooSTLRefCountList.*$', RooSTLRefCountListPrinter)
   pp.add_printer('RooPrintable', '^RooPrintable$', NonePrinter)

   return pp


gdb.printing.register_pretty_printer(gdb.current_objfile(),
    build_pretty_printer())
