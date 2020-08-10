ROOTIO FILES

   A ROOTIO file consists of one "file header", one or more "data
records," and zero or more "free segments".  The file header is always
at the beginning of the file, while the data records and free segments
may in principle appear in any order.

   The file header is fixed length (64 bytes in the current
release.)  It's detailed format is given in "header.txt".

   A free segment is of variable length.  One free segment is a set
of contiguous bytes that are unused, and are available for ROOTIO to use
for new or resized data records.  The first four bytes of a a free
segment contain the negative of the number of bytes in the segment.  The
contents of the remainder of the free segment are irrelevant.

   A data record represents either user data or data used
internally by ROOTIO.  All data records have two portions, a "key"
portion and a "data" portion.  The key portion precedes the data
portion.  The format of the key portion is the same for all data.
(The key portion corresponds to a class TKey object).  The object name
and they key cycle are together sufficient to uniquely determine the
record within the file.  The file "dobject.txt" describes the format
of the data portion of a record for an object that uses the default
streamer.

DATA RECORD TYPES

"core" record types

   There are several types of data records used internally by
ROOTIO to support the storage of byte sequences.  These record types
are "TFile", "TDirectory", "KeysList", and "FreeSegments".  These types
can be considered to be in the "core" layer of ROOTIO.

   A file always contains exactly one "TFile" data record, which
(nearly?) always immediately follows the file header.  The TFile record
consists of either data pertaining to the file as a whole, or data
pertaining to the root "directory" of records in the file.  Its detailed
format is given in "tfile.txt".

   A file contains zero or more "TDirectory" data records, each
representing a subdirectory in the directory tree that has the "TFile"
record at its root.  The detailed format is given in "tdirectory.txt".

   A file contains one or more "KeysList" data records.  There is
one corresponding to the root directory (represented by the TFile
record), and one corresponding to each (non-empty) subdirectory in the
tree (each represented by a TDirectory record).  The data portion of
each KeysList record consists of the sequential keys of those data
records in that directory.  The detailed format is given in
"keyslist.txt".  Note that keys for "TFile", "KeysList", "FreeSegments",
and "StreamerInfo" data records never appear in the data portion of
a KeysList data record.

   A file always contains exactly one "FreeSegments" data record,
which keeps track of the free segments in the file.  Its detailed format
is given in "freesegments.txt".  Note that the list of free segments
contains one additional free segment that is not in the file itself,
because it represents the free space after the end of the file.

"streamer" layer record types

There is an additional data record type ("StreamerInfo") needed
internally to support the storage of self-identifying objects.  Its
detailed format is given in "streamerinfo.txt".  Note that the
StreamerInfo data record itself and the "core" data records described
above are not self-identifying objects.  A ROOTIO file contains
exactly one StreamerInfo record.  The use of the "StreamerInfo" record
is described under the "STREAMERINFO" heading below.

"pointer to persistent object" object types

There are three object types ("TProcessID", "TRef", and "TRefArray") used
internally to support pointers to persistent objects.  Their formats are
given in "tprocessid.txt", "tref.txt", and "trefarray.txt" respectively.
Of these three objects, only TProcessID objects necessarily comprise
a complete data record (a "TProcessID" record).  TRef and TRefArray
objects typically are data members of larger objects, and therefore are
only a part of the data portion of a record.  In addition, objects that
are referenced by such a pointer have an additional field in the base TObject.
See "tobject.txt".  A description of how these pointers work is given under
the "POINTERS TO PERSISTENT OBJECTS" heading below.

"application" layer record types

These are either user defined record types, or record types supplied
by ROOT that are not needed by ROOTIO. The format of such an object that
uses the default streamer is shown in "dobject.txt".

DATA COMPRESSION

The user can set the data compression level for new or modified data records
when creating or opening a file.  When an existing file is opened for update,
the compression level selected need not match that used previously.  The
compression level of existing records is not modified unless the record itself
is modified.

There are ten compression levels, 0-9, ranging from 0 (no compression) to 9
(maximum compression), with level 1 being the default.  The level chosen is
a tradeoff between disk space and compression performance.  The decompression
speed is independent of level.  Currently, in release 3.2.6, level 2 is not used.
If level 2 is selected, level 1 is used with no notification to the user.

The chosen compression level is not applied to the entire file.  The following
portions of the file are not compressed, regardless of the compression level
selected:
   1) the file header
   2) the KeysList data record
   3) the FreeSegments data record
   4) any data record (outside of a TTree) where the uncompressed size of
      the data portion is 256 bytes or less.
   5) the key portion of any data record
Furthermore, the data portion of the StreamerInfo data record is always
compressed at level 1 (if over 256 bytes uncompressed), regardless of the
compression level selected (even if no compression is selected).

The compression algorithm used is an in memory ZIP compression written for the
DELPHI collaboration at CERN.  Its author is E. Chernyaev (IHEP/Protvino).
The source code is internal to ROOTIO.

STREAMERINFO

The "StreamerInfo" data record is used by ROOTIO to support the storage of
self-identifying objects.  Its detailed format is given in "streamerinfo.txt".
A ROOTIO file contains exactly one StreamerInfo record, which is written to disk
automatically when a new or modified file is closed.

The StreamerInfo record is a list (ROOTIO class TList) of "StreamerInfo" objects
(ROOTIO class TStreamerInfo).  There is one StreamerInfo object in the list for
every class used in the file in a data record, other than a core layer record.
There is no streamerinfo object for a class used in a core layer record unless the
class is also used elsewhere in a data record.  When reading a self-identifying
object from a file, the system uses the StreamerInfo list to decompose the object
recursively into its simple data members.

Each streamerinfo object is an array of "streamer element" objects, each of which
describes a base class of the object or a (non-static and non-transient) data member
of the object.  If the base class or data member is itself a class, then there will
also be a streamerinfo object in the record for that class.  In this way, each
class is recursively decomposed into its atomic elements, each of which is a simple
type (e.g. "int").  A "long" or "unsigned long" member is always written
as an 8 byte quantity, even if it occupies only 4 bytes in memory.

A data member of a class is marked transient on the line of its declaration by a
comment beginning with "//!".  Such members are not written to disk, nor is there
any streamerinfo for such a member.

A data member that is a C++ pointer (not to be confused with "pointers to persistent
objects" described below) is never written to disk as a pointer value.  If it is a
pointer to an object, the object itself (or 0 (4 bytes) if the pointer value is NULL)
is written.  If the declaration line has a comment beginning with "//->", this indicates
that the pointer value will never be null, which allows a performance optimization.
Another optimization is that if two or more pointers pointing to the same object are
streamed in the same I/O operation, the object is written only once.  The remaining
pointers reference the object through a unique object identifier.  This saves space
and avoids the infinite loop that might otherwise arise if the directed graph of object
instance pointer references contains a cycle.

If a data member is a pointer to a simple type, the Streamer presumes it is an array,
with the dimension defined in a comment of the form "//[<length>]", where length is
either an integer constant or a variable that is an integer data member of the class.
If a variable is used, it must be defined ahead of its use or in a base class.

The above describes the function of the StreamerInfo record in decomposing a
self-identifying object if the user uses the streamer generated by "rootcint".
There are two reasons why a user may need to write a specialized streamer for a class.
One reason is that it may be necessary to execute some code before or after data is read
or written, for example, to initialize some non-persistent data members after the
persistent data is read.  In this case, the custom streamer can use the StreamerInfo record
to decompose a self-identifying object in the exact same manner as the generated
streamer would have done.  An example is given (for the Event class) in the Root User's
Guide (URL below) (Input/Output chapter, Streamers subchapter).  On the other hand, if
the user needs to write a streamer for a class that ROOT cannot handle, the user may need
to explicitly code the decomposition and composition of the object to its members.
In this case, the StreamerInfo for that class might not be used.  In any case, if the
composition/decomposition of the class is explicitly coded, the user should include
the byte count, class information, and version number of the class before the data on
disk as shown in "dobject.txt".

The special method used for streaming a TClonesArray is described in the TClonesArray
section below.

More information on the StreamerInfo record and its use is found in the Input/Output
chapter of the Root Users Guide:  http://root.cern/root/RootDoc.html

NOTE:  Some of the classes used internally in ROOTIO (e.g. TObject, TRef, TRefArray)
have explicitly coded (de)compositions, and do not use the information in the
StreamerInfo record to do the (de)composition.  In this case, the StreamerInfo for
the class may still be present in the StreamerInfo record, but may not match what is
actually written to disk for those objects.

POINTERS TO PERSISTENT OBJECTS

Information on how these work in memory can be found at:
http://root.cern/root/htmldoc/examples/Version302.news.html
These were introduced in release 3.02, so there is not yet a description in the current
Root Users Guide, which is for a version release 3.1.  Here we discuss only the information
on disk.

A ROOT file contains zero or more TProcessID records.  Each such record contains a globally
unique ID defining a given ROOT job that wrote a referenced object (see tprocessid.txt).
Each referenced object contains a "pidf" field referencing the corresponding TProcessID
record and an "fUniqueID" field uniquely identifying the referenced object among those
written by that process (see tobject.txt).  Similarly, every persistent reference to that
object (a TRef Object, see tref.txt) also contains "pidf" and "fUniqueID" fields with the
same value, thereby uniquely determining the referenced object (which need not even be in the
same file).  In the case of an array of references (a TRefArray object, see "trefarray.txt"),
there is one "pidf" value for the entire array, and a separate "fUniqueID" value for each
reference.  For further information, see the above URL.

SOME USEFUL CONTAINER CLASSES

1) TObjArray and TClonesArray

The TObjArray class can be used to support an array of objects.  The objects need not be of the
same type, but each object must be of a class type that inherits from TObject.  We have already
seen a specific example of the use of TObjArray, in the StreamerInfo record, where it is used
to hold an array of TStreamerElement objects, each of which is of a class inheriting from
TStreamerElement, which in turn inherits from TObject.

The TClonesArray class is a specialization of the TObjArray class for holding an array
of objects that are all of the same type.  The format of a TClonesArray object
is given in "tclonesarray.txt".

There are two great advantages in the use of TClonesArray over TObjArray when the objects
all will be of the same class:
  1) Memory for the objects will be allocated only once for the entire array, rather
     than the per-object allocation for TObjArray.  This can be done because all the
     objects are the same size.
  2) In the case of TObjArray, the stored objects are written sequentially. However,
     in a TClonesArray, by default, each object is split one level deep into its base
     class(es) and data members, and each of these members is written sequentially for
     all objects in the array before the next member is written.  This has two advantages:
      a) Greater compression can be achieved when similar data is consecutive.
      b) The object's data members can easily be split into different TTree branches
         (TTrees are discussed below).

2) TTree

A TTree is a highly specialized container class for efficient storage and retrieval of user data.
The use of TTrees is discussed in detail in the Trees chapter of the Root Users Guide:
http://root.cern/root/RootDoc.html
Here we discuss in particular how a TTree is stored in a ROOTIO file.

A TTree object is split into one or more branches (class TBranch), each of which may have its own
(sub)branches, recursively to any depth.  Each TBranch contains an array of zero or more leaves
(class TLeaf), each corresponding to a basic variable type or a class object that has not been split.
The TLeaf object does not actually contain variable values, only information about the variables.
The actual data on each branch is physically stored in basket objects (class TBasket).  The user
can set the basket size on a per TBranch basis.  The default basket size is 32000 bytes.
This should be viewed as an approximate number.

There is one TTree data record per file for each tree in the file, corresponding to a TTree
class object.  The TTree class object recursively contains TBranch objects, each of which
contains an array of TBasket objects to hold its data.

However, the TTree data record does not necessarily contain the entire TTree object.  For each
branch, exactly one TBasket object is contained in the TTree data record.  If the data on a
given branch fits in one basket, then all the data for that branch will be in the TTree record
itself.  Otherwise, there will be a separate TBasket data record for each additional basket used on
the branch, each containing a TBasket object containing user data.

By default, the additional TBasket data records are stored in the same file as that of the
corresponding TTree data record.  However, the user may specify a separate file for a given
branch.  If the data for that branch fits into one basket, this option has no effect.  Otherwise,
the additional TBasket records are written into the specified file, rather than the file containing
the TTree data record itself.  In this case, a TBranch data record for the specified branch is also
written to the specified file, containing the TBranch object for the specified branch.

The file "ttree.txt" shows the streamer information for the TTree, TBranch, TLeaf, and some related
classes, together with some additional commentary.  For writing to a ROOTIO file, the streamers
for these three classes act exactly as those of default generated streamers, except that, if the
user has specified a separate file for a branch, the TBranch streamer also writes the TBranch object
as a keyed data record to the specified file.

There is no streamer information for the TBasket class.  The custom written TBasket streamer
internally handles the packing of data into fixed size TBasket objects.
