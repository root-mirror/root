// @(#)root/mpi / Author: Omar.Zapata@cern.ch 2016-2017 http://oproject.org
#ifndef ROOT_Mpi_TCommunicator
#define ROOT_Mpi_TCommunicator

#ifndef ROOT_Mpi_Globals
#include<Mpi/Globals.h>
#endif

#ifndef ROOT_Mpi_TMpiMessage
#include <Mpi/TMpiMessage.h>
#endif

#ifndef ROOT_Mpi_TStatus
#include <Mpi/TStatus.h>
#endif

#ifndef ROOT_Mpi_TRequest
#include<Mpi/TRequest.h>
#endif


#include<memory>

/**
 * @namespace ROOT::Mpi
 * namespace associated RMpi package for ROOT.
 * @defgroup Mpi Message Passing Interface
 */


namespace ROOT {

   namespace Mpi {

      class TMpiMessage;

      /**
       *      \class TCommunicator
       *         Class for communicator, with this class you can to communicate the processes using messages,
       *      the messages can be any serializable object supported by ROOT like object from standart c++ libraries or
       *      objects that inherits from TObject.
       *
       *      You can to create your own classes and communicate it just creating its dictionaries
       *         \ingroup Mpi
       */

      class TCommunicator: public TObject {
      private:
         MPI_Comm fComm;           //! Raw communicator
         Int_t fMainProcess;    // Rank used like a main process
      public:
         /**
          *         Copy constructor for communicator
          *              \param comm other TCommunicator object
          */
         TCommunicator(const TCommunicator &comm);
         TCommunicator(const MPI_Comm &comm = MPI_COMM_WORLD);
         ~TCommunicator();

         TCommunicator &operator=(const MPI_Comm &comm)
         {
            fComm = comm;
            return *this;
         }


         /**
          *         Method to get the current rank or process id
          *              \return integer with the rank value
          */
         inline Int_t GetRank() const
         {
            Int_t rank;
            MPI_Comm_rank(fComm, &rank);
            return rank;
         }

         /**
          *         Method to get the total number of ranks or processes
          *              \return integer with the number of processes
          */
         inline Int_t GetSize() const
         {
            Int_t size;
            MPI_Comm_size(fComm, &size);
            return size;
         }

         /**
          *         Method to know if the current rank us the main process
          *              \return boolean true if it is the main rank
          */
         inline Bool_t IsMainProcess() const
         {
            return GetRank() == fMainProcess;
         }

         /**
          *         Method to set the main process rank
          *              \param Int_t main process rank number
          */
         inline void SetMainProcess(Int_t p)
         {
            fMainProcess = p;
         }

         /**
          *         Method to get the main process id
          *              \return integer with the main rank
          */
         inline Int_t GetMainProcess() const
         {
            return fMainProcess;
         }

         /**
          *         Method to abort  processes
          *              \param integer with error code
          */
         inline void Abort(Int_t err) const
         {
            MPI_Abort(fComm, err);
         }

         /**
          *         Method for synchronization between MPI processes in a communicator
          */
         virtual void Barrier() const;

         /**
          *         Method for synchronization between MPI nonblocking processes in a communicator
          *         \param req request object
          */
         virtual void IBarrier(TRequest &req) const;

         /**
          *            Nonblocking test for a message. Operations  allow checking of incoming messages without actual receipt of them.
          *              \param source Source rank or ROOT::Mpi::ANY_SOURCE (integer).
          *              \param tag Tag value or ROOT::Mpi::ANY_TAG (integer).
          *              \param status TStatus object with extra information.
          *              \return boolean true if the probe if ok
          */
         virtual Bool_t IProbe(Int_t source, Int_t tag, TStatus &status) const;

         /**
          *            Nonblocking test for a message. Operations  allow checking of incoming messages without actual receipt of them.
          *              \param source Source rank or ROOT::Mpi::ANY_SOURCE (integer).
          *              \param tag Tag value or ROOT::Mpi::ANY_TAG (integer).
          *              \return boolean true if the probe if ok
          */
         virtual Bool_t IProbe(Int_t source, Int_t tag) const;

         /**
          *            Test for a message. Operations  allow checking of incoming messages without actual receipt of them.
          *              \param source Source rank or ROOT::Mpi::ANY_SOURCE (integer).
          *              \param tag Tag value or ROOT::Mpi::ANY_TAG (integer).
          *              \param status TStatus object with extra information.
          *              \return boolean true if the probe if ok
          */
         virtual void Probe(Int_t source, Int_t tag, TStatus &status) const;

         /**
          *            Test for a message. Operations  allow checking of incoming messages without actual receipt of them.
          *              \param source Source rank or ROOT::Mpi::ANY_SOURCE (integer).
          *              \param tag Tag value or ROOT::Mpi::ANY_TAG (integer).
          *              \return boolean true if the probe if ok
          */
         virtual void Probe(Int_t source, Int_t tag) const;

         /**
          *         Method to send a message for p2p communication
          *              \param var any selializable object
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> void Send(const Type &var, Int_t dest, Int_t tag) const;

         /**
          *         Method to receive a message for p2p communication
          *              \param var any selializable object reference to receive the message
          *              \param source id with the origin(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type>  void Recv(Type &var, Int_t source, Int_t tag) const; //must be changed by ROOOT::Mpi::TStatus& Recv(...)

         /**
          *            Starts a standard-mode, nonblocking send.
          *              \param var any selializable object
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> TRequest ISend(const Type &var, Int_t dest, Int_t tag);

         /**
          *         Starts a nonblocking synchronous send
          *              \param var any selializable object
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> TRequest ISsend(const Type &var, Int_t dest, Int_t tag);
         /**
          *         Starts a ready-mode nonblocking send.
          *              \param var any selializable object
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> TRequest IRsend(const Type &var, Int_t dest, Int_t tag);

         /**
          *         Method to receive a message from nonblocking send (ISend, ISsend, IRsend)
          *         to receive the object you need to call the methods Complete() and Wait()
          *         TGrequest req=comm.IRecv(..);
          *         req.Complete();
          *         req.Wait();
          *
          *              \param var any selializable object reference to receive the message
          *              \param source id with the origin(Rank/Process) of the message
          *              \param tag id of the message
          *              \return TGrequest object.
          */
         template<class Type> TRequest IRecv(Type &var, Int_t source, Int_t tag) const;


         /**
          *          Broadcasts a message from the process with rank root to all other processes of the group.
          *              \param var any selializable object reference to send/receive the message
          *              \param root id of the main message where message was sent
          */
         template<class Type> void Bcast(Type &var, Int_t root) const;

         /**
          *          Broadcasts a message from the process with rank root to all other processes of the group.
          *              \param var any selializable object reference to send/receive the message
          *              \param root id of the main message where message was sent
          *              \return TGrequest obj
          */
         template<class Type> TRequest IBcast(Type &var, Int_t root) const;

         /**
          *          Sends data from one task to all tasks in a group.
          *              \param in_vars any selializable object vector reference to send the message
          *              \param incount Number of elements in receive in \p in_vars
          *              \param out_var any selializable object vector reference to receive the message
          *              \param outcount Number of elements in receive in \p out_vars
          *              \param root id of the main message where message was sent
          *              \return TGrequest obj
          */
         template<class Type> void Scatter(const Type *in_vars, Int_t incount, Type *out_vars, Int_t outcount, Int_t root) const;

         /**
          *          Each process (root process included) sends the contents of its send buffer to the root process.
         *          The root process receives the messages and stores them in rank order.
         *          The outcome is as if each of the n processes in the group (including the root process)
          *              \param in_vars any selializable object vector reference to send the message
          *              \param incount Number of elements in receive in \p in_vars
          *              \param out_var any selializable object vector reference to receive the message
          *              \param outcount Number of elements in receive in \p out_vars
          *              \param root id of the main message where message was sent
          *              \return TGrequest obj
          */
         template<class Type> void Gather(const Type *in_vars, Int_t incount, Type *out_vars, Int_t outcount, Int_t root) const;

         /**
          *         Method to apply reduce operation using binary tree reduction.
          *                        \param in_var variable to eval in the reduce operation
          *                        \param out_var variable to receive the variable operation
          *                        \param op function the perform operation
          *                        \param root id of the main process where the result was received
          */
         template<class Type> void Reduce(const Type &in_var, Type &out_var, Op<Type> (*opf)(), Int_t root) const;


         /////////////////////////////////
         //methods with arrary arguments//
         /////////////////////////////////

         /**
          *         Method to send a message for p2p communication
          *              \param vars array of any selializable objects
          *              \param count number of elements in array \p vars
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> void Send(const Type *vars, Int_t count, Int_t dest, Int_t tag) const;

         /**
          *         Method to receive a message for p2p communication
          *              \param vars array of any selializable objects
          *              \param count number of elements in array \p vars
          *              \param source id with the origin(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type>  void Recv(Type *vars, Int_t count, Int_t source, Int_t tag) const;

         /**
          *            Starts a standard-mode, nonblocking send.
          *              \param vars any selializable object
          *              \param count number of elements in array \p vars
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> TRequest ISend(const Type *vars, Int_t count, Int_t dest, Int_t tag);

         /**
          *         Starts a nonblocking synchronous send
          *              \param vars any selializable object
          *              \param count number of elements in array \p vars
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> TRequest ISsend(const Type *var, Int_t count, Int_t dest, Int_t tag);

         /**
          *         Starts a ready-mode nonblocking send.
          *              \param vars any selializable object
          *              \param count number of elements in array \p vars
          *              \param dest id with the destination(Rank/Process) of the message
          *              \param tag id of the message
          */
         template<class Type> TRequest IRsend(const Type *vars, Int_t count, Int_t dest, Int_t tag);

         /**
          *         Method to receive a message from nonblocking send (ISend, ISsend, IRsend)
          *         to receive the object you need to call the methods Complete() and Wait()
          *         TGrequest req=comm.IRecv(..);
          *         req.Complete();
          *         req.Wait();
          *
          *              \param vars any selializable object reference to receive the message
          *              \param count number of elements in array \p vars
          *              \param source id with the origin(Rank/Process) of the message
          *              \param tag id of the message
          *              \return TGrequest object.
          */
         template<class Type> TRequest IRecv(Type *vars, Int_t count, Int_t source, Int_t tag) const;

         /**
          *          Broadcasts a message from the process with rank root to all other processes of the group.
          *              \param vars any selializable object reference to send/receive the message
          *              \param count number of elements in array \p vars
          *              \param root id of the main message where message was sent
          *              \return TGrequest obj
          */
         template<class Type> TRequest IBcast(Type *vars, Int_t count, Int_t root) const;

         /**
          *          Broadcasts a message from the process with rank root to all other processes of the group.
          *              \param vars any selializable objects pointer to send/receive the message
          *              \param count Number of elements to broadcast in \p in_vars
          *              \param root id of the main message where message was sent
          */
         template<class Type> void Bcast(Type *vars, Int_t count, Int_t root) const;

         /**
          *         Method to apply reduce operation over and array of elements using binary tree reduction.
          *                        \param in_vars variable to eval in the reduce operation
          *                        \param out_vars variable to receive the variable operation
          *                        \param count Number of elements to reduce in \p in_vars and \p out_vars
          *                        \param op function the perform operation
          *                        \param root id of the main process where the result was received
          */
         template<class Type> void Reduce(const Type *in_vars, Type *out_vars, Int_t count, Op<Type> (*opf)(), Int_t root) const;

         /**
          *         static method to serialize objects. used in the multiple communication schemas.
          *         \param buffer double pointer to Char_t to save the serialized data
          *         \param size   reference to Int_t with the size of the buffer with serialized data
          *         \param vars   any selializable object
          *         \param count  number of elements to serialize in \p in_vars
          *         \param comm   communicator object
          *         \param dest   (optional) destination of the serialized information, must be the same unserializing
          *         \param source (optional) source of the serialized information, must be the same unserializing
          *         \param root   (optional) root of collective operation, must be the same unserializing
          */
         template<class T> static void Serialize(Char_t **buffer, Int_t &size, const T *vars, Int_t count, const TCommunicator *comm, Int_t dest = 0, Int_t source = 0, Int_t tag = 0, Int_t root = 0)
         {
            std::vector<TMpiMessageInfo> msgis(count);
            for (auto i = 0; i < count; i++) {
               TMpiMessage msg;
               msg.WriteObject(vars[i]);
               auto mbuffer = msg.Buffer();
               auto msize   = msg.BufferSize();
               if (mbuffer == NULL) {
                  comm->Error(__FUNCTION__, "Error serializing object type %s \n", ROOT_MPI_TYPE_NAME(T));
                  comm->Abort(ERR_BUFFER);
               }
               TMpiMessageInfo msgi(mbuffer, msize);
               msgi.SetSource(comm->GetRank());
               msgi.SetDestination(dest);
               msgi.SetSource(source);
               msgi.SetRoot(root);
               msgi.SetTag(tag);
               msgi.SetDataTypeName(ROOT_MPI_TYPE_NAME(T));
               msgis[i] = msgi;
            }
            TMpiMessage msg;
            msg.WriteObject(msgis);
            auto ibuffer = msg.Buffer();
            size = msg.BufferSize();
            *buffer = new Char_t[size];
            if (ibuffer == NULL) {
               comm->Error(__FUNCTION__, "Error serializing object type %s \n", ROOT_MPI_TYPE_NAME(msgis));
               comm->Abort(ERR_BUFFER);
            }
            memcpy(*buffer, ibuffer, size);
         }

         /**
          *         static method to unserialize objects. used in the multiple communication schemas.
          *         \param buffer pointer to Char_t to read the serialized data
          *         \param size   size of the buffer with serialized data
          *         \param vars   any selializable object
          *         \param count  number of elements to serialize in \p in_vars
          *         \param comm   communicator object
          *         \param dest   (optional) destination of the serialized information, must be the same serializing
          *         \param source (optional) source of the serialized information, must be the same serializing
          *         \param root   (optional) root of collective operation, must be the same serializing
          */
         template<class T> static  void Unserialize(Char_t *buffer, Int_t size, T *vars, Int_t count, const TCommunicator *comm, Int_t dest = 0, Int_t source = 0, Int_t tag = 0, Int_t root = 0)
         {
            TMpiMessage msg(buffer, size);
            auto cl = gROOT->GetClass(typeid(std::vector<TMpiMessageInfo>));
            auto msgis = (std::vector<TMpiMessageInfo> *)msg.ReadObjectAny(cl);
            if (msgis == NULL) {
               comm->Error(__FUNCTION__, "Error unserializing object type %s \n", cl->GetName());
               comm->Abort(ERR_BUFFER);
            }

            if (msgis->data()->GetDataTypeName() != ROOT_MPI_TYPE_NAME(T)) {
               comm->Error(__FUNCTION__, "Error unserializing objects type %s where objects are %s \n", ROOT_MPI_TYPE_NAME(T), msgis->data()->GetDataTypeName().Data());
               comm->Abort(ERR_TYPE);
            }

            ROOT_MPI_ASSERT(msgis->data()->GetDestination() == dest, comm)
            ROOT_MPI_ASSERT(msgis->data()->GetSource() == source, comm)
            ROOT_MPI_ASSERT(msgis->data()->GetRoot() == root, comm)
            ROOT_MPI_ASSERT(msgis->data()->GetTag() == tag, comm)

            for (auto i = 0; i < count; i++) {
               //passing information from TMpiMessageInfo to TMpiMessage
               auto isize = msgis->data()[i].GetBufferSize();
               Char_t *ibuffer = new Char_t[isize];//this memory dies when the unserialized object dies
               memcpy(ibuffer, msgis->data()[i].GetBuffer(), isize);
               TMpiMessage vmsg(ibuffer, isize);
               auto vcl = gROOT->GetClass(typeid(T));
               auto vobj_tmp = vmsg.ReadObjectAny(vcl);
               if (vobj_tmp == NULL) {
                  comm->Error(__FUNCTION__, "Error unserializing objects type %s \n", vcl->GetName());
                  comm->Abort(ERR_BUFFER);
               }
               memmove((void *)&vars[i], vobj_tmp, sizeof(T));
            }
         }

         ClassDef(TCommunicator, 1)
      };

      //______________________________________________________________________________
      template<class Type> void TCommunicator::Send(const Type &var, Int_t dest, Int_t tag) const
      {
         Send(&var, 1, dest, tag);
      }

      //______________________________________________________________________________
      template<class Type> void TCommunicator::Send(const Type *vars, Int_t count, Int_t dest, Int_t tag) const
      {
         if (std::is_class<Type>::value) {
            Char_t *buffer;
            Int_t size;
            Serialize(&buffer, size, vars, count, this, dest, GetRank(), tag);
            MPI_Send(buffer, size, MPI_CHAR, dest, tag, fComm);
            delete buffer;
         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Send((void *)vars, count, GetDataType<Type>(), dest, tag, fComm);
         }
      }


      //______________________________________________________________________________
      template<class Type>  void TCommunicator::Recv(Type &var, Int_t source, Int_t tag) const
      {
         Recv(&var, 1, source, tag);
      }

      //______________________________________________________________________________
      template<class Type>  void TCommunicator::Recv(Type *vars, Int_t count, Int_t source, Int_t tag) const
      {
         if (std::is_class<Type>::value) {
            Int_t size;
            TStatus s;
            Probe(source, tag, s);

            MPI_Get_elements(&s.fStatus, MPI_CHAR, &size);

            Char_t *buffer = new Char_t[size];
            MPI_Recv(buffer, size, MPI_CHAR, source, tag, fComm, MPI_STATUS_IGNORE);
            Unserialize<Type>(buffer, size, vars, count, this, GetRank(), source, tag, 0);

         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            //TODO: added status argument to this method
            MPI_Recv((void *)vars, count, GetDataType<Type>(), source, tag, fComm, MPI_STATUS_IGNORE);
         }
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::ISend(const Type &var, Int_t dest, Int_t tag)
      {
         return ISend(&var, 1, dest, tag);
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::ISend(const Type *vars, Int_t count, Int_t dest, Int_t tag)
      {
         TRequest req;
         if (std::is_class<Type>::value) {
            Char_t *buffer;
            Int_t size;
            Serialize(&buffer, size, vars, count, this, dest, GetRank(), tag);
            MPI_Isend(buffer, size, MPI_CHAR, dest, tag, fComm, &req.fRequest);
            req.fCallback = [buffer]()mutable { //use to clean memory after wait
               if (buffer) delete buffer;
               buffer = NULL;
            };
         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Isend((void *)vars, count, GetDataType<Type>(), dest, tag, fComm, &req.fRequest);
         }
         return req;
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::ISsend(const Type &var, Int_t dest, Int_t tag)
      {
         return ISsend(&var, 1, dest, tag);
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::ISsend(const Type *vars, Int_t count, Int_t dest, Int_t tag)
      {
         TRequest req;
         if (std::is_class<Type>::value) {
            Char_t *buffer;
            Int_t size;
            Serialize(&buffer, size, vars, count, this, dest, GetRank(), tag);
            MPI_Issend(buffer, size, MPI_CHAR, dest, tag, fComm, &req.fRequest);
            req.fCallback = [buffer]()mutable { //use to clean memory after wait
               if (buffer) delete buffer;
               buffer = NULL;
            };
         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Issend((void *)vars, count, GetDataType<Type>(), dest, tag, fComm, &req.fRequest);
         }
         return req;
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::IRsend(const Type &var, Int_t dest, Int_t tag)
      {
         return IRsend(&var, 1, dest, tag);
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::IRsend(const Type *vars, Int_t count, Int_t dest, Int_t tag)
      {
         TRequest req;
         if (std::is_class<Type>::value) {
            Char_t *buffer;
            Int_t size;
            Serialize(&buffer, size, vars, count, this, dest, GetRank(), tag);
            MPI_Irsend(buffer, size, MPI_CHAR, dest, tag, fComm, &req.fRequest);
            req.fCallback = [buffer]()mutable { //use to clean memory after wait
               if (buffer) delete buffer;
               buffer = NULL;
            };
         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Irsend((void *)vars, count, GetDataType<Type>(), dest, tag, fComm, &req.fRequest);
         }
         return req;
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::IRecv(Type &var, Int_t source, Int_t tag) const
      {
         return IRecv(&var, 1, source, tag);
      }
      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::IRecv(Type *vars, Int_t count, Int_t source, Int_t tag) const
      {
         TRequest req;
         if (std::is_class<Type>::value) {
            Int_t size;
            TStatus s;
            while (!IProbe(source, tag, s)) {
               gSystem->Sleep(100);
            }
            MPI_Get_elements(&s.fStatus, MPI_CHAR, &size);

            Char_t *buffer = new Char_t[size];
            MPI_Irecv(buffer, size, MPI_CHAR, source, tag, fComm, &req.fRequest);

            req.fCallback = std::bind(Unserialize<Type>, buffer, size, vars, count, this, GetRank(), source, tag, 0);

         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Irecv((void *)vars, count, GetDataType<Type>(), source, tag, fComm, &req.fRequest);
         }
         return req;
      }

      //______________________________________________________________________________
      template<class Type> void TCommunicator::Bcast(Type &var, Int_t root) const
      {
         Bcast(&var, 1, root);
      }

      //______________________________________________________________________________
      template<class Type> void TCommunicator::Bcast(Type *vars, Int_t count, Int_t root) const
      {
         if (std::is_class<Type>::value) {
            Int_t size;
            Char_t *buffer;

            if (GetRank() == root) Serialize(&buffer, size, vars, count, this, 0, 0, 0, root);

            Bcast(size, root);

            if (GetRank() != root) {
               buffer = new Char_t[size];
            }

            Bcast(buffer, size, root);

            Unserialize(buffer, size, vars, count, this, 0, 0, 0, root);

         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Bcast((void *)vars, count, GetDataType<Type>(), root, fComm);
         }

      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::IBcast(Type *vars, Int_t count, Int_t root) const
      {
         //NOTE: may is good idea to consider to implement tree broadcast algorithm,
         //because I am sending just one integer with the size firts
         TRequest req;
         if (std::is_class<Type>::value) {
            TRequest prereq;
            Int_t size;
            Char_t *buffer;

            if (GetRank() == root) {
               Serialize(&buffer, size, vars, count, this, 0, 0, 0, root);
            }

            prereq = IBcast(size, root);
            prereq.Wait();

            if (GetRank() != root) buffer = new Char_t[size];
            req = IBcast(buffer, size, root);
            req.fCallback = std::bind(Unserialize<Type>, buffer, size, vars, count, this, 0, 0, 0, root);

         } else {
            ROOT_MPI_CHECK_DATATYPE(Type);
            MPI_Ibcast((void *)vars, count, GetDataType<Type>(), root, fComm, &req.fRequest);
         }
         return req;
      }

      //______________________________________________________________________________
      template<class Type> TRequest TCommunicator::IBcast(Type &var, Int_t root) const
      {
         return IBcast(&var, 1, root);
      }

      //______________________________________________________________________________
      template<class Type> void TCommunicator::Scatter(const Type *in_vars, Int_t incount, Type *out_vars, Int_t outcount, Int_t root) const
      {
         if (GetRank() == root) {
            if (incount % (GetSize()*outcount) != 0) {
               Fatal("TCommunicator::Scatter", "Number of elements sent and elements in receive are not divisible. Can't no split to scatter message");
               Abort(ERR_COUNT);
            }
            for (auto i = 0 ; i < GetSize(); i++) {
               if (i == root) continue;
               auto stride = outcount * i;
               Send(&in_vars[stride], outcount, i, MPI_TAG_UB);
            }
            auto stride = outcount * root;
            memcpy((void *)out_vars, (void *)&in_vars[stride], sizeof(Type)*outcount);
         } else {
            Recv(out_vars, outcount, root, MPI_TAG_UB);
         }
      }

      //______________________________________________________________________________
      template<class Type> void TCommunicator::Gather(const Type *in_vars, Int_t incount, Type *out_vars, Int_t outcount, Int_t root) const
      {
         if (GetRank() == root) {
            //TODO: check special cases to improved this error handling
            if ((GetSize()*incount) % outcount   != 0) {
               Fatal("TCommunicator::Gather", "Number of elements sent can't be fitted in gather message");
               Abort(ERR_COUNT);
            }
            for (auto i = 0 ; i < GetSize(); i++) {
               if (i == root) continue;
               auto stride = incount * i;
               Recv(&out_vars[stride], incount, i, MPI_TAG_UB);
            }
            //NOTE: copy memory with memmove because memcpy() with overlapping areas produces undefined behavior
            //In scatter is not same because out_vars have not overlapping, I mean I just need to fill the entire vector not a region
            auto stride = incount * root;
            memmove((void *)&out_vars[stride], (void *)in_vars, sizeof(Type)*incount);
         } else {
            Send(in_vars, incount, root, MPI_TAG_UB);
         }
      }


      template<class Type> void TCommunicator::Reduce(const Type &in_var, Type &out_var, Op<Type> (*opf)(), Int_t root) const
      {
         Reduce(&in_var, &out_var, 1, opf, root);
      }

      template<class Type> void TCommunicator::Reduce(const Type *in_var, Type *out_var, Int_t count, Op<Type> (*opf)(), Int_t root) const
      {
         auto op = opf();

         if (!std::is_class<Type>::value) memmove((void *)out_var, (void *)in_var, sizeof(Type)*count);
         else {
            for (auto i = 0; i < count; i++) {
               TMpiMessage msgi;
               msgi.WriteObject(in_var[i]);
               Char_t *buffer = new Char_t[msgi.BufferSize()]; //this pointer can't be freed, it will be free when the object dies
               memcpy((void *)buffer, (void *)msgi.Buffer(), sizeof(Char_t)*msgi.BufferSize());
               TMpiMessage msgo(buffer, msgi.BufferSize()); //using serialization to copy memory without copy constructor
               auto cl = gROOT->GetClass(typeid(Type));
               auto obj_tmp = msgo.ReadObjectAny(cl);
               memmove((void *)&out_var[i], obj_tmp, sizeof(Type));
            }
         }

         auto size = GetSize();
         auto lastpower = 1 << (Int_t)log2(size);

         for (Int_t i = lastpower; i < size; i++)
            if (GetRank() == i)
               Send(in_var, count, i - lastpower, MPI_TAG_UB);
         for (Int_t i = 0; i < size - lastpower; i++)
            if (GetRank() == i) {
               Type recvbuffer[count];
               Recv(recvbuffer, count, i + lastpower, MPI_TAG_UB);
               for (Int_t j = 0; j < count; j++) out_var[j] = op(in_var[j], recvbuffer[j]);
            }

         for (Int_t d = 0; d < (Int_t)log2(lastpower); d++)
            for (Int_t k = 0; k < lastpower; k += 1 << (d + 1)) {
               auto receiver = k;
               auto sender = k + (1 << d);
               if (GetRank() == receiver) {
                  Type recvbuffer[count];
                  Recv(recvbuffer, count, sender, MPI_TAG_UB);
                  for (Int_t j = 0; j < count; j++) out_var[j]  = op(out_var[j], recvbuffer[j]);
               } else if (GetRank() == sender)
                  Send(out_var, count, receiver, MPI_TAG_UB);
            }
         if (root != 0 && GetRank() == 0) Send(out_var, count, root, MPI_TAG_UB);
         if (root == GetRank() && GetRank() != 0) Recv(out_var, count, 0, MPI_TAG_UB);
      }

      //______________________________________________________________________________
      template<> void TCommunicator::Serialize<TMpiMessage>(Char_t **buffer, Int_t &size, const TMpiMessage *vars, Int_t count, const TCommunicator *comm, Int_t dest, Int_t source, Int_t tag, Int_t root);

      //______________________________________________________________________________
      template<> void TCommunicator::Unserialize<TMpiMessage>(Char_t *ibuffer, Int_t isize, TMpiMessage *vars, Int_t count, const TCommunicator *comm, Int_t dest , Int_t source , Int_t tag , Int_t root);


   }
}
R__EXTERN ROOT::Mpi::TCommunicator *gComm;
R__EXTERN ROOT::Mpi::TCommunicator COMM_WORLD;


#endif
