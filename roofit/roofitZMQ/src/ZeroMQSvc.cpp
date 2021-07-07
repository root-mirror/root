#include "RooFit_ZMQ/ZeroMQSvc.h"

#include <functional> // std::ref
#include <ROOT/RMakeUnique.hxx>

ZeroMQSvc &zmqSvc()
{
   static std::unique_ptr<ZeroMQSvc> svc;
   if (!svc) {
      svc = std::make_unique<ZeroMQSvc>();
   }
   return *svc;
}

ZeroMQSvc::Encoding ZeroMQSvc::encoding() const
{
   return m_enc;
}

void ZeroMQSvc::setEncoding(const ZeroMQSvc::Encoding &e)
{
   m_enc = e;
}

zmq::context_t &ZeroMQSvc::context() const
{
   if (!m_context) {
      try {
         m_context = new zmq::context_t;
      } catch (zmq::error_t &e) {
         std::cerr << "ERROR: Creating ZeroMQ context failed. This only happens when PGM initialization failed or when "
                      "a nullptr was returned from zmq_ctx_new because the created context was invalid. Contact ZMQ "
                      "experts when this happens, because it shouldn't.\n";
         throw;
      }
   }
   return *m_context;
}

zmq::socket_t ZeroMQSvc::socket(int type) const
{
   try {
      // the actual work this function should do, all the rest is error handling:
      return zmq::socket_t{context(), type};
   } catch (zmq::error_t &e) {
      // all zmq errors not recoverable from here, only at call site
      std::cerr << "ERROR in ZeroMQSvc::socket: " << e.what() << " (errno: " << e.num() << ")\n";
      throw;
   }
}

zmq::socket_t *ZeroMQSvc::socket_ptr(int type) const
{
   try {
      // the actual work this function should do, all the rest is error handling:
      return new zmq::socket_t(context(), type);
   } catch (zmq::error_t &e) {
      // all zmq errors not recoverable from here, only at call site
      std::cerr << "ERROR in ZeroMQSvc::socket_ptr: " << e.what() << " (errno: " << e.num() << ")\n";
      throw;
   }
}

void ZeroMQSvc::close_context() const
{
   if (m_context) {
      delete m_context;
      m_context = nullptr;
   }
}

zmq::message_t ZeroMQSvc::encode(const char *item) const
{
   std::function<size_t(const char &t)> fun = ZMQ::stringLength;
   return encode(*item, fun);
}

zmq::message_t ZeroMQSvc::encode(const std::string &item) const
{
   return encode(item.c_str());
}

// zmq::message_t ZeroMQSvc::encode(const TObject& item) const {
//  auto deleteBuffer = []( void* data, void* /* hint */ ) -> void {
//    delete [] (char*)data;
//  };
//
//  TBufferFile buffer(TBuffer::kWrite);
//  buffer.WriteObject(&item);
//
//  // Create message and detach buffer
//  // This is the only ZMQ thing that can throw, and only when memory ran out
//  // (errno ENOMEM), and that is something only the caller can fix, so we don't
//  // catch it here:
//  zmq::message_t message(buffer.Buffer(), buffer.Length(), deleteBuffer);
//  buffer.DetachBuffer();
//
//  return message;
//}

bool ZeroMQSvc::send(zmq::socket_t &socket, const char *item, int flags) const
{
   return retry_send(socket, 2, encode(item), flags);
}

bool ZeroMQSvc::send(zmq::socket_t &socket, zmq::message_t &msg, int flags) const
{
   return retry_send(socket, 2, std::ref(msg), flags);
}

bool ZeroMQSvc::send(zmq::socket_t &socket, zmq::message_t &&msg, int flags) const
{
   return retry_send(socket, 2, std::move(msg), flags);
}
