#ifndef ZEROMQPOLLER_H
#define ZEROMQPOLLER_H 1
#include <iostream>
#include <vector>
#include <deque>
#include <exception>
#include <unordered_map>

#include "IZeroMQSvc.h"
#include "functions.h"

class ZeroMQPoller {
public:

   using entry_t = std::tuple<size_t, zmq::PollType, const zmq::socket_t*>;
   // The key is what zmq::socket_t stores inside, and what goes into
   // pollitem_t through zmq::socket_t's conversion to void* operator
   using sockets_t = std::unordered_map<void*, entry_t>;

   using fd_entry_t = std::tuple<size_t, zmq::PollType>;
   using fds_t = std::unordered_map<int, fd_entry_t>;

   using free_t = std::deque<int>;

   ZeroMQPoller() = default;

   std::vector<std::pair<size_t, int>> poll(int timeo = -1)
   {
      std::vector<std::pair<size_t, int>> r;
      if (m_items.empty()) {
         throw std::runtime_error("No sockets registered");
      }
      while (true) {
         try {
            auto n = zmq::poll(&m_items[0], m_items.size(), timeo);
            if (n == 0) return r;
            break;
         } catch (const zmq::error_t& e) {
            if (e.num() != EINTR) {
               std::cerr << e.what() << std::endl;
               throw;
            }
         }
      }
      // TODO: replace this with ranges::v3::zip
      for (size_t i = 0; i < m_items.size(); ++i) {
         void* socket = m_items[i].socket;
         size_t index = 0;
         int flags = 0;
         if (socket == nullptr) {
            // an fd was registered
            std::tie(index, flags) = m_fds[m_items[i].fd];
         } else {
            // a socket was registered
            const zmq::socket_t* s;
            std::tie(index, flags, s) = m_sockets[socket];
         }
         if (m_items[i].revents & short(flags)) {
            r.emplace_back(index, flags);
         }
      }
      return r;
   }

   size_t size() const
   {
      return m_items.size();
   }

   size_t register_socket(zmq::socket_t& socket, zmq::PollType type)
   {
      zmq::socket_t* s = &socket;
      auto it = m_sockets.find(s);
      if (it != m_sockets.end()) {
         return std::get<0>(it->second);
      }
      size_t index = m_free.empty() ? m_items.size() : m_free.front();
      if (!m_free.empty()) m_free.pop_front();
      // NOTE: tis uses the conversion-to-void* operator of
      // zmq::socket_t, which returns the wrapped object
      m_items.push_back({socket, 0, type, 0});

      // We need to lookup by the pointer to the object wrapped by zmq::socket_t
      m_sockets.emplace(m_items.back().socket, std::make_tuple(index, type, s));
      return index;
   }

   size_t register_socket(int fd, zmq::PollType type)
   {
      auto it = m_fds.find(fd);
      if (it != m_fds.end()) {
         return std::get<0>(it->second);
      }
      size_t index = m_free.empty() ? m_items.size() : m_free.front();
      if (!m_free.empty()) m_free.pop_front();
      // NOTE: tis uses the conversion-to-void* operator of
      // zmq::socket_t, which returns the wrapped object
      m_items.push_back({nullptr, fd, type, 0});

      // We need to lookup by the pointer to the object wrapped by zmq::socket_t
      m_fds.emplace(fd, std::make_tuple(index, type));
      return index;
   }


   size_t unregister_socket(zmq::socket_t& socket)
   {
      if (!m_sockets.count(socket.operator void*())) {
         throw std::out_of_range("Socket is not registered");
      }
      // Remove from m_sockets
      // Can't search by the key of m_sockets, as that is the wrapped
      // object, but have to use the pointer to the wrapper
      // (zmq::socket_t)
      auto it = std::find_if(begin(m_sockets), end(m_sockets),
                             [&socket](const decltype(m_sockets)::value_type& entry) {
                                return &socket == std::get<2>(entry.second);
                             });
      auto index = std::get<0>(it->second);
      m_free.push_back(index);
      m_sockets.erase(it);

      // Remove from m_items
      auto iit = std::find_if(begin(m_items), end(m_items), [&it](const zmq::pollitem_t& item) {
            return it->first == item.socket;
         });
      assert(iit != end(m_items));
      m_items.erase(iit);

      return index;
   }

   size_t unregister_socket(int fd)
   {
      if (!m_fds.count(fd)) {
            throw std::out_of_range("fileno is not registered");
      }
      // Remove from m_fds
      auto it = m_fds.find(fd);
      auto index = std::get<0>(it->second);
      m_free.push_back(index);
      m_fds.erase(it);

      // Remove from m_items
      auto iit = std::find_if(begin(m_items), end(m_items), [&it](const zmq::pollitem_t& item) {
            return it->first == item.fd;
         });
      assert(iit != end(m_items));
      m_items.erase(iit);

      return index;
   }


private:

   // Vector of (socket, flags)
   std::vector<zmq::pollitem_t> m_items;
   sockets_t m_sockets;
   fds_t m_fds;

   // free slots in items
   free_t m_free;
};

#endif // ZEROMQPOLLER_H
