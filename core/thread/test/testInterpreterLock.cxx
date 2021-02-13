#include "RConfig.h"
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif
#include "TVirtualRWMutex.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include <iostream>
#include <string>
#include <vector>
#include <thread>

#include "gtest/gtest.h"

static bool gDeclarator =
   (ROOT::EnableThreadSafety(),
    gInterpreter->Declare("int f(int n) { int x = 0; for (auto i = 0u; i < n; ++i) ++x; return x; }"));

constexpr bool gDebugOrder = false;
constexpr unsigned int nThreads = 16;

TEST(InterpreterLock, ConcurrentCalc)
{
   ASSERT_TRUE(gDeclarator && gGlobalMutex);
   bool flag = true;
   auto func = [&](int) {
      gInterpreter->Calc("f(10)");
      R__LOCKGUARD(gROOTMutex);
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : start" << std::endl;

      EXPECT_TRUE(flag) << "Broken flag update on thread " << std::this_thread::get_id();
      flag = false;
      gInterpreter->Calc("f(100000)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : end" << std::endl;
      flag = true;
   };
#ifdef R__USE_IMT
   ROOT::TThreadExecutor pool(32);
   std::vector<int> id;
   for (auto i = 0u; i < nThreads; i++)
      id.push_back(i);
   pool.Foreach(func, id);
#else
   std::vector<std::thread> threads;
   for (unsigned int i = 0; i < nThreads; ++i) {
      threads.emplace_back([&, i] { func(i); });
   }
   std::for_each(threads.begin(), threads.end(), [](std::thread &thr) { thr.join(); });
#endif
}

TEST(InterpreterLock, ReadLocks)
{
   ASSERT_TRUE(gDeclarator && gGlobalMutex);
   auto func = [&](int) {
      R__READ_LOCKGUARD(ROOT::gCoreMutex);
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : start" << std::endl;
      R__READ_LOCKGUARD(ROOT::gCoreMutex);
      gInterpreter->Calc("f(10)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : middle" << std::endl;

      R__LOCKGUARD(gROOTMutex);
      gInterpreter->Calc("f(100000)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : end" << std::endl;
   };
#ifdef R__USE_IMT
   ROOT::TThreadExecutor pool(32);
   std::vector<int> id;
   for (auto i = 0u; i < nThreads; i++)
      id.push_back(i);
   pool.Foreach(func, id);
#else
   std::vector<std::thread> threads;
   for (unsigned int i = 0; i < nThreads; ++i) {
      threads.emplace_back([&, i] { func(i); });
   }
   std::for_each(threads.begin(), threads.end(), [](std::thread &thr) { thr.join(); });
#endif
}

TEST(InterpreterLock, BalancedUserReadLock)
{
   ASSERT_TRUE(gDeclarator && gGlobalMutex);
   auto func = [&](int) {
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : start" << std::endl;
      gInterpreter->Calc("R__READ_LOCKGUARD(ROOT::gCoreMutex); f(10)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : middle" << std::endl;

      R__LOCKGUARD(gROOTMutex);
      gInterpreter->Calc("R__READ_LOCKGUARD(ROOT::gCoreMutex); f(100000)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : end" << std::endl;
   };
#ifdef R__USE_IMT
   ROOT::TThreadExecutor pool(32);
   std::vector<int> id;
   for (auto i = 0u; i < nThreads; i++)
      id.push_back(i);
   pool.Foreach(func, id);
#else
   std::vector<std::thread> threads;
   for (unsigned int i = 0; i < nThreads; ++i) {
      threads.emplace_back([&, i] { func(i); });
   }
   std::for_each(threads.begin(), threads.end(), [](std::thread &thr) { thr.join(); });
#endif
}

TEST(InterpreterLock, BalancedUserWriteLock)
{
   ASSERT_TRUE(gDeclarator && gGlobalMutex);
   auto func = [&](int) {
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : start" << std::endl;
      gInterpreter->Calc("R__WRITE_LOCKGUARD(ROOT::gCoreMutex); f(10)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : middle" << std::endl;

      R__LOCKGUARD(gROOTMutex);
      gInterpreter->Calc("R__WRITE_LOCKGUARD(ROOT::gCoreMutex); f(100000)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : end" << std::endl;
   };
#ifdef R__USE_IMT
   ROOT::TThreadExecutor pool(32);
   std::vector<int> id;
   for (auto i = 0u; i < nThreads; i++)
      id.push_back(i);
   pool.Foreach(func, id);
#else
   std::vector<std::thread> threads;
   for (unsigned int i = 0; i < nThreads; ++i) {
      threads.emplace_back([&, i] { func(i); });
   }
   std::for_each(threads.begin(), threads.end(), [](std::thread &thr) { thr.join(); });
#endif
}

TEST(InterpreterLock, UnBalancedUserReadLock)
{
   ASSERT_TRUE(gDeclarator && gGlobalMutex);
   auto func = [&](int) {
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : start" << std::endl;
      gInterpreter->Calc("ROOT::gCoreMutex->ReadLock(); f(10)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : middle" << std::endl;

      {
         R__LOCKGUARD(gROOTMutex);
         gInterpreter->Calc("R__READ_LOCKGUARD(ROOT::gCoreMutex); f(100000)");
         if (gDebugOrder)
            std::cerr << std::this_thread::get_id() << " : end" << std::endl;
      }

      ROOT::gCoreMutex->ReadUnLock(nullptr);
   };
#ifdef R__USE_IMT
   ROOT::TThreadExecutor pool(32);
   std::vector<int> id;
   for (auto i = 0u; i < nThreads; i++)
      id.push_back(i);
   pool.Foreach(func, id);
#else
   std::vector<std::thread> threads;
   for (unsigned int i = 0; i < nThreads; ++i) {
      threads.emplace_back([&, i] { func(i); });
   }
   std::for_each(threads.begin(), threads.end(), [](std::thread &thr) { thr.join(); });
#endif
}

TEST(InterpreterLock, UnBalancedUserWriteLock)
{
   ASSERT_TRUE(gDeclarator && gGlobalMutex);
   auto func = [&](int) {
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : start" << std::endl;
      gInterpreter->Calc("ROOT::gCoreMutex->ReadLock(); f(10)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : middle" << std::endl;

      gInterpreter->Calc("ROOT::gCoreMutex->WriteLock(); f(100000)");
      if (gDebugOrder)
         std::cerr << std::this_thread::get_id() << " : end" << std::endl;

      ROOT::gCoreMutex->WriteUnLock(nullptr);
      ROOT::gCoreMutex->ReadUnLock(nullptr);
   };
#ifdef R__USE_IMT
   ROOT::TThreadExecutor pool(32);
   std::vector<int> id;
   for (auto i = 0u; i < nThreads; i++)
      id.push_back(i);
   pool.Foreach(func, id);
#else
   std::vector<std::thread> threads;
   for (unsigned int i = 0; i < nThreads; ++i) {
      threads.emplace_back([&, i] { func(i); });
   }
   std::for_each(threads.begin(), threads.end(), [](std::thread &thr) { thr.join(); });
#endif
}
