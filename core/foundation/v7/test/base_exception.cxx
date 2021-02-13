#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <ROOT/RError.hxx>

#include <memory>
#include <stdexcept>

using RException = ROOT::Experimental::RException;
template<typename T>
using RResult = ROOT::Experimental::RResult<T>;

namespace {

/// Used to verify that wrapped return values are not unnecessarily copied
struct ComplexReturnType {
   static int gNCopies;
   ComplexReturnType() { gNCopies++; }
   ComplexReturnType(const ComplexReturnType &) { gNCopies++; }
   ComplexReturnType(ComplexReturnType &&other) = default;
   ComplexReturnType &operator= (const ComplexReturnType &) { return *this; }
   ComplexReturnType &operator= (ComplexReturnType &&other) = default;
};
int ComplexReturnType::gNCopies = 0;

static ROOT::Experimental::RResult<void> TestFailure()
{
   return R__FAIL("test failure");
}

static ROOT::Experimental::RResult<void> TestSuccess()
{
   return ROOT::Experimental::RResult<void>::Success();
}

static ROOT::Experimental::RResult<void> TestSuccessOrFailure(bool succeed)
{
   if (succeed)
      return R__FORWARD_RESULT(TestSuccess());
   return R__FORWARD_RESULT(TestFailure());
}

static ROOT::Experimental::RResult<int> TestSyscall(bool succeed)
{
   if (succeed)
      return 42;
   return R__FAIL("failure");
}

static ROOT::Experimental::RResult<int> TestChain(bool succeed)
{
   auto rv = TestSyscall(succeed);
   return R__FORWARD_RESULT(rv);
}

static ROOT::Experimental::RResult<int> TestChainMultiTypes(bool succeed)
{
   auto rv = TestSuccessOrFailure(succeed);
   if (!rv)
      return R__FORWARD_ERROR(rv);
   return 0;
}

static ROOT::Experimental::RResult<ComplexReturnType> TestComplex()
{
   return ComplexReturnType();
}

class ExceptionX : public std::runtime_error {
public:
   explicit ExceptionX(const std::string &what) : std::runtime_error(what) {}
};

} // anonymous namespace


TEST(Exception, Report)
{
   try {
      TestChain(false);
      EXPECT_TRUE(false) << "Above line should have thrown!";
   } catch (const RException& e) {
      ASSERT_EQ(2U, e.GetError().GetStackTrace().size());
      EXPECT_THAT(e.GetError().GetStackTrace().at(0).fFunction, ::testing::HasSubstr("TestSyscall(bool)"));
      EXPECT_THAT(e.GetError().GetStackTrace().at(1).fFunction, ::testing::HasSubstr("TestChain(bool)"));
   }
}


TEST(Exception, ForwardResult)
{
   auto res = TestChain(true);
   ASSERT_TRUE(static_cast<bool>(res));
   EXPECT_EQ(42, res.Inspect());
}


TEST(Exception, ForwardError)
{
   EXPECT_THROW(TestSuccessOrFailure(false), RException);
   EXPECT_NO_THROW(TestSuccessOrFailure(true));

   auto res = TestChainMultiTypes(true);
   ASSERT_TRUE(static_cast<bool>(res));
   EXPECT_EQ(0, res.Inspect());

   EXPECT_THROW(TestChainMultiTypes(false), RException);
}


TEST(Exception, DiscardReturnValue)
{
   EXPECT_THROW(TestFailure(), RException);
   EXPECT_NO_THROW(TestSuccess());
}


TEST(Exception, CheckReturnValue)
{
   auto rv = TestFailure();
   EXPECT_FALSE(rv);
   // No exception / crash when the scope closes
}


TEST(Exception, DoubleThrow)
{
   try {
      auto rv = TestFailure();
      // Throwing ExceptionX will destruct rv along the way. Since rv carries an error state, it would normally
      // throw an exception itself. In this test, we verify that rv surpresses throwing an exception if another
      // exception is currently active.
      throw ExceptionX("something else went wrong");
   } catch (const ExceptionX&) {
      // This will only catch ExceptionX but not RException. In case rv mistakenly throws an exception,
      // we would notice the test failure by a crash of the unit test.
   }
}

TEST(Exception, VoidThrowOnError)
{
   // no throw on success
   TestSuccess().ThrowOnError();
   // throw on failure
   EXPECT_THROW(TestFailure().ThrowOnError(), RException);
}

TEST(Exception, Syscall)
{
   auto fd = TestSyscall(true);
   if (!fd) {
      // In production code, we would expect error handling code other than throw
      EXPECT_THROW(fd.Throw(), RException);
   }
   EXPECT_EQ(42, fd.Inspect());

   EXPECT_THROW(TestSyscall(false).Inspect(), RException);
}

TEST(Exception, ComplexReturnType)
{
   auto res = TestComplex();
   EXPECT_EQ(1, ComplexReturnType::gNCopies);
}

TEST(Exception, MoveOnlyReturnType)
{
   auto TestMoveOnly = []() -> RResult<std::unique_ptr<int>> {
      return std::make_unique<int>(1);
   };
   auto res = TestMoveOnly();

   // Using Inspect to make a copy won't compile
   // auto copy_inner = res.Inspect();

   // This will compile, but we only have read-only access
   const auto& copy_inner = res.Inspect();
   EXPECT_EQ(1, *copy_inner);

   // Instead, Unwrap is required to get ownership of the move-only type
   auto move_inner = res.Unwrap();
   EXPECT_EQ(1, *move_inner);
   move_inner.reset();
   move_inner = std::make_unique<int>(2);
   EXPECT_EQ(2, *move_inner);
}
