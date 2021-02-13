
#include "ROOT/TIOFeatures.hxx"

#include "TBasket.h"

#include "gtest/gtest.h"

#include <vector>
#include "TBasket.h"

TEST(TIOFeatures, IOBits)
{
   EXPECT_EQ(static_cast<Int_t>(ROOT::EIOFeatures::kSupported) |
                static_cast<Int_t>(ROOT::Experimental::EIOFeatures::kSupported) |
                static_cast<Int_t>(ROOT::Experimental::EIOUnsupportedFeatures::kUnsupported),
             (1 << static_cast<Int_t>(TBasket::kIOBitCount)) - 1);

   EXPECT_EQ(static_cast<Int_t>(ROOT::EIOFeatures::kSupported) &
                static_cast<Int_t>(ROOT::Experimental::EIOUnsupportedFeatures::kUnsupported),
             0);

   EXPECT_EQ(static_cast<Int_t>(ROOT::EIOFeatures::kSupported) &
                static_cast<Int_t>(ROOT::Experimental::EIOFeatures::kSupported),
             0);

   EXPECT_EQ(static_cast<Int_t>(ROOT::Experimental::EIOFeatures::kSupported) &
                static_cast<Int_t>(ROOT::Experimental::EIOUnsupportedFeatures::kUnsupported),
             0);

   // These are currently defined but empty.
   EXPECT_EQ(static_cast<Int_t>(ROOT::Experimental::EIOUnsupportedFeatures::kUnsupported), 0);
   EXPECT_EQ(static_cast<Int_t>(ROOT::EIOFeatures::kSupported), 0);

   // Currently, the experimental features are identical to TBasket::EIOBits
   EXPECT_EQ(static_cast<Int_t>(ROOT::Experimental::EIOFeatures::kSupported),
             static_cast<Int_t>(TBasket::EIOBits::kSupported));
}
