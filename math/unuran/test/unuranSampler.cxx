// test TUnuran Sampler class using t he DistSampler interface 
#include "gtest/gtest.h"

#include "Math/DistSampler.h"
#include "Math/DistSamplerOptions.h"
#include "Math/Factory.h"
#include "Math/DistFunc.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"

using namespace ROOT::Math; 

std::unique_ptr<TH1D> FillHisto1D( DistSampler & s) {
    // create histogram using automatic binning
    auto h1 = std::make_unique<TH1D>("h1","h1",100, 1, 0);  
    const int nevt = 10000; 
    for (int i = 0; i < nevt; i++) {
        h1->Fill( s.Sample1D());
    }
    return h1;
}

// test using Unuran string API
TEST(OneDim, StringAPI)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    bool ret = sampler->Init("distr = normal(3.,0.75); domain = (0,6) & method = tdr; c = 0");
    EXPECT_EQ(ret, true);
    if (!ret) return;

    sampler->SetMode(0.75);
  
    auto h1 = FillHisto1D(*sampler);
    EXPECT_NEAR(h1->GetMean(), 3., 5 * h1->GetMeanError());
    EXPECT_NEAR(h1->GetRMS(), 0.75, 5*h1->GetRMSError());
}

// test passing a PDF distribution object
TEST(OneDim, DistPdfAPI)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    auto pdf = [](double x){ return ROOT::Math::normal_pdf(x,0.75,3);};
    Functor1D pdfDist(pdf);
    sampler->SetFunction(pdfDist);

    bool ret = sampler->Init("tdr");
    EXPECT_EQ(ret, true);
    if (!ret) return;

    sampler->SetRange(0,6.);
  
    auto h1 = FillHisto1D(*sampler);
    EXPECT_NEAR(h1->GetMean(), 3., 5 * h1->GetMeanError());
    EXPECT_NEAR(h1->GetRMS(), 0.75, 5*h1->GetRMSError());
}

// test passing a CDF distribution object
TEST(OneDim, DistCdfAPI)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));
    auto cdf = [](double x){ return ROOT::Math::normal_cdf(x,0.75,3);};
    Functor1D cdfDist(cdf);
    sampler->SetCdf(cdfDist);

    bool ret = sampler->Init("ninv");
    EXPECT_EQ(ret, true);
    if (!ret) return;

    sampler->SetRange(0,6.);
  
    auto h1 = FillHisto1D(*sampler);
    EXPECT_NEAR(h1->GetMean(), 3., 5 * h1->GetMeanError());
    EXPECT_NEAR(h1->GetRMS(), 0.75, 5*h1->GetRMSError());
}

// test using DistSampler options
TEST(OneDim, DistSamplerOptAPI){

    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    auto pdf = [](double x){ return ROOT::Math::normal_pdf(x,0.75,3);};
    Functor1D pdfDist(pdf);
    sampler->SetFunction(pdfDist);
    DistSamplerOptions opt; 
    opt.SetAlgorithm("tdr");
    // set specific options for the tdr method
    opt.SetAlgoOption("c",0.);
    opt.SetAlgoOption("cpoints",50);
    opt.SetAlgoOption("usedars","true");
    opt.SetAlgoOption("variant_gw","");
    //opt.Print();
    
    bool ret = sampler->Init(opt);
    EXPECT_EQ(ret, true);
    if (!ret) return;

    sampler->SetRange(0,6.);
  
    auto h1 = FillHisto1D(*sampler);
    EXPECT_NEAR(h1->GetMean(), 3., 5 * h1->GetMeanError());
    EXPECT_NEAR(h1->GetRMS(), 0.75, 5*h1->GetRMSError());
}
// test using a discreate PDF
TEST(Discrete, DistAPI)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    auto pdf = [](double x){ return ROOT::Math::poisson_pdf(x,10);};
    Functor1D pdfDist(pdf);
    sampler->SetFunction(pdfDist);
    sampler->SetMode(10);
    sampler->SetArea(1.);

    bool ret = sampler->Init("dari");
    EXPECT_EQ(ret, true);
    if (!ret) return;
  
    auto h1 = FillHisto1D(*sampler);
    EXPECT_NEAR(h1->GetMean(), 10, 5 * h1->GetMeanError());
    EXPECT_NEAR(h1->GetRMS(), sqrt(10), 5*h1->GetRMSError());
}

TEST(Discrete, StringAPI)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    // test using string API from a probability vector
    bool ret = sampler->Init("distr = discr; pv = (0.5,0.2,0.3) & method = dau");
    EXPECT_EQ(ret, true);
    if (!ret) return;
  
    TH1D h1("h1","h1",3,0,3);
    const int nevt = 10000;
    for (int i = 0; i < nevt; i++) { h1.Fill(sampler->Sample1D());}

    double prob[] = {0.5,0.2,0.3};
    for (int ibin = 1; ibin <=3; ibin++) {
        EXPECT_NEAR(h1.GetBinContent(ibin), prob[ibin-1]*nevt, 5*h1.GetBinError(ibin));
    }
}

std::unique_ptr<TH2D> FillHisto2D( DistSampler & s) {
    // create histogram using automatic binning
    auto h2 = std::make_unique<TH2D>("h2","h2",100, 1, 0, 100, 1, 0);  
    const int nevt = 100000;
    for (int i = 0; i < nevt; i++) {
        auto x = s.Sample();
        h2->Fill(x[0],x[1]);
    }
    return h2;
}

// test using a Multidim pdf
TEST(MultiDim, DistAPI)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    auto pdf = [](const double * x){ return ROOT::Math::bigaussian_pdf(x[0],x[1],1.,2.,0.7,0.,5.);};
    Functor pdfDist(pdf,2);
    sampler->SetFunction(pdfDist);

    bool ret = sampler->Init("vnrou");
    EXPECT_EQ(ret, true);
    if (!ret) return;

    auto h1 = FillHisto2D(*sampler); 

    EXPECT_NEAR(h1->GetMean(1), 0, 5 * h1->GetMeanError(1));
    EXPECT_NEAR(h1->GetRMS(1), 1, 5*h1->GetRMSError(1));
    EXPECT_NEAR(h1->GetMean(2), 5, 5 * h1->GetMeanError(2));
    EXPECT_NEAR(h1->GetRMS(2), 2, 5*h1->GetRMSError(2));
    EXPECT_NEAR(h1->GetCorrelationFactor(1,2), 0.7, 0.05);
    
}

// test using a Multidim pdf using also mode
TEST(MultiDim, DistAPIWithOpt)
{
    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    auto pdf = [](const double * x){ return ROOT::Math::bigaussian_pdf(x[0],x[1],1.,2.,0.7,0.,5.);};
    Functor pdfDist(pdf,2);
    sampler->SetFunction(pdfDist);
    sampler->SetMode({0,5});
    sampler->SetRange({-5,-5},{5,15});

    bool ret = sampler->Init("hitro");
    EXPECT_EQ(ret, true);
    if (!ret) return;

    auto h1 = FillHisto2D(*sampler);
    // relax test tolelrance. Method is Markov-Chain method
    EXPECT_NEAR(h1->GetMean(1), 0, 10 * h1->GetMeanError(1));
    EXPECT_NEAR(h1->GetRMS(1), 1, 10*h1->GetRMSError(1));
    EXPECT_NEAR(h1->GetMean(2), 5, 10 * h1->GetMeanError(2));
    EXPECT_NEAR(h1->GetRMS(2), 2, 10*h1->GetRMSError(2));
    EXPECT_NEAR(h1->GetCorrelationFactor(1,2), 0.7, 0.1);
}
// test using DistSampler options
TEST(MultiDim, DistSamplerOptAPI){

    std::unique_ptr<DistSampler> sampler(Factory::CreateDistSampler("Unuran"));

    auto pdf = [](const double * x){ return ROOT::Math::bigaussian_pdf(x[0],x[1],1.,2.,0.7,0.,5.);};
    Functor pdfDist(pdf,2);
    sampler->SetFunction(pdfDist);
    sampler->SetMode({0,5});
    sampler->SetRange({-5,-5},{5,15});

    DistSamplerOptions opt; 
    opt.SetAlgorithm("hitro");
    // set specific options for the tdr method
    opt.SetAlgoOption("burnin",200);
    //opt.SetAlgoOption("startingpoint","(0,5)");
    //opt.Print();

    bool ret = sampler->Init(opt);
    EXPECT_EQ(ret, true);
    if (!ret) return;

    auto h1 = FillHisto2D(*sampler);
    // relax test tolelrance. Method is Markov-Chain method
    EXPECT_NEAR(h1->GetMean(1), 0, 10 * h1->GetMeanError(1));
    EXPECT_NEAR(h1->GetRMS(1), 1, 10*h1->GetRMSError(1));
    EXPECT_NEAR(h1->GetMean(2), 5, 10 * h1->GetMeanError(2));
    EXPECT_NEAR(h1->GetRMS(2), 2, 10*h1->GetRMSError(2));
    EXPECT_NEAR(h1->GetCorrelationFactor(1,2), 0.7, 0.1);
    
}