#include "TROOT.h"
#include <ROOT/RVec.hxx>
#include "TH1D.h"
#include <TLorentzVector.h>
#include <TH2.h>
#include <TFitResult.h>
#include <bits/stdc++.h>
#include <TRandom3.h>
#include <ROOT/RDF/RInterface.hxx>
#include <TF1.h>

using namespace ROOT;

// Pair two vectors according to ids
auto id_pair(RVec<float> v, RVec<int> idx) {
    RVec<float> result;
    int size = Max(idx)+1;
    result.resize(size);

    for (int i=0; i<idx.size(); i++) {
        if (idx[i] >= 0) {
            result[idx[i]] = v[i];
        }
    }

    return result;
}

// Return array that has only the indices included in idx
auto id_filter(RVec<float> v, RVec<int> idx) {
    RVec<float> result;
    int size = Max(idx)+1, idx_size = idx.size();
    result.resize(size);

    for (int i=0; i<idx_size; i++){
        if (idx[i] >= 0) {
            result[idx[i]] = v[idx[i]];
        }
    }

    return result;
}

// Match num elements with denom elements based on idx and return the ratio of the elements
auto vector_division(RVec<float> num, RVec<float> denom) {
    RVec<float> result;
    int size = denom.size();
    result.resize(size);
    
    for (int i = 0; i < size; i++) {
        result[i] = num[i] / denom[i];
    }
    
    return result;
}

// Compare flavours of arr to given type, return arr elements accordingly
auto flavour_match(RVec<float> arr, RVec<int> flavs, int type) {
    RVec<float> result;
    int size = flavs.size();
    for (int i=0; i < size; i++) {
        if (flavs[i] == type) {
            result.push_back(arr[i]);
        }
    }
    return result;
}

// Compare flavours of arr to multiple given types, return arr elements accordingly
auto flavour_match(RVec<float> arr, RVec<int> flavs, RVec<int> type) {
    RVec<float> result;
    int size = flavs.size(), tsize = type.size();
    for (int i=0; i < size; i++) {
        for (int j=0; j < tsize; j++) {
            if (type[j] == flavs[i]) {
                result.push_back(arr[i]);
                break;
            }
        }
    }
    return result;
}

// Return a vector containing element wise deltaR of two four vectors 
auto dR_filt(RVec<float> eta1, RVec<float> phi1, RVec<float> eta2, RVec<float> phi2) {
    RVec<float> result;
    TLorentzVector v1, v2;
    int esize = eta1.size();
    result.resize(esize);

    for (int i=0; i < esize; i++) {
        v1 = TLorentzVector(0.0, eta2[i], phi2[i], 0.0);
        v1.SetPtEtaPhiM(0,eta1[i], phi1[i],0);
        v2 = TLorentzVector(0.0, eta2[i], phi2[i], 0.0);
        v2.SetPtEtaPhiM(0,eta2[i], phi2[i],0);
        result[i] = v1.DeltaR(v2);
    }

    return result;
}

// Return two largest values of a vector
auto leading_two(RVec<float> v) {
    auto indices = Argsort(v);
    return Take(v, {indices.at(v.size() - 1), indices.at(v.size() - 2)});
}

// Return indices of two largest values of a vector RVec<float> v)
auto leading_twoarg(RVec<float> arr) {
    RVec<int> result;
    auto indices = Argsort(arr);
    if (arr.size()==0) {
        result.resize(1);
    } else if (arr.size()==1) {
        result.resize(1);
        result.at(0) = indices.at(arr.size() - 1);
    } else{
        result.resize(2);
        result.at(0) = indices.at(arr.size() - 1);
        result.at(1) = indices.at(arr.size() - 2);
    }
    
    return result;
}

// Median of a histogram
bool interpolateMedian = true;
double Median(const TH1D * h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   if (!interpolateMedian) {
     return TMath::Median(n, &x[0], &y[1]); 
   }
   else {
     double xmid = TMath::Median(n, &x[0], &y[1]);
     int ix = h1->GetXaxis()->FindBin(xmid);
     double nlow = h1->Integral(0,ix-1);
     //double nhigh = h1->Integral(ix+1,h1->GetNbinsX()+1);
     //double nmid = h1->Integral(ix,ix);
     double nhigh = h1->Integral(0,ix);
     //double ntot = nlow+nhigh+nmid; //h1->Integral();
     double ntot = h1->Integral();
     assert(nlow<=0.5*ntot);
     //assert(nhigh<=0.5*ntot);
     assert(nhigh>=0.5*ntot);
     double w = (nhigh>nlow ? (0.5*ntot-nlow) / (nhigh-nlow) : 0);
     return (h1->GetBinLowEdge(ix) + w * h1->GetBinWidth(ix));
   }
}

// Fill a nxbins histogram with 2D histograms y-projections medians
RVec<double> fill_median(TH2 *h, int nxbins) {
    double med;
    RVec<double> result;
    result.resize(nxbins);

    for (int i=1; i<nxbins+1; i++){
        result[i-1] = Median(h->ProjectionY(Form("histo_%d",i), i, i));
    }
    return result;
}

// Fill a nxbins histogram with 2D histograms y-projections Gaussian means
RVec<double> fill_gmean(TH2 *h, int nxbins, double Nsigma) {
    double mean, RMS, sigma, meanGauss, sigma_err;
    TH1D *projection;
    TFitResultPtr gaussFit;
    RVec<double> result;
    result.resize(nxbins);

    for (int i=1; i<nxbins+1; i++){

        // Initial fit based on the mean and RMS of the histogram
        projection = h->ProjectionY(Form("histo_%d",i), i, i);
        RMS = projection->GetRMS();
        mean = projection->GetMean();

        // if structure to avoid empty hists
        if (mean != 0) {
            // Initial fit
            gaussFit = projection -> Fit("gaus","S Q","",mean-RMS, mean+RMS);

            // Final fit parameters
            meanGauss = gaussFit -> Parameter(1);
            sigma = gaussFit -> Parameter(2);

            // result[i-1] = meanGauss;
            // Final fit and collect result
            gaussFit = projection -> Fit("gaus","S Q","", meanGauss-Nsigma*sigma, meanGauss+Nsigma*sigma);
            result[i-1] = gaussFit -> Parameter(1);
        } else {
            result[i-1] = 0;
        }
    }
    return result;
}

// Fill a nxbins histogram with 2D histograms y-projections Gaussian errors
RVec<double> fill_gmean_err(TH2 *h, int nxbins, double Nsigma) {
    double mean, RMS, sigma, meanGauss, sigma_err;
    TH1D *projection;
    TFitResultPtr gaussFit;
    RVec<double> result;
    result.resize(nxbins);

    for (int i=1; i<nxbins+1; i++){
        projection = h->ProjectionY(Form("histo_%d",i), i, i);
        RMS = projection->GetRMS();
        mean = projection->GetMean();

        if (mean != 0) {
            gaussFit = projection -> Fit("gaus","S Q","",mean-RMS, mean+RMS);

            meanGauss = gaussFit -> Parameter(1);
            sigma = gaussFit -> Parameter(2);

            gaussFit = projection -> Fit("gaus","S Q","", meanGauss-Nsigma*sigma, meanGauss+Nsigma*sigma);

            result[i-1] = gaussFit -> ParError(1);// Parameter(1);
        } else {
            result[i-1] = 0;
        }
    }
    return result;
}

// Create a 2D vector that includes the borders of the veto area given in TH2 with format v(x_lower, x_upper, y_lower, y_upper)
// TODO: optimize
auto veto_map(TH2 *h, int pass) {
    RVec<RVec<float>> result(4);
    bool in = false;
    int s = 0;

    TAxis *xAxis = h->GetXaxis(); TAxis *yAxis = h->GetYaxis();

    // Find areas to veto
    // loop goes through columns, if a veto area begins or ends the borders are saved
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        for (int j = 1; j <= h->GetNbinsY(); j++) {
            if (h->GetBinContent(i, j) > pass) {
                if (!in) {
                    result[2].push_back(yAxis->GetBinLowEdge(j));
                    result[0].push_back(xAxis->GetBinLowEdge(i));
                    result[1].push_back(xAxis->GetBinUpEdge(i));
                    in = true;           
                } 
                
            } else if (h->GetBinContent(i, j) < pass) {
                if (in) {
                    result[3].push_back(yAxis->GetBinLowEdge(j));
                    in = false;                       
                }
            }
        }
        in = false;
    }

    // Connect the borders of areas, this could probably be done inside the upper loop
    while(true) {
        if ((result[0][s] == result[1][s-1]) && (result[2][s] == result[2][s-1]) && (result[2][s] == result[2][s-1])) {
            result[1][s-1] = result[1][s];
            result[0].erase(result[0].begin() + s);
            result[1].erase(result[1].begin() + s);
            result[2].erase(result[2].begin() + s);
            result[3].erase(result[3].begin() + s);
        }
        else {
            s = s+1;
        }
        if (s==result[0].size()-1) {
            break;
        }
    }

    return result;
}

// Filter RVec based on veto map limits in the form limits(eta_lower, eta_upper, phi_lower, phi_upper)
auto veto_map_filt(RVec<RVec<float>> limits, RVec<float> eta, RVec<float> phi, RVec<float> input) {
    RVec<float> result;
    bool in = false;
    int esize = eta.size(), lsize = limits[0].size();

    // Go through elements of eta, if the element isn't inside a veto area return it
    for (int k = 0; k < esize; k++) {
        for (int i = 0; i < lsize; i++) {
            if (((eta[k] >= limits[0][i]) && (eta[k] < limits[1][i]) && (phi[k] >= limits[2][i]) && (phi[k] < limits[3][i]))) {
                in = true;
            }
        }
        if (!in) {
            result.push_back(input[k]);
        }
        in = false;
    }

    return result;
}

// Return 1 or 0 randomly
int rand_oneZero() {
    TRandom3 gen;
    gen.SetSeed(1556);
    double r = gen.Rndm();
    if (r > 0.5) {
        return 1;
    } else {
        return 0;
    }
}


float get_weight(float pT, UChar_t nConstituents, RVec<float> bins, RVec<TF1> fits) {
    int idx, bin_size = bins.size();
    float result;
    
    for (int i=1; i < bin_size; i++) {
        if (pT > bins[i-1] && pT < bins[i]) {
            idx = i;
            break;
        }
    }
    
    result = fits[idx].Eval(nConstituents);
    
    return result;
}

ROOT::RDataFrame create_weights(ROOT::RDataFrame df_DT, ROOT::RDataFrame df_MC, RVec<float> bins, int cnt) {
    auto MC_hist = df_MC.Histo2D({'mcHist', 'mcHistTitle', cnt, bins, 100u, 0., 100.}, 'pTtag', 'Jet_nConstituents');
    auto DT_hist = df_DT.Histo2D({'dtHist', 'dtHistTitle', cnt, bins, 100u, 0., 100.}, 'pTtag', 'Jet_nConstituents');
    
    RVec<TF1> fits;
    
    for (int i=1; i <= cnt; i++) {
        auto h1 = MC_hist -> ProjectionY('mcProj', i, i);
        auto h2 = DT_hist -> ProjectionY('dtProj', i, i);
        
        if ((h1->Integral() != 0.0) && (h2->Integral() != 0.0)) {
            h1->Scale(1/(h1 -> Integral()));
            h2->Scale(1/(h2 -> Integral()));
            
            h1 -> Divide(h2);
            h1 -> Fit('chebyshev4', 'S');
            
            fits.push_back(h1->GetFunction("chebyshev4"));
        }
        else {
            auto func = TF1(); 
            fits.push_back(func);
        }
    }
    
    ROOT::RDataFrame df = df_MC.Define('dtWeights', get_weight, '{pTtag, Jet_nConstituents, bins, fits}');
    
    return df;
    
}
