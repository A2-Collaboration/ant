#include "Clustering.h"
#include "detail/Clustering_NextGen.h"

#include "base/Detector_t.h"

#include "tree/TCluster.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

bool check_TClusterHit(const TClusterHit& hit, const ClusterDetector_t& clusterdetector) {
    if(hit.IsSane())
        return true;
    if(clusterdetector.HasElementFlags(hit.Channel, Detector_t::ElementFlag_t::BadTDC)
       && isfinite(hit.Energy))
        return true;
    return false;
}

void Clustering_NextGen::Build(const ClusterDetector_t& clusterdetector,
                               const TClusterHitList& clusterhits,
                               TClusterList& clusters) const
{
    // clustering detector, so we need additional information
    // to build the crystals_t
    list<clustering::crystal_t> crystals;
    for(const TClusterHit& hit : clusterhits) {
        // try to include as many hits as possible
        if(!check_TClusterHit(hit, clusterdetector)) {
            continue;
        }
        crystals.emplace_back(
                    hit.Energy,
                    clusterdetector.GetClusterElement(hit.Channel),
                    addressof(hit)
                    );
    }

    // do the clustering (calls detail/Clustering_NextGen.h code)
    vector< clustering::cluster_t > crystal_clusters;
    clustering::do_clustering(crystals, crystal_clusters);

    // now calculate some cluster properties,
    // and create TCluster out of it

    for(const clustering::cluster_t& cluster : crystal_clusters) {
        const double cluster_energy = clustering::calc_total_energy(cluster);

        clusters.emplace_back(
                    vec3(0,0,0),
                    cluster_energy,
                    std_ext::NaN,
                    clusterdetector.Type,
                    0
                    );
        auto& the_cluster = clusters.back();

        auto& clusterhits = the_cluster.Hits;
        clusterhits.reserve(cluster.size());

        double weightedSum = 0;
        double cluster_maxenergy   = 0;
        bool crystalTouchesHole = false;
        for(const clustering::crystal_t& crystal : cluster) {

            clusterhits.emplace_back(*crystal.Hit);

            double wgtE = clustering::calc_energy_weight(crystal.Energy, cluster_energy);
            the_cluster.Position += crystal.Element->Position * wgtE;
            weightedSum += wgtE;

            crystalTouchesHole |= crystal.Element->TouchesHole;

            // search for crystal with maximum energy
            // which is defined as the central element
            if(crystal.Energy >= cluster_maxenergy) {
                cluster_maxenergy = crystal.Energy;

                the_cluster.SetFlag(TCluster::Flags_t::TouchesHoleCentral, crystal.Element->TouchesHole);
                the_cluster.Time = crystal.Hit->Time;
                the_cluster.CentralElement = crystal.Element->Channel;
                // search for short energy
                for(const TClusterHit::Datum& datum : crystal.Hit->Data) {
                    if(datum.Type == Channel_t::Type_t::IntegralShort) {
                        the_cluster.ShortEnergy = datum.Value.Calibrated;
                        break;
                    }
                }
            }
        }
        the_cluster.Position *= 1.0/weightedSum;

        if(cluster.Split)
            the_cluster.SetFlag(TCluster::Flags_t::Split);

        if(crystalTouchesHole)
            the_cluster.SetFlag(TCluster::Flags_t::TouchesHoleCrystal);
    }
}

// Clustering_Sergey as PIMPL, a lot of copy-pasted code


#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/TAPS.h"

#include "TVector3.h"
#include "TMath.h"

struct Clustering_Sergey::Impl {
    enum {
        ENullHit = 0xFFFFFFFF,       // undefined hit index (end of hit buffer)
        EBufferEnd = 0xFFFFFFFF,     // end of file marker
    };
    enum { EFalse, ETrue };         // Logic...should use kTRUE, kFALSE

    struct TA2ClusterDetector;

    struct HitCluster_t {
        Double_t fEnergy =0;        // Total energy deposited in cluster
        Double_t fTheta =0;         // Cluster's theta
        Double_t fPhi =0;           // Cluster's phi
        UInt_t fNhits =0;          // # of hits in cluster
        UInt_t fIndex =0;           // index of central element
        Double_t fSqrtEtot =0;      // Sum of sqrt( energy[i] )
        Double_t fTime =0;          // Time of max-energy element of cluster
        Double_t fSqrtEtUp =0;
        Double_t fSqrtEtDn =0;
        TVector3 *fMeanPosUp =0;    // energy-weighted mean pos. of clust.
        TVector3 *fMeanPosDn =0;    // energy-weighted mean pos. of clust.
        TVector3 *fMeanPosition =0; // energy-weighted mean pos. of clust.
        UInt_t *fHits =0;          // indices of hit elements
        UInt_t fNNeighbour =0;     // # neighbour elements in array
        Double_t fCentralFrac =0;   // Fractional energy in central crystal
        UInt_t *fNeighbour =0;      // indices of neighbouring elements
        Double_t fRadius =0;       // effective radius of cluster
        Int_t fMaxHits =0;         // size of hits array
        UInt_t fNNearNeighbour =0; // # nearest neighbours

        Double_t GetEnergy() { return fEnergy; }
        Double_t GetTheta() { return fTheta; }
        Double_t GetPhi() { return fPhi; }
        UInt_t GetNhits() { return fNhits; }
        Double_t GetTime() { return fTime; }
        UInt_t GetIndex() { return fIndex; }
        TVector3 *GetMeanPosition() { return fMeanPosition; }
        UInt_t *GetHits() { return fHits; }

        HitCluster_t(const Char_t *line, UInt_t index, Int_t sizefactor = 1);
        Double_t ClusterRadius(TA2ClusterDetector *cldet);
        void ClusterDetermine(TA2ClusterDetector *cldet);
        Bool_t ClusterDetermine2(TA2ClusterDetector *cldet);
        void Cleanup();
    };


    struct TA2ClusterDetector {

        UInt_t *fTempHits =0;          // Element-Hit store
        Int_t* fHits =0;                          // array elements which fired in event
        UInt_t fNhits =0;                         // No. detector hits in event
        UInt_t fNCluster =0;           // # of clusters
        UInt_t fNelement =0;                      // # elements in array
        UInt_t *fTryHits =0;           // Element-Hit store
        UInt_t fMaxCluster =0;         // Max # of clusters
        Double_t* fEnergy =0;                     // stored hit energies if any
        UInt_t *fTempHits2 =0;         // Element-Hit store
        HitCluster_t **fCluster =0;    // Clusters of hits
        Double_t *fClEnergyOR =0;      // OR of cluster energies
        Double_t *fTheta =0;           // theta of cluster hit
        Double_t *fPhi =0;             // phi of cluster hit
        UInt_t *fClustHit =0;          // Cluster indices
        UInt_t *fNClustHitOR =0;       // OR of #hits in individuyal clusters
        Double_t *fClRadius =0;        // cluster radius
        Double_t* fTime =0;                       // stored hit times
        TVector3** fPosition =0;                  // stored hit positions
        Int_t fClustSizeFactor;     // enlarge factor, hit cluster buffers

        Double_t* GetEnergy(){ return fEnergy; }          // ptr to energy array
        Double_t* GetTime(){ return fTime; }              // ptr to time array
        UInt_t *GetTempHits2() { return fTempHits2; }
        TVector3** GetPosition(){ return fPosition; }     // ptr position array
        UInt_t GetNhits(){ return fNhits; }              // No. hits in event
        UInt_t GetNelement(){ return fNelement; }        // max detector elements

        TA2ClusterDetector(const ClusterDetector_t& ant_det);
        void DecodeCluster();
        void Cleanup();
    };

    // internal dispatch to Acqu classes
    // by emulating the two detectors internally
    TA2ClusterDetector cb;
    TA2ClusterDetector taps;

    Impl(const ClusterDetector_t& ant_cb, const ClusterDetector_t& ant_taps) :
        cb(ant_cb),
        taps(ant_taps)
    {

    }

    // non-const as TA2Detector is not const-correct at all
    void Build(const ClusterDetector_t& clusterdetector,
               const TClusterHitList& clusterhits,
               TClusterList& clusters)
    {
        auto& det = clusterdetector.Type == Detector_t::Type_t::TAPS ? taps : cb;

        // fill the hits, similar to TA2Detector::DecodeBasic()
        det.fNhits = 0;
        for(auto& clusterhit : clusterhits) {
            // try to include as many hits as possible
            if(!check_TClusterHit(clusterhit, clusterdetector)) {
                continue;
            }
            const auto ch = clusterhit.Channel;
            det.fHits[det.fNhits] = ch;
            det.fTime[ch] = clusterhit.Time;
            det.fEnergy[ch] = clusterhit.Energy;
            det.fNhits++;
        }

        // decode the clusters
        det.DecodeCluster();

        // "readout" the clusters
        for(auto cl=0u;cl<det.fNCluster;cl++) {
            auto acqu_hitcluster = det.fCluster[det.fClustHit[cl]];
            clusters.emplace_back(
                        *acqu_hitcluster->GetMeanPosition(),
                        acqu_hitcluster->GetEnergy(),
                        acqu_hitcluster->GetTime(), // timing
                        clusterdetector.Type,
                        acqu_hitcluster->GetIndex() // central element
                        );
            // copy the hits
            auto& the_cluster = clusters.back();
            auto nHits = acqu_hitcluster->GetNhits();
            auto Hits = acqu_hitcluster->GetHits();
            for(auto i=0u;i<nHits;i++) {
                auto ch = Hits[i];
                // search the initial hits,
                // that might be slow...
                auto it_hit = find_if(clusterhits.begin(), clusterhits.end(),
                                      [ch] (const TClusterHit& hit) {
                    return ch == hit.Channel;
                });
                if(it_hit != clusterhits.end())
                    the_cluster.Hits.emplace_back(*it_hit);
                else
                    throw runtime_error("Did not find TClusterHit for ch="+to_string(ch));
            }
        }

        // clean the detector
        det.Cleanup();
    }
};

// implement stuff

Clustering_Sergey::Impl::TA2ClusterDetector::TA2ClusterDetector(const ClusterDetector_t& ant_det)
{
    fNelement = ant_det.GetNChannels();
    fMaxCluster = 30;
    fClustHit = new UInt_t[fMaxCluster];
    fTheta = new Double_t[fNelement];
    fPhi = new Double_t[fNelement];
    fNClustHitOR = new UInt_t[fNelement];
    fClEnergyOR = new Double_t[fNelement];
    fClRadius = new Double_t[fNelement];
    fEnergy = new Double_t[fNelement];
    fTime = new Double_t[fNelement];
    fHits = new Int_t[fNelement];
    fTempHits = new UInt_t[fNelement];
    fTryHits = new UInt_t[fNelement];
    fTempHits2 = new UInt_t[fNelement];
    fCluster = new HitCluster_t *[fNelement];
    fClustSizeFactor = 1;
    fPosition = new TVector3*[fNelement];

    for(auto i=0u;i<fNelement;i++) {
        // fake AcquRoot Next-Neighbour: config line
        // first number of neighbours (including itself, thus +1 below)
        auto clusterelem = ant_det.GetClusterElement(i);
        stringstream configline;
        configline << (clusterelem->Neighbours.size()+1) << ' ' << i;
        for(auto n : clusterelem->Neighbours)
            configline << ' ' << n;
        fCluster[i] = new HitCluster_t(configline.str().c_str(),i,fClustSizeFactor);

        // copy position of each element
        fPosition[i] = new TVector3(clusterelem->Position);
    }

    // clean the detector
    Cleanup();
}

void Clustering_Sergey::Impl::TA2ClusterDetector::Cleanup() {
    UInt_t i;
    for (i = 0; i < fNelement; i++) {
        fEnergy[i] = (Double_t)ENullHit;
        fTime[i] = (Double_t)ENullHit;
    }
    //  would also clean up any cluster info here
    for (UInt_t i = 0; i < fNCluster; i++)
      fCluster[fClustHit[i]]->Cleanup();
}

void Clustering_Sergey::Impl::TA2ClusterDetector::DecodeCluster() {

    // Determine clusters of hits
    // Search around peak energies absorbed in individual crystals
    // Make copy of hits array as the copy will be altered

    memcpy(fTempHits, fHits, sizeof(UInt_t) * fNhits); // temp copy

    const Double_t fEthcrs_CB = 4., fEthcls_CB = 12.;
    const Double_t fEthcrs_TAPS = 4., fEthcls_TAPS = 12.;
    const Double_t fEthcls2_CB = 70.; // version 2
    const Double_t fEthcls2_TAPS = 50.;
    Double_t opangl;

    const Double_t opangl1 = 24., opangl2 = 8.; // ver 5
    const Double_t difmax = 24.;
    Double_t fEthcrs, fEthclsi, fEthcls, fEthcls2, thet, phi, Ecl, Ecli, oang;
    //-
    Double_t maxenergy;
    UInt_t i, j, k, kmax=0, jmax=0, m, ntaken, nc;
    UInt_t nomit = 0;
    TVector3 vcl, vcln, vdif;

    fNCluster = 0;
    if (fNhits > 250 || fNhits < 1)
        goto OUT;

    if (fNelement == 720) {
        fEthcrs = fEthcrs_CB;
        fEthcls = fEthcls_CB;
        fEthcls2 = fEthcls2_CB;
    } else {
        fEthcrs = fEthcrs_TAPS;
        fEthcls = fEthcls_TAPS;
        fEthcls2 = fEthcls2_TAPS;
    }
    // Find hit with maximum energy

    for (m = 0; m < fNhits; m++)
        fTryHits[m] = fTempHits[m];

    for (i = 0; i < fMaxCluster;) {
        maxenergy = 0;
        for (j = 0; j < fNhits; j++) {
            if ((k = fTryHits[j]) == ENullHit)
                continue;
            if ((k = fTempHits[j]) == ENullHit)
                continue;
            if (maxenergy < fEnergy[k] && fEnergy[k] > fEthcrs) {
                maxenergy = fEnergy[k];
                kmax = k;
                jmax = j;
            }
        }
        if (maxenergy == 0)
            break; // no more isolated hits
        if (kmax < fNelement) {
            for (m = 0; m < fNhits; m++)
                fTempHits2[m] = fTempHits[m];
            fCluster[kmax]->ClusterDetermine(this); // determine the cluster
            Ecl = fCluster[kmax]->GetEnergy();
            if (Ecl >= fEthcls) {
                // printf("Ecl* %lf\n",Ecl);
                if (Ecl < fEthcls2 && fNCluster > 0) {
                    for (j = 0; j < fNCluster; j++) {
                        Ecli = fClEnergyOR[j];
                        fEthclsi = 26. + 24. * Ecli / 1000.;
                        thet = fTheta[j] * TMath::DegToRad();
                        phi = fPhi[j] * TMath::DegToRad();
                        vcl.SetMagThetaPhi(146., thet, phi);
                        thet = fCluster[kmax]->GetTheta() * TMath::DegToRad();
                        phi = fCluster[kmax]->GetPhi() * TMath::DegToRad();
                        vcln.SetMagThetaPhi(146., thet, phi);
                        if (fNelement == 720) {
                            oang = vcl.Angle(vcln) * TMath::RadToDeg();
                            opangl = opangl1 + opangl2 * Ecli / 1000.; // ver 2
                            if (oang < opangl && Ecl < fEthclsi) {
                                nomit++;
                                fTryHits[jmax] = ENullHit;
                                goto NEXTCR;
                            }
                        } else {
                            vdif = vcl - vcln;
                            if (vdif.Mag() < difmax && Ecl < fEthclsi) {
                                nomit++;
                                fTryHits[jmax] = ENullHit;
                                goto NEXTCR;
                            }
                        }
                    }
                }
                for (m = 0; m < fNhits; m++)
                    fTempHits[m] = fTempHits2[m];
                fClustHit[i] = kmax;
                fTheta[i] = fCluster[kmax]->GetTheta();
                fPhi[i] = fCluster[kmax]->GetPhi();
                fNClustHitOR[i] = fCluster[kmax]->GetNhits();
                fClEnergyOR[i] = Ecl;
                fClRadius[i] = fCluster[kmax]->ClusterRadius(this);
                i++;
                fNCluster = i;
            } else
                fTryHits[jmax] = ENullHit;
        }
        // If you reach here then there is an error in the decode
        // possible bad detector ID
        else
            fTryHits[jmax] = ENullHit;
NEXTCR:
        continue;
    }
OUT:
    fClustHit[fNCluster] = EBufferEnd;
    fTheta[fNCluster] = EBufferEnd;
    fPhi[fNCluster] = EBufferEnd;
    fNClustHitOR[fNCluster] = EBufferEnd;
    fClEnergyOR[fNCluster] = EBufferEnd;
    fClRadius[fNCluster] = EBufferEnd;

    if (fNCluster == 0)
        return;

    ntaken = 0;
    for (m = 0; m < fNhits; m++) {
        fTempHits2[m] = fTempHits[m];
        if (fTempHits2[m] == ENullHit)
            ntaken++;
    }
    if (ntaken == fNhits)
        return;

    for (j = 0; j < fNCluster; j++) {
        kmax = fClustHit[j];
        if (fCluster[kmax]->ClusterDetermine2(this)) { // the wider cluster
            fTheta[j] = fCluster[kmax]->GetTheta();
            fPhi[j] = fCluster[kmax]->GetPhi();
            fNClustHitOR[j] = fCluster[kmax]->GetNhits();
            fClEnergyOR[j] = fCluster[kmax]->GetEnergy();
            fClRadius[j] = fCluster[kmax]->ClusterRadius(this);
        }
    }

    if (nomit == 0)
        return;
    ntaken = 0;
    for (m = 0; m < fNhits; m++) {
        if (fTempHits2[m] == ENullHit)
            ntaken++;
        fTryHits[m] = fTempHits[m] = fTempHits2[m];
    }
    if (ntaken == fNhits)
        return;

    nc = fNCluster;
    for (i = nc; i < fMaxCluster;) {
        maxenergy = 0;
        for (j = 0; j < fNhits; j++) {
            if ((k = fTryHits[j]) == ENullHit)
                continue;
            if ((k = fTempHits[j]) == ENullHit)
                continue;
            if (maxenergy < fEnergy[k] && fEnergy[k] > fEthcrs) {
                maxenergy = fEnergy[k];
                kmax = k;
                jmax = j;
            }
        }
        if (maxenergy == 0)
            break; // no more isolated hits
        if (kmax < fNelement) {
            for (m = 0; m < fNhits; m++)
                fTempHits2[m] = fTempHits[m];
            fCluster[kmax]->ClusterDetermine(this); // determine the cluster
            Ecl = fCluster[kmax]->GetEnergy();
            if (Ecl >= fEthcls) {
                for (m = 0; m < fNhits; m++)
                    fTempHits[m] = fTempHits2[m];
                fClustHit[i] = kmax;
                fTheta[i] = fCluster[kmax]->GetTheta();
                fPhi[i] = fCluster[kmax]->GetPhi();
                fNClustHitOR[i] = fCluster[kmax]->GetNhits();
                fClEnergyOR[i] = Ecl;
                fClRadius[i] = fCluster[kmax]->ClusterRadius(this);
                i++;
                fNCluster = i;
            } else
                fTryHits[jmax] = ENullHit;
        }
        // If you reach here then there is an error in the decode
        // possible bad detector ID
        else
            fTryHits[jmax] = ENullHit;
    }
    fClustHit[fNCluster] = EBufferEnd;
    fTheta[fNCluster] = EBufferEnd;
    fPhi[fNCluster] = EBufferEnd;
    fNClustHitOR[fNCluster] = EBufferEnd;
    fClEnergyOR[fNCluster] = EBufferEnd;
    fClRadius[fNCluster] = EBufferEnd;

    if (nc == fNCluster)
        return;
    ntaken = 0;
    for (m = 0; m < fNhits; m++)
        if (fTempHits2[m] == ENullHit)
            ntaken++;
    if (ntaken == fNhits)
        return;

    for (j = nc; j < fNCluster; j++) {
        kmax = fClustHit[j];
        if (fCluster[kmax]->ClusterDetermine2(this)) { // the wider cluster
            fTheta[j] = fCluster[kmax]->GetTheta();
            fPhi[j] = fCluster[kmax]->GetPhi();
            fNClustHitOR[j] = fCluster[kmax]->GetNhits();
            fClEnergyOR[j] = fCluster[kmax]->GetEnergy();
            fClRadius[j] = fCluster[kmax]->ClusterRadius(this);
        }
    }
}

Clustering_Sergey::Impl::HitCluster_t::HitCluster_t(const Char_t* line, UInt_t index, Int_t sizefactor) {
    // store input parameters
    // # inner nearest neighbours (outer calculated from total read)
    // coordinates of center of front face of central element
    // List of nearest neighbours inner & outer
    UInt_t hit[64], nw;
    fIndex = index;
    UInt_t n = sscanf(
        line, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%"
              "d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
        &fNNearNeighbour, hit, hit + 1, hit + 2, hit + 3, hit + 4, hit + 5,
        hit + 6, hit + 7, hit + 8, hit + 9, hit + 10, hit + 11, hit + 12,
        hit + 13, hit + 14, hit + 15, hit + 16, hit + 17, hit + 18, hit + 19,
        hit + 20, hit + 21, hit + 22, hit + 23, hit + 24, hit + 25, hit + 26,
        hit + 27, hit + 28, hit + 29, hit + 30, hit + 31, hit + 32, hit + 33,
        hit + 34, hit + 35, hit + 36, hit + 37, hit + 38, hit + 39, hit + 40,
        hit + 41, hit + 42, hit + 43, hit + 44, hit + 45, hit + 46, hit + 47,
        hit + 48, hit + 49, hit + 50, hit + 51, hit + 52, hit + 53, hit + 54,
        hit + 55, hit + 56, hit + 57, hit + 58, hit + 59, hit + 60, hit + 61,
        hit + 62, hit + 63);

    // Consistency check...1st hit must be the index
    if ((n < (fNNearNeighbour + 1)) || (index != *hit)) {
      printf(" Error in nearest neighbour input at line:\n %s\n", line);
      return;
    }
    n -= 2; // # neighbours around central element
    fNNeighbour = n;
    fNeighbour = new UInt_t[n];
    if (n > 7)
      nw = 22 + 230;
    else
      nw = 19 + 230;
    fHits = new UInt_t[nw];
    fMaxHits = n * sizefactor;
    // fHits = new UInt_t[ fMaxHits ];
    fHits[0] = ENullHit;
    fNhits = 0;
    fEnergy = (Double_t)ENullHit;
    fSqrtEtot = (Double_t)ENullHit;
    fSqrtEtUp = 0.;
    fSqrtEtDn = 0.;
    for (UInt_t i = 0; i < n; i++)
      fNeighbour[i] = hit[i + 1];
    fMeanPosition = new TVector3(0.0, 0.0, 0.0);
    fMeanPosUp = new TVector3(0.0, 0.0, 0.0);
    fMeanPosDn = new TVector3(0.0, 0.0, 0.0);
    fTheta = (Double_t)ENullHit;
    fPhi = (Double_t)ENullHit;
}

void Clustering_Sergey::Impl::HitCluster_t::Cleanup() {
    // End-of-Event cleanup
    fNhits = 0;
    *fHits = ENullHit;
    fEnergy = (Double_t)ENullHit;
    fSqrtEtot = (Double_t)ENullHit;
    fSqrtEtUp = 0.;
    fSqrtEtDn = 0.;
    fRadius = (Double_t)ENullHit;
}

void Clustering_Sergey::Impl::HitCluster_t::ClusterDetermine(TA2ClusterDetector *cldet) {
    // Determine the boundary of the cluster the local total energy
    // and the sqrt(energy)-weighted centre-of-gravity vector

    const Double_t Peng = 2. / 3.;
    UInt_t i, k, m, icl;
    Double_t energyi, wtime;
    Double_t *energy = cldet->GetEnergy();
    Double_t *time = cldet->GetTime();
    UInt_t *hits = cldet->GetTempHits2();
    TVector3 **pos = cldet->GetPosition();
    UInt_t nhits = cldet->GetNhits();
    UInt_t nelem = cldet->GetNelement();
    TVector3 vcr, vcl, vdif;

    fEnergy = energy[fIndex]; // energy in "central" element
    if (fEnergy > 2000.)
        printf("fEnergy = %lf, %d, %d\n", fEnergy, fIndex, nhits);
    Double_t sqrtE = pow(fEnergy, Peng);
    fSqrtEtot = sqrtE;
    if (nelem == 720)
        wtime = energy[fIndex];
    else
        wtime = sqrtE;
    fTime = time[fIndex] * wtime; // time in central element
    for (m = 0; m < nhits; m++)
        if (fIndex == hits[m])
            hits[m] = ENullHit;
    vcr = *(pos[fIndex]);
    if (nelem == 720) {
        if (vcr.Y() > 0.) {
            fSqrtEtUp = sqrtE;
            *fMeanPosUp = vcr * sqrtE;
        } else {
            fSqrtEtDn = sqrtE;
            *fMeanPosDn = vcr * sqrtE;
        }
    }
    *fMeanPosition = vcr * sqrtE; // position = "centre"
    fHits[0] = fIndex;
    k = 1;

    // Accumulate weighted mean position
    UInt_t nneib = fNNeighbour;
    for (m = 0; m < nneib; m++) {
        icl = fNeighbour[m];

        for (i = 0; i < nhits; i++) {
            if (icl == hits[i]) {
                hits[i] = ENullHit;
                fHits[k] = icl; // add to cluster hits collection
                energyi = energy[icl];
                sqrtE = pow(energyi, Peng);
                if (energyi > 2000.)
                    printf("energy[j] = %lf, %d, %d\n", energy[icl], icl, nhits);
                fEnergy += energyi;
                if (nelem == 720)
                    wtime = energyi;
                else
                    wtime = sqrtE;
                fTime += time[icl] * wtime;
                fSqrtEtot += sqrtE;
                vcr = *(pos[icl]);
                if (nelem == 720) {
                    if (vcr.Y() > 0.) {
                        fSqrtEtUp += sqrtE;
                        *fMeanPosUp += vcr * sqrtE;
                    } else {
                        fSqrtEtDn += sqrtE;
                        *fMeanPosDn += vcr * sqrtE;
                    }
                }
                *fMeanPosition += vcr * sqrtE;
                k++;
            }
        }
    }

    const Double_t difmax1 = 18., difmax2 = 8.;
    const Double_t opangl1 = 13.;
    Double_t oang, difmax, opangl;
    for (i = 0; i < nhits; i++) {
        if ((icl = hits[i]) == ENullHit)
            continue; // was previously counted
        vcl = (*fMeanPosition) * (1. / fSqrtEtot);
        vcr = *(pos[icl]);
        energyi = energy[icl];
        if (nelem == 720) {
            oang = vcl.Angle(vcr) * TMath::RadToDeg();
            opangl = opangl1;
            if (oang > opangl) {
                continue;
            }
        } else {
            difmax = difmax1 + difmax2 * fEnergy / 1000.;
            vdif = vcl - vcr;
            if (vdif.Mag() > difmax) {
                continue;
            }
        }
        hits[i] = ENullHit; // so its not double counted
        fHits[k] = icl;     // add to cluster hits collection
        sqrtE = pow(energyi, Peng);
        fEnergy += energyi;
        if (nelem == 720)
            wtime = energyi;
        else
            wtime = sqrtE;
        fTime += time[icl] * wtime;
        fSqrtEtot += sqrtE;
        *fMeanPosition += vcr * sqrtE;
        k++;
    }

    fNhits = k;
    fHits[k] = EBufferEnd;
    *fMeanPosition =
            (*fMeanPosition) * (1. / fSqrtEtot); // normalise weighted mean
    if (nelem == 720) {
        if (fSqrtEtUp > 0.)
            *fMeanPosUp = (*fMeanPosUp) * (1. / fSqrtEtUp);
        if (fSqrtEtDn > 0.)
            *fMeanPosDn = (*fMeanPosDn) * (1. / fSqrtEtDn);
    }
    fTheta = TMath::RadToDeg() * fMeanPosition->Theta();
    fPhi = TMath::RadToDeg() * fMeanPosition->Phi();
    fCentralFrac = energy[fIndex] / fEnergy;
    if (nelem == 720)
        fTime /= fEnergy;
    else
        fTime /= fSqrtEtot;
}

Double_t Clustering_Sergey::Impl::HitCluster_t::ClusterRadius(TA2ClusterDetector *cldet) {
    // Determine the boundary of the cluster the local total energy
    // and the sqrt(energy)-weighted centre-of-gravity vector

    UInt_t i, ind;
    TVector3 vcl, vcr, vdif;

    Double_t *energy = cldet->GetEnergy();
    TVector3 **pos = cldet->GetPosition();
    UInt_t nelem = cldet->GetNelement();
    fRadius = 0.;
    Double_t thet = fTheta * TMath::DegToRad(), phi = fPhi * TMath::DegToRad();
    Double_t vdev, vmag = 1.;
    if (nelem != 720) {
        vcr = *(pos[1]);
        vmag = vcr.Z();
    }
    vcl.SetMagThetaPhi(vmag / cos(thet), thet, phi);
    for (i = 0; i < fNhits; i++) {
        ind = fHits[i];
        vcr = *(pos[ind]);
        if (nelem == 720) {
            vdev = vcl.Angle(vcr) * TMath::RadToDeg();
            if (vdev > 0.)
                fRadius += energy[ind] * vdev * vdev;
        } else {
            vdif = vcl - vcr;
            vdev = vdif.Mag();
            if (vdev > 0.)
                fRadius += energy[ind] * vdev * vdev;
        }
    }
    if (fRadius > 0.)
        fRadius = sqrt(fRadius / fEnergy);
    return fRadius;
}

Bool_t Clustering_Sergey::Impl::HitCluster_t::ClusterDetermine2(TA2ClusterDetector *cldet) {
    // Determine the boundary of the cluster the local total energy
    // and the sqrt(energy)-weighted centre-of-gravity vector

    const Double_t Peng = 2. / 3.;
    const Double_t difmax1 = 24., difmax2 = 10.;
    const Double_t opangl1 = 30., opangl2 = 7.;
    Double_t difmax, opangl;
    UInt_t k, ind, l;
    Double_t sqrtE;

    Double_t energyi, oang, wtime;
    Double_t *energy = cldet->GetEnergy();
    Double_t *time = cldet->GetTime();
    UInt_t *hits = cldet->GetTempHits2();
    TVector3 **pos = cldet->GetPosition();
    UInt_t nhits = cldet->GetNhits();
    UInt_t nelem = cldet->GetNelement();
    TVector3 vcr, vcl, vdif;

    UInt_t nhitold = fNhits;

    *fMeanPosition = (*fMeanPosition) * fSqrtEtot; // unnormalise weighted mean
    if (nelem == 720)
        fTime *= fEnergy;
    else
        fTime *= fSqrtEtot;

    if (nelem == 720) {
        if (fSqrtEtUp > 0.)
            *fMeanPosUp = (*fMeanPosUp) * fSqrtEtUp;
        if (fSqrtEtDn > 0.)
            *fMeanPosDn = (*fMeanPosDn) * fSqrtEtDn;
    }
    k = fNhits;

    for (l = 0; l < nhits; l++) {
        if ((ind = hits[l]) == ENullHit)
            continue; // was previously counted
        vcl = (*fMeanPosition) * (1. / fSqrtEtot);
        vcr = *(pos[ind]);
        energyi = energy[ind];
        if (nelem == 720) {
            oang = vcl.Angle(vcr) * TMath::RadToDeg();
            opangl = opangl1 + opangl2 * fEnergy / 1000.;
            if (oang > opangl) {
                continue;
            }
        } else {
            difmax = difmax1 + difmax2 * fEnergy / 1000.;
            vdif = vcl - vcr;
            if (vdif.Mag() > difmax) {
                continue;
            }
        }
        hits[l] = ENullHit; // so its not double counted
        fHits[k] = ind;     // add to cluster hits collection
        sqrtE = pow(energyi, Peng);
        fEnergy += energyi;
        if (nelem == 720)
            wtime = energyi;
        else
            wtime = sqrtE;
        fTime += time[ind] * wtime;
        fSqrtEtot += sqrtE;
        *fMeanPosition += vcr * sqrtE;
        k++;
    }
    fNhits = k;
    fHits[k] = EBufferEnd;
    *fMeanPosition =
            (*fMeanPosition) * (1. / fSqrtEtot); // normalise weighted mean
    if (nelem == 720) {
        if (fSqrtEtUp > 0.)
            *fMeanPosUp = (*fMeanPosUp) * (1. / fSqrtEtUp);
        if (fSqrtEtDn > 0.)
            *fMeanPosDn = (*fMeanPosDn) * (1. / fSqrtEtDn);
    }
    fCentralFrac = energy[fIndex] / fEnergy;
    if (nelem == 720)
        fTime /= fEnergy;
    else
        fTime /= fSqrtEtot;
    fTheta = TMath::RadToDeg() * fMeanPosition->Theta();
    fPhi = TMath::RadToDeg() * fMeanPosition->Phi();

    if (fNhits == nhitold)
        return EFalse;
    else
        return ETrue;
}

// dispatch PIMPL

Clustering_Sergey::Clustering_Sergey() :
    impl(std_ext::make_unique<Impl>(
             *ExpConfig::Setup::GetDetector<expconfig::detector::CB>(),
             *ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>()
             ))
{}

void Clustering_Sergey::Build(const ClusterDetector_t& clusterdetector,
                              const TClusterHitList& clusterhits,
                              TClusterList& clusters) const
{
    impl->Build(clusterdetector, clusterhits, clusters);
}

Clustering_Sergey::~Clustering_Sergey()
{}
