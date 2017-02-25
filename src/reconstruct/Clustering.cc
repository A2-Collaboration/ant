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

#include "TVector3.h"
#include "TMath.h"

struct Clustering_Sergey::Impl {
    enum {
        ENullHit = 0xFFFFFFFF,       // undefined hit index (end of hit buffer)
        EBufferEnd = 0xFFFFFFFF,     // end of file marker
    };

    struct TA2ClusterDetector;

    struct HitCluster_t {
        Double_t fEnergy;        // Total energy deposited in cluster
        Double_t fTheta;         // Cluster's theta
        Double_t fPhi;           // Cluster's phi
        UInt_t fNhits;          // # of hits in cluster

        Double_t GetEnergy() { return fEnergy; }
        Double_t GetTheta() { return fTheta; }
        Double_t GetPhi() { return fPhi; }
        UInt_t GetNhits() { return fNhits; }

        Double_t ClusterRadius(TA2ClusterDetector *cldet);
        void ClusterDetermine(TA2ClusterDetector *cldet);
        Bool_t ClusterDetermine2(TA2ClusterDetector *cldet);
    };


    struct TA2ClusterDetector {

        UInt_t *fTempHits;          // Element-Hit store
        Int_t* fHits;                          // array elements which fired in event
        UInt_t fNhits;                         // No. detector hits in event
        UInt_t fNCluster;           // # of clusters
        UInt_t fNelement;                      // # elements in array
        UInt_t *fTryHits;           // Element-Hit store
        UInt_t fMaxCluster;         // Max # of clusters
        Double_t* fEnergy;                     // stored hit energies if any
        UInt_t *fTempHits2;         // Element-Hit store
        HitCluster_t **fCluster;    // Clusters of hits
        Double_t *fClEnergyOR;      // OR of cluster energies
        Double_t *fTheta;           // theta of cluster hit
        Double_t *fPhi;             // phi of cluster hit
        UInt_t *fClustHit;          // Cluster indices
        UInt_t *fNClustHitOR;       // OR of #hits in individuyal clusters
        Double_t *fClRadius;        // cluster radius


        void DecodeCluster() {
            // Determine clusters of hits
            // Search around peak energies absorbed in individual crystals
            // Make copy of hits array as the copy will be altered

            memcpy(fTempHits, fHits, sizeof(UInt_t) * fNhits); // temp copy
            //  fNCluster = 0;
            //+
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
            UInt_t i, j, k, kmax, jmax, m, ntaken, nc;
            UInt_t nomit = 0;
            TVector3 vcl, vcln, vdif;

            //- hist
            // printf("Nelement Nhits %d %d\n",fNelement,fNhits);
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
                        // printf("%lf\n",fElement[k]->GetA1Value());
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
                                // fEthclsi = 25.+ 25. * Ecli / 1000.; // ver 1
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
                                        // printf("CB-omit %lf %lf\n",oang,fEthclsi);
                                        nomit++;
                                        fTryHits[jmax] = ENullHit;
                                        goto NEXTCR;
                                    }
                                } else {
                                    vdif = vcl - vcln;
                                    if (vdif.Mag() < difmax && Ecl < fEthclsi) {
                                        // printf("TAPS-omit %lf %lf\n",vdif.Mag(),fEthclsi);
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

            // printf("Nelement, Ncl %d %d\n",fNelement,fNCluster);
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
                    // printf("Ecl** %lf\n",fClEnergyOR[j]);
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
            // printf("Nelement, Ncl %d %d\n",fNelement,fNCluster);
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
                        // printf("Ecl*** %lf\n",Ecl);
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

            // printf("Nelement, Ncl %d %d\n",fNelement,fNCluster);

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
                    // printf("Ecl**** %lf\n",fClEnergyOR[j]);
                }
            }
        }



    };

    // implement stuff from "forward" declared HitCluster_t



    // internal dispatch to Acqu classes

    void Build(const ClusterDetector_t& clusterdetector,
               const TClusterHitList& clusterhits,
               TClusterList& clusters) const {

    }
};



// dispatch PIMPL

Clustering_Sergey::Clustering_Sergey() :
    impl(std_ext::make_unique<Impl>())
{}

void Clustering_Sergey::Build(const ClusterDetector_t& clusterdetector,
                              const TClusterHitList& clusterhits,
                              TClusterList& clusters) const
{
    impl->Build(clusterdetector, clusterhits, clusters);
}

Clustering_Sergey::~Clustering_Sergey()
{}
