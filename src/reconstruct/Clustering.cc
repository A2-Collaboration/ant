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
    enum { EFalse, ETrue };         // Logic...should use kTRUE, kFALSE

    struct TA2ClusterDetector;

    struct HitCluster_t {
        Double_t fEnergy;        // Total energy deposited in cluster
        Double_t fTheta;         // Cluster's theta
        Double_t fPhi;           // Cluster's phi
        UInt_t fNhits;          // # of hits in cluster
        UInt_t fIndex;           // index of central element
        Double_t fSqrtEtot;      // Sum of sqrt( energy[i] )
        Double_t fTime;          // Time of max-energy element of cluster
        Double_t fSqrtEtUp;
        Double_t fSqrtEtDn;
        TVector3 *fMeanPosUp;    // energy-weighted mean pos. of clust.
        TVector3 *fMeanPosDn;    // energy-weighted mean pos. of clust.
        TVector3 *fMeanPosition; // energy-weighted mean pos. of clust.
        UInt_t *fHits;          // indices of hit elements
        UInt_t fNNeighbour;     // # neighbour elements in array
        Double_t fCentralFrac;   // Fractional energy in central crystal
        UInt_t *fNeighbour;      // indices of neighbouring elements
        Double_t fRadius;       // effective radius of cluster

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
        Double_t* fTime;                       // stored hit times
        TVector3** fPosition;                  // stored hit positions

        Double_t* GetEnergy(){ return fEnergy; }          // ptr to energy array
        Double_t* GetTime(){ return fTime; }              // ptr to time array
        UInt_t *GetTempHits2() { return fTempHits2; }
        TVector3** GetPosition(){ return fPosition; }     // ptr position array
        UInt_t GetNhits(){ return fNhits; }              // No. hits in event
        UInt_t GetNelement(){ return fNelement; }        // max detector elements

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

void Clustering_Sergey::Impl::HitCluster_t::ClusterDetermine(TA2ClusterDetector *cldet) {
  // Determine the boundary of the cluster the local total energy
  // and the sqrt(energy)-weighted centre-of-gravity vector

  const Double_t Peng = 2. / 3.;
  UInt_t i, j, k, m, icl;
  Double_t energyi, wtime;
  Double_t *energy = cldet->GetEnergy();
  Double_t *time = cldet->GetTime();
  UInt_t *hits = cldet->GetTempHits2();
  TVector3 **pos = cldet->GetPosition();
  UInt_t nhits = cldet->GetNhits();
  UInt_t nelem = cldet->GetNelement();
  TVector3 vcr, vcl, vdif;
  //  for(m=0;m<fNNeighbour;m++) printf("%d ",fNeighbour[m]);

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
  // printf("fEnergy,fTime %d %lf %lf\n",nelem,fEnergy,fTime);
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
        // printf("%lf %lf\n",energyi,time[icl]);
        sqrtE = pow(energyi, Peng);
        if (energyi > 2000.)
          printf("energy[j] = %lf, %d, %d\n", energy[icl], icl, nhits);
        fEnergy += energyi;
        if (nelem == 720)
          wtime = energyi;
        else
          wtime = sqrtE;
        fTime += time[icl] * wtime;
        // if ( nelem!=720 ) printf("time %f\n",time[icl]);
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
  // const Double_t opangl1 = 14., opangl2 = 16.; // version 1
  // const Double_t opangl1 = 14., opangl2 = 0.;
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
      // opangl = opangl1 + opangl2 * fEnergy / 1000.;
      opangl = opangl1;
      if (oang > opangl) {
        // printf("CB %lf %lf %lf %lf\n",fEnergy,energyi,oang,opangl);
        continue;
      }
    } else {
      difmax = difmax1 + difmax2 * fEnergy / 1000.;
      vdif = vcl - vcr;
      if (vdif.Mag() > difmax) {
        // printf("TAPS %lf %lf %lf %lf\n",fEnergy,energyi,vdif.Mag(),difmax);
        continue;
      }
    }
    hits[i] = ENullHit; // so its not double counted
    fHits[k] = icl;     // add to cluster hits collection
    sqrtE = pow(energyi, Peng);
    fEnergy += energyi;
    if (nelem == 720)
      wtime = energyi;
    // if ( nelem==720 ) wtime = sqrtE;
    // else wtime = 1.;
    // else wtime = energyi;
    else
      wtime = sqrtE;
    fTime += time[icl] * wtime;
    fSqrtEtot += sqrtE;
    *fMeanPosition += vcr * sqrtE;
    k++;
    // printf("%d %lf %lf\n",k,energyi,fEnergy);
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
  // if ( nelem==720 ) fTime /= fSqrtEtot;
  // else fTime /= (Double_t)fNhits;
  // else fTime /= fEnergy;
  else
    fTime /= fSqrtEtot;
}

Double_t Clustering_Sergey::Impl::HitCluster_t::ClusterRadius(TA2ClusterDetector *cldet) {
  // Determine the boundary of the cluster the local total energy
  // and the sqrt(energy)-weighted centre-of-gravity vector

  UInt_t i, m, ind;
  // TVector3 vcl, vcr, vcr1, vdif, vycorr(0.,0.476,0.);
  // TVector3 vcl, vcr, vdif, vycorr(0.,0.576,0.);
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
  // printf("%lf %lf %lf\n",vmag,thet,cos(thet));
  vcl.SetMagThetaPhi(vmag / cos(thet), thet, phi);
  for (i = 0; i < fNhits; i++) {
    ind = fHits[i];
    vcr = *(pos[ind]);
    // printf("R calc: ind, energ %d %lf\n",ind,energy[ind]);
    if (nelem == 720) {
      vdev = vcl.Angle(vcr) * TMath::RadToDeg();
      if (vdev > 0.)
        fRadius += energy[ind] * vdev * vdev;
      // printf("dangl %lf %lf\n", energy[ind],
      // vcl.Angle(vcr)*TMath::RadToDeg());
    } else {
      vdif = vcl - vcr;
      // printf("taps cl %lf %lf %lf\n",vcl.X(),vcl.Y(),vcl.Z());
      // printf("taps cr %lf %lf %lf\n",vcr.X(),vcr.Y(),vcr.Z());
      vdev = vdif.Mag();
      if (vdev > 0.)
        fRadius += energy[ind] * vdev * vdev;
    }
  }
  if (fRadius > 0.)
    fRadius = sqrt(fRadius / fEnergy);
  // printf("R, E, Nh = %lf, %lf, %d\n",fRadius,fEnergy,fNhits);
  return fRadius;
}

Bool_t Clustering_Sergey::Impl::HitCluster_t::ClusterDetermine2(TA2ClusterDetector *cldet) {
  // Determine the boundary of the cluster the local total energy
  // and the sqrt(energy)-weighted centre-of-gravity vector

  const Double_t Peng = 2. / 3.;
  const Double_t difmax1 = 24., difmax2 = 10.;
  const Double_t opangl1 = 30., opangl2 = 7.;
  // const Double_t difmax=21.;
  // const Double_t opangl = 30.;
  Double_t difmax, opangl;
  UInt_t k, ind, l;
  Double_t sqrtE;

  // static UInt_t nhitma =0;
  //
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
  // if ( nelem==720 ) fTime *= fSqrtEtot;
  // else fTime *= (Double_t)fNhits;
  // else fTime *= fEnergy;
  else
    fTime *= fSqrtEtot;

  if (nelem == 720) {
    if (fSqrtEtUp > 0.)
      *fMeanPosUp = (*fMeanPosUp) * fSqrtEtUp;
    if (fSqrtEtDn > 0.)
      *fMeanPosDn = (*fMeanPosDn) * fSqrtEtDn;
  }
  //   printf("\n fIndex, fNeighbour = %d, %d\n",fIndex, fNeighbour[0]);
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
        // printf("CB2 %lf %lf %lf %lf\n",fEnergy,energyi,oang,opangl);
        continue;
      }
    } else {
      difmax = difmax1 + difmax2 * fEnergy / 1000.;
      vdif = vcl - vcr;
      if (vdif.Mag() > difmax) {
        // printf("TAPS2 %lf %lf %lf %lf\n",fEnergy,energyi,vdif.Mag(),difmax);
        continue;
      }
    }
    hits[l] = ENullHit; // so its not double counted
    fHits[k] = ind;     // add to cluster hits collection
    sqrtE = pow(energyi, Peng);
    fEnergy += energyi;
    if (nelem == 720)
      wtime = energyi;
    // if ( nelem==720 ) wtime = sqrtE;
    // else wtime = 1.;
    // else wtime = energyi;
    else
      wtime = sqrtE;
    fTime += time[ind] * wtime;
    fSqrtEtot += sqrtE;
    *fMeanPosition += vcr * sqrtE;
    k++;
  }
  fNhits = k;
  // if ( fNhits>nhitma ) { nhitma=fNhits ; printf("%d\n",nhitma);}
  //  printf("\n fNhits, nhits = %d, %d\n",fNhits,nhits);
  // for(m=0;m<nhits;m++) printf("%d ",hits[m]);
  fHits[k] = EBufferEnd;
  *fMeanPosition =
      (*fMeanPosition) * (1. / fSqrtEtot); // normalise weighted mean
  if (nelem == 720) {
    if (fSqrtEtUp > 0.)
      *fMeanPosUp = (*fMeanPosUp) * (1. / fSqrtEtUp);
    if (fSqrtEtDn > 0.)
      *fMeanPosDn = (*fMeanPosDn) * (1. / fSqrtEtDn);
  }
  //  if( fNhits != nhitold )  printf("Wider Mean X= %lf\n",fMeanPosition->X());
  fCentralFrac = energy[fIndex] / fEnergy;
  if (nelem == 720)
    fTime /= fEnergy;
  // if ( nelem==720 ) fTime /= fSqrtEtot;
  // else fTime /= (Double_t)fNhits;
  // else fTime /= fEnergy;
  else
    fTime /= fSqrtEtot;
  fTheta = TMath::RadToDeg() * fMeanPosition->Theta();
  fPhi = TMath::RadToDeg() * fMeanPosition->Phi();

  //  if( fNhits != nhitold )  printf("nhitold,fNhits = %d,
  //  %d\n",nhitold,fNhits);
  //  printf("nhitold,fNhits = %d, %d\n",nhitold,fNhits);
  if (fNhits == nhitold)
    return EFalse;
  else
    return ETrue;
}

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
