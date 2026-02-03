#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TChain.h>
#include <vector>
#include <TVector3.h>
#include <fstream>
#include <string>
#include <iostream>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <variant>
#include <stdexcept>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TMath.h>
#include <algorithm>


/*Function to read file.root*/
std::vector<std::string> readFileList(const std::string& filename) {
    std::vector<std::string> files;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line)) {
        if (!line.empty())
            files.push_back(line);
    }

    return files;
}

/*Function to split a TString into vector<TString> using ';' as delimiter*/
std::vector<TString> SplitTString(const TString &str) {
    std::vector<TString> parts;
    TObjArray* tokens = str.Tokenize(";");

    for (int i = 0; i < tokens->GetEntries(); ++i) {
        TObjString* objStr = dynamic_cast<TObjString*>(tokens->At(i));
        if (objStr) {
            parts.emplace_back(objStr->GetString());
        }
    }

    delete tokens;
    return parts;
}

/*Function to classify the process type*/
TString FilterProcs(const std::vector<TString>& input) {
    for (const auto& element : input) {
        if (element.BeginsWith("proc:")) {
            return element; 
        }
    }
    return "";
}

/*Function to classify the target type*/
TString FilterTGT(const std::vector<TString>& input) {
    for (const auto& element : input) {
        if (element.BeginsWith("tgt:")) {
            return element; 
        }
    }
    return "";
}

/*funtion to classify the events*/   
int Evt_cat(Int_t nu_type , TString target, TString st_proc_type, TString Where_int) {

    if (nu_type == -14 && target == "tgt:1000010010" && st_proc_type == "proc:Weak[CC],QES" && Where_int == "C3H6") // anti mu_nu on H in plastic
        return 1;

    else if (nu_type == -14 && target == "tgt:1000060120" && st_proc_type == "proc:Weak[CC],QES" && Where_int == "C3H6" ) // anti mu_nu on C in plastic
        return 2;

    else if (nu_type == -14 && target == "tgt:1000010010" && st_proc_type != "proc:Weak[CC],QES" && Where_int == "C3H6" ) // anti mu_nu on H in plastic non QES
        return 3;

    else if (nu_type == -14 && target == "tgt:1000060120" && st_proc_type == "proc:Weak[CC],QES" && Where_int == "Graphite" ) // anti mu_nu on C in graphite
        return 4;
     
     else if (nu_type == -14 && target == "tgt:1000010010" && st_proc_type != "proc:Weak[CC],QES" && Where_int != "Graphite" && Where_int != "C3H6") // on H everywhere else
        return 5;

     else return 0;
} 

/*function to obtain the PDG of the primaries*/
ROOT::VecOps::RVec<Int_t> Primaries_PDG(const ROOT::VecOps::RVec<Int_t>& prim, const ROOT::VecOps::RVec<Int_t>& pdg) {
    ROOT::VecOps::RVec<Int_t> OUT;

    for (size_t i = 0; i < prim.size(); ++i) {
        if (prim[i] == 1) {
            OUT.push_back(pdg[i]);
        }
    }

    return OUT;
}

//function to obtain the momentum of the primaries
ROOT::VecOps::RVec<Double_t> Primaries_momentum(const ROOT::VecOps::RVec<Int_t>& primaries, const ROOT::VecOps::RVec<Double_t>& mom) {
    //ROOT::VecOps::RVec<Int_t> output_vec(primaries.size());
    ROOT::VecOps::RVec<Double_t> output;
    for (size_t i = 0; i < primaries.size(); ++i) {
        if (primaries[i] == 1) {
            output.push_back(mom[i]);  
        }
    }  
    return output;
}

//function to obtain the momentum of the muon
ROOT::VecOps::RVec<Double_t> mu_momentum(const ROOT::VecOps::RVec<Int_t>& pdg, const ROOT::VecOps::RVec<Double_t>& mom) {
    ROOT::VecOps::RVec<Double_t> matchedMom;

    for (size_t i = 0; i < pdg.size(); ++i) {
        if (pdg[i] == -13) {
            matchedMom.push_back(mom[i]);
        }
    }

    return matchedMom;  
}

//function to obtain the momentum of the neutron
ROOT::VecOps::RVec<Double_t> n_momentum(const ROOT::VecOps::RVec<Int_t>& pdg, const ROOT::VecOps::RVec<Double_t>& mom) {
    ROOT::VecOps::RVec<Double_t> matchedMom;

    for (size_t i = 0; i < pdg.size(); ++i) {
        if (pdg[i] == 2112) {
            matchedMom.push_back(mom[i]);
        }
    }

    return matchedMom;  
}


//funtion to know in which part of the detector the interaction takes place 
TString Where_int(double x, double y, double z) {
    if (!gGeoManager) return "NO_GEOMETRY";

    TGeoNode* node = gGeoManager->FindNode(x, y, z);
    if (!node) return "OUTSIDE";

    TGeoVolume* volume = node->GetVolume();
    if (!volume) return "NO_VOLUME";
    //return node->GetName();
    
    TGeoMedium* medium = volume->GetMedium();
    if (!medium) return "NO_MEDIUM";

    TGeoMaterial* material = medium->GetMaterial();
    if (!material) return "NO_MATERIAL";

    return material->GetName();
    
}

//function to access cell Energies
std::vector<std::vector<double>> ExtractEnergies(const ROOT::VecOps::RVec<cluster>& clusters) {
    std::vector<std::vector<double>> energies_per_cluster;

    for (const auto& clu : clusters) {
         if (clu.reco_cells.empty()) continue;  // Skip empty clusters
        std::vector<double> energies;
        energies.reserve(clu.reco_cells.size());  // Pre-allocate memory
        for (const auto& cell : clu.reco_cells) {
            energies.push_back(cell.e); 
        }
        energies_per_cluster.push_back(std::move(energies));
    }
    return energies_per_cluster;
}

//functions to access cell position  
std::vector<std::vector<double>> Extract_X_cell(const ROOT::VecOps::RVec<cluster>& clusters) {
    std::vector<std::vector<double>> position_per_cluster;

    for (const auto& clu : clusters) {
         if (clu.reco_cells.empty()) continue;  // Skip empty clusters
        std::vector<double> position;
        position.reserve(clu.reco_cells.size());  // Pre-allocate memory
        for (const auto& cell : clu.reco_cells) {
            position.push_back(cell.x); 
        }
        position_per_cluster.push_back(std::move(position));
    }
    return position_per_cluster;
}

std::vector<std::vector<double>> Extract_Y_cell(const ROOT::VecOps::RVec<cluster>& clusters) {
    std::vector<std::vector<double>> position_per_cluster;

    for (const auto& clu : clusters) {
         if (clu.reco_cells.empty()) continue;  // Skip empty clusters
        std::vector<double> position;
        position.reserve(clu.reco_cells.size());  // Pre-allocate memory
        for (const auto& cell : clu.reco_cells) {
            position.push_back(cell.y); 
        }
        position_per_cluster.push_back(std::move(position));
    }
    return position_per_cluster;
}

std::vector<std::vector<double>> Extract_Z_cell(const ROOT::VecOps::RVec<cluster>& clusters) {
    std::vector<std::vector<double>> position_per_cluster;

    for (const auto& clu : clusters) {
         if (clu.reco_cells.empty()) continue;  // Skip empty clusters
        std::vector<double> position;
        position.reserve(clu.reco_cells.size());  // Pre-allocate memory
        for (const auto& cell : clu.reco_cells) {
            position.push_back(cell.z); 
        }
        position_per_cluster.push_back(std::move(position));
    }
    return position_per_cluster;
}

//funtion to extract the cell id referred to the cluster 
std::vector<int> CellID(const ROOT::VecOps::RVec<cluster>& clusters) {
    std::vector<int> cell_ID;

    for (std::size_t i = 0; i < clusters.size(); ++i) {
        const auto& clu = clusters[i];
        std::size_t n_cells = clu.reco_cells.size();

        // Append the cluster index 'i', repeated n_cells times
        cell_ID.insert(cell_ID.end(), n_cells, static_cast<int>(i));
    }

    return cell_ID;
}

/*function to get the particles with a track*/
std::vector<int> Particles_with_track(const ROOT::VecOps::RVec<particle>& particles){

    std::vector<int> particles_with_track;

        for (const auto& p : particles) {
        if (p.primary == 1 && p.has_track) {  
            particles_with_track.push_back(p.pdg);
        }
    }

    return particles_with_track;
}

/*function to get the nuber of points of the track*/
// std::vector<int> N_points(ROOT::VecOps::RVec<track>& track){
//     std::vector<int> N_points;

//     for (std::size_t i = 0; i < track.size(); ++i){
//          N_points.push_back(track[i].n_points);
//     }

//     return N_points;
// }
std::vector<int> N_points(const ROOT::VecOps::RVec<particle>& particles){
    std::vector<int> N_points;

        for (const auto& p : particles) {
        if (p.primary == 1 && p.has_track) {  
            N_points.push_back(p.tr.n_points);
        }
    }
    return N_points;
}

/*function to get the chi2 of the track*/
// std::vector<double> Chi2(ROOT::VecOps::RVec<track>& track){
//     std::vector<double> chi2;

//     for (std::size_t i = 0; i < track.size(); ++i){
//          chi2.push_back(track[i].chi2_cr);
//     }

//     return chi2;
// }
std::vector<double> Chi2(const ROOT::VecOps::RVec<particle>& particles){
    std::vector<double> chi2;

        for (const auto& p : particles) {
        if (p.primary == 1 && p.has_track) {  
            chi2.push_back(p.tr.chi2_cr);
        }
    }
    return chi2;
}

/*function to get the neutron cluster*/
std::vector<double> Neutron_cluster(const ROOT::VecOps::RVec<particle>& particles){
    std::vector<double> n_cl(5,0.0);

        for (const auto& p : particles) {
        if (p.primary == 1 && p.has_cluster && p.pdg == 2112) {  
            n_cl[0] = p.cl.x;
            n_cl[1] = p.cl.y;
            n_cl[2] = p.cl.z;
            n_cl[3] = p.cl.t;
            n_cl[4] = p.cl.e;
            break; //to take just the firt 
        }
    }
    return n_cl;
}



/*function to reconstruct Neutrino-energy*/
double NeutrinoEnergy(ROOT::VecOps::RVec<double> E_mu, ROOT::VecOps::RVec<double> px_mu, ROOT::VecOps::RVec<double> py_mu, ROOT::VecOps::RVec<double> pz_mu,
                      Double_t px_nu , Double_t py_nu, Double_t pz_nu) {

    const double M_n = 939.5654205;  // Neutron mass (MeV/c^2)
    const double M_p = 938.27208816;  // Proton mass (MeV/c^2)
    const double m_mu = 105.6583755; // Muon mass (MeV/c^2)

    double mu_mom = TMath::Sqrt(px_mu[0]*px_mu[0] + py_mu[0]*py_mu[0] + pz_mu[0]*pz_mu[0]);
    double nu_mom = TMath::Sqrt(px_nu*px_nu + py_nu*py_nu+ pz_nu*pz_nu);
    double mu_versor_x = px_mu[0] / mu_mom;
    double mu_versor_y = py_mu[0] / mu_mom;
    double mu_versor_z = pz_mu[0] / mu_mom;
    double nu_versor_x = px_nu / nu_mom;
    double nu_versor_y = py_nu / nu_mom;
    double nu_versor_z = pz_nu / nu_mom;


    double numerator = (M_n * M_n) - (m_mu * m_mu) - (M_p * M_p) + (2.0 * E_mu[0] * M_p);
    double denominator = 2.0 * (M_p - E_mu[0] + mu_mom * (mu_versor_x*nu_versor_x + mu_versor_y*nu_versor_y + mu_versor_z*nu_versor_z));

    return numerator/denominator;
}

 
/*function to get the neutrino reco momentum*/
ROOT::VecOps::RVec<Double_t> NeutrinoMomentum(const double nu_Ereco) {
    ROOT::VecOps::RVec<Double_t> nu_Preco(3);

    // Beam direction
    const double theta = -5.3 * TMath::Pi() / 180.0;
                          
    double versor_x = 0.0;
    double versor_y = TMath::Sin(theta);
    double versor_z = TMath::Cos(theta);

    // Neutrino energy (magnitude)
    double mod = nu_Ereco;

    // // Generatore random locale (puoi fissare il seed per riproducibilità)
    // static TRandom3 randGen(0);

    // // Distribuzione gaussiana centrata a 0 con sigma = 20/3 ≈ 6.67 MeV
    // // → 99.7 dei valori cade tra -20 e +20 MeV (±3sigma)
    // double px_rand = randGen.Gaus(0.0, 15.0 / 3.0);

    // Imposta le componenti del vettore
    // nu_Preco[0] = px_rand;     
    nu_Preco[0] = mod * versor_x;
    nu_Preco[1] = mod * versor_y;      
    nu_Preco[2] = mod * versor_z;      

    return nu_Preco;
}

/*function to get neutron predicted momentum direction*/
Double_t Neutron_predicted_momentumX(ROOT::VecOps::RVec<Double_t> nu_Preco, ROOT::VecOps::RVec<Double_t> mu_preco){
    Double_t n_px_pred;
    n_px_pred = nu_Preco[0] - mu_preco[0];
    return n_px_pred;
}
Double_t Neutron_predicted_momentumY(ROOT::VecOps::RVec<Double_t> nu_Preco, ROOT::VecOps::RVec<Double_t> mu_preco){
    Double_t n_py_pred;
    n_py_pred = nu_Preco[1] - mu_preco[0];
    return n_py_pred;
}
Double_t Neutron_predicted_momentumZ(ROOT::VecOps::RVec<Double_t> nu_Preco, ROOT::VecOps::RVec<Double_t> mu_preco){
    Double_t n_pz_pred;
    n_pz_pred = nu_Preco[2] - mu_preco[0];
    return n_pz_pred;
}

/*function to get neutron predicted energy*/
double Neutron_predicted_energy( const double nu_Ereco, ROOT::VecOps::RVec<double> E_mu){ 

    const double M_p = 938.27208816;
    double n_E_pred = nu_Ereco + M_p - E_mu[0];

    return n_E_pred;
}


/*function to predict the hits in the ECAl barrel + endcaps*/
ROOT::VecOps::RVec<double> Predicted_hit(double px, double py, double pz, double x0, double y0, double z0){

    // ECAL cylinder dimensions (mm)
    const double r_out = 2230.0;
    const double r_in  = 2000.0;
    const double x_min = -2150.0;
    const double x_max =  2150.0;

    // Cylinder center
    const double cx = 0.0;
    const double cy = -2384.73;
    const double cz = 23910.0;

    // Normalize direction
    double norm = std::sqrt(px*px + py*py + pz*pz);
    if (norm == 0) return {0,0,0,0,0,0};

    double nx = px / norm;
    double ny = py / norm;
    double nz = pz / norm;

                          
    //cilinder intersection
    auto intersect = [&](double R)
    {
        // Shift ray origin into cylinder frame
        double Y = y0 - cy;
        double Z = z0 - cz;

        // Quadratic coefficients for (y^2 + z^2 = R^2)
        double A = ny*ny + nz*nz;
        double B = 2*(Y*ny + Z*nz);
        double C = Y*Y + Z*Z - R*R;

        double D = B*B - 4*A*C;
        if (D < 0 || A == 0)
            return std::vector<double>{};

        double sqrtD = std::sqrt(D);
        double t1 = (-B + sqrtD)/(2*A);
        double t2 = (-B - sqrtD)/(2*A);

        double t = 1e99;
        if (t1 > 0) t = std::min(t, t1);
        if (t2 > 0) t = std::min(t, t2);

        if (t == 1e99)
            return std::vector<double>{};

        // Intersection point
        double xi = x0 + t*nx;
        double yi = y0 + t*ny;
        double zi = z0 + t*nz;
        return std::vector<double>{xi, yi, zi};
    };
                          

    // plane intersection
    auto intersect_plane = [&](double x_plane) {
                          
        if (nx == 0) return std::vector<double>{0.0,0.0,0.0}; 
        double t = (x_plane - x0) / nx;
        if (t <= 0) return std::vector<double>{0.0,0.0,0.0};  // plane behind start

        double xi = x_plane;
        double yi = y0 + t*ny;
        double zi = z0 + t*nz;
        return std::vector<double>{xi, yi, zi};
    };

                      
    // Compute inner & outer intersections with cylinder
    auto in  = intersect(r_in);
    auto out = intersect(r_out);

    // If intersection missing
    if (in.empty())  in  = {0,0,0};
    if (out.empty()) out = {0,0,0};

    ROOT::VecOps::RVec<double> result;
    result.insert(result.end(), in.begin(), in.end());
    result.insert(result.end(), out.begin(), out.end());
    
    //Endcpas
        if (result[0] < x_min) {
        double xA = -1690;
        double xB = -1690 - 460.0;

        auto pA = intersect_plane(xA);
        auto pB = intersect_plane(xB);

        ROOT::VecOps::RVec<double> res_planes;
        res_planes.insert(res_planes.end(), pA.begin(), pA.end());
        res_planes.insert(res_planes.end(), pB.begin(), pB.end());
        return res_planes;  
    }

        if (result[0] > x_max) {
        double xA = 1690;
        double xB = 1690 + 460.0;

        auto pA = intersect_plane(xA);
        auto pB = intersect_plane(xB);

        ROOT::VecOps::RVec<double> res_planes;
        res_planes.insert(res_planes.end(), pA.begin(), pA.end());
        res_planes.insert(res_planes.end(), pB.begin(), pB.end());
        return res_planes;  
    }                                 
                                           
    return result;
                          
}


/*function to obtain the predicted time (time of flight from the vertex to the center of the cluster)*/
ROOT::VecOps::RVec<double> Predicted_time(const Double_t n_E_pred, const Double_t n_px_pred, const Double_t n_py_pred, const Double_t n_pz_pred, 
                            const double t0, const double x0, const double y0, const double z0,
                            const ROOT::VecOps::RVec<double> x_ecal, const ROOT::VecOps::RVec<double> y_ecal, const ROOT::VecOps::RVec<double> z_ecal){
                                
                                ROOT::VecOps::RVec<double> time(x_ecal.size());                                
                                double beta = TMath::Sqrt(n_px_pred*n_px_pred + n_py_pred*n_py_pred + n_pz_pred*n_pz_pred) / n_E_pred;
                                const double c = 299.7924580;//mm/ns

                                for (size_t i = 0; i < x_ecal.size(); ++i){
                                    double s = std::sqrt((x_ecal[i] - x0)*(x_ecal[i] - x0) + (y_ecal[i] - y0)*(y_ecal[i] - y0) + (z_ecal[i] - z0)*(z_ecal[i] - z0));
                                    time[i] = s/(beta * c );
                                }
                                return time;
}


/*function to get the true neutron hit*/
// ROOT::VecOps::RVec<TLorentzVector> true_n_start(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const ROOT::VecOps::RVec<std::vector<TG4PrimaryParticle>>& particles)
// {
//     ROOT::VecOps::RVec<TLorentzVector> starts;

//     for (const auto& primaryVec : particles) {
//         for (const auto& particle : primaryVec) {

//             if (particle.GetPDGCode() == 2112) {
//                 int neutronTrackId = particle.GetTrackId();

//                 for (const auto& [name, hitSegments] : SegmentDetectors) {

//                     if (name == "EMCalSci") {

//                         for (const auto& hitSegment : hitSegments) {
//                             if (hitSegment.GetPrimaryId() == neutronTrackId) {
//                                 starts.push_back(hitSegment.GetStart());
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return starts;
//}

// ROOT::VecOps::RVec<double> true_n_start(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
// {
//     std::cout<<"sono_qui"<<std::endl;

//     ROOT::VecOps::RVec<double> starts;
//     //ROOT::VecOps::RVec<int> neutron_Id;

//         std::cout<<"sono dentro al primo for"<<std::endl;

//         for (const auto& particle : particles) {

//             if (particle.GetPDGCode() == 2112) {
//                 int neutronTrackId = particle.GetTrackId();

//                 for (const auto& [name, hitSegments] : SegmentDetectors) {

//                     if (name == "EMCalSci") {

//                         for (const auto& hitSegment : hitSegments) {
//                             if (hitSegment.GetPrimaryId() == neutronTrackId) {

//                                 const TLorentzVector& start = hitSegment.GetStart();
//                                 //starts.push_back(hitSegment.GetStart());
//                             starts.push_back(start.X());
//                             starts.push_back(start.Y());
//                             starts.push_back(start.Z());
//                             starts.push_back(start.T());
//                             }
//                         }
//                     }
//                 }
//             }
//         }


//     return starts;
// }
ROOT::VecOps::RVec<double> true_n_startX(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> starts;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& start = hitSegment.GetStart();
                                //starts.push_back(hitSegment.GetStart());
                            starts.push_back(start.X());
                            }
                        }
                    }
                }
            }
        }


    return starts;
}
ROOT::VecOps::RVec<double> true_n_startY(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> starts;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& start = hitSegment.GetStart();
                                //starts.push_back(hitSegment.GetStart());
                            starts.push_back(start.Y());
                            }
                        }
                    }
                }
            }
        }


    return starts;
}
ROOT::VecOps::RVec<double> true_n_startZ(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> starts;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& start = hitSegment.GetStart();
                                //starts.push_back(hitSegment.GetStart());
                            starts.push_back(start.Z());
                            }
                        }
                    }
                }
            }
        }


    return starts;
}
ROOT::VecOps::RVec<double> true_n_startT(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> starts;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& start = hitSegment.GetStart();
                                //starts.push_back(hitSegment.GetStart());
                            starts.push_back(start.T());
                            }
                        }
                    }
                }
            }
        }


    return starts;
}

//___________________________________________________________________
ROOT::VecOps::RVec<double> true_n_stopX(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> stops;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& stop = hitSegment.GetStop();
                                //starts.push_back(hitSegment.GetStart());
                            stops.push_back(stop.X());
                            }
                        }
                    }
                }
            }
        }


    return stops;
}

ROOT::VecOps::RVec<double> true_n_stopY(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> stops;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& stop = hitSegment.GetStop();
                                //starts.push_back(hitSegment.GetStart());
                            stops.push_back(stop.Y());
                            }
                        }
                    }
                }
            }
        }


    return stops;
}

ROOT::VecOps::RVec<double> true_n_stopZ(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> stops;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& stop = hitSegment.GetStop();
                                //starts.push_back(hitSegment.GetStart());
                            stops.push_back(stop.Z());
                            }
                        }
                    }
                }
            }
        }


    return stops;
}


ROOT::VecOps::RVec<double> true_n_stopT(const std::map<std::string, std::vector<TG4HitSegment>>& SegmentDetectors, const std::vector<TG4PrimaryParticle>& particles)
{
    ROOT::VecOps::RVec<double> stops;

        for (const auto& particle : particles) {

            if (particle.GetPDGCode() == 2112) {
                int neutronTrackId = particle.GetTrackId();

                for (const auto& [name, hitSegments] : SegmentDetectors) {

                    if (name == "EMCalSci") {

                        for (const auto& hitSegment : hitSegments) {
                            if (hitSegment.GetPrimaryId() == neutronTrackId) {

                                const TLorentzVector& stop = hitSegment.GetStop();
                                //starts.push_back(hitSegment.GetStart());
                            stops.push_back(stop.T());
                            }
                        }
                    }
                }
            }
        }
    return stops;
}

/*Function to get neutron trajectory points positions*/
ROOT::VecOps::RVec<double> true_x_pos(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
    
    ROOT::VecOps::RVec<double> true_x_pos;
    //ROOT::VecOps::RVec<TG4TrajectoryPoint> points;

        for (const auto& trajectory : trajectories) {

            if (trajectory.GetPDGCode() == 2112 && trajectory.GetParentId() == -1) {

                const auto& points = trajectory.Points;

                for (const auto& point : points) {
                    
                    const TLorentzVector& pos = point.GetPosition();
                    true_x_pos.push_back(pos.X());
                }
            }
        }
    return true_x_pos;
}
ROOT::VecOps::RVec<double> true_y_pos(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
    
    ROOT::VecOps::RVec<double> true_y_pos;
    //ROOT::VecOps::RVec<TG4TrajectoryPoint> points;

        for (const auto& trajectory : trajectories) {

            if (trajectory.GetPDGCode() == 2112 && trajectory.GetParentId() == -1) {

                const auto& points = trajectory.Points;

                for (const auto& point : points) {
                    
                    const TLorentzVector& pos = point.GetPosition();
                    true_y_pos.push_back(pos.Y());
                }
            }
        }
    return true_y_pos;
}
ROOT::VecOps::RVec<double> true_z_pos(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
    
    ROOT::VecOps::RVec<double> true_z_pos;
    //ROOT::VecOps::RVec<TG4TrajectoryPoint> points;

        for (const auto& trajectory : trajectories) {

            if (trajectory.GetPDGCode() == 2112 && trajectory.GetParentId() == -1) {

                const auto& points = trajectory.Points;

                for (const auto& point : points) {
                    
                    const TLorentzVector& pos = point.GetPosition();
                    true_z_pos.push_back(pos.Z());
                }
            }
        }
    return true_z_pos;
}
ROOT::VecOps::RVec<double> true_t_pos(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
    
    ROOT::VecOps::RVec<double> true_t;
    //ROOT::VecOps::RVec<TG4TrajectoryPoint> points;

        for (const auto& trajectory : trajectories) {

            if (trajectory.GetPDGCode() == 2112 && trajectory.GetParentId() == -1 ) {

                const auto& points = trajectory.Points;

                for (const auto& point : points) {
                    
                    const TLorentzVector& pos = point.GetPosition();
                    true_t.push_back(pos.T());
                }
            }
        }
    return true_t;
}



// ROOT::VecOps::RVec<TLorentzVector> true_trajectory_points(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
//     ROOT::VecOps::RVec<TLorentzVector> positions;
//     for (const auto& trajectory : trajectories) {
//         if (trajectory.GetPDGCode() == 2112) {
//             const auto& points = trajectory.GetPoints(); 
//             for (const auto& point : points) {
//                 positions.push_back(point.GetPosition());
//             }
//         }
//     }

//     return positions;
// }



void test_snap(){

    //leoading geometry files
    // TFile* file = TFile::Open("/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/SAND/SAND_0/sand-events.0.edep.root");
    TFile* file = TFile::Open("/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/SAND_innerVOL_10/SAND_innerVOL_10_0/sand-events.0.edep.root");
    TGeoManager* gGeoManager = dynamic_cast<TGeoManager*>(file->Get("EDepSimGeometry"));

    //read the correct files
    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/GSIM_reco_list.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/GSIM_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/GSIM_edep_list.txt");

    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_reco_list.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_edep_list.txt");

    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/EVT1_reco_list.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/EVT1_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/EVT1_edep_list.txt");
    
    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/tot_EVT_reco_list.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/tot_EVT_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/tot_EVT_edep_list.txt");
     
    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_reco_list2.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_gtrac_list2.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_edep_list2.txt");

    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_reco_listGS.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_gtrac_listGS.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_edep_listGS.txt");

    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_reco_listGS_D.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_gtrac_listGS_D.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/inner_edep_listGS_D.txt");

    /*PRODUZIONE GS DA 1M DI EVENTI FATTA CON 1M DI EVENTI*/
    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/GS_inner_reco_list.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/GS_inner_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/GS_inner_edep_list.txt");

    /*PRODUZIONE GS con pochi eventi, 9000*/
    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/MINI_reco_list.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/MINI_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/MINI_edep_list.txt");
    
    /*PRODUZIONE FATTA CON 2M DI EVENTI GS*/
    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/FINAL_inner_reco_list.txt"); CORROTTI
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/FINAL_inner_gtrac_list.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/FINAL_inner_edep_list.txt");

    // std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/battisti/camilla_tests/list/FINAL_inner_reco_list_cleaned.txt");
    // std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/battisti/camilla_tests/list/FINAL_inner_gtrac_list_cleaned.txt");
    // std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/battisti/camilla_tests/list/FINAL_inner_edep_list_cleaned.txt");


    std::vector<std::string> reco_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/AGG_inner_reco_list.txt");
    std::vector<std::string> gtrac_files  = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/AGG_inner_gtrac_list.txt");
    std::vector<std::string> edep_files   = readFileList("/storage/gpfs_data/neutrino/users/croselli/Chek_file/AGG_inner_edep_list.txt");


    //crating chains
    TChain* ch_tEvent = new TChain("tEvent");
    TChain* ch_gRooTracker = new TChain("gRooTracker");
    TChain* ch_EDepSimEvents = new TChain("EDepSimEvents");
    TChain* ch_tReco = new TChain("tReco");

    // Function to add files from a vector
    auto add_files = [](TChain* chain, const std::vector<std::string>& files) {
        for (const auto& f : files)
            chain->Add(f.c_str());
    };

    add_files(ch_tEvent, reco_files);
    add_files(ch_tReco, reco_files);
    add_files(ch_gRooTracker, gtrac_files);
    add_files(ch_EDepSimEvents, edep_files);


    // ch_tEvent->Add("/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/SAND_AN/SAND_AN_*/sand-events.*.reco.root");
   // ch_gRooTracker->Add("/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/SAND_AN/SAND_AN_*/sand-events.*.gtrac.root");
   // ch_EDepSimEvents->Add("/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/SAND_AN/SAND_AN_*/sand-events.*.edep.root");
   // ch_tReco->Add("/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/SAND_AN/SAND_AN_*/sand-events.*.reco.root");

// ch_tEvent->Add("@reco_list.txt");
// ch_gRooTracker->Add("@gtrac_list.txt");
// ch_EDepSimEvents->Add("@edep_list.txt");
// ch_tReco->Add("@reco_list.txt");

ch_tEvent->AddFriend(ch_gRooTracker, "gRooTracker");
ch_tEvent->AddFriend(ch_EDepSimEvents, "EDepSimEvents");
ch_tEvent->AddFriend(ch_tReco, "tReco");

//create the RDF
ROOT::RDataFrame df(*ch_tEvent);

/*
auto column_names = df.GetColumnNames();
    for (size_t i = 0; i < column_names.size(); ++i) {
        std::cout << i << ". " << column_names[i] << std::endl;
    }
*/

auto tRough = df.Define("vertex_x", [](const ROOT::VecOps::RVec<double>& pos) { return pos[0] * 1000; }, {"gRooTracker.EvtVtx"}) //mm to match the geometry
                .Define("vertex_y", [](const ROOT::VecOps::RVec<double>& pos) { return pos[1] * 1000 ; }, {"gRooTracker.EvtVtx"}) 
                .Define("vertex_z", [](const ROOT::VecOps::RVec<double>& pos) { return pos[2] * 1000; }, {"gRooTracker.EvtVtx"})
                .Define("vertex_t", [](const ROOT::VecOps::RVec<double>& pos) { return pos[3] ; }, {"gRooTracker.EvtVtx"})
                .Define("nu_type", [](const ROOT::VecOps::RVec<Int_t>& p) {return p[0]; }, {"gRooTracker.StdHepPdg"})
                // .Define("nu_E", [](const ROOT::VecOps::RVec<double>& mom4) {return mom4[3]; }, {"gRooTracker.StdHepP4"})
                //.Define("target", [](const ROOT::VecOps::RVec<Int_t>& p) {return p[1]; }, {"StdHepPdg"})
                //.Define("primaries_PDG", [](const ROOT::VecOps::RVec<Int_t>& p) {return ROOT::VecOps::RVec<Int_t>(p.begin() + 2, p.end());}, {"StdHepPdg"})
                .Define("SplitParts", SplitTString , {"gRooTracker.EvtCode.fString"})
                .Define("st_proc_type", FilterProcs, {"SplitParts"})
                .Define("target", FilterTGT, {"SplitParts"})                          
                .Define("Where_int", Where_int, {"vertex_x", "vertex_y","vertex_z"})
                .Define("Evt_cat", Evt_cat, {"nu_type", "target", "st_proc_type","Where_int"})
                .Define("primaries_PDG", Primaries_PDG, {"event.particles.primary","event.particles.pdg"})
                    .Define("primaries_px", Primaries_momentum, {"event.particles.primary","event.particles.pxtrue"})
                    .Define("primaries_py", Primaries_momentum, {"event.particles.primary","event.particles.pytrue"})
                    .Define("primaries_pz", Primaries_momentum, {"event.particles.primary","event.particles.pztrue"})
                    .Define("primaries_pxreco", Primaries_momentum, {"event.particles.primary","event.particles.pxreco"})
                    .Define("primaries_pyreco", Primaries_momentum, {"event.particles.primary","event.particles.pyreco"})
                    .Define("primaries_pzreco", Primaries_momentum, {"event.particles.primary","event.particles.pzreco"})
                    .Define("primaries_E", Primaries_momentum, {"event.particles.primary","event.particles.Etrue"})
                    .Define("primaries_Ereco", Primaries_momentum, {"event.particles.primary","event.particles.Ereco"})
                    .Define("mu_px", mu_momentum, {"primaries_PDG","primaries_px"})
                    .Define("mu_py", mu_momentum, {"primaries_PDG","primaries_py"})
                    .Define("mu_pz", mu_momentum, {"primaries_PDG","primaries_pz"})
                    .Define("mu_pxreco", mu_momentum, {"primaries_PDG","primaries_pxreco"})
                    .Define("mu_pyreco", mu_momentum, {"primaries_PDG","primaries_pyreco"})
                    .Define("mu_pzreco", mu_momentum, {"primaries_PDG","primaries_pzreco"})
                    .Define("mu_E", mu_momentum, {"primaries_PDG","primaries_E"})
                    .Define("mu_Ereco", mu_momentum, {"primaries_PDG","primaries_Ereco"})
                    .Define("mu_P", "sqrt(mu_px[0]*mu_px[0] + mu_py[0]*mu_py[0] + mu_pz[0]*mu_pz[0])")
                    .Define("mu_Preco", "sqrt(mu_pxreco[0]*mu_pxreco[0] + mu_pyreco[0]*mu_pyreco[0] + mu_pzreco[0]*mu_pzreco[0])")
                    .Define("n_px", n_momentum, {"primaries_PDG","primaries_px"})
                    .Define("n_py", n_momentum, {"primaries_PDG","primaries_py"})
                    .Define("n_pz", n_momentum, {"primaries_PDG","primaries_pz"})
                    .Define("n_P", "sqrt(n_px[0]*n_px[0] + n_py[0]*n_py[0] + n_pz[0]*n_pz[0])")
                    // .Define("n_pxreco", n_momentum, {"primaries_PDG","primaries_pxreco"})
                    // .Define("n_pyreco", n_momentum, {"primaries_PDG","primaries_pyreco"})
                    // .Define("n_pzreco", n_momentum, {"primaries_PDG","primaries_pzreco"})
                    .Define("n_E", n_momentum, {"primaries_PDG","primaries_E"})
                    .Define("n_Ereco", n_momentum, {"primaries_PDG","primaries_Ereco"})
                    .Define("nu_Ereco", NeutrinoEnergy, {"mu_Ereco", "mu_pxreco", "mu_pyreco", "mu_pzreco", "pxnu","pynu","pznu"})
                    //.Define("nu_Ereco", NeutrinoEnergy, {"mu_Ereco", "mu_pxreco", "mu_pyreco", "mu_pzreco", "pxnu","pynu","pznu","primaries_PDG","Enureco"})
                    .Define("n_E_pred", Neutron_predicted_energy, {"nu_Ereco", "mu_Ereco"})
                    // .Define("nu_Preco", NeutrinoMomentum, {"nu_Ereco","pxnu","pynu","pznu"})
                    .Define("nu_Preco", NeutrinoMomentum, {"nu_Ereco"})
                    // .Define("n_P_pred", Neutron_predicted_momentum, {"nu_Preco", "mu_pxreco", "mu_pyreco", "mu_pzreco"})
                    .Define("n_px_pred", Neutron_predicted_momentumX, {"nu_Preco", "mu_pxreco"})
                    .Define("n_py_pred", Neutron_predicted_momentumY, {"nu_Preco", "mu_pyreco"})
                    .Define("n_pz_pred", Neutron_predicted_momentumZ, {"nu_Preco", "mu_pzreco"})
                    .Define("n_Ppred", "sqrt(n_px_pred*n_px_pred + n_py_pred*n_py_pred + n_pz_pred*n_pz_pred)")
                    .Define("primaries_particles", "EDepSimEvents.Primaries[0].Particles")
                    .Define("true_n_startX",true_n_startX, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_startY",true_n_startY, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_startZ",true_n_startZ, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_startT",true_n_startT, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_stopX",true_n_stopX, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_stopY",true_n_stopY, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_stopZ",true_n_stopZ, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_n_stopT",true_n_stopT, {"EDepSimEvents.SegmentDetectors","primaries_particles"})
                    .Define("true_x_pos", true_x_pos,{"EDepSimEvents.Trajectories"})
                    .Define("true_y_pos", true_y_pos,{"EDepSimEvents.Trajectories"})
                    .Define("true_z_pos", true_z_pos,{"EDepSimEvents.Trajectories"})
                    .Define("true_t", true_t_pos,{"EDepSimEvents.Trajectories"})
                    .Define("hit_pred", Predicted_hit, {"n_px_pred","n_py_pred", "n_pz_pred","vertex_x","vertex_y","vertex_z"})
                    .Define("Beta", "TMath::Sqrt(n_px_pred*n_px_pred + n_py_pred*n_py_pred + n_pz_pred*n_pz_pred) / n_E_pred")
                    //.Define("t_pred", Predicted_time, {"Beta","vertex_t","hit_pred","vertex_x","vertex_y","vertex_z"})
                          .Define("x_ecal", [](const ROOT::VecOps::RVec<Double_t>& x) {return x;}, {"tReco.cluster.x"})
                          .Define("y_ecal", [](const ROOT::VecOps::RVec<Double_t>& y) {return y;}, {"tReco.cluster.y"})
                          .Define("z_ecal", [](const ROOT::VecOps::RVec<Double_t>& z) {return z;}, {"tReco.cluster.z"})
                          .Define("t_ecal", [](const ROOT::VecOps::RVec<Double_t>& t) {return t;}, {"tReco.cluster.t"})
                          .Define("e_ecal", [](const ROOT::VecOps::RVec<Double_t>& e) {return e;}, {"tReco.cluster.e"})
                          .Define("cell_energies", ExtractEnergies, {"tReco.cluster"})
                          .Define("cell_x", Extract_X_cell ,{"tReco.cluster"})
                          .Define("cell_y", Extract_Y_cell ,{"tReco.cluster"})
                          .Define("cell_z", Extract_Z_cell ,{"tReco.cluster"})
                          .Define("cell_ID", CellID, {"tReco.cluster"})
                          .Define("N_tracks", [](const ROOT::VecOps::RVec<track>& t){ return (int)t.size();}, {"tReco.track"})
                          .Define("tracks_ID", Particles_with_track, {"particles"})
                          .Define("N_points", N_points, {"particles"})
                          .Define("chi2", Chi2, {"particles"})
                          .Define("neutron_cluster", Neutron_cluster, {"particles"})
                          .Define("t_pred", Predicted_time,{"n_E_pred", "n_px_pred","n_py_pred","n_pz_pred","vertex_t","vertex_x","vertex_y","vertex_z","x_ecal","y_ecal","z_ecal"});

// auto tFiltered = tRough.Filter("nu_type == -14");
//                         .Filter("Evt_cat != 0");
//Filter("st_proc_type == \"proc:Weak[CC],QES\"")
                    //    .Filter("nu_type == -14");



//  tFiltered.Snapshot("Filtered_tree", "B_Filtered_tree_15_09.root", {"EDepSimEvents.EventId", "vertex_x", "vertex_y","vertex_z","vertex_t", "nu_type","st_proc_type","target","Where_int", "Evt_cat",
//                                                        "Enu","pxnu","pynu","pznu","nu_Ereco","mu_transverse_mom","n_transverse_mom", "primaries_PDG","mu_px","mu_py","mu_pz", "mu_pxreco","mu_pyreco","mu_pzreco","mu_E","mu_Ereco","n_px","n_py","n_pz","n_E",
//                                                         "x_ecal","y_ecal","z_ecal", "t_ecal", "e_ecal", "cell_energies","cell_y","cell_z","cell_ID","N_Tracks"});

 tRough.Snapshot("Rough_tree", "/storage/gpfs_data/neutrino/users/croselli/root_macros/Rough_tree_CONTROLLO.root", {"EDepSimEvents.EventId", "vertex_x", "vertex_y","vertex_z","vertex_t", "nu_type","st_proc_type","target","Where_int", "Evt_cat",
                                                       "Enu","pxnu","pynu","pznu","nu_Ereco", "nu_Preco", 
                                                       "primaries_PDG", "primaries_px","primaries_py", "primaries_pz",
                                                       "mu_px","mu_py","mu_pz", "mu_pxreco","mu_pyreco","mu_pzreco","mu_P","mu_Preco","mu_E","mu_Ereco",
                                                       "n_px","n_py","n_pz","n_P","n_E", "n_E_pred", "n_px_pred", "n_py_pred", "n_pz_pred","n_Ppred", "hit_pred", "t_pred",
                                                       "true_n_startX", "true_n_startY", "true_n_startZ","true_n_startT",
                                                       "true_n_stopX", "true_n_stopY", "true_n_stopZ","true_n_stopT",
                                                       "true_x_pos","true_y_pos","true_z_pos","true_t",
                                                       "x_ecal","y_ecal","z_ecal", "t_ecal", "e_ecal", "cell_energies","cell_y","cell_z","cell_ID",
                                                       "N_tracks","tracks_ID","N_points", "chi2","neutron_cluster"});
                                                        
//  tRough.Snapshot("NoN_EVT_Filtered_tree", "NON_EVT_Filtered_tree_11_09.root", {"EDepSimEvents.EventId", "vertex_x", "vertex_y","vertex_z","vertex_t", "nu_type","st_proc_type","target","Where_int", "Evt_cat",
//                                                        "Enu","pxnu","pynu","pznu","nu_Ereco", "primaries_PDG","mu_px","mu_py","mu_pz", "mu_pxreco","mu_pyreco","mu_pzreco","mu_E","mu_Ereco","n_px","n_py","n_pz","n_E",
//                                                         "x_ecal","y_ecal","z_ecal", "t_ecal", "e_ecal", "cell_energies","cell_y","cell_z","cell_ID","N_Tracks"});

}