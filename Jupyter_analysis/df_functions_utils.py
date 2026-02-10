import ROOT
from ROOT import TFile
from ROOT import TChain
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import RDataFrame
from ROOT import RVec
from ROOT import kRed, kBlue, gPad
import numpy as np
import array

ROOT.gSystem.Load("libROOTDataFrame")

#select the TRUE events with a antimuon and a nuetron in the final state 
ROOT.gInterpreter.Declare("""
bool IsMuonNeutron(const ROOT::VecOps::RVec<int>& pdgs) {
    if (pdgs.size() != 2) return false;
    return ((pdgs[0] == -13 && pdgs[1] == 2112) || (pdgs[0] == 2112 && pdgs[1] == -13));
}
""")

#select the events with a successfull reconstruction
ROOT.gInterpreter.Declare("""
    bool Reco_ok( const Int_t& N_tracks, ROOT::VecOps::RVec<double> mu_pxreco, ROOT::VecOps::RVec<double> mu_pyreco, ROOT::VecOps::RVec<double> mu_pzreco, double mu_Preco ) {
                          if (N_tracks == 0) return false; 
                          if (mu_pxreco.empty()) return false;
                          if (mu_pzreco[0] < 0.0) return false;
                          //if (mu_Preco == 1.000000 || mu_Preco == 0.000000) return false;
                          if (std::abs(mu_Preco - 1.0) < 1e-9 || std::abs(mu_Preco - 0.0) < 1e-9) return false;
                          if (std::abs(mu_Preco - 0.0) < 1e-9 || std::abs(mu_Preco - 0.0) < 1e-9) return false;
                          return true; 
                          }
""")


#Filter on neutron predicetd position
ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
                          
RVec<bool> neutron_correspondence(RVec<double> hit_pred, RVec<double> x_ecal, RVec<double> y_ecal, RVec<double> z_ecal ){

double r = 100.0; //(mm) radius of the cilinder almost half of the ECAL trapezoidal module
//double r = 200.0;
//double r = 180.0;
//double r = 150.0;
// double r = 140.0;                        
//double r = 120.0;
                          

//cilinder axis
double x1 = hit_pred[0], y1 = hit_pred[1], z1 = hit_pred[2];
double x2 = hit_pred[3], y2 = hit_pred[4], z2 = hit_pred[5];

double vx = x2 - x1;
double vy = y2 - y1;
double vz = z2 - z1;
double vlen2 = vx*vx + vy*vy + vz*vz; 

//return vector
size_t n = x_ecal.size();
RVec<bool> inside(n, false);
                          
for (size_t i = 0; i < n; ++i) {
    
    //cluster point
    double qx = x_ecal[i] - x1;
    double qy = y_ecal[i] - y1;
    double qz = z_ecal[i] - z1;
                          
    // Parametro t della proiezione
    double t = (qx*vx + qy*vy + qz*vz) / vlen2;

    if (t < 0.0 || t > 1.0) {
        inside[i] = false;
        continue;
    }
                          
    // Punto proiettato sull'asse
    double px = x1 + t * vx;
    double py = y1 + t * vy;
    double pz = z1 + t * vz;
                          
    // Distanza punto-asse
    double dx = x_ecal[i] - px;
    double dy = y_ecal[i] - py;
    double dz = z_ecal[i] - pz;
    double dist2 = dx*dx + dy*dy + dz*dz;
                          
    inside[i] = (dist2 <= r*r);
}

return inside; 
}                          
""")

#function to find the pair of closer points
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<double> ComparisonPoint( const ROOT::VecOps::RVec<double>& x_ecal, const ROOT::VecOps::RVec<double>& y_ecal, const ROOT::VecOps::RVec<double>& z_ecal, const ROOT::VecOps::RVec<double>& t_ecal,
                                       const ROOT::VecOps::RVec<double>& true_x_pos, const ROOT::VecOps::RVec<double>& true_y_pos, const ROOT::VecOps::RVec<double>& true_z_pos, const ROOT::VecOps::RVec<double>& true_t){
                          
    ROOT::VecOps::RVec<double> selected(8, -1e10);
    double min_dist = 1e9;
    int best_i = -1;
    int best_j = -1;

    //comparison between cluster points and hit_Start points                      
    for (size_t i = 0; i < x_ecal.size(); ++i) {
        for (size_t j = 0; j < true_x_pos.size(); ++j) {
            double dx = x_ecal[i] - true_x_pos[j];
            double dy = y_ecal[i] - true_y_pos[j];
            double dz = z_ecal[i] - true_z_pos[j];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (dist < min_dist) {
                min_dist = dist;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Se ha trovato un match valido, salva i valori corrispondenti
    if (best_i >= 0 && best_j >= 0) {
        selected[0] = x_ecal[best_i];
        selected[1] = y_ecal[best_i];
        selected[2] = z_ecal[best_i];
        selected[3] = t_ecal[best_i];
        selected[4] = true_x_pos[best_j];
        selected[5] = true_y_pos[best_j];
        selected[6] = true_z_pos[best_j];
        selected[7] = true_t[best_j];
    }

    return selected;                          
                          
                          }
                                                    """)




#function to find the nearest point to the cilinder axis 
ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
                          
RVec<double> Nearest_point(RVec<double> hit_pred, RVec<double> true_x_pos, RVec<double> true_y_pos, RVec<double> true_z_pos, RVec<double> true_t){

double r = 10000.0; //(mm) radius of the cilinder
     

//cilinder axis
double x1 = hit_pred[0], y1 = hit_pred[1], z1 = hit_pred[2];
double x2 = hit_pred[3], y2 = hit_pred[4], z2 = hit_pred[5];

double vx = x2 - x1;
double vy = y2 - y1;
double vz = z2 - z1;
double vlen2 = vx*vx + vy*vy + vz*vz; 

//return vector (x,y,z,t,dist)
RVec<double> nearest(5, -1e9);

double min_dist2 = 1e99; // starting min dis 
                          
for (size_t i = 0; i < true_x_pos.size(); ++i) {
    
    //true point
        double qx = true_x_pos[i] - x1;
        double qy = true_y_pos[i] - y1;
        double qz = true_z_pos[i] - z1;
                          
    // projection
    double t = (qx*vx + qy*vy + qz*vz) / vlen2;

    if (t < 0.0 || t > 1.0) continue; //if the point is not inside C, continue
    
                          
    // projection of the point on the axis
    double px = x1 + t * vx;
    double py = y1 + t * vy;
    double pz = z1 + t * vz;
                          
    // distance point-axis
    double dx = true_x_pos[i] - px;
    double dy = true_y_pos[i] - py;
    double dz = true_z_pos[i] - pz;
    double dist2 = dx*dx + dy*dy + dz*dz;
                          
    // Check if the point is inside C
        if (dist2 <= r*r) {
        //check if the point is the nearest         
            if (dist2 < min_dist2) {
                min_dist2 = dist2;
                nearest[0] = true_x_pos[i];
                nearest[1] = true_y_pos[i];
                nearest[2] = true_z_pos[i];
                nearest[3] = true_t[i];
                nearest[4] = std::sqrt(dist2);
            }
        }
    }

return nearest; 
}                          
""")

ROOT.gInterpreter.Declare("""
double Predicted_time(const Double_t n_E_pred, const Double_t n_px_pred, const Double_t n_py_pred, const Double_t n_pz_pred, 
                            const double t0, const double x0, const double y0, const double z0,
                            const ROOT::VecOps::RVec<double> nearest_point, 
                            const ROOT::VecOps::RVec<double> true_t){
                                
                          double time;
                                //ROOT::VecOps::RVec<double> time(x_ecal.size());                                
                                double beta = TMath::Sqrt(n_px_pred*n_px_pred + n_py_pred*n_py_pred + n_pz_pred*n_pz_pred) / n_E_pred;
                                const double c = 299.7924580;//mm/ns

                                
                                double s = std::sqrt((nearest_point[0] - x0)*(nearest_point[0] - x0) + (nearest_point[1] - y0)*(nearest_point[1] - y0) + (nearest_point[2] - z0)*(nearest_point[2] - z0));
                                time = true_t[0] + s/(beta * c );
                                
                                return time;
                                
                          }
""")




ROOT.gInterpreter.Declare("""
int Cat(TString target, TString st_proc_type, TString Where_int) {

    if (target == "tgt:1000010010" && st_proc_type == "proc:Weak[CC],QES" && Where_int == "C3H6") // anti mu_nu on H in plastic
        return 1;

    else if (target == "tgt:1000060120" && Where_int == "C3H6" ) // anti mu_nu on C in plastic
        return 2;

    else if (target == "tgt:1000010010" && st_proc_type != "proc:Weak[CC],QES" && Where_int == "C3H6" ) // anti mu_nu on H in plastic non QES
        return 3;

    else if (target == "tgt:1000060120" && Where_int == "Graphite" ) // anti mu_nu on C in graphite
        return 4;
     
     else if (Where_int != "Graphite" && Where_int != "C3H6") // on H everywhere else
        return 5;

     else return 0;
} 
                          """)

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
                          
RVec<bool> neutron_correspondenceTRUE(RVec<double> hit_pred, RVec<double> x_ecal, RVec<double> y_ecal, RVec<double> z_ecal ){

double r = 100.0; //(mm) 
//double r = 200.0;
//double r = 180.0;
//double r = 150.0;
// double r = 140.0;                        
//double r = 120.0;
                          

//cilinder axis
double x1 = hit_pred[0], y1 = hit_pred[1], z1 = hit_pred[2];
double x2 = hit_pred[3], y2 = hit_pred[4], z2 = hit_pred[5];

double vx = x2 - x1;
double vy = y2 - y1;
double vz = z2 - z1;
double vlen2 = vx*vx + vy*vy + vz*vz; 

//return vector
size_t n = x_ecal.size();
RVec<bool> inside(n, false);
                          
for (size_t i = 0; i < n; ++i) {
    
    //cluster point
    double qx = x_ecal[i] - x1;
    double qy = y_ecal[i] - y1;
    double qz = z_ecal[i] - z1;
                          
    // Parametro t della proiezione
    double t = (qx*vx + qy*vy + qz*vz) / vlen2;

    if (t < 0.0 || t > 1.0) {
        inside[i] = false;
        continue;
    }
                          
    // Punto proiettato sull'asse
    double px = x1 + t * vx;
    double py = y1 + t * vy;
    double pz = z1 + t * vz;
                          
    // Distanza punto-asse
    double dx = x_ecal[i] - px;
    double dy = y_ecal[i] - py;
    double dz = z_ecal[i] - pz;
    double dist2 = dx*dx + dy*dy + dz*dz;
                          
    inside[i] = (dist2 <= r*r);
}

return inside; 
}                          
""")

ROOT.gInterpreter.Declare("""
Double_t transverse_momentum( ROOT::VecOps::RVec<double> primaries_px, ROOT::VecOps::RVec<double> primaries_py, ROOT::VecOps::RVec<double> primaries_pz,
                              Double_t px_nu , Double_t py_nu, Double_t pz_nu) {

    TVector3 v1( px_nu, py_nu,pz_nu);
    TVector3 v2(Sum(primaries_px), Sum(primaries_py), Sum(primaries_pz));

    Double_t transvers_mom = v2.Perp(v1);

    return transvers_mom;
}
                          
""")


def df_extend(df):
    df_extended = (df.Define("nearest_segpoint", "Nearest_point(hit_pred, true_n_startX, true_n_startY, true_n_startZ, true_n_startT)")
                    .Define("nearest_ecalpoint", "Nearest_point(hit_pred, x_ecal, y_ecal, z_ecal, t_ecal)")
                    .Define("distance_segpoint", "nearest_segpoint[4]")
                    .Define("distance_ecalpoint", "nearest_ecalpoint[4]")
                    .Define("t_pred_seg", "Predicted_time(n_E_pred,n_px_pred, n_py_pred, n_pz_pred, vertex_t,vertex_x, vertex_y,vertex_z, nearest_segpoint, true_t)")
                    .Define("t_pred_ecal", "Predicted_time(n_E_pred,n_px_pred, n_py_pred, n_pz_pred, vertex_t,vertex_x, vertex_y,vertex_z, nearest_ecalpoint, true_t)")
                    #.Define("inside2", "neutron_correspondence(hit_pred, x_ecal, y_ecal, z_ecal)")
                    .Define("inside", "neutron_correspondence(hit_pred, x_ecal, y_ecal, z_ecal)")
                    .Define("n_inside", "Sum(inside)")
                    .Define("insideTRUE", "neutron_correspondenceTRUE(hit_pred, true_n_startX,true_n_startY, true_n_startZ)")
                    .Define("n_insideTRUE", "Sum(insideTRUE)")
                    .Define("t_res_seg", "t_pred_seg - nearest_segpoint[3]")
                    .Define("t_res_ecal", "t_pred_ecal - nearest_ecalpoint[3]")
                    .Define("Cat", "Cat(target, st_proc_type, Where_int)")
                    )
    return df_extended

def df_addflags(df):
    true_filter = 'IsMuonNeutron(primaries_PDG) && target == "tgt:1000010010" && st_proc_type == "proc:Weak[CC],QES"'
    volume_filter = "vertex_z > 22850"
    reco_filter = "Reco_ok(N_tracks, mu_pxreco, mu_pyreco, mu_pzreco, mu_Preco )"
    track_filter = "N_tracks == 1"
    space_filter = "n_inside > 0"
    spaceTRUE_filter = "n_insideTRUE > 0"
    time_filter = "t_res_ecal > -20 && t_res_ecal < 20"
    timeTRUE_filter = "t_res_seg > -20 && t_res_seg < 20"
    df_flag = (df.Define("true_flag", f"({true_filter}) ? 1 : 0")
                      .Define("volume_flag", f"({volume_filter}) ? 1 : 0")
                      .Define("reco_flag", f"({reco_filter}) ? 1 : 0")
                      .Define("track_flag", f"({track_filter}) ? 1 : 0")
                      .Define("space_flag", f"({space_filter}) ? 1 : 0")
                      .Define("spaceTRUE_flag", f"({spaceTRUE_filter}) ? 1 : 0")
                      .Define("time_flag", f"({time_filter}) ? 1 : 0")
                      .Define("timeTRUE_flag", f"({timeTRUE_filter}) ? 1 : 0")
    ) 
    return df_flag

def df_addtransversemomentum(df):
    df_tr = (df.Define("transverse_mom", "transverse_momentum(primaries_px, primaries_py, primaries_pz, pxnu, pynu, pznu)")
      .Define("transverse_mom_GeV", "transverse_mom/1000.0")
      )
    
    return df_tr