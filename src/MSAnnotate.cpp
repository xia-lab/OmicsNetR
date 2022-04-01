#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This script is developed to accerate NetID annotation process
// Author: Zhiqiang Pang [zhiqiang.pang@mail.mcgill.ca]
// ref: NetID paper (https://www.nature.com/articles/s41592-021-01303-3)

int find_int_pos (int num, IntegerVector numVec) {
  for (int i = 0; i < numVec.length(); i++) {
    if(numVec[i] == num) {
      return(i);
    }
  }
  return(-1);
}

IntegerVector find_int_poss (int num, IntegerVector numVec) {
  IntegerVector res;
  for (int i = 0; i < numVec.length(); i++) {
    if(numVec[i] == num) {
      res.push_back(i);
    }
  }
  return(res);
}


void coreFUN (int i, 
              StringVector &res,
              IntegerVector &solRes, 
              IntegerVector &ilpNodeIdVec,
              CharacterVector &colNMs,
              CharacterVector &rowNMs,
              IntegerVector &core_ilpNIdVec,
              StringVector &core_annoVec,
              NumericMatrix &dis_mat,
              List &g_annotation,
              Function &f,
              IntegerVector &IP_edge_From,
              IntegerVector &IP_edge_To,
              StringVector &classVec,
              StringVector &IP_edge_cat,
              StringVector &IP_edge_linktp,
              StringVector &IP_edge_fm1,
              StringVector &IP_edge_fm2,
              IntegerVector &IP_edge_dir,
              string igraphmode) {
  
  int IlpNodeID = ilpNodeIdVec[i]; // this equals query_ilp_id

  if(std::find(colNMs.begin(), colNMs.end(), IlpNodeID) == colNMs.end()) {
    res[i] = "not exited in dm";
    return;
  }
  //std::to_string()
  if(std::find(rowNMs.begin(), rowNMs.end(), IlpNodeID) != rowNMs.end()){
    int core_pos = find_int_pos(IlpNodeID, core_ilpNIdVec);
    res[i] = core_annoVec[core_pos];
    return;
  }
  
  int tmp_pos = std::find(colNMs.begin(), colNMs.end(), IlpNodeID).index();
  NumericVector v = dis_mat(_, tmp_pos);
  double dm_min = min(v);
  double inf = std::numeric_limits<double>::infinity();
  
  if(dm_min == inf) {
    res[i] = "No edge existed!";
    return;
  }
  int np = which_min(v);
  String parent_selected = rowNMs[np];
  
  List path_nodes = f(g_annotation, 
                      Named("from", parent_selected),
                      Named("to", std::to_string(IlpNodeID)),
                      Named("mode", igraphmode),
                      Named("output", "vpath"));
  List tmp1 = as<List>(path_nodes[0]);
  IntegerVector tmp2 = as<IntegerVector>(tmp1[0]);
  StringVector tmp3 = tmp2.names();
  IntegerVector ilp_node_path(tmp3.size());
  for(int j = 0; j < tmp3.size(); j++){
    String tmp4 = tmp3[j];
    ilp_node_path[j] = stoi(tmp4.get_cstring());
  }
  
  int corePos = find_int_pos(ilp_node_path[0], core_ilpNIdVec);
  
  String temp_core = core_annoVec[corePos];
  String transform_path;
  int dir;
  
  for(int k =0; k < ilp_node_path.size()-1; k++) {
    IntegerVector resFrom = find_int_poss(ilp_node_path[k], IP_edge_From);
    IntegerVector resTo = find_int_poss(ilp_node_path[k+1], IP_edge_To);
    IntegerVector resposs = intersect(resFrom, resTo);
    resposs = resposs.sort();
    
    if(resposs.size() > 0) {
      dir = IP_edge_dir[resposs[0]];
    } else if(resposs.size() == 0) {
      resFrom = find_int_poss(ilp_node_path[k], IP_edge_To);
      resTo = find_int_poss(ilp_node_path[k+1], IP_edge_From);
      resposs = intersect(resFrom, resTo);
      dir = -1;
    }
    
    string temp_sign;
    if(classVec[i] == "Artifact") {
      if(IP_edge_cat[resposs[0]] == "Oligomer") {
        temp_sign = "*";
      } else if (IP_edge_cat[resposs[0]] == "Multicharge") {
        temp_sign = "/";
      } else if (IP_edge_cat[resposs[0]] == "Heterodimer") {
        temp_sign = "+ Peak";
      } else if (dir == 1) {
        temp_sign = "+";
      } else if (dir == -1){
        temp_sign = "-";
      }
    } else {
      temp_sign = (dir == 1) ? "+" : "-";
    }
    
    String temp_linktype = IP_edge_linktp[resposs[0]];
    String temp_formula = (dir == 1) ? IP_edge_fm2[resposs[0]] : IP_edge_fm1[resposs[0]];
    
    string seperator = " ";
    transform_path =
      transform_path.get_cstring() + seperator +
      temp_sign + seperator +
      temp_linktype.get_cstring() + seperator + "-> " + 
      temp_formula.get_cstring();
  }
  string reso = transform_path.get_cstring();
  string resx= temp_core.get_cstring() + reso;
  res[i] = resx;
  
}

// [[Rcpp::export]]
StringVector path_annotate(DataFrame ilp_nodes, 
                                DataFrame canu_met, //core_annotation_unique_*
                                DataFrame ilp_edges_anno_met, //ilp_edges_annotate_*
                                NumericMatrix dis_mat_met,
                                List g_annotation,
                                DataFrame canu_nonmet, //core_annotation_unique_*
                                DataFrame ilp_edges_anno_nonmet, //ilp_edges_annotate_*
                                NumericMatrix dis_mat_nonmet,
                                List g_anno_non) {
  
  IntegerVector solRes = as<IntegerVector>(ilp_nodes[27]);
  IntegerVector ilpNodeIdVec = as<IntegerVector>(ilp_nodes[0]);
  StringVector classVec = as<StringVector>(ilp_nodes[12]);
  StringVector res(solRes.length());
  
  // for met
  IntegerVector core_ilpNIdVec1 = as<IntegerVector>(canu_met[0]);
  StringVector core_annoVec1 = as<StringVector>(canu_met[2]);
  
  CharacterVector colNMs1 = colnames(dis_mat_met);
  CharacterVector rowNMs1 = rownames(dis_mat_met);
  
  IntegerVector IP_edge_From1 = as<IntegerVector>(ilp_edges_anno_met[0]);
  IntegerVector IP_edge_To1 = as<IntegerVector>(ilp_edges_anno_met[1]);
  IntegerVector IP_edge_dir1 = as<IntegerVector>(ilp_edges_anno_met[8]);
  StringVector IP_edge_linktp1 = as<StringVector>(ilp_edges_anno_met[7]);
  StringVector IP_edge_fm1_1 = as<StringVector>(ilp_edges_anno_met[10]);
  StringVector IP_edge_fm2_1 = as<StringVector>(ilp_edges_anno_met[11]);
  StringVector IP_edge_cat1 = as<StringVector>(ilp_edges_anno_met[6]);
  
  // for nonmet
  IntegerVector core_ilpNIdVec2 = as<IntegerVector>(canu_nonmet[0]);
  StringVector core_annoVec2 = as<StringVector>(canu_nonmet[2]);
  
  CharacterVector colNMs2 = colnames(dis_mat_nonmet);
  CharacterVector rowNMs2 = rownames(dis_mat_nonmet);
  
  IntegerVector IP_edge_From2 = as<IntegerVector>(ilp_edges_anno_nonmet[0]);
  IntegerVector IP_edge_To2 = as<IntegerVector>(ilp_edges_anno_nonmet[1]);
  IntegerVector IP_edge_dir2 = as<IntegerVector>(ilp_edges_anno_nonmet[8]);
  StringVector IP_edge_linktp2 = as<StringVector>(ilp_edges_anno_nonmet[7]);
  StringVector IP_edge_fm1_2 = as<StringVector>(ilp_edges_anno_nonmet[10]);
  StringVector IP_edge_fm2_2 = as<StringVector>(ilp_edges_anno_nonmet[11]);
  StringVector IP_edge_cat2 = as<StringVector>(ilp_edges_anno_nonmet[6]);
  
  // igraph FUN
  Environment pkg = Environment::namespace_env("igraph");
  Function f = pkg["shortest_paths"];
  
  for(int i = 0; i < solRes.size(); i++) {
    if(solRes[i] == 0) {continue;}
    if(classVec[i] == "Metabolite" || classVec[i] == "Putative metabolite") {
      coreFUN(i, 
              res,
              solRes, 
              ilpNodeIdVec,
              colNMs1,
              rowNMs1,
              core_ilpNIdVec1,
              core_annoVec1,
              dis_mat_met,
              g_annotation,
              f,
              IP_edge_From1,
              IP_edge_To1,
              classVec,
              IP_edge_cat1,
              IP_edge_linktp1,
              IP_edge_fm1_1,
              IP_edge_fm2_1,
              IP_edge_dir1,
              "all");
    } else if (classVec[i] == "Artifact") {
      coreFUN(i, 
              res,
              solRes, 
              ilpNodeIdVec,
              colNMs2,
              rowNMs2,
              core_ilpNIdVec2,
              core_annoVec2,
              dis_mat_nonmet,
              g_anno_non,
              f,
              IP_edge_From2,
              IP_edge_To2,
              classVec,
              IP_edge_cat2,
              IP_edge_linktp2,
              IP_edge_fm1_2,
              IP_edge_fm2_2,
              IP_edge_dir2,
              "out");
    } else {
      res[i] = "Unknown";
    }
  }
  
  return res;
}


// [[Rcpp::export]]
StringVector matchHMDB(DataFrame annotation, 
                       DataFrame Restable,
                       int peakidx,
                       int formularidx) {
  IntegerVector PeakIDs = as<IntegerVector>(Restable[peakidx]);
  StringVector formulasIDs = as<StringVector>(Restable[formularidx]);
  StringVector res(PeakIDs.size());
  IntegerVector NodeIDs = as<IntegerVector>(annotation[3]);
  StringVector formulas = as<StringVector>(annotation[4]);
  StringVector HMDBIDs = as<StringVector>(annotation[23]);
  for(int i = 0; i < PeakIDs.size(); i++) {
    for(int j = 0; j < NodeIDs.size(); j++) {
      if(PeakIDs[i] == NodeIDs[j] && formulasIDs[i] == formulas[j]) {
        res[i] = HMDBIDs[j];
        break;
      }
    }
  }
  return res;
}


// [[Rcpp::export]]
StringVector convert2KEGG(StringVector HMDBIDs, DataFrame database) {
  StringVector res(HMDBIDs.size());
  StringVector HMDBDBids = as<StringVector>(database[0]);
  StringVector KEGGDBids = as<StringVector>(database[9]);
  for(int i=0; i< HMDBIDs.size(); i++){
    if(HMDBIDs[i] != "") {
      for(int j=0; j < HMDBDBids.size(); j++) {
        if(HMDBIDs[i] == HMDBDBids[j]) {
          res[i] = KEGGDBids[j];
          break;
        }
      }
    }
  }
  
  return res;
}
