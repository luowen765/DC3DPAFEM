
#include "parahandler.h"
#include "em.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

ParaHandler::ParaHandler(char *model_file) {
  this->read_model_info(model_file);
}

ParaHandler::~ParaHandler() {}

void ParaHandler::skip_comments(std::istream &in,
                                std::vector<std::string> &para_str_vec,
                                std::string comment_str) {
  std::string line;
  while (std::getline(in, line)) {
    for (char &c : line) // loop string in C++11 grammar
    {
      if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n') {
        c = ' ';
      }
    }

    line.erase(0, line.find_first_not_of(
                      " ")); // delete the space at the beginning of the line
    line.erase(line.find_last_not_of(" ") +
               1); // delete the space at the ending of the line

    int n_comment_start = line.find_first_of(comment_str);
    if (n_comment_start != std::string::npos) {
      line.erase(n_comment_start); // delete comments
    }

    if (line.empty())
      continue;

    para_str_vec.push_back(line);
  }
}

void ParaHandler::read_model_info(char *model_file) {
  std::ifstream in_stream(model_file);
  myassert(in_stream.good());
  int para_lines = 18; // number of parameters' lines
  std::vector<std::string> para_str_vec;
  skip_comments(in_stream, para_str_vec);
  // std::cout << "para_str_vec.size(): " << para_str_vec.size() << std::endl;
  // std::cout << " para_lines: " << para_lines << std::endl;
  myassert(para_str_vec.size() == para_lines);
  in_stream.close();
  std::stringstream ss;

  // survey mode file
  ss << para_str_vec[0];
  ss >> survey_input;
  ss.clear();

  // Model parameters
  ss << para_str_vec[1];
  ss >> model_parameters_file;
  ss.clear();

  ss << para_str_vec[2];
  ss >> N_uniformRefine;
  ss.clear();
  myassert(N_uniformRefine >= 0);

  ss << para_str_vec[3];
  ss >> maxit;
  ss.clear();
  myassert(maxit >= 1);

  ss << para_str_vec[4];
  ss >> beta;
  ss.clear();
  myassert(beta >= 0.0 && beta <= 1);

  ss << para_str_vec[5];
  ss >> max_dofs;
  ss.clear();
  myassert(max_dofs >= 10);

  ss << para_str_vec[6];
  ss >> linear_solver;
  ss.clear();
  myassert(linear_solver == "amg");

  // Mesh parameters
  ss << para_str_vec[7];
  ss >> save_amr_mesh;
  ss.clear();

  ss << para_str_vec[8];
  ss >> mesh_file;
  ss.clear();

  ss << para_str_vec[9];
  ss >> find_points_by;
  ss.clear();
  myassert(find_points_by == "FindPoints");

  ss << para_str_vec[10];
  ss >> pcg_maxit;
  ss.clear();
  myassert(pcg_maxit > 5);

  ss << para_str_vec[11];
  ss >> pcg_primal_tol;
  ss.clear();
  myassert(pcg_primal_tol > 0 && pcg_primal_tol < 1);

  ss << para_str_vec[12];
  ss >> pcg_dual_tol;
  ss.clear();
  myassert(pcg_dual_tol > 0 && pcg_dual_tol < 1 &&
           pcg_dual_tol >= pcg_primal_tol);

  ss << para_str_vec[13];
  ss >> pcg_print_level;
  ss.clear();
  myassert(pcg_print_level >= 0 && pcg_print_level < 10);

  ss << para_str_vec[14];
  ss >> amg_print_level;
  ss.clear();
  myassert(amg_print_level >= 0 && amg_print_level < 10);

  ss << para_str_vec[15];
  ss >> if_solve_dual;
  ss.clear();
  myassert(if_solve_dual == "true" || if_solve_dual == "false");

  ss << para_str_vec[16];
  ss >> error_esti_way;
  ss.clear();
  myassert(error_esti_way == "jn");

  ss << para_str_vec[17];
  ss >> if_compute_E_J;
  ss.clear();
  myassert(if_compute_E_J == "true" || if_compute_E_J == "false");

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  sorted_survey = "sorted_" + survey_input;
  unsigned int n_data = 0;
  int mode = 0;
  if (myid == 0) {
    // check survey_input file
    std::ifstream survey_in(survey_input);
    myassert(survey_in.good());
    std::ofstream survey_out(sorted_survey);
    myassert(survey_out.good());
    survey_in >> n_data // nunmber of measurments
        >> mode;        // measurment mode
    myassert(n_data > 0);
    myassert(mode == 11 || mode == 12 || mode == 21 || mode == 22 ||
             mode == 41);
    if (mode == 11) { // pole-pole configuration
      std::vector<std::vector<double>> lines;
      for (int i = 0;; i++) {
        std::vector<double> templine;
        templine.resize(6); // two points
        for (int j = 0; j < 6; j++) {
          survey_in >> templine[j];
        }
        if (survey_in.eof()) {
          break;
        }
        lines.push_back(templine);
      }
      myassert(lines.size() == n_data);
      std::sort(lines.begin(), lines.end(), cmp2);
      // generate sorted survey file
      survey_out << n_data << "\t" << mode << "\n";
      for (int i = 0; i < n_data; i++) {
        for (int j = 0; j < 6; j++) {
          survey_out << lines[i][j] << "\t";
        }
        survey_out << "\n";
      }
    } else {
      std::cout << "Unsopported mode temporarily!" << std::endl;
      std::abort();
    }
    survey_in.close();
    survey_out.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  // get sources and potentions by sorted survey file
  std::ifstream survey_in_stream(sorted_survey);
  myassert(survey_in_stream.good());
  survey_in_stream >> n_data // nunmber of measurments
      >> mode;               // measurment mode
  Vertex last_ts(-1, -1, -1);
  Vertex ts(-1, -1, -1);
  Vertex tm(-1, -1, -1);
  Vertex ts2(-1, -1, -1);
  Vertex ts3(-1, -1, -1);
  Vertex ts4(-1, -1, -1);
  std::vector<Vertex> temp_sites;
  for (int i = 0; i < n_data; i++) {
    if (mode == 11) {
      survey_in_stream >> ts(0) >> ts(1) >> ts(2);
      survey_in_stream >> tm(0) >> tm(1) >> tm(2);
      if (i == 0) {
        last_ts = ts;
        sources.push_back(ts);
        temp_sites.push_back(tm);
      } else if (i == n_data - 1) {
        temp_sites.push_back(tm);
        sites.push_back(temp_sites);
        temp_sites.clear();
      } else {
        if (ts(0) == last_ts(0) && ts(1) == last_ts(1) && ts(2) == last_ts(2)) {
          temp_sites.push_back(tm);
        } else {
          last_ts = ts;
          sources.push_back(ts);
          sites.push_back(temp_sites);
          temp_sites.clear();
          temp_sites.push_back(tm); 
        }
      }
    } else {
      std::cout << "Unsopported mode temporarily!\n";
      std::abort();
    }
  }
  survey_in_stream.close();
  source_number = sources.size();
  myassert(source_number > 0);

  // Read conductivity model
  std::ifstream cond_in_stream(model_parameters_file.c_str());
  myassert(cond_in_stream.good());
  cond_in_stream >> n_regions;
  for (int i = 0; i < n_regions; i++) {
    int marker;
    Vector main_cond(3);
    main_cond = 0.0;
    Vector angle(3);
    angle = 0.0;
    // double permeability;
    cond_in_stream >> marker >> main_cond[0] >> main_cond[1] >> main_cond[2] >>
        angle[0] >> angle[1] >> angle[2]; // >> permeability;
    marker_vec.push_back(marker);
    region2conductivity[marker] = cal_conductivity(main_cond, angle);

  }
  cond_in_stream.close();
  // check mesh file stream
  std::ifstream msh_stream(mesh_file);
  myassert(msh_stream.good());
  msh_stream.close();
}

DenseMatrix ParaHandler::cal_conductivity(Vector &main_cond, Vector &_angle) {
  Vector angle(3);
  angle = 0.0;
  angle[0] = _angle[0] * (EM::PI / 180);
  angle[1] = _angle[1] * (EM::PI / 180);
  angle[2] = _angle[2] * (EM::PI / 180);
  DenseMatrix sigma(3);
  sigma = 0.0;
  sigma(0, 0) = main_cond[0];
  sigma(1, 1) = main_cond[1];
  sigma(2, 2) = main_cond[2];

  DenseMatrix R_1(3), R_2(3), R_3(3), R(3);
  R_1 = 0.0;
  R_2 = 0.0;
  R_3 = 0.0;
  R_1(0, 0) = cos(angle[0]);
  R_1(0, 1) = -sin(angle[0]);
  R_1(0, 2) = 0;
  R_1(1, 0) = sin(angle[0]);
  R_1(1, 1) = cos(angle[0]);
  R_1(1, 2) = 0;
  R_1(2, 0) = 0;
  R_1(2, 1) = 0;
  R_1(2, 2) = 1;
  R_2(0, 0) = 1;
  R_2(0, 1) = 0;
  R_2(0, 2) = 0;
  R_2(1, 0) = 0;
  R_2(1, 1) = cos(angle[1]);
  R_2(1, 2) = -sin(angle[1]);
  R_2(2, 0) = 0;
  R_2(2, 1) = sin(angle[1]);
  R_2(2, 2) = cos(angle[1]);
  R_3(0, 0) = cos(angle[2]);
  R_3(0, 1) = -sin(angle[2]);
  R_3(0, 2) = 0;
  R_3(1, 0) = sin(angle[2]);
  R_3(1, 1) = cos(angle[2]);
  R_3(1, 2) = 0;
  R_3(2, 0) = 0;
  R_3(2, 1) = 0;
  R_3(2, 2) = 1;

  DenseMatrix temp1(3);
  temp1 = 0.0;
  Mult_DenseMatrix3(R_1, R_2, temp1);
  Mult_DenseMatrix3(temp1, R_3, R);

  DenseMatrix R_sigma_RT(3), temp2(3),temp3(3);
  R_sigma_RT = 0.0;
  temp2 = 0.0;
  temp3 = R;
  temp3.Transpose();
  Mult_DenseMatrix3(R, sigma, temp2);
  Mult_DenseMatrix3(temp2, temp3, R_sigma_RT);
  return R_sigma_RT;
}

DenseMatrix ParaHandler::get_elem_conductivity(int marker) {
  myassert(!region2conductivity.empty());
  std::map<int, DenseMatrix>::iterator it = region2conductivity.find(marker);
  if (it == region2conductivity.end()) {
    std::cout << "not found marker: " << marker << std::endl;
  }
  myassert(it != region2conductivity.end());
  return (*it).second;
}


Vector ParaHandler::get_sources_center() {
  Vector gpoint(3);
  gpoint = 0.0;
  for (int i = 0; i < sources.size(); i++) {
    for (int j = 0; j < 3; j++) {
      gpoint[j] += sources[i](j);
    }
  }
  gpoint /= sources.size();
  return gpoint;
}

void ParaHandler::preProcess() {
  for (int i = 0; i < sites.size(); i++) {
    for (int j = 0; j < sites[i].size(); j++) {
      s_plus_m.push_back(sites[i][j]);
    }
  }
  remove_duplicate_vertex(s_plus_m);
}



void ParaHandler::find_point_tets(ParMesh *pmesh, std::vector<Vertex> &points,
                                  Array<int> &find_tets,
                                  std::string find_points_by) {
  // Find the potential tets which contain sites
  int npts = points.size();
  DenseMatrix point_mat(3, npts);
  for (int j = 0; j < npts; j++) {
    point_mat(0, j) = points[j](0);
    point_mat(1, j) = points[j](1);
    point_mat(2, j) = points[j](2);
  }
  Array<int> temp_tet_id;
  Array<IntegrationPoint> ips;
  pmesh->FindPoints(point_mat, temp_tet_id, ips);
  for (int j = 0; j < temp_tet_id.Size(); j++) {
    int tid = temp_tet_id[j];
    if (tid == -1) {
      std::cout << "measurements point not found! please check your parameters";
    } else {
      find_tets.Append(tid); 
    }
  }
  myassert(find_tets.Size() == npts);
}

void ParaHandler::printForwardParameters() {
  std::cout << survey_input << std::endl;
  std::cout << model_parameters_file<< std::endl;
  std::cout << N_uniformRefine << "\t--N_uniformRefine" << std::endl;
  std::cout << maxit << "\t--maxit" << std::endl;
  std::cout << beta << "\t--beta" << std::endl;
  std::cout << max_dofs << "\t--max_dofs" << std::endl;
  std::cout << linear_solver << "\t--linear_solver" << std::endl;
  std::cout << save_amr_mesh << "\t--save_amr_mesh" << std::endl;
  std::cout << mesh_file<< std::endl;
  std::cout << find_points_by << "\t--find_points_by" << std::endl;
  std::cout << pcg_maxit << "\t--pcg_maxit" << std::endl;
  std::cout << pcg_primal_tol << "\t--pcg_primal_tol" << std::endl;
  std::cout << pcg_dual_tol << "\t--pcg_dual_tol" << std::endl;
  std::cout << pcg_print_level << "\t--pcg_print_level_tol" << std::endl;
  std::cout << amg_print_level << "\t--amg_print_level" << std::endl;
  std::cout << if_solve_dual << "\t--if_solve_dual" << std::endl;
  std::cout << error_esti_way << "\t--error_esti_way" << std::endl;
  std::cout << if_compute_E_J << "\t--if_compute_U_E_J" << std::endl;
}
