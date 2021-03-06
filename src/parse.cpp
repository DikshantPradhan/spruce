#include "utils.h"
#include "solutionset.h"
#include <fstream>
#include <iostream>

std::vector<std::string> split(const std::string& s, char delimiter){
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)){
    tokens.push_back(token);
  }
  return tokens;
}

std::string writeOutputFile(const std::string inputFile, std::string add){
  std::vector<std::string> input_vec = split(inputFile, '/');

  int size = input_vec.size();
  std::string outputFile = split(input_vec[size-1], '.')[0];
  outputFile = outputFile + add;
  // std::cout << outputFile;

  return outputFile;
}

int getReadNum(const std::string inputFile){
  std::vector<std::string> input_vec = split(inputFile, '/');

  int size = input_vec.size();
  std::string outputFile = split(input_vec[size-1], '.')[0];
  outputFile = split(outputFile, 'c')[1];
  // std::cout << outputFile;

  return stoi(outputFile);
}

int writeParamFile(int reads, std::string outputFile){

  std::filebuf fb;
  fb.open (outputFile,std::ios::out);
  std::ostream os(&fb);

  os << "reads" << "\n";
  os << reads << "\n";

  fb.close();

  return 0;
}

std::vector<std::string> parseLabel(std::string label){

  // char delim = ')';
  std::string string = split(label, ')')[0];
  string = split(string, '(')[1];
  std::vector<std::string> vector = split(string, ',');

  return vector;
}

double** returnStates(gm::Solution sol, int vectorIdx){
  //
  int numVertices = sol.S()[vectorIdx].numVertices();
  double** states = 0;
  states = new double*[3];

  states[0] = new double[numVertices];
  states[1] = new double[numVertices];
  states[2] = new double[numVertices];

  for (int i = 0; i < numVertices; i++){
    std::string label = sol.S()[vectorIdx].label(i);
    std::vector<std::string> labels = parseLabel(label);
    // std::string a = labels[0];
    // double c = atof(a.c_str());
    states[0][i] = atof(labels[0].c_str());
    states[1][i] = atof(labels[1].c_str());
    states[2][i] = atof(labels[2].c_str());

    std::cout << states[0][i] << " " << states[1][i] << " " << states[2][i] << " | ";
  }

  return states;//states; // change this to return array
}

// use this for observed and/or inferred
int returnVAF(gm::Solution sol, std::string outputFile){

  gm::RealTensor F = sol.inferredF();

  int numVertices = sol.S()[0].numVertices();
  // GET VAF

  int k = F.k(); //
  int m = F.m(); //
  int n = F.n(); //

  // std::cout << "count comparison " << numVertices << " " << k << " " << m << " " << n << '\n';
  if (numVertices != k){std::cout << "vaf size mismatch";}

  // std::cout << k;
  // get states (need to put extract into another function)

  // multiply states by frequency

  double calculatedVAF [m][n];
  // double calculatedR [m][n];

  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      double total_r = 0;
      double total_v = 0;

      for (int h = 0; h < k; h++){
        double freq = F(h, i, j);
        //double** states = returnStates(sol, n);
        std::string label = sol.S()[j].label(h);
        std::vector<std::string> labels = parseLabel(label);

        double x = atof(labels[0].c_str());
        double y = atof(labels[1].c_str());
        double z = atof(labels[2].c_str());

        total_r = total_r + ((x + y)*freq);
        total_v = total_v + (z*freq);
      }

      calculatedVAF[i][j] = total_v/total_r;
    }
  }


  // std::cout << "calculated: " << calculatedVAF << '\n';

  std::filebuf fb;
  // outputFile = "data/" + outputFile;

  fb.open (outputFile,std::ios::out);
  //const char *path="/home/user/file.txt";
  std::ostream os(&fb);

  // sols.solution(1).S()[1].writeEdgeList(os);

  // os << "Test sentence\n";
  // OUTPUT COLUMN NAMES
  for (int i = 0; i < n; i++){
    std::string colName = "SNV_" + std::to_string(i+1);
    os << colName << ",";
  }
  os << "\n";
  // OUTPUT NUMBERS
  for (int i = 0; i < m; i++){
    for (int j = n-1; j > -1; j--){

      os << calculatedVAF[i][j] << ",";
    }
    os << "\n";
  }


  fb.close();

  // return tensor of VAFs


  return 0;
}

int writeSolutionEdgeList(gm::Solution sol, std::string outputFile){


  std::filebuf fb;
  fb.open (outputFile,std::ios::out);
  std::ostream os(&fb);

  int vecSize = sol.S().size();
  for (int i = 0; i < vecSize; i++){
    os << "SNV: " << i << "\n";
    sol.S()[i].writeEdgeList(os);
    os << "\n";
  }

  fb.close();

  return 0;
  // sols.solution(1).S()[1].writeEdgeList(os);
}

int comparePhyloMatrix(gm::PerfectPhyloMatrix mtx){

  // std::cout << mtx << '\n';

  int n = mtx.n();
  int k = mtx.k();

  int** A = 0;
  A = new int*[n];

  int** trueA = 0;
  trueA = new int*[n];

  for (int i = 0; i < n; i++){
    A[i] = new int[k];
    trueA[i] = new int[k];
    for (int j = 0; j < k; j++){
      A[i][j] = 0;
      trueA[i][j] = 0;
    }
  }

  // set up true matrix
  trueA[0][1] = 1;
  trueA[0][4] = 1;
  trueA[0][5] = 1;
  trueA[0][7] = 1;
  trueA[1][2] = 1;
  trueA[2][1] = 1;
  trueA[2][3] = 1;
  trueA[2][8] = 1;
  trueA[3][4] = 1;
  trueA[3][5] = 1;
  trueA[4][5] = 1;
  trueA[1][6] = 2;
  trueA[2][4] = 2;
  trueA[2][5] = 2;
  trueA[2][7] = 2;
  trueA[3][8] = 2;

  for (int l = 0; l < n; l++){
    for (int i = 0; i < n; i ++){
      for (int j = 0; j < k; j++){
        int val = mtx(i,j, l);
        if (val == 1){
          A[i][j] = 1;
        }
        if (val > 1){
          A[i][j] = 2;
        }
        // std::cout << mtx(i,j, l) << ' ';
      }
      // std::cout << '\n';
    }
    // std::cout << '\n';
  }

  // COMPARISON...NEED TO FIX THIS
  int total_rows = 0;
  for (int jt = 0; jt < k; jt++){
    for (int j = 0; j < k; j++){
      bool flag = 1;
      for (int i = 0; i < n; i++){
        if (A[i][j] != trueA[i][jt]){flag = 0;}
      }
      if (flag){total_rows++;}
    }
    // std::cout << '\n';
  }
  return total_rows;
}

int compareStateTree(gm::Solution sol){

  int correct = 0;

  int trueStateTree[5][5] = {{-1, 0, -2, -2, -2}, {-1, 0, 1, -2, -2}, {-1, 0, -2, 1, -2}, {-1, 0, -2, -2, 0}, {-1, 0, -2, -2, -2}};

  int nState = sol.S()[0].numVertices();
  int nSNV = sol.S().size();

  for (int i = 0; i < nSNV; i++){
    for (int j = 0; j < nState; j++){
      if (sol.S()[i].parent(j) == trueStateTree[i][j]){
        correct++;
      }
    }
  }

  return correct;
}

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <SPRUCE_OUTPUT>" << std::endl;
    return 1;
  }

  std::ifstream inFile(argv[1]);

  if (!inFile.good())
  {
    std::cerr << "Error: could not open '" << argv[1] << "' for reading" << std::endl;
    return 1;
  }

  gm::SolutionSet sols;
  inFile >> sols;
  inFile.close();

  int solCount = sols.solutionCount();
  std::cout << "num of sols: " << solCount  << '\n';

  // std::string edgeCompar = "trees/test_comparison.csv";
  // std::filebuf tf;
  // tf.open (edgeCompar,std::ios::out);
  // std::ostream os(&tf);
  //
  // std::string stateCompar = "trees/test_state_comparison.csv";
  // std::filebuf sf;
  // sf.open (stateCompar,std::ios::out);
  // std::ostream ot(&sf);

  std::string distCompar = "trees/" + writeOutputFile(argv[1], "_dist_comparison.csv");
  std::filebuf uf;
  uf.open (distCompar,std::ios::out);
  std::ostream ou(&uf);

  sols.assignDistancesByOccurenceCounts();

  for (int i = 0; i < solCount; i++){
    // std::cout << "solnum: " << i << '\n';
    int sol_idx = i;
    std::string file_add = "_sol_" + std::to_string(sol_idx) + "_parsed.csv";
    // outputFile = outputFile + add;
    std::string parsedFile = "data/" +  writeOutputFile(argv[1], file_add);
    std::string readFile = "data/" +  writeOutputFile(argv[1], "_reads.csv");

    returnVAF(sols.solution(sol_idx), parsedFile);
    writeParamFile(getReadNum(argv[1]), readFile);
    gm::PerfectPhyloMatrix mtx = sols.solution(sol_idx).A();

    // os << sol_idx << "," << comparePhyloMatrix(mtx) << "\n";
    // ot << sol_idx << "," << compareStateTree(sols.solution(sol_idx)) << "\n";
    ou << sol_idx << "," << sols.solution(sol_idx).distance() << "\n";

  }
  // tf.close();
  // sf.close();
  uf.close();
  // int sol_idx = 179308;
  // std::cout << "Sol: " << sol_idx << '\n';
  // writeSolutionEdgeList(sols.solution(sol_idx), "trees/test_edge_list.txt");
  // gm::PerfectPhyloMatrix mtx = sols.solution(sol_idx).A();
  //
  // std::cout << comparePhyloMatrix(mtx) << '\n';

  std::cout << '\n';
  // int n =
  // std::cout << sols.solution(sol_idx).U() << '\n';
  // std::cout << '\n';

  // returnVAF(sols.solution(2).inferredF(), sols.solution(2), 1);

  // char a = ')';
  // std::vector<std::string> test_vector = split("(test,vector)", a);
  // std::cout << test_vector[0];//sols.solution(1).U();

  std::cout << '\n';

  return 0;
}
