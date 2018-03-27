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

double returnStates(gm::Solution sol, int vectorIdx){
  //
  int numVertices = sol.S()[vectorIdx].numVertices();
  double states [3][numVertices];

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

  return 0;//states; // change this to return array
}

// use this for observed and/or inferred
int returnVAF(gm::RealTensor F, gm::Solution sol, int vectorIdx, std::string outputFile){

  // GET STATES
  int numVertices = sol.S()[vectorIdx].numVertices();
  double states [3][numVertices];

  for (int i = 0; i < numVertices; i++){
    std::string label = sol.S()[vectorIdx].label(i);
    std::vector<std::string> labels = parseLabel(label);
    // std::string a = labels[0];
    // double c = atof(a.c_str());
    states[0][i] = atof(labels[0].c_str());
    states[1][i] = atof(labels[1].c_str());
    states[2][i] = atof(labels[2].c_str());

    // std::cout << states[0][i] << " " << states[1][i] << " " << states[2][i] << " | ";
  }

  // GET VAF

  int k = F.k(); //
  int m = F.m(); //
  int n = F.n(); //

  std::cout << "count comparison " << numVertices << " " << k << " " << m << " " << n << '\n';
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
        total_r = total_r + ((states[0][h] + states[1][h])*freq);
        total_v = total_v + (states[2][h]*freq);
      }

      calculatedVAF[i][j] = total_v/total_r;
    }
  }


  // std::cout << "calculated: " << calculatedVAF << '\n';

  std::filebuf fb;
  fb.open (outputFile,std::ios::out);
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
    for (int j = 0; j < n; j++){

      os << calculatedVAF[i][j] << ",";
    }
    os << "\n";
  }


  fb.close();

  // return tensor of VAFs


  return 0;
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

  // std::filebuf fb;
  // fb.open ("test.txt",std::ios::out);
  // std::ostream os(&fb);
  //
  // sols.solution(1).S()[1].writeEdgeList(os);
  //
  // os << "Test sentence\n";
  // fb.close();

  // returnStates(sols.solution(1), 1);

  std::cout << '\n';

  std::string parsedFile = writeOutputFile(argv[1], "_parsed.csv");
  std::string readFile = writeOutputFile(argv[1], "_reads.csv");

  returnVAF(sols.solution(2).observedF(), sols.solution(2), 1, parsedFile);
  writeParamFile(getReadNum(argv[1]), readFile);

  std::cout << '\n';

  // returnVAF(sols.solution(2).inferredF(), sols.solution(2), 1);

  // char a = ')';
  // std::vector<std::string> test_vector = split("(test,vector)", a);
  // std::cout << test_vector[0];//sols.solution(1).U();

  std::cout << '\n';

  return 0;
}
