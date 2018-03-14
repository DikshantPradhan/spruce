#include "utils.h"
#include "solutionset.h"
#include <fstream>

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
  
  std::cout << sols;
  
  return 0;
}
