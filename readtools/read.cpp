#include <hps/src/hps.h>
#include <vector>
#include <fstream>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cout << "Usage: [executable_name] indices#.dat values#.dat\n";
    std::exit(EXIT_FAILURE);
  }
  
  {
    std::ifstream indices_file(argv[1], std::ifstream::binary);
    auto indices = hps::from_stream<std::vector<long int>>(indices_file);
    for (const auto& index: indices) std::cout<<index<<" ";
    std::cout<<"\n";
  }
  {
    std::ifstream values_file(argv[2], std::ifstream::binary);
    auto values = hps::from_stream<std::vector<double>>(values_file);
    for (const auto& value: values) std::cout<<value<<" ";
    std::cout<<"\n";
  }
  
}
