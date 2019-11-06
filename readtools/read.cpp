#include <hps/src/hps.h>
#include <vector>
#include <fstream>
#include <cmath>

int main(int argc, char *argv[]) {

  if (argc == 2) {
    std::ifstream values_file(argv[1], std::ifstream::binary);
    std::cout<<"reading "<<argv[1]<<std::endl;
    auto values = hps::from_stream<std::vector<double>>(values_file);

    std::cout<<"No. nonzero entries: "<<values.size()<<std::endl;
    double normalization = 0., max_elem = 0., min_elem = 2.;
    size_t count1 = 0;
    double normalization1 =0.;
    for (const double value: values) {
      double abs_value = std::abs(value);
      if (abs_value > max_elem) max_elem = abs_value;
      else if (abs_value < min_elem) min_elem = abs_value;
      normalization += value * value;
      if (abs_value > 1e-4) {
        count1++;
        normalization1 += value * value;
      }
    }
    std::cout<<"Max magnitude element: "<<max_elem<<std::endl;
    std::cout<<"Min magnitude element: "<<min_elem<<std::endl;
    std::cout<<"Normalization: "<<normalization<<std::endl;
    std::cout<<"Above 1e-4: No. entries: "<<count1<<" normalization: "<<normalization1<<std::endl;
    std::cout<<"\n";

  } else {
    if (argc != 3) {
      std::cout << "Usage: [executable_name] indices#.dat values#.dat OR [excutable_name] values#.dat\n";
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
}
