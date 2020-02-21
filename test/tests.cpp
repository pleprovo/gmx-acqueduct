
// #include <iostream>
// #include <fstream>
// #include <cmath>
// #include <future>

// #include "cgal.hpp"
// #include "boost/version.hpp"

// #include <chrono>

// std::vector<int> sub (std::vector<int> &nums, int start, int stop)
// {
//     std::vector<int> list;
//     list.reserve(stop-start+1);
//     for ( int i = start; i < stop; ++i )
//     {
// 	list.push_back(i*i);
//     }
//     return list;
// }
 

int main (int argc, char *argv[])
{
    // using namespace boost;
    // std::cout << "Using BOOST : " << ((BOOST_VERSION / 100) % 1000) << std::endl;
    // std::cout << " --- TESTS --- " << std::endl;
    // std::cout << " --- Test Distance Switch : " << std::flush;
    
    // // Testing radius switch
    // std::ofstream oss1;
    // oss1.open("switch_radius.dat");
    // for (double r = 2.0; r < 7.0; r+=0.1) {
    // 	oss1 << r << " ";	
    // }
    // oss1 << "\n";
    // for (double r = 2.0; r < 7.0; r+=0.1) {
    // 	oss1 << cgal::switch_function_radius(r, 5.5, 6.5) << " ";
    // }
    // oss1.close();
    // std::cout << " DONE --- " << std::endl;
    
    // std::cout << " --- Test Angle Switch : " << std::flush;
    
    // // Testing angle switch
    // std::ofstream oss2;
    // oss2.open("switch_angle.dat");
    // for (double cos = 0.0; cos < 1.0; cos+=0.001) {
    // 	oss2 << cos << " ";	
    // }
    // oss2 << "\n";
    // for (double cos = 0.0; cos < 1.0; cos+=0.001) {
    // 	oss2 << 1-cgal::switch_function_angle(cos, 0.0301, 0.25) << " ";	
    // }
    // oss2.close();
    // std::cout << " DONE --- " << std::endl;
    
    // // Testing Hydrogen Bond potential
    // std::cout << " --- Test Hydrogen Bond Potential : " << std::flush;
    // std::ofstream oss;
    // oss.open("potential.dat");
    // for (double r = 2.5; r < 7.0; r+=0.001) {
    // 	for (double theta = 0.0; theta < 1.0; theta+=0.01) {
    // 	    oss << cgal::computeEnergy(r, theta) << " " ;
    // 	}
    //     oss << "\n";
    // }
    // oss.close();
    // std::cout << " DONE --- " << std::endl;

    // std::cout << " --- Test future : " << std::flush;
    // std::vector<int> nums;
    // for (int i = 0; i < 1000000; ++i)
    // {
    // 	nums.push_back(i);
    // }

    // std::chrono::system_clock::time_point start0 = std::chrono::system_clock::now();
    
    // std::vector<int> resultsingle = sub(nums, 0, nums.size());
    
    // auto stop0 = std::chrono::system_clock::now();
    
    // std::cout << "\n" << "Single done in "
    // 	      << std::chrono::duration_cast<std::chrono::milliseconds>(stop0 - start0).count()
    // 	      << "ms" <<std::flush;

    // auto start1 = std::chrono::system_clock::now();
    
    // std::future<std::vector<int>> fut = std::async(std::launch::async,
    // 						   sub,
    // 						   std::ref(nums), 0, nums.size()/2);
    // std::vector<int> chunk2 = sub(nums, nums.size()/2, nums.size());
    // fut.wait();
    // std::vector<int> chunk1 = fut.get();
    // std::vector<int> results;
    // results.reserve(chunk1.size()+chunk2.size());
    // results.insert(results.end(), chunk1.begin(), chunk1.end());
    // results.insert(results.end(), chunk2.begin(), chunk2.end());
    
    // auto stop1 = std::chrono::system_clock::now();

    // std::cout << "\n" << "Parallel done in "
    // 	      << std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1).count()
    // 	      << "ms" <<std::flush;
    
    // std::cout << " DONE --- " << std::endl;
    // std::cout << 101 / 4 << " " << 101 % 4 << std::endl;
    return 0;
}
