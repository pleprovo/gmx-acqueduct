
#include "brute-find-pairs.hpp"

#include <future>


std::vector<std::pair<int, int>>
BruteFindPairs::find(std::vector<Point>& points)
{
    std::vector<std::pair<int, int>> pairs;
    
    /* Make chunks */
    // int n = std::thread::hardware_concurrency();
    // int chunk_size = sites.size() / n;
    
    /* Make Futures */
    // std::vector<std::future<std::vector<HydrogenBond>>> futures;
    
    float cutoff = 6.5;
    cutoff *= cutoff;
    for (unsigned int i = 0; i < points.size(); ++i)
    {
	for (unsigned int j = i+1; j < points.size(); ++j)
	{
	    if (10.0*CGAL::squared_distance(points.at(i), points.at(j)) < cutoff)
	    {
		pairs.push_back(std::make_pair<int, int>(i, j));
	    }
	}
    }

    return pairs;
//     // /* Make chunks */
//     // int n = std::thread::hardware_concurrency();
//     // int chunk_size = sites.size() / n;

//     // /* Make Futures */
//     // std::vector<std::future<std::vector<HydrogenBond>>> futures;

//     // /* Launch Futures */
//     // for ( int k = 0; k < n-1; k++ )
//     // {
//     // 	futures.push_back(std::async(std::launch::async,
//     // 				     [],
//     // 				     std::ref(sites),
//     // 				     chunk_size*k,
//     // 				     chunk_size*(k+1)));
//     // }
    
//     // std::vector<HydrogenBond> hbsVector;
//     // hbsVector = findImpl(sites, chunk_size*(n-1), sites.size());

//     // /* Await results */
//     // for ( auto &fut : futures )
//     // {
//     // 	fut.wait();
//     // }

//     // /* Get result */
//     // for ( auto &fut : futures )
//     // {
//     // 	std::vector<HydrogenBond> chunk = fut.get();
//     // 	hbsVector.insert(hbsVector.end(),
//     // 			 std::make_move_iterator(chunk.begin()),
//     // 			 std::make_move_iterator(chunk.end()));
//     // }
//     std::vector<HydrogenBond> hbsVector;

//     for (unsigned int i = 0; i < sites.size(); ++i)
//     {
//     	for ( unsigned int j = 1; j < i; ++j )
//     	{
//     	    bool isedge = false;

//     	    Vector_3 OO = sites.at(j).first - sites.at(i).first;
//     	    float distance = CGAL::sqrt(OO*OO);
//     	    OO /= distance;
//     	    distance *= 10.0;
    
//     	    if ( sites.at(i).second->info->type != sites.at(j).second->info->type )
//     	    {
//     		isedge = true;
//     	    }
//     	    if ( sites.at(i).second->info->type == WATER and
//     		 sites.at(j).second->info->type == WATER )
//     	    {
//     		isedge = true;
//     	    }
//     	    if ( distance > distanceCutoff_ ) { isedge = false; } 
//     	    if ( sites.at(i).second->hydrogens.empty() and
//     		 sites.at(j).second->hydrogens.empty() )
//     	    {
//     		isedge = false;
//     	    } 

//     	    if ( isedge )
//     	    {
//     		std::vector<float> coss;
		
//     		if ( !sites.at(i).second->hydrogens.empty() )
//     		{
//     		    for ( Point_ptr& hydrogen : sites.at(i).second->hydrogens )
//     		    {
//     			Vector_3 DH = *hydrogen - sites.at(i).first;
//     			DH /= CGAL::sqrt(DH*DH);
//     			float proj = DH*OO;
//     			coss.push_back(proj);
//     		    }
//     		}
		
//     		if ( !sites.at(j).second->hydrogens.empty() )
//     		{
//     		    for ( Point_ptr& hydrogen : sites.at(j).second->hydrogens )
//     		    {
//     			Vector_3 DH = *hydrogen - sites.at(j).first;
//     			DH /= CGAL::sqrt(DH*DH);
//     			float proj = DH*(-1.0*OO);
//     			coss.push_back(proj);
//     		    }
//     		}
		
//     		int pos = -1;
//     		float cosmax = angleCutoff_;
//     		for ( unsigned int i = 0; i < coss.size(); ++i )
//     		{
//     		    if ( coss.at(i) > cosmax )
// 		    {
// 			cosmax = coss.at(i);
// 			pos = i;
// 		    }
// 		}

// 		if ( pos != -1 )
// 		{
// 		    hbsVector.push_back(HydrogenBond{
// 			    std::make_pair(sites.at(i).second->info,
// 					   sites.at(j).second->info),
// 				FORWARD,
// 				distance,
// 				coss.at(pos)});
// 		}
// 	    }
// 	}
//     }
//     return hbsVector;
}

