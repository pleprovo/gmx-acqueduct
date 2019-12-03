
#include "brute-find-hydrogen-bonds.hpp"

#include <future>

std::vector<HydrogenBond> find(const std::vector<Site>& sites)
{

    /* Make chunks */
    int n = std::thread::hardware_concurrency();
    int chunk_size = nodes.size() / n;

    /* Make Futures */
    std::vector<std::future<std::vector<HydrogenBond>>> futures;
	    
    /* Launch Futures */
    for ( int k = 0; k < n-1; k++ )
    {
	futures.push_back(std::async(std::launch::async,
				     this->findImplementation,
				     std::ref(nodes),
				     chunk_size*k,
				     chunk_size*(k+1)));
    }
    std::vector<HydrogenBond> hbsVector;
    hbsVector = this->findImplementation(nodes, chunk_size*(n-1), nodes.size());

    /* Await results */
    for ( auto &fut : futures )
    {
	fut.wait();
    }

    /* Get result */
    for ( auto &fut : futures )
    {
	std::vector<HydrogenBond> chunk = fut.get();
	hbsVector.insert(hbsVector.end(),
			 std::make_move_iterator(chunk.begin()),
			 std::make_move_iterator(chunk.end()));
    }
	    
    return hbsVector;
}

std::vector<HydrogenBond> findImplementation(const std::vector<Site>& sites,
					     unsigned int start, unsigned int stop) 
{
    std::vector<HydrogenBond> hbsVector;

    for (unsigned int i = start; i < stop; ++i)
    {
	for ( unsigned int j = 1; j < i; ++j )
	{
	    bool isedge = false;

	    cgal::Vector_3 OO = *nodes.at(j).point - *nodes.at(i).point;
	    float distance = CGAL::sqrt(OO*OO);
	    OO /= distance;
	    distance *= 10.0;
    
	    if ( nodes.at(i).site->type != nodes.at(j).site->type ) { isedge = true; }
	    if ( nodes.at(i).site->type == WATER and nodes.at(j).site->type == WATER )
	    {
		isedge = true;
	    }
	    if ( distance > mDistanceCutoff ) { isedge = false; } 
	    if ( nodes.at(i).hydrogens.empty() and nodes.at(j).hydrogens.empty() )
	    {
		isedge = false;
	    } 

	    if ( isedge )
	    {
		std::vector<float> energies;
		std::vector<float> coss;
		
		if ( !nodes.at(i).hydrogens.empty() )
		{
		    for ( cgal::Point_ptr& hydrogen : nodes.at(i).hydrogens )
		    {
			cgal::Vector_3 DH = *hydrogen - *nodes.at(i).point;
			DH /= CGAL::sqrt(DH*DH);
			float proj = DH*OO;
			energies.push_back(cgal::computeEnergy(distance, proj));
			coss.push_back(proj);
		    }
		}
		
		if ( !nodes.at(j).hydrogens.empty() )
		{
		    for ( cgal::Point_ptr& hydrogen : nodes.at(j).hydrogens )
		    {
			cgal::Vector_3 DH = *hydrogen - *nodes.at(j).point;
			DH /= CGAL::sqrt(DH*DH);
			float proj = DH*(-1.0*OO);
			energies.push_back(cgal::computeEnergy(distance, proj));
			coss.push_back(proj);
		    }
		}
		
		int pos = -1;
		for ( unsigned int i = 0; i < coss.size(); ++i )
		{
		    if ( coss.at(i) > mAngleCutoff )
		    {
			cosmax = coss.at(i);
			pos = i;
		    }
		}

		if ( pos != -1 )
		{
		    results.push_back(HydrogenBond{
			    std::make_pair(nodes.at(i).site, nodes.at(j).site),
				FORWARD,
				distance,
				coss.at(pos),
				energies.at(pos)});
		}
	    }
	}
    }
}
