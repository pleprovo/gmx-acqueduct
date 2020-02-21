
#include "alpha-shape-surface.hpp"

#include <future>

void AlphaShapeSurface::make(std::vector<Point>& points)
{
    CGAL::spatial_sort(points.begin(), points.end());
    dt_ = std::make_unique<Alpha_shape_3>(points.begin(), points.end());
    dt_->set_alpha(alpha_);
}

float AlphaShapeSurface::area()
{
    return 0.0;
}

float AlphaShapeSurface::volume()
{
    float volume = 0.0;

    for (Alpha_shape_3::Finite_cells_iterator it = dt_->finite_cells_begin();
	 it != dt_->finite_cells_end(); ++it)
    {
	if (dt_->classify(it) == Alpha_shape_3::INTERIOR)
	{
	    volume += dt_->tetrahedron(it).volume();
	}
    }
    
    return volume;
}

std::vector<int> locate_async(int start,
			      int stop,
			      const std::vector<Point>& positions,
			      const std::unique_ptr<Alpha_shape_3>& dt)
{
    std::vector<int> locatedPositions;
    for ( int i = start; i < stop; i++ )
    {
        if (dt->classify(positions.at(i)) == Alpha_shape_3::INTERIOR)
    	{
    	    locatedPositions.push_back(i);
    	}
    }
    return locatedPositions;
}

int AlphaShapeSurface::locate(const std::vector<Point>& positions,
			      std::vector<int>& locatedPositions)
{
    
    /* Make chunks */
    int n = std::thread::hardware_concurrency();
    int chunk_size = positions.size() / n;

    /* Make Futures */
    std::vector<std::future<std::vector<int>>> futures;

    /* Launch Futures */
    for ( int k = 0; k < n-1; k++ )
    {
    	futures.push_back(std::async(std::launch::async,
    				     locate_async,
    				     k*chunk_size,
				     (k+1)*chunk_size,
				     std::cref(positions),
				     std::cref(dt_)));
    }
    
    locatedPositions = locate_async(chunk_size*(n-1),
				    positions.size(), positions, dt_);

    /* Await results */
    for ( std::future<std::vector<int>>& fut : futures )
    {
    	fut.wait();
    }

    /* Get result */
    for ( std::future<std::vector<int>>& fut : futures )
    {
	fut.wait();
    	std::vector<int>chunk = fut.get();
    	locatedPositions.insert(locatedPositions.end(),
				std::make_move_iterator(chunk.begin()),
				std::make_move_iterator(chunk.end()));
    }
    
    
    // for ( int i = 0; i < positions.size(); i++)
    // {
    // 	if (dt_->classify(positions.at(i)) == Alpha_shape_3::INTERIOR)
    // 	{
    // 	    locatedPositions.push_back(i);
    // 	}
    // }
    return locatedPositions.size();
}

void AlphaShapeSurface::setAlphaValue(float alpha)
{
    alpha_ = alpha;
}
