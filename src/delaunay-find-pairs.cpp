

#include "delaunay-find-pairs.hpp"

std::vector<std::pair<int, int>> DelaunayFindPairs::find(std::vector<Point>& points)
{
    std::vector<std::size_t> indices;
    indices.reserve(points.size());
    
    std::copy(boost::counting_iterator<std::size_t>(0),
	      boost::counting_iterator<std::size_t>(points.size()),
	      std::back_inserter(indices));
    
    CGAL::spatial_sort( indices.begin(),
    			indices.end(),
    			Search_traits_3(CGAL::make_property_map(points)) );
    DelaunayWithInfo di(boost::make_zip_iterator(boost::make_tuple(
						     points.begin(), indices.begin() )),
			boost::make_zip_iterator(boost::make_tuple(
						     points.end(), indices.end())));
    CGAL_assertion( di.number_of_vertices() == points.size() );

    std::vector<std::pair<int, int>> pairs;
    for (DelaunayWithInfo::Finite_edges_iterator eit = di.finite_edges_begin();
	 eit != di.finite_edges_end(); ++eit)
    {
	pairs.push_back(std::make_pair<int, int>(
			    eit->first->vertex(eit->second)->info(),
			    eit->first->vertex(eit->third)->info()));
    }
    return pairs;
}
