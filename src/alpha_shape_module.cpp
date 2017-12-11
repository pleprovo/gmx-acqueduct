/*
 *  Class for finding point inside the Alpha Shape Surface
 *  Pierre Leprovost, University of Oulu, Biocenter Oulu
 *  1/10/2015
 */

#include "alpha_shape_module.hpp"


AlphaShapeModule::AlphaShapeModule() {}

AlphaShapeModule::~AlphaShapeModule() {}

AlphaShapeModule::AlphaShapeModule(const std::vector<Point_3> &pointVector, 
				   const float &alpha) 
{
    this->build(pointVector, alpha);
}

void AlphaShapeModule::build(const std::vector<Point_3> &pointVector,
				       const float &alpha)
{
    // std::shared_ptr<Delaunay_hierarchy> delaunayTesselation =
    // 	std::make_shared<Delaunay_hierarchy>(pointVector.begin(), pointVector.end());
    // alphaShape_ = std::make_shared<Alpha_shape_3>(*delaunayTesselation);
    alphaShape_ = std::make_shared<Alpha_shape_3>(pointVector.begin(), pointVector.end());
    alphaShape_->set_alpha(alpha);
}

void AlphaShapeModule::update(const std::vector<Point_3> &pointVector)
{
    // // Iterate over alpha shape vertex and update their position
    // int count=0;
    // for (Vb::Alpha_shape_vertices_iterator
    // 	     avit = alphaShape_->alpha_shape_vertices_begin(),
    // 	     avit_end = alphaShape->alpha_shape_vertices_end();
    // 	 avit!=avit_end; ++avit)
    // {
    // 	count++;
    // }
}

std::vector<int> AlphaShapeModule::locate(const std::vector<Point_3> &pointVector, 
					  const bool &location) const
{
    std::vector<int> pointLocated;
    Classification locationClassifier;

    if (location) {
        locationClassifier = Alpha_shape_3::INTERIOR;
    }
    else {
        locationClassifier = Alpha_shape_3::EXTERIOR;
    }
    
    for (unsigned int i = 0; i < pointVector.size(); i++) {
        if (alphaShape_->classify(pointVector.at(i)) == locationClassifier) {
            pointLocated.push_back(i);
        }        
    }

    return pointLocated;
}

float AlphaShapeModule::volume()
{
    float volume = 0;

    for(Alpha_shape_3::Finite_cells_iterator it = alphaShape_->finite_cells_begin(); 
	it != alphaShape_->finite_cells_end(); it++) {
        if(alphaShape_->classify(it)==Alpha_shape_3::INTERIOR) { 
	    Tetrahedron_3 tetr = alphaShape_->tetrahedron(it);
	    volume += tetr.volume();
	}
    }

    return volume;
}

void AlphaShapeModule::writeOFF(std::string &outputString) 
{
    std::vector<Facet> facets;
    std::stringstream outputStreamPositions;
    std::stringstream outputStreamIndices;

    alphaShape_->get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    alphaShape_->get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::SINGULAR);	
    std::size_t facetsNumber = facets.size();

    outputStreamPositions << "OFF "<< 3*facetsNumber << " " << 
	facetsNumber << "\n";

    for (std::size_t i=0;i<facetsNumber;i++) { 
	if ( alphaShape_->classify(facets[i].first)!=Alpha_shape_3::EXTERIOR) { 
	    facets[i]=alphaShape_->mirror_facet(facets[i]);
	}
        CGAL_assertion(alphaShape_->classify(facets[i].first)==Alpha_shape_3::EXTERIOR);
	
        int indices[3]={(facets[i].second+1)%4,(facets[i].second+2)%4,(facets[i].second+3)%4};

        if (facets[i].second%2==0) {
	    std::swap(indices[0], indices[1]);
        }

        outputStreamPositions << facets[i].first->vertex(indices[0])->point() << "\n" <<
	    facets[i].first->vertex(indices[1])->point() << "\n" <<
	    facets[i].first->vertex(indices[2])->point() << "\n"; 
        outputStreamIndices << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
    }
    
    outputString += outputStreamPositions.str();
    outputString += outputStreamIndices.str();
}
