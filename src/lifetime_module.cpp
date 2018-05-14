

#include "lifetime_module.hpp"


LifetimeModule::LifetimeModule () {

}

void LifetimeModule::initialise(int nb_frames, int nb_elements)
{
    this->presenceMatrix_.resize(nb_frames, std::vector<bool>(nb_elements, false));
    
}

void LifetimeModule::analyseFrame(int frame_nb, std::vector<int> selected)
{
    for (auto &id : selected) {
	this->presenceMatrix_.at(frame_nb).at(id) = true;
    }
}
