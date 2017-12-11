
#include "dipole_module.hpp"
#include <iostream>

DipoleModule::DipoleModule() {

}


void DipoleModule::initialise(const gmx::ConstArrayRef<rvec>& waterCoordinates) {
    
    this->initDipoles_ = std::vector<rvec>(waterCoordinates.size()/3);
    for (unsigned int i = 0; i < waterCoordinates.size()-3; i += 3) {
	rvec OH1, OH2, mu;
	rvec_sub(waterCoordinates.at(i+1), waterCoordinates.at(i), OH1);
	rvec_sub(waterCoordinates.at(i+2), waterCoordinates.at(i), OH2);
	rvec_add(OH1, OH2, mu);
        copy_rvec(mu, initDipoles_.at(i/3));
    }
    
}


void DipoleModule::analyseFrame(const gmx::ConstArrayRef<rvec>& waterCoordinates) {
    std::vector<real> currentCos(waterCoordinates.size()/3);
    for (unsigned int i = 0; i < waterCoordinates.size()-3; i += 3) {
	rvec OH1, OH2, mu;
	rvec_sub(waterCoordinates.at(i+1), waterCoordinates.at(i), OH1);
	rvec_sub(waterCoordinates.at(i+2), waterCoordinates.at(i), OH2);
	rvec_add(OH1, OH2, mu);
	currentCos.at(i/3) = cos_angle(this->initDipoles_.at(i/3), mu);
    }
    this->dipolesFrames_.push_back(currentCos);
}


real DipoleModule::averageFrame(const std::vector<int>& selection = std::vector<int>()) {
    if (selection.empty()) {
	real output = 0.0;
	std::vector<real> frame = this->dipolesFrames_.back();
	for (auto &cos : frame) {
	    output += cos;
	}
	return output/frame.size();
    }
    else {
	real output = 0.0;
	std::vector<real> frame = this->dipolesFrames_.back();
	for (auto &id : selection) {
	    output += frame.at(id/3);
	}
	return output/selection.size();
    }
}


std::vector<real> DipoleModule::average(const std::vector<int>& selection = std::vector<int>()) {
    if (selection.empty()) {
	std::vector<real> output;
	for (auto &frame : dipolesFrames_) {
	    real corr = 0.0;
	    for (auto &cos : frame) {
		corr += cos;
	    }
	    output.push_back(corr/selection.size());
	}
	return output;
    }
    else {
	std::vector<real> output;
	for (auto &frame : dipolesFrames_) {
	    real corr = 0.0;
	    for (auto &id : selection) {
		corr += frame.at(id/3);
	    }
	    output.push_back(corr/selection.size());
	}
	return output;
    }

}
