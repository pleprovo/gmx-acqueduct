/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef WATERNETWORK_HPP
#define WATERNETWORK_HPP

#include <string>
#include <vector>
#include <fstream>

#include "cgal.hpp"

#include <gromacs/trajectoryanalysis.h>



/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
struct Site;

class WaterNetwork : public gmx::TrajectoryAnalysisModule
{
public:
    WaterNetwork();

    virtual void initOptions(gmx::IOptionsContainer          *options,
			     gmx::TrajectoryAnalysisSettings *settings) override;
    void optionsFinished(gmx::TrajectoryAnalysisSettings *settings) override;
    virtual void initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			      const gmx::TopologyInformation        &top) override;

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			      gmx::TrajectoryAnalysisModuleData *pdata) override;
    virtual void finishAnalysis(int nframes) override;
    virtual void writeOutput() override;

private:
    class ModuleData;
    
    std::string fnFilter_;
    std::string fnGraph_;
    
    double alphaValue_;
    
    double lengthOn_;
    double lengthOff_;
    
    double angleOn_;
    double angleOff_;
    
    gmx::Selection solventSel_;
    gmx::Selection alphaSel_;
    gmx::Selection proteinSel_;
    
    gmx::AnalysisData filterData_; 
    gmx::AnalysisData graphData_;

    std::unique_ptr<FindHydrogenBonds> computator_;
    
    std::vector<const Site> solventSites_;
    std::vector<const Site> proteinSites_;

    std::ofstream outputStream_;
    
    //! Topology
    const gmx::TopologyInformation *top_;
};

#endif



