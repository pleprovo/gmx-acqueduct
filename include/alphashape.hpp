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


#ifndef ALPHASHAPE_HPP
#define ALPHASHAPE_HPP

#include <string>
#include <vector>

#include <gromacs/analysisdata/modules/lifetime.h> 
#include <gromacs/trajectoryanalysis.h>


/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */

class AlphaShape : public gmx::TrajectoryAnalysisModule
{
public:
    AlphaShape();

    virtual void initOptions(gmx::Options                    *options,
			     gmx::TrajectoryAnalysisSettings *settings);
    virtual void initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			      const gmx::TopologyInformation        &top);

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			      gmx::TrajectoryAnalysisModuleData *pdata);

    virtual void finishAnalysis(int nframes);
    virtual void writeOutput();

private:
    //class ModuleData;
    
    gmx::SelectionList               selectionListAlpha_;
    gmx::Selection                   selectionWater_;
    
    std::string                      filenameWater_;
    std::string                      filenameVolume_;
    std::string                      filenameLifetime_;
    std::string                      filenameStats_;
    
    double                           alphaValue_;
    double                           numFrameValue_;
    
    gmx::AnalysisData                waterData_;
    gmx::AnalysisData                volumeData_;

    gmx::AnalysisData                lifetimeData_;
    gmx::AnalysisDataLifetimeModulePointer  lifetimeModule_;
};

#endif 



