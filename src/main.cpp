

#include "waternetwork.hpp"

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WaterNetwork>(argc, argv);
}
