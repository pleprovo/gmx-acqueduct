

#include <string>
#include <utility>
#include <memory>

enum SiteType {DONOR, ACCEPTOR, WATER};

struct SiteInfo
{
    SiteType type;
    std::string name;
    unsigned id;
    unsigned int atmIndex;
    unsigned int resIndex;
    std::string resName;
    short unsigned int nbHydrogen;
};

using SiteInfo_ptr = std::shared_ptr<SiteInfo>;
using SiteInfoPair = std::pair<SiteInfo_ptr, SiteInfo_ptr>;

enum Direction {FORWARD, BACKWARD};

struct HydrogenBond
{
    SiteInfoPair sites;
    Direction direction;
    double length;
    double angle;
    double energy;
};
