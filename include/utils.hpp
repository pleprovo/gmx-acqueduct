

enum SiteType {DONOR, ACCEPTOR, WATER};

struct Site
{
    SiteType type;
    std::string name;
    unsigned id;
    unsigned int atmIndex;
    unsigned int resIndex;
    std::string resName;
    short unsigned int nbHydrogen;
};

using Site_ptr = std::shared_ptr<Site>;
using SitePair = std::pair<Site_ptr, Site_ptr>;

enum Direction {FORWARD, BACKWARD};

struct HydrogenBond
{
    SitePair sites;
    Direction direction;
    double length;
    double angle;
    double energy;
};
