#ifndef SITE_HPP
#define SITE_HPP

#include "cgal.hpp"


class Site
{
public:
    Site(unsigned id,
	 unsigned atmId,
	 unsigned resId,
	 SiteType type,
	 std::string atmName,
	 std::string resName);

    int& getId() const
	{
	    return mId;
	};

    int& getAtmId() const
	{
	    return mAtmId;
	};

    int& getResId() const
	{
	    return mResId;
	};

    SiteType& getType() const
	{
	    return mType;
	};

    std::string& getAtmName() const
	{
	    return mdAtmName;
	};

    std::string& getResName() const
	{
	    return mdResName;
	};

    void setSitePosition(const Point& position)
	{
	    mSitePosition = std::make_shared(position);
	};

    Point_ptr& getSitePosition() const
	{
	    return mSitePosition;
	};
    
    void addHydrogenPosition(const Point_ptr& position)
	{
	    mHydrogenPositions.push_back(std::make_shared(position));
	};

    std::vector<Point_ptr>& getHydrogenPositions() const
	{
	    return mHydrogenPosition;
	};
    
private:
    unsigned mId;
    unsigned mAtmId;
    unsigned mResId;
    SiteType mType;
    std::string mAtmName;
    std::string mResName;

    Point_ptr mSitePosition;
    std::vector<Point_ptr> mHydrogenPositions;
};

#endif
