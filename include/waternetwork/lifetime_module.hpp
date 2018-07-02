

#ifndef LIFETIME_MODULE_HPP
#define LIFETIME_MODULE_HPP

#include <vector>

class LifetimeModule
{
public:
    LifetimeModule();
    void initialise(int nb_frames, int nb_elements);
    void analyseFrame(int frame_nb, std::vector<int> selected);

private:
    std::vector<std::vector<bool> > presenceMatrix_;
};

#endif
