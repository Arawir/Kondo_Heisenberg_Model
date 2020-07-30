#include "experiment1.h"


int main(int argc, char *argv[])
{
    Experiment1::readExperimentParameters(argc,argv);
    Experiment1::writeExperimentParameters();
    Experiment1::run();

    return 0;
}
