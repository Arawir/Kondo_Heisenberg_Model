#include "experiment1.h"


int main(int argc, char *argv[])
{
    readExperimentParameters(argc,argv);
    writeExperimentParameters();
    Experiment1::run();

    return 0;
}
