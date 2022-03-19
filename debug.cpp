//
// Created by ubuntu on 11/3/21.
//
#include "debug.hpp"
#define STRLEN 256
void DEBUG_INFO(std::string info){
    std::cout << "\x1b[34m" << info << "\x1b[0m" << std::endl;
}

void MEM_USAGE()
{
    // extremely non-portable
//    std::cout << "\x1b[36m";
//    unsigned int pid = static_cast<unsigned>(getpid());
//    char statFileName[STRLEN];
//    memset(statFileName, 0, STRLEN*sizeof(char));
//    sprintf(statFileName, "/proc/%u/status", pid);
//
//    std::ifstream statFile(statFileName, std::ios::in);
//    if (!statFile.is_open()) {
//        std::cerr << "ERROR: DDgrp::MEM_USAGE: "
//             << " couldn't open file " << statFileName << "for read" << std::endl;
//        return;
//    }
//
//    double vmRss  = -1.;
//    double vmSize = -1.;
//    double vmPeak = -1.;
//
//    std::string buffer;
//    unsigned int pos = 0;
//
//    while (!statFile.eof()) {
//        pos = statFile.tellg();
//        getline(statFile, buffer);
//        if (buffer.find("VmPeak:") != std::string::npos) {
//            statFile.seekg(pos);
//            statFile >> buffer >> vmPeak >> buffer;
//        } else if (buffer.find("VmSize:") != std::string::npos) {
//            statFile.seekg(pos);
//            statFile >> buffer >> vmSize >> buffer;
//        } else if (buffer.find("VmRSS:") != std::string::npos) {
//            statFile.seekg(pos);
//            statFile >> buffer >> vmRss >> buffer;
//        }
//    }
//    statFile.close();
//
//    vmPeak /= 1024.;
//    vmSize /= 1024.;
//    vmRss  /= 1024.;
//
//    std::cout.setf(std::ios::fixed, std::ios::floatfield);
//    std::cout.precision(3);
//    std::cout << "Memory Usage: VmRSS = " << vmRss << "; VmSize = " << vmSize;
//    if (vmPeak > 0.) {
//        std::cout << "; VmPeak = " << vmPeak;
//    }
//    std::cout << " (MB)" << std::endl;
//    std::cout << "\x1b[0m";
}
