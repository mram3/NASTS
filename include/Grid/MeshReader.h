#pragma once

#include "Grid/Grid.h"

#include <string>

namespace MeshReader
{
    //for testing purposes
    Grid readGmshAscii22(const std::string& filePath);

    //Prime reader for main
    Grid readGmshAscii22();
}