//
//  main.cpp
//  Canny Edge Detector
//
//  Created by Hasan Akgün on 21/03/14.
//  Copyright (c) 2014 Hasan Akgün. All rights reserved.
//

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "canny.h"

using namespace std;
using namespace cimg_library;

int main()
{
    const char* filePath = "./test/lena.jpg"; //Filepath of input image
    canny cny(filePath);
        
    return 0;
}

