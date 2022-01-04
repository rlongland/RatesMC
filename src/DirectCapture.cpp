/* ======================================================================
   DirectCapture.cpp
   Author: R. Longland
   Date: 2018-11-13
   
   Description: Direct Capture description
   ======================================================================
*/
#include <iostream>

#include "DirectCapture.h"


DirectCapture::DirectCapture(){
  std::cout << "This is direct capture!" << std::endl;
}
DirectCapture::~DirectCapture(){
  std::cout << "DirectCapture destructor!" << std::endl;
}
