/*
 * Copyright (C) 2005-2018 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "otbObiaStreamUtils.h"

#include <vector>

int otbObiaStreamUtils()
{
  std::vector<unsigned long> myVector;
  myVector.push_back(0);
  myVector.push_back(1);
  myVector.push_back(2);
  myVector.push_back(3);

  size_t vectorOffset = otb::obia::stream_offset(myVector);

  double myDouble = 100.;

  size_t doubleOffset = otb::obia::stream_offset(myDouble);

  unsigned short myShort = 12;

  size_t shortOffset = otb::obia::stream_offset(myShort);

  std::vector<char> myStream(vectorOffset+doubleOffset+shortOffset);
  
  size_t offset = 0;

  otb::obia::to_stream(myStream,myDouble,offset);
  otb::obia::to_stream(myStream,myVector,offset);
  otb::obia::to_stream(myStream,myShort,offset);

  offset = 0;
  
  double decodedDouble;
  std::vector<unsigned long> decodedVector;
  unsigned short decodedShort;

  otb::obia::from_stream(myStream,decodedDouble,offset);
  otb::obia::from_stream(myStream,decodedVector,offset);
  otb::obia::from_stream(myStream,decodedShort,offset);

  return EXIT_SUCCESS;
}
